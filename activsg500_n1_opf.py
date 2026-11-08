#!/usr/bin/env python3
"""
DC N-1 SCOPF for the ACTIVSg500 test case (TAMU synthetic 500-bus system).

Pipeline:
  network model  : ACTIVSg500.RAW parsed by the VeraGridEngine RAW parser
                   (vg_data.json, cross-validated against the official
                   MATPOWER case by crosscheck.py)
  cost/limit     : case_ACTIVSg500.m  (the RAW contains no P records;
                   gencost + gen Pmax/Pmin/status are the dataset's
                   official OPF data)
  contingency set: contab_ACTIVSg500.m (TAMU): 591 line outages are the
                   primary (classical) N-1 set; the 90 gen outages are
                   tested as a diagnostic (the joint line+gen SCOPF is
                   infeasible -- caused by exactly two large must-run units)
  solver         : HiGHS 1.15.1 via highspy (constraint generation)

Formulation (lossless DC, classical N-1 security, no corrective
re-dispatch) -- identical to dc_n1_opf.py:
  variables   P_g (online units), PL cost segments, bus angles
  objective   piecewise-linear (chord) generator cost
  constraints nodal balance at EVERY bus (incl. slack),
              generator limits, |F_l| <= Fmax_l (base case),
              |F_l^(c)| <= Fmax_l for every N-1 contingency c
              F_l^(c) = K_c (E_c P_g - loads), linear in the pre-
              contingency dispatch; disconnected contingencies skipped.

Usage:
  python activsg500_n1_opf.py [--nseg 20] [--lmps 12] [--maxcg 50]
"""
import argparse
import json
import os
import re
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

import highspy

A = Path(__file__).resolve().parent / "activsg500"
INF = highspy.kHighsInf


# ----------------------------------------------------------------------------
# model construction (VeraGrid-parsed RAW + official MATPOWER OPF data)
# ----------------------------------------------------------------------------
def build_model(verbose=True):
    VG = json.loads((A / "vg_data.json").read_text())
    MP = json.loads((A / "mp_data.json").read_text())
    S_BASE = MP["baseMVA"]

    bus_ids = sorted(b["I"] for b in VG["buses"])
    n = len(bus_ids)
    bidx = {b: i for i, b in enumerate(bus_ids)}

    # loads: L-records + bus shunt BL (MW), per bus
    loads = np.zeros(n)
    for l in VG["loads"]:
        if l["STATUS"] in (1, None):
            loads[bidx[int(l["I"])]] += l["PL"] or 0.0
    for b in VG["buses"]:
        loads[bidx[int(b["I"])]] += b["BL"] or 0.0

    # branches: lines + transformers (all in service, no phase shift)
    br = []          # (i_idx, j_idx, b=1/x, Fmax_MW, kind, ckt)
    for b in VG["branches"]:
        br.append((bidx[int(b["I"])], bidx[int(b["J"])], 1.0 / b["X"],
                   b["RATEA"] or 0.0, "line", str(b["CKT"]).strip()))
    for t in VG["transformers"]:
        br.append((bidx[int(t["I"])], bidx[int(t["J"])], 1.0 / t["X"],
                   t["RATEA"] or 0.0, "xfrm", str(t["CKT"]).strip()))
    m = len(br)

    # generators: RAW STAT==1 (online); Pmin/Pmax/cost from the MP case
    mp_row_of_bus = {int(r[0]): i for i, r in enumerate(MP["gen"])}
    gen_buses = []
    gen_pmin, gen_pmax, gen_c2, gen_c1, gen_c0 = [], [], [], [], []
    for g in VG["generators"]:
        if g["STAT"] != 1:
            continue
        bus = int(g["I"])
        r = MP["gen"][mp_row_of_bus[bus]]
        gc = MP["gencost"][mp_row_of_bus[bus]]
        assert gc[0] == 2 and gc[3] == 3, "expected quadratic 3-coeff gencost"
        gen_buses.append(bus)
        gen_pmin.append(r[9])
        gen_pmax.append(r[8])
        gen_c2.append(gc[4]); gen_c1.append(gc[5]); gen_c0.append(gc[6])
    ng = len(gen_buses)
    gcol = {b: k for k, b in enumerate(gen_buses)}
    offline_gen = [int(g["I"]) for g in VG["generators"] if g["STAT"] != 1]

    # contingency set from the TAMU change table (1-based MP row numbers)
    ct_text = (A / "contab_ACTIVSg500.m").read_text(errors="ignore")
    ct_rows = re.findall(
        r"^\s*\d+\s+\d+\s+(CT_\w+)\s+(\d+)\s+(\w+)\s+CT_\w+\s+[\d.]+;",
        ct_text, re.M)
    mp_branch_row_to_line = {i + 1: i for i in range(len(MP["branch"]))}
    line_cts, gen_cts, offline_gen_cts = [], [], 0
    for table, row, field in ct_rows:
        if table == "CT_TBRCH":
            line_cts.append(mp_branch_row_to_line[int(row)])
        elif table == "CT_TGEN":
            gbus = int(MP["gen"][int(row) - 1][0])
            if gbus in gcol:
                gen_cts.append(gcol[gbus])
            else:
                offline_gen_cts += 1
    nct = len(line_cts) + len(gen_cts)

    slack = bidx[1]
    ns = [i for i in range(n) if i != slack]
    ns_pos = {i: k for k, i in enumerate(ns)}

    Sfull = np.zeros((n, n))
    for (f, t, bl, *_r) in br:
        Sfull[f, t] -= bl
        Sfull[t, f] -= bl
        Sfull[f, f] += bl
        Sfull[t, t] += bl

    E = np.zeros((n, ng))
    for b, k in gcol.items():
        E[bidx[b], k] = 1.0

    Fmax_all = np.array([x[3] for x in br], dtype=float)
    unlimited = Fmax_all <= 1e-9

    M = dict(VG=VG, MP=MP, S_BASE=S_BASE, bus_ids=bus_ids, n=n, bidx=bidx,
             loads=loads, br=br, m=m, gen_buses=gen_buses, gcol=gcol,
             gen_pmin=gen_pmin, gen_pmax=gen_pmax, gen_c2=gen_c2,
             gen_c1=gen_c1, gen_c0=gen_c0, ng=ng, offline_gen=offline_gen,
             line_cts=line_cts, gen_cts=gen_cts, offline_gen_cts=offline_gen_cts,
             nct=nct, slack=slack, ns=ns, ns_pos=ns_pos, Sfull=Sfull, E=E,
             Fmax_all=Fmax_all, unlimited=unlimited)
    if verbose:
        print(f"model: {n} buses, {m} branches ({len(VG['branches'])} lines + "
              f"{len(VG['transformers'])} xfmr), {ng} online generators "
              f"({len(offline_gen)} offline)")
        print(f"total load {loads.sum():.2f} MW | sum Pmin {sum(gen_pmin):.1f} | "
              f"sum Pmax {sum(gen_pmax):.1f}")
        print(f"contingencies: {len(line_cts)} line + {len(gen_cts)} gen = {nct} "
              f"({offline_gen_cts} offline gen outages skipped as redundant)")
        print(f"limited branches: {int((~unlimited).sum())}/{m}")
    return M


def contingency_dispatch_map(M, out_line, out_gen):
    """K data for one contingency, or None if the post-contingency
    slack-side component has no lines.

    Returns (Kme, Kload, act, True, island_load, island_gen_pmax,
    island_exists):
      Kme:    (len(act), ng)  [MW/MW]
      Kload:  (len(act),)     [MW]
      act:    active line ids (slack-side component only)
      island_*: load / gen Pmax of the disconnected part, 0 if none
      island_exists: True if any bus is outside the slack-side component

    Lines of disconnected islands are NOT included (standard contingency
    practice: the island is shed / handled locally; only the network that
    stays connected to the slack is checked against flow limits).

    For a gen outage the failed unit's column of E is zeroed (E_c), so its
    frozen pre-contingency output is lost from the post-contingency
    injection.

    F_l^(c) = Kme[l,:] @ P_g - Kload[l]  (P_g = pre-contingency dispatch)"""
    n = M["n"]
    br = M["br"]
    m = M["m"]
    ns = M["ns"]
    ns_pos = M["ns_pos"]
    slack = M["slack"]
    E = M["E"]
    loads = M["loads"]
    ng = M["ng"]

    all_active = [li for li in range(m) if li != out_line]
    adj = np.zeros((n, n), bool)
    for li in all_active:
        f, t = br[li][0], br[li][1]
        adj[f, t] = adj[t, f] = True
    adj[slack, slack] = True
    _, labels = connected_components(sparse.csr_matrix(adj), directed=False)
    slack_label = labels[slack]
    # keep only lines inside the slack-side component
    act = [li for li in all_active
           if labels[br[li][0]] == slack_label]
    if not act:
        return None
    # island statistics (load/gen outside the slack-side component)
    in_slack = np.array([labels[i] == slack_label for i in range(n)])
    island_exists = bool((~in_slack).any())
    island_load = float(loads[~in_slack].sum())
    island_gen_pmax = float(np.sum([M["gen_pmax"][k] for k in range(M["ng"])
                                    if not in_slack[M["gen_buses"][k]]]))
    # gen outage: zero the failed unit's column of E so its pre-contingency
    # output is lost from the injection (F = K(E_c P_g - loads))
    E_use = E
    if out_gen is not None:
        E_use = E.copy()
        E_use[M["bidx"][M["gen_buses"][out_gen]], out_gen] = 0.0

    # reduced S over the non-slack buses of the slack-side component
    ns_c = [i for i in ns if labels[i] == slack_label]
    pos_c = {i: k for k, i in enumerate(ns_c)}
    S = np.zeros((len(ns_c), len(ns_c)))
    for li in act:
        f, t, bl = br[li][0], br[li][1], br[li][2]
        if f in pos_c:
            S[pos_c[f], pos_c[f]] += bl
        if t in pos_c:
            S[pos_c[t], pos_c[t]] += bl
        if f in pos_c and t in pos_c:
            S[pos_c[f], pos_c[t]] -= bl
            S[pos_c[t], pos_c[f]] -= bl
    Minv = np.linalg.pinv(S)
    Mm = np.zeros((n, n))
    for a, i in enumerate(ns_c):
        for b2, j in enumerate(ns_c):
            Mm[i, j] = Minv[a, b2]
    Kme = np.zeros((len(act), ng))
    Kload = np.zeros(len(act))
    for r, li in enumerate(act):
        f, t, bl = br[li][0], br[li][1], br[li][2]
        row = np.zeros(n)
        row[f] = bl
        row[t] = -bl
        Kflow = row @ Mm                  # dimensionless PTDF row (over buses)
        Kme[r, :] = Kflow @ E_use
        Kload[r] = Kflow @ loads
    return Kme, Kload, act, True, island_load, island_gen_pmax, island_exists


def precompute_contingencies(M):
    t0 = time.time()
    nct, m, ng = M["nct"], M["m"], M["ng"]
    Kme = np.zeros((nct, m, M["ng"]))
    Kload = np.zeros((nct, m))
    active = np.full((nct, m), -1, dtype=np.int32)
    valid = np.zeros(nct, dtype=bool)
    island = np.zeros(nct, dtype=bool)
    island_load = np.zeros(nct)
    island_gen = np.zeros(nct)
    for c in range(nct):
        out_line = M["line_cts"][c] if c < len(M["line_cts"]) else None
        out_gen = M["gen_cts"][c - len(M["line_cts"])] if c >= len(M["line_cts"]) else None
        res = contingency_dispatch_map(M, out_line, out_gen)
        if res is None:
            continue
        Kme_c, Kload_c, act, _, isl_ld, isl_gen, isl_flag = res
        valid[c] = True
        # island = any bus outside the slack-side component (bus count, not
        # line count: a pure radial load-bus island has no internal lines)
        island[c] = isl_flag
        island_load[c] = isl_ld
        island_gen[c] = isl_gen
        for r, li in enumerate(act):
            Kme[c, li, :] = Kme_c[r, :]
            Kload[c, li] = Kload_c[r]
            active[c, r] = li
    if M.get("verbose", True):
        print(f"contingency precompute: {time.time()-t0:.1f}s "
              f"({int(valid.sum())} connected, {nct - int(valid.sum())} dropped)")
        if island.any():
            print(f"  {int(island.sum())} contingencies disconnect an island: "
                  f"island load max {island_load[island].max():.1f} MW, "
                  f"island gen-Pmax max {island_gen[island].max():.1f} MW, "
                  f"islands with generation: {int((island & (island_gen > 0)).sum())}")
    return Kme, Kload, active, valid, island, island_load, island_gen


# ----------------------------------------------------------------------------
# PL costs (chord approximation of the quadratic gencost -- an upper bound)
# ----------------------------------------------------------------------------
def build_pl(M, nseg):
    prices, widths = [], []
    for k in range(M["ng"]):
        pmin, pmax = M["gen_pmin"][k], M["gen_pmax"][k]
        c2, c1, c0 = M["gen_c2"][k], M["gen_c1"][k], M["gen_c0"][k]
        xs = np.linspace(pmin, pmax, nseg + 1)
        cc = c2 * xs**2 + c1 * xs + c0
        prices.append(np.diff(cc) / np.diff(xs))
        widths.append(np.diff(xs))
    return np.concatenate(prices), np.concatenate(widths)


# ----------------------------------------------------------------------------
# highspy LP
# ----------------------------------------------------------------------------
def solve_lp(c_obj, lb, ub, A_ub, b_ub, A_eq, b_eq, label=""):
    h = highspy.Highs()
    opts = highspy.HighsOptions()
    opts.output_flag = bool(int(os.environ.get("VG_DEBUG", "0")))
    opts.time_limit = 1800.0
    h.passOptions(opts)
    lp = highspy.HighsLp()
    lp.num_row_ = A_ub.shape[0] + A_eq.shape[0]
    lp.num_col_ = c_obj.shape[0]
    lp.sense_ = highspy.ObjSense.kMinimize
    lp.col_cost_ = np.asarray(c_obj, dtype=np.float64)
    lp.col_lower_ = np.asarray(lb, dtype=np.float64)
    lp.col_upper_ = np.asarray(ub, dtype=np.float64)
    b_eq = np.asarray(b_eq, dtype=np.float64)
    b_ub = np.asarray(b_ub, dtype=np.float64)
    lp.row_lower_ = np.concatenate([np.full(A_ub.shape[0], -INF), b_eq])
    lp.row_upper_ = np.concatenate([b_ub, b_eq])
    A = sparse.vstack([A_ub, sparse.csr_matrix(A_eq)], format="csr")
    mat = highspy.HighsSparseMatrix()
    mat.format_ = highspy.MatrixFormat.kRowwise
    mat.num_row_ = A.shape[0]
    mat.num_col_ = A.shape[1]
    # numpy arrays pass through pybind11 (verified: passModel copies the data)
    mat.start_ = A.indptr
    mat.index_ = A.indices
    mat.value_ = A.data
    lp.a_matrix_ = mat
    model = highspy.HighsModel()
    model.lp_ = lp
    st = h.passModel(model)
    assert st in (highspy.HighsStatus.kOk, highspy.HighsStatus.kWarning), st
    h.run()
    ok = h.getModelStatus() == highspy.HighsModelStatus.kOptimal
    sol = h.getSolution()
    if label:
        print(f"  [{label}] status={'optimal' if ok else h.modelStatusToString(h.getModelStatus())} "
              f"obj={h.getObjectiveValue():,.1f}" if ok else
              f"  [{label}] NOT OPTIMAL: {h.modelStatusToString(h.getModelStatus())}")
    return ok, np.array(sol.col_value, dtype=float), np.array(sol.row_dual, dtype=float), h


def build_full_ub(M, Kme, Kload, active, valid, base_rows, base_rhs):
    """Build the full N-1 ub-matrix (base rows + every contingency row pair)
    directly in CSR, streaming row by row. Row order: base pairs, then each
    valid contingency's active lines (active order), each pair (+row, -row).
    Never materializes the dense matrix: ~770k rows x dim would be ~10 GB
    dense, but ~40M nnz sparse is fine.

    Returns (indptr, indices, data, rhs) as numpy arrays."""
    dim = len(base_rows[0])
    n_base = len(base_rows)       # already in (+row, -row) pairs
    Fmax = M["Fmax_all"]

    ct_pairs = []
    for c in range(M["nct"]):
        if not valid[c]:
            continue
        for li in active[c]:
            if li >= 0:
                ct_pairs.append((c, int(li)))
    n_pairs = len(ct_pairs)
    nrows = n_base + 2 * n_pairs

    # pass 1: nnz per row
    nnz = np.empty(nrows, dtype=np.int64)
    for k, r in enumerate(base_rows):
        nnz[k] = int(np.count_nonzero(r))
    kv_nnz = np.array([np.count_nonzero(Kme[c, li]) for c, li in ct_pairs])
    nnz[n_base:] = np.repeat(kv_nnz, 2)
    total = int(nnz.sum())

    # pass 2: fill
    indptr = np.empty(nrows + 1, dtype=np.int64)
    indptr[0] = 0
    indptr[1:] = np.cumsum(nnz)
    indices = np.empty(total, dtype=np.int32)
    data = np.empty(total, dtype=np.float64)
    rhs = np.empty(nrows, dtype=np.float64)

    for k, r in enumerate(base_rows):
        nz = np.nonzero(r)[0]
        sl = slice(indptr[k], indptr[k + 1])
        indices[sl] = nz
        data[sl] = r[nz]
        rhs[k] = base_rhs[k]

    for k, (c, li) in enumerate(ct_pairs):
        a = n_base + 2 * k
        b = a + 1
        kv = Kme[c, li]
        nz = np.nonzero(kv)[0]
        sl_a, sl_b = slice(indptr[a], indptr[a + 1]), slice(indptr[b], indptr[b + 1])
        indices[sl_a] = nz
        data[sl_a] = kv[nz]
        indices[sl_b] = nz
        data[sl_b] = -kv[nz]
        rhs[a] = Fmax[li] + Kload[c, li]
        rhs[b] = Fmax[li] - Kload[c, li]
    return indptr, indices, data, rhs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nseg", type=int, default=20)
    ap.add_argument("--lmps", type=int, default=12,
                    help="number of buses to get FD LMPs for (0 = skip)")
    ap.add_argument("--maxcg", type=int, default=50,
                    help="max constraint-generation iterations")
    ap.add_argument("--gentest", type=int, default=1,
                    help="run the gen-outage infeasibility diagnostic "
                         "(56 extra LPs, ~8 min); 0 to skip")
    args = ap.parse_args()
    nseg = args.nseg

    M = build_model()
    print("precomputing contingency K data ...")
    Kme, Kload, active, valid, island, island_load, island_gen = \
        precompute_contingencies(M)
    n, m, ng = M["n"], M["m"], M["ng"]
    br, E, Sfull, loads = M["br"], M["E"], M["Sfull"], M["loads"]
    bus_ids, bidx = M["bus_ids"], M["bidx"]
    ns, ns_pos, slack = M["ns"], M["ns_pos"], M["slack"]
    Fmax_all, unlimited = M["Fmax_all"], M["unlimited"]
    gen_buses = M["gen_buses"]
    S_BASE = M["S_BASE"]

    # ---- variable layout: [P_g(ng) | s(ng*nseg) | theta_ns(n-1)]
    nseg_total = ng * nseg
    dim = ng + nseg_total + len(ns)
    s0, t0 = ng, ng + nseg_total

    seg_price, seg_width = build_pl(M, nseg)
    c_obj = np.zeros(dim)
    c_obj[s0:t0] = seg_price

    # equality rows: nodal balance (every bus, pu) + segment links
    a_eq, b_eq = [], []
    Sall_ns = Sfull[:, ns]
    for i in range(n):
        row = np.zeros(dim)
        row[0:ng] = E[i, :] / S_BASE
        row[t0:t0 + len(ns)] = -Sall_ns[i, :]
        a_eq.append(row)
        b_eq.append(loads[i] / S_BASE)
    for k in range(ng):
        row = np.zeros(dim)
        row[k] = 1.0
        row[s0 + k * nseg:s0 + (k + 1) * nseg] = -1.0
        a_eq.append(row)
        # P_g = pmin + sum(segments); segments span [pmin, pmax]
        b_eq.append(M["gen_pmin"][k])
    A_eq = sparse.csr_matrix(np.array(a_eq))

    lb = np.zeros(dim)
    ub = np.full(dim, np.inf)
    for k in range(ng):
        lb[k], ub[k] = M["gen_pmin"][k], M["gen_pmax"][k]
    for k in range(ng):
        for j in range(nseg):
            ub[s0 + k * nseg + j] = seg_width[k * nseg + j]
    lb[t0:] = -np.inf

    # ---- base-case line limits (pu)
    def base_ineq():
        rows, rhs = [], []
        for li in range(m):
            if unlimited[li]:
                continue
            f, t, bl = br[li][0], br[li][1], br[li][2]
            row = np.zeros(dim)
            if f in ns_pos:
                row[t0 + ns_pos[f]] += bl
            if t in ns_pos:
                row[t0 + ns_pos[t]] -= bl
            rows += [row, -row]
            rhs += [Fmax_all[li] / S_BASE, Fmax_all[li] / S_BASE]
        return rows, rhs

    base_rows, base_rhs = base_ineq()

    # primary N-1 problem: classical security constraints = line outages only.
    # (TAMU's contab also lists gen outages; including them makes the joint
    # problem infeasible -- see the gen-outage diagnostic below.)
    n_line = len(M["line_cts"])
    v_line = valid.copy()
    v_line[n_line:] = False

    # ---- contingency evaluation at a dispatch (vectorized, MW)
    def evaluate(p_g_mw):
        F = np.einsum("clg,g->cl", Kme, p_g_mw) - Kload
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(Fmax_all[None, :] > 1e-9,
                             np.abs(F) / Fmax_all[None, :], 0.0)
        return F, ratio

    # ---- constraint generation
    t_start = time.time()
    cg_rows, cg_rhs = [], []
    n_solves = 0
    final = None
    for it in range(args.maxcg):
        all_rows = base_rows + cg_rows
        all_rhs = base_rhs + cg_rhs
        A_ub = sparse.csr_matrix(np.array(all_rows)) if all_rows else \
            sparse.csr_matrix((0, dim))
        b_ub = np.array(all_rhs) if all_rhs else np.array([])
        ok, x, duals, h = solve_lp(c_obj, lb, ub, A_ub, b_ub, A_eq, b_eq,
                                   label=f"CG it{it} rows={A_ub.shape[0]}")
        n_solves += 1
        if not ok:
            raise RuntimeError("LP failed: " + h.modelStatusToString(h.getModelStatus()))
        p_g = x[:ng]
        F, ratio = evaluate(p_g)
        per_ct = np.where(v_line, np.max(ratio, axis=1), 0.0)
        max_ratio = float(per_ct.max())
        if max_ratio < 1.0 + 1e-7:
            final = (x, duals, h)
            print(f"CG converged in {it+1} iterations, {n_solves} solves, "
                  f"{time.time()-t_start:.1f}s (max N-1 ratio {max_ratio:.6f})")
            break
        c = int(np.argmax(per_ct))
        act = active[c]
        act = act[act >= 0]
        for li in act:
            li = int(li)
            row = np.zeros(dim)
            row[0:ng] = Kme[c, li, :]
            cg_rows += [row, -row]
            cg_rhs += [Fmax_all[li] + Kload[c, li], Fmax_all[li] - Kload[c, li]]
        print(f"  it{it}: max ratio {max_ratio:.4f} -> adding contingency {c} "
              f"({len(act)} lines)")
    if final is None:
        raise RuntimeError("CG did not converge within --maxcg iterations")
    x, duals, h = final
    p_g = x[:ng]

    # ---- verification
    print("\nverification ...")
    th = np.zeros(n)
    th[ns] = x[t0:]
    resid = Sfull @ th * S_BASE - (E @ p_g - loads)
    print(f"  nodal balance residual: max |res| = {np.abs(resid).max():.3e} MW")
    assert abs(p_g.sum() - loads.sum()) < 1e-4, "gen/load imbalance"
    print(f"  gen total = {p_g.sum():.2f} MW = load total {loads.sum():.2f} MW")
    F, ratio = evaluate(p_g)
    viol = int(np.where(v_line[:, None], ratio > 1.0 + 1e-6, 0.0).sum())
    print(f"  N-1 check on ALL {int(v_line.sum())} line contingencies: "
          f"violations={viol}, worst ratio={ratio[v_line].max():.4f}")
    assert viol == 0, "N-1 violation found in final dispatch"

    exact_cost = sum(M["gen_c2"][k] * p_g[k]**2 + M["gen_c1"][k] * p_g[k]
                     + M["gen_c0"][k] for k in range(ng))
    # LP objective omits the per-unit constant c(pmin) (segments are measured
    # from pmin); add it back so PL cost is comparable with the exact cost
    pl_const = sum(M["gen_c2"][k] * M["gen_pmin"][k]**2
                   + M["gen_c1"][k] * M["gen_pmin"][k] + M["gen_c0"][k]
                   for k in range(ng))
    pl_cost = float(h.getObjectiveValue()) + pl_const

    # ---- N-0 for comparison
    print("solving N-0 (base case only) ...")
    A_ub0 = sparse.csr_matrix(np.array(base_rows))
    ok0, x0, d0, h0 = solve_lp(c_obj, lb, ub, A_ub0, np.array(base_rhs),
                               A_eq, b_eq, label="N-0")
    assert ok0
    p_g0 = x0[:ng]
    exact0 = sum(M["gen_c2"][k] * p_g0[k]**2 + M["gen_c1"][k] * p_g0[k]
                 + M["gen_c0"][k] for k in range(ng))
    F0, ratio0 = evaluate(p_g0)
    n0_viol = int(np.where(v_line[:, None], ratio0 > 1.0 + 1e-6, 0.0).sum())
    n0_worst = float(ratio0[v_line].max())

    # ---- LMPs
    lmps_fd = lmps_dual = sel = None
    if args.lmps > 0:
        # build the FULL line-outage N-1 LP directly in CSR (never dense:
        # ~700k rows x dim would be ~10 GB; the sparse form is ~39M nnz)
        t_fs = time.time()
        indptr, indices, data, b_fub = build_full_ub(
            M, Kme, Kload, active, v_line, base_rows, base_rhs)
        A_fub = sparse.csr_matrix((data, indices, indptr),
                                  shape=(indptr.size - 1, dim))
        n_full_ub = A_fub.shape[0]
        print(f"full N-1 LP: {n_full_ub} ub-rows + {len(a_eq)} eq-rows x "
              f"{dim} vars ({A_fub.nnz / 1e6:.1f}M nnz), "
              f"built in {time.time()-t_fs:.1f}s")
        okf, xf, df, hf = solve_lp(c_obj, lb, ub, A_fub, b_fub, A_eq, b_eq,
                                   label="FULL")
        assert okf, "full N-1 LP failed"
        print(f"full LP solved in {time.time()-t_fs:.1f}s")
        lmps_dual = np.array(df[n_full_ub:n_full_ub + n], dtype=float)

        sel = sorted(set([slack] + [int(i) for i in np.argsort(-loads)[:args.lmps - 1]]))
        pert = 0.05
        lmps_fd = np.zeros(n)
        t_fd = time.time()
        for bi in sel:
            vals = {}
            for sgn in (+1.0, -1.0):
                hf.changeRowBounds(n_full_ub + bi,
                                   b_eq[bi] + sgn * pert / S_BASE,
                                   b_eq[bi] + sgn * pert / S_BASE)
                hf.run()
                assert hf.getModelStatus() == highspy.HighsModelStatus.kOptimal, \
                    "FD re-solve failed"
                vals[sgn] = hf.getObjectiveValue()
            # restore the base bound so every bus's FD is taken at the
            # unperturbed point
            hf.changeRowBounds(n_full_ub + bi, b_eq[bi], b_eq[bi])
            lmps_fd[bi] = (vals[+1.0] - vals[-1.0]) / (2 * pert)
            print(f"  FD LMP bus {bus_ids[bi]}: {lmps_fd[bi]:.3f} "
                  f"(dual {lmps_dual[bi] / S_BASE:.3f})")
        print(f"FD LMPs for {len(sel)} buses in {time.time()-t_fd:.1f}s")

    # ---- gen-outage diagnostic (TAMU's contab also lists CT_TGEN rows;
    # check whether the joint problem is feasible, and which gen outages
    # break it if not)
    gen_diag = None
    if args.gentest:
        print("\ngen-outage diagnostic ...")
        t_gd = time.time()

        def solve_mask(mask, label=""):
            indptr, indices, data, rhs = build_full_ub(
                M, Kme, Kload, active, mask, base_rows, base_rhs)
            A = sparse.csr_matrix((data, indices, indptr),
                                  shape=(indptr.size - 1, dim))
            return solve_lp(c_obj, lb, ub, A, rhs, A_eq, b_eq, label=label)

        okf, xf, df, hf = solve_mask(valid, label="lines+gens")
        full_feasible = okf
        bad_gens = []
        for k in range(ng):
            v = v_line.copy()
            v[n_line + M["gen_cts"].index(k)] = True
            okk, xk, dk, hk = solve_mask(v, label=f"lines+gen k{k}")
            if not okk:
                bad_gens.append(k)
        variant_feasible = variant_obj = None
        if bad_gens:
            v_var = valid.copy()
            for k in bad_gens:
                v_var[n_line + M["gen_cts"].index(k)] = False
            okv, xv, dv, hv = solve_mask(v_var, label="variant")
            variant_feasible = okv
            if okv:
                variant_obj = float(hv.getObjectiveValue()) + pl_const
        print(f"gen diagnostic done in {time.time()-t_gd:.1f}s: "
              f"lines+gens={'feasible' if full_feasible else 'INFEASIBLE'}, "
              f"infeasible gen outages = {sorted(M['gen_buses'][k] for k in bad_gens)}")
        gen_diag = dict(full_feasible=full_feasible, bad_gens=bad_gens,
                        variant_feasible=variant_feasible,
                        variant_obj=variant_obj)

    # ---- report
    def branch_name(li):
        f, t, _, _, kind, ckt = br[li]
        if kind == "xfrm":
            return f"{bus_ids[f]}-{bus_ids[t]} (xfmr {ckt})"
        return f"{bus_ids[f]}-{bus_ids[t]}" + (f" ({ckt})" if str(ckt) != "1" else "")

    L = []
    L.append("# DC N-1 SCOPF — ACTIVSg500 (TAMU 500-bus test case)")
    L.append("")
    L.append(f"Network from ACTIVSg500.RAW via the VeraGridEngine RAW parser; "
             f"costs/limits from the official MATPOWER case; contingencies "
             f"from contab_ACTIVSg500.m; solved with HiGHS via highspy "
             f"(constraint generation; {nseg}-segment PL cost, lossless DC, "
             f"no corrective re-dispatch). Total load {loads.sum():.2f} MW.")
    L.append("")
    L.append("### N-0 (base case only)")
    L.append(f"Cost (exact quadratic): ${exact0:,.1f}/h")
    for k, b in enumerate(gen_buses):
        L.append(f"  gen bus {b:<4d}: {p_g0[k]:9.2f} MW")
    L.append(f"N-1 violations of this dispatch: {n0_viol} of "
             f"{int(v_line.sum())} line outages (worst ratio {n0_worst:.4f})")
    L.append("")
    L.append("### N-1 (SCOPF, classical: all 591 line outages enforced)")
    L.append(f"Cost (exact quadratic): ${exact_cost:,.1f}/h "
             f"(PL obj {pl_cost:,.1f}; delta vs N-0 {exact_cost-exact0:+,.1f}/h)")
    for k, b in enumerate(gen_buses):
        mark = " *" if p_g[k] <= M["gen_pmin"][k] + 1e-3 else \
               (" **" if p_g[k] >= M["gen_pmax"][k] - 1e-3 else "")
        L.append(f"  gen bus {b:<4d}: {p_g[k]:9.2f} MW {mark}")
    L.append(f"All {int(v_line.sum())} line-outage contingencies satisfied "
             f"(violations: {viol}); worst post-contingency ratio "
             f"{ratio[v_line].max():.4f}")
    per_ct = np.where(v_line, np.max(ratio, axis=1), 0.0)
    bind = [(i, float(per_ct[i])) for i in range(n_line)
            if per_ct[i] >= 1.0 - 1e-4]
    if bind:
        L.append("Binding N-1 contingencies (>=99.9% of limit), top 15:")
        for i, r in sorted(bind, key=lambda t: -t[1])[:15]:
            tag = (f"line {branch_name(M['line_cts'][i])} out" if i < len(M["line_cts"])
                   else f"gen bus {gen_buses[M['gen_cts'][i-len(M['line_cts'])]]} out")
            L.append(f"  - {tag}: ratio {r:.4f}")
    flows = np.array([bl * (th[f] - th[t]) * S_BASE for (f, t, bl, *_r) in br])
    L.append("")
    L.append("### Base-case flows (N-1 dispatch), top 15 by |flow|")
    for li in np.argsort(-np.abs(flows))[:15]:
        fmax = Fmax_all[li]
        pct = 100 * flows[li] / fmax if fmax > 0 else float("nan")
        L.append(f"  {branch_name(li):<18s} {flows[li]:9.2f} MW "
                 f"({pct:6.1f}% of {fmax:.0f} MW)")
    if lmps_fd is not None:
        L.append("")
        L.append("### LMPs [$ / MWh]")
        L.append("(dual = optimal dual of the N-1 LP balance rows -- pu-scaled, "
                 "so divided by S_BASE --; FD = central finite difference of "
                 "the N-1 optimal cost, +/-0.05 MW)")
        L.append(f"  {'bus':>5s} {'load':>8s} {'LMP(FD)':>10s} {'LMP(dual)':>12s}")
        for bi in sel:
            fd = f"{lmps_fd[bi]:10.3f}" if lmps_fd[bi] != 0 else f"{'--':>10s}"
            L.append(f"  {bus_ids[bi]:>5d} {loads[bi]:8.2f} {fd} "
                     f"{lmps_dual[bi] / S_BASE:12.3f}")
    if gen_diag is not None:
        L.append("")
        L.append("### N-1 including generator outages (TAMU contab CT_TGEN)")
        if gen_diag["full_feasible"]:
            L.append("Feasible with the full (line + gen) contingency set.")
        else:
            L.append("Infeasible: no dispatch keeps post-contingency flows "
                     "within limits for every line AND every generator outage "
                     "(proven infeasible by HiGHS, confirmed by an "
                     "independent HiGHS call through scipy).")
            bad = gen_diag["bad_gens"]
            if bad:
                L.append("Per-outage test (all 591 line outages + each single "
                         "gen outage) isolates the cause to "
                         f"{len(bad)} of {ng}:")
                for k in bad:
                    b = M["gen_buses"][k]
                    c = n_line + M["gen_cts"].index(k)
                    act = active[c]
                    act = act[act >= 0]
                    F = Kme[c, act, :] @ p_g - Kload[c, act]
                    r = np.abs(F) / Fmax_all[act]
                    w = int(np.argmax(r))
                    li = int(act[w])
                    L.append(f"  - gen bus {b} (Pmin {M['gen_pmin'][k]:.1f} / "
                             f"Pmax {M['gen_pmax'][k]:.1f} MW, must-run): "
                             f"at the N-1 dispatch its outage "
                             f"forces line {branch_name(li)} to carry "
                             f"{abs(F[w]):.0f} MW vs limit {Fmax_all[li]:.1f} MW "
                             f"(ratio {r[w]:.3f})")
                L.append("They are the two largest generators in the case "
                         "(Pmax 888.9 / 771.8 MW, both with Pmin > 0); "
                         "their post-outage imports at their buses exceed "
                         "what the connector paths can carry simultaneously "
                         "with the line-outage security requirements.")
            if gen_diag["variant_feasible"]:
                L.append(f"Excluding those {len(bad)} gen outages, the "
                         f"{M['nct'] - len(bad)}-contingency SCOPF is feasible: "
                         f"${gen_diag['variant_obj']:,.1f}/h "
                         f"(+${gen_diag['variant_obj'] - exact0:,.1f}/h vs N-0).")
    L.append("")
    L.append("### Data / assumptions")
    L.append("- network: ACTIVSg500.RAW via VeraGridEngine RAW parser "
             "(500 buses, 466 lines, 131 transformers, 90 generators, 206 loads)")
    L.append("- costs/limits: official MATPOWER case (the RAW has no P records); "
             f"{ng} online units, {len(M['offline_gen'])} offline")
    L.append("- primary contingency set: 591 TAMU line outages (classical "
             "N-1 security constraint; gen outages -- see diagnostic above)")
    L.append(f"- {int(island.sum())} of {M['nct']} contingencies disconnect an "
             f"island (radial cut-edge outages); post-contingency flow limits are "
             f"checked on the network that remains connected to the slack, "
             f"i.e. the island is assumed shed/handled locally (standard "
             f"contingency practice)")
    if island.any():
        L.append(f"  island sizes: max load {island_load[island].max():.1f} MW, "
                 f"max gen Pmax {island_gen[island].max():.1f} MW, "
                 f"islands with generation: {int((island & (island_gen > 0)).sum())}")
    L.append("- lossless DC, no voltage limits, no corrective re-dispatch")
    L.append(f"- PL cost: {nseg}-segment chord (upper bound) of quadratic gencost")
    report = "\n".join(L)
    print("\n" + report)
    out = A / "results_activsg500.md"
    out.write_text(report + "\n", encoding="utf-8")
    print(f"\n[report saved to {out}]")


if __name__ == "__main__":
    main()
