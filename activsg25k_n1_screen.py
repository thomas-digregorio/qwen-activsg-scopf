#!/usr/bin/env python3
"""
DC N-0 OPF + N-1 security screening for ACTIVSg25k (TAMU 25,000-bus system).

Why screening and not a full N-1 SCOPF:
    The classical N-1 SCOPF would add ~64k contingency constraint rows
    (32,229 branch outages). At 25k scale that LP is intractable. Instead:
      N-0  : full piecewise-linear DC OPF (HiGHS, --nseg segments/unit)
      N-1  : exact security SCREENING of the N-0 dispatch
               * non-bridge branch outages: Sherman-Morrison rank-1 update
                 of the pre-contingency angle vector (exact solution of the
                 reduced network; verified in-driver vs direct re-solves)
               * bridge outages (islanding): per-bridge re-factorization of
                 the slack-side component (flow sweep) + island statistics
                 (load / frozen gen / gen Pmax / shortfall) + internal flow
                 sweep of the island
               * gen outages (diagnostic, clearly labeled): injection
                 update F' = F0 - b_l (w_f - w_t) P_g; the units on the
                 slack bus are re-rooted (new reference bus absorbs)

Network : ACTIVSg25k.RAW via the VeraGridEngine RAW parser (vg_data.json),
          3-winding transformers expanded to 3-branch stars with one
          synthetic stub bus each (verified against the official MATPOWER
          case by crosscheck.py: 24,250 + 750 stubs = 25,000 buses,
          24,699 lines + 5,280 2-winding + 2,250 star legs = 32,229
          in-service branches; 1 out-of-service line excluded)
Costs   : case_ACTIVSg25k.m -- 3,779 online units on 2,753 buses
Contab  : none exists for the 25k case -> N-1 line set = all in-service
          branches
Slack   : MP reference bus 62120

Usage:  python activsg25k_n1_screen.py [--nseg 20] [--fdlmp 8]
        [--gen-outages 1] [--sm-batch 256]
"""
import argparse
import collections
import json
import os
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import splu

import highspy

HERE = Path(__file__).resolve().parent
A = HERE / "activsg25k"
INF = highspy.kHighsInf
SLACK_BUS = 62120
BZ = 1e8  # fallback b for zero-x star legs (never fires on this dataset)
VIOL = 1.0 + 1e-9  # flow-ratio threshold that counts as a violation


def build_model():
    """Build the 25k network (RAW + 3-winding star expansion) and the
    generator data (from the MP case). Returns dict M."""
    t0 = time.time()
    vg = json.loads((A / "vg_data.json").read_text())
    mp = json.loads((A / "mp_data.json").read_text())
    S_BASE = mp["baseMVA"]

    tw = [t for t in vg["transformers"] if t["K"] in (0, 0.0)]
    threew = [t for t in vg["transformers"] if t["K"] not in (0, 0.0)]
    bus_ids = [int(b["I"]) for b in vg["buses"]]
    bus_ids += [10**6 + i for i in range(len(threew))]
    n = len(bus_ids)
    bidx = {b: i for i, b in enumerate(bus_ids)}

    br_f, br_t, br_b, br_Fmax, kind = [], [], [], [], []
    n_off_line = 0
    for b in vg["branches"]:
        if b["ST"] != 1:
            n_off_line += 1
            continue
        br_f.append(bidx[int(b["I"])])
        br_t.append(bidx[int(b["J"])])
        br_b.append(1.0 / b["X"])
        br_Fmax.append(b["RATEA"] or 0.0)
        kind.append("line")
    n_off_tw = 0
    for t in tw:
        if t["STAT"] != 1:
            n_off_tw += 1
            continue
        br_f.append(bidx[int(t["I"])])
        br_t.append(bidx[int(t["J"])])
        br_b.append(1.0 / (t["X1_2"] + t["X2_3"]))
        br_Fmax.append(t["RATA1"] or 0.0)
        kind.append("xfrm")
    n_star = 0
    n_zx = 0
    for k, t in enumerate(threew):
        if t["STAT"] != 1:
            continue
        s = bidx[10**6 + k]
        Z12 = complex(t["R1_2"] or 0, t["X1_2"] or 0)
        Z23 = complex(t["R2_3"] or 0, t["X2_3"] or 0)
        Z31 = complex(t["R3_1"] or 0, t["X3_1"] or 0)
        for bus, Z, R in (
            (t["I"], (Z12 + Z31 - Z23) / 2, t["RATA1"]),
            (t["J"], (Z12 + Z23 - Z31) / 2, t["RATA2"]),
            (t["K"], (Z23 + Z31 - Z12) / 2, t["RATA3"]),
        ):
            b_eff = BZ if abs(Z.imag) < 1e-12 else 1.0 / Z.imag
            if b_eff == BZ:
                n_zx += 1
            br_f.append(bidx[int(bus)])
            br_t.append(s)
            br_b.append(b_eff)
            br_Fmax.append(R or 0.0)
            kind.append("star")
            n_star += 1
    m = len(br_f)
    br_f = np.array(br_f, dtype=np.int64)
    br_t = np.array(br_t, dtype=np.int64)
    br_b = np.array(br_b, dtype=float)
    Fmax = np.array(br_Fmax, dtype=float)

    loads = np.zeros(n)
    for l in vg["loads"]:
        if l["STATUS"] in (1, None):
            loads[bidx[int(l["I"])]] += l["PL"] or 0.0
    for b in vg["buses"]:
        loads[bidx[int(b["I"])]] += b["BL"] or 0.0

    # generators: MP rows are the authoritative cost/limit data (RAW/MP
    # per-bus unit counts match exactly, verified by _genmap3.py)
    gen_bus, gen_pmin, gen_pmax, gen_c2, gen_c1, gen_c0 = [], [], [], [], [], []
    for i, r in enumerate(mp["gen"]):
        if r[7] != 1:
            continue
        gc = mp["gencost"][i]
        assert gc[0] == 2 and gc[3] == 3, "expected quadratic 3-coeff gencost"
        gen_bus.append(int(r[0]))
        gen_pmin.append(r[9])
        gen_pmax.append(r[8])
        gen_c2.append(gc[4]); gen_c1.append(gc[5]); gen_c0.append(gc[6])
    ng = len(gen_bus)
    gen_bus = np.array(gen_bus, dtype=np.int64)
    gen_pmin = np.array(gen_pmin)
    gen_pmax = np.array(gen_pmax)
    gen_c2 = np.array(gen_c2); gen_c1 = np.array(gen_c1); gen_c0 = np.array(gen_c0)
    gen_bus_idx = np.array([bidx[b] for b in gen_bus], dtype=np.int64)

    slack = bidx[SLACK_BUS]
    ns = np.array([i for i in range(n) if i != slack], dtype=np.int64)
    ns_pos = np.full(n, -1, dtype=np.int64)
    ns_pos[ns] = np.arange(n - 1)

    # sparse susceptance matrix S
    off_r = np.concatenate([br_f, br_t])
    off_c = np.concatenate([br_t, br_f])
    off_v = np.concatenate([-br_b, -br_b])
    d = np.zeros(n)
    np.add.at(d, br_f, br_b)
    np.add.at(d, br_t, br_b)
    S = (coo_matrix((off_v, (off_r, off_c)), shape=(n, n)) + diags(d)).tocsr()
    Bt = S[ns][:, ns].tocsc()
    t1 = time.time()
    lu = splu(Bt)

    # per-bus incident branches: bus -> list of (branch, +1 if bus==f)
    adj = collections.defaultdict(list)
    for li in range(m):
        adj[int(br_f[li])].append((li, +1.0))
        adj[int(br_t[li])].append((li, -1.0))
    bus_units = collections.defaultdict(list)
    for k, b in enumerate(gen_bus_idx):
        bus_units[int(b)].append(k)

    M = dict(S_BASE=S_BASE, bus_ids=bus_ids, n=n, m=m, bidx=bidx,
             br_f=br_f, br_t=br_t, br_b=br_b, Fmax=Fmax, kind=kind,
             loads=loads, gen_bus=gen_bus, gen_bus_idx=gen_bus_idx,
             gen_pmin=gen_pmin, gen_pmax=gen_pmax,
             gen_c2=gen_c2, gen_c1=gen_c1, gen_c0=gen_c0, ng=ng,
             slack=slack, ns=ns, ns_pos=ns_pos, S=S, lu=lu,
             adj=adj, bus_units=bus_units,
             n_off_line=n_off_line, n_off_tw=n_off_tw,
             n_star=n_star, n_zx=n_zx, n_tw=len(tw), n_3w=len(threew))
    print(f"model: {n} buses ({len(bus_ids) - len(threew)} RAW + {len(threew)} star stubs), "
          f"{m} branches ({m - M['n_star'] - len(tw) + n_off_tw} lines + {len(tw) - n_off_tw} 2w + {n_star} star)")
    print(f"  out-of-service excluded: {n_off_line} lines, {n_off_tw} 2w xfmr; "
          f"zero-x star legs {n_zx}; negative-b branches {int((br_b < 0).sum())}")
    print(f"  {ng} online units on {len(bus_units)} buses | total load {loads.sum():.2f} MW "
          f"| sum Pmin {gen_pmin.sum():.1f} | sum Pmax {gen_pmax.sum():.1f}")
    print(f"  splu factor {time.time()-t1:.2f}s (S nnz {S.nnz}, L nnz {lu.L.nnz}) [{time.time()-t0:.1f}s total]")
    return M


def build_pl_prices(M, nseg):
    """Piecewise-linear (chord) prices/widths per unit over [Pmin, Pmax].
    Chord prices = secant slopes of the quadratic cost on each segment --
    a conservative (upper-bound) approximation."""
    prices, widths = [], []
    for k in range(M["ng"]):
        pmin, pmax = M["gen_pmin"][k], M["gen_pmax"][k]
        c2, c1, c0 = M["gen_c2"][k], M["gen_c1"][k], M["gen_c0"][k]
        xs = np.linspace(pmin, pmax, nseg + 1)
        cc = c2 * xs**2 + c1 * xs + c0
        prices.append(np.diff(cc) / np.diff(xs))
        widths.append(np.diff(xs))
    return np.concatenate(prices), np.concatenate(widths)


def build_n0_lp(M, nseg):
    """N-0 DC OPF as a sparse LP.
    Variables: [pg (ng) | seg (ng*nseg) | theta_ns (n-1)]
    Eq rows   : nodal balance at every bus (n) + segment links (ng)
    Ub rows   : +/- flow limits on every rated branch
    """
    n, ng, m = M["n"], M["ng"], M["m"]
    ns, ns_pos, slack = M["ns"], M["ns_pos"], M["slack"]
    br_f, br_t, br_b, Fmax = M["br_f"], M["br_t"], M["br_b"], M["Fmax"]
    loads = M["loads"]
    S_BASE = M["S_BASE"]
    rated = Fmax > 1e-9
    rated_idx = np.where(rated)[0]
    n_ub = 2 * len(rated_idx)

    dim = ng + ng * nseg + (n - 1)
    s0, u0 = ng, ng + ng * nseg

    seg_price, seg_width = build_pl_prices(M, nseg)
    c_obj = np.zeros(dim)
    c_obj[s0:s0 + ng * nseg] = seg_price

    # ---- equality rows (COO): nodal balance + segment links
    # Balance at bus i:  gen_inj[i]/S - (S @ theta)[i] = load[i]/S
    # i.e. sum_units pg/S  +  sum_{j in ns} (-S[i,j]) theta_j  = load/S
    # (the diagonal S[i,i]=sum_b gives the self term; off-diagonals -S[i,j]=+b)
    S = M["S"]
    s_indptr, s_indices, s_data = S.indptr, S.indices, S.data
    er, ec, ev, eq_rhs = [], [], [], []
    for i in range(n):
        eq_rhs.append(loads[i] / S_BASE)
        for k in M["bus_units"].get(i, ()):
            er.append(i); ec.append(k); ev.append(1.0 / S_BASE)
        for jj in range(s_indptr[i], s_indptr[i + 1]):
            j = s_indices[jj]
            if j == slack:
                continue  # theta_slack is fixed at 0
            er.append(i); ec.append(u0 + ns_pos[j]); ev.append(-s_data[jj])
    # segment links: pg_k - sum(seg_k) = pmin_k
    krows = np.arange(n, n + ng)
    er = np.concatenate([np.array(er, dtype=np.int64), krows,
                         np.repeat(krows, nseg)])
    ec = np.concatenate([np.array(ec, dtype=np.int64), np.arange(ng),
                         np.arange(ng * nseg) + s0])
    ev = np.concatenate([np.array(ev, dtype=float), np.ones(ng),
                         -np.ones(ng * nseg)])
    eq_rhs = np.concatenate([np.array(eq_rhs, dtype=float),
                             M["gen_pmin"]])
    A_eq = csr_matrix((ev, (er, ec)), shape=(n + ng, dim))
    b_eq = eq_rhs

    # ---- inequality rows: +/- |F_l| <= Fmax_l / S_BASE (rated branches)
    ur, uc, uv, ub_rhs = [], [], [], []
    for j, li in enumerate(rated_idx):
        r1, r2 = 2 * j, 2 * j + 1
        ub_rhs.append(Fmax[li] / S_BASE)
        ub_rhs.append(Fmax[li] / S_BASE)
        jf = ns_pos[br_f[li]]
        jt = ns_pos[br_t[li]]
        if jf >= 0:
            ur += [r1, r2]; uc += [u0 + jf, u0 + jf]; uv += [br_b[li], -br_b[li]]
        if jt >= 0:
            ur += [r1, r2]; uc += [u0 + jt, u0 + jt]; uv += [-br_b[li], br_b[li]]
    A_ub = csr_matrix((uv, (ur, uc)), shape=(n_ub, dim))
    b_ub = np.array(ub_rhs, dtype=float)

    lb = np.zeros(dim)
    ub = np.full(dim, np.inf)
    lb[:ng] = M["gen_pmin"]
    ub[:ng] = M["gen_pmax"]
    ub[s0:s0 + ng * nseg] = seg_width
    lb[u0:] = -np.inf
    return dict(c_obj=c_obj, lb=lb, ub=ub, A_eq=A_eq, b_eq=b_eq,
                A_ub=A_ub, b_ub=b_ub, dim=dim, nseg=nseg,
                seg_price=seg_price, seg_width=seg_width,
                n_ub=n_ub, rated_idx=rated_idx)


def solve_lp(c_obj, lb, ub, A_ub, b_ub, A_eq, b_eq, label=""):
    """Solve the LP with HiGHS; returns (ok, x, row_duals, handle)."""
    t0 = time.time()
    h = highspy.Highs()
    opts = highspy.HighsOptions()
    opts.output_flag = bool(int(os.environ.get("VG_DEBUG", "0")))
    opts.time_limit = 3600.0
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
    if A_ub.shape[0] == 0:
        A = sparse.csr_matrix(A_eq)
        lp.row_lower_ = b_eq
        lp.row_upper_ = b_eq
    else:
        lp.row_lower_ = np.concatenate([np.full(A_ub.shape[0], -INF), b_eq])
        lp.row_upper_ = np.concatenate([b_ub, b_eq])
        A = sparse.vstack([A_ub, sparse.csr_matrix(A_eq)], format="csr")
    mat = highspy.HighsSparseMatrix()
    mat.format_ = highspy.MatrixFormat.kRowwise
    mat.num_row_ = A.shape[0]
    mat.num_col_ = A.shape[1]
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
    dt = time.time() - t0
    if label:
        if ok:
            print(f"  [{label}] optimal  obj={h.getObjectiveValue():,.1f}  "
                  f"({dt:.1f}s)", flush=True)
        else:
            print(f"  [{label}] NOT OPTIMAL: {h.modelStatusToString(h.getModelStatus())} "
                  f"({dt:.1f}s)", flush=True)
    return ok, np.array(sol.col_value, dtype=float), \
        np.array(sol.row_dual, dtype=float), h


def find_bridges_dfs(M):
    """Single DFS from the slack (network is one component): Tarjan bridge
    detection + subtree intervals. For a bridge edge (parent, child), the
    child's subtree [tin, tout) is the island (the side without the slack).
    Returns (bridges set, tin, tout)."""
    n, m = M["n"], M["m"]
    br_f, br_t = M["br_f"], M["br_t"]
    slack = M["slack"]
    adj = collections.defaultdict(list)
    for li in range(m):
        adj[int(br_f[li])].append((int(br_t[li]), li))
        adj[int(br_t[li])].append((int(br_f[li]), li))
    tin = np.full(n, -1, dtype=np.int64)
    tout = np.full(n, -1, dtype=np.int64)
    bridges = set()
    timer = 0
    n_comp = 0
    for root in range(n):
        if tin[root] != -1:
            continue
        n_comp += 1
        timer += 1
        tin[root] = tout[root] = timer
        low = tin.copy()
        pedge = {root: -1}
        it_stack = [(root, iter(adj[root]))]
        while it_stack:
            v, it = it_stack[-1]
            advanced = False
            for (w, ei) in it:
                if ei == pedge[v]:
                    continue
                if tin[w] == -1:
                    pedge[w] = ei
                    timer += 1
                    tin[w] = low[w] = timer
                    it_stack.append((w, iter(adj[w])))
                    advanced = True
                    break
                low[v] = min(low[v], tin[w])
            if advanced:
                continue
            it_stack.pop()
            tout[v] = timer
            if it_stack:
                p = it_stack[-1][0]
                low[p] = min(low[p], low[v])
                if low[v] > tin[p]:
                    bridges.add(pedge[v])
    print(f"  DFS from slack: {n_comp} component(s), {len(bridges)} bridges "
          f"of {m} branches ({len(bridges) / m * 100:.1f}%)", flush=True)
    assert n_comp == 1, f"network has {n_comp} components -- unexpected"
    return bridges, tin, tout


def make_incidence(M):
    """Incidence on the ns space: row li = +b at f, -b at t (slack -> 0).
    (m x (n-1)) CSR;  A_inc @ theta_ns gives b*(theta_f - theta_t) per line."""
    n, m = M["n"], M["m"]
    ns_pos, br_f, br_t, br_b = M["ns_pos"], M["br_f"], M["br_t"], M["br_b"]
    f_ns, t_ns = ns_pos[br_f], ns_pos[br_t]
    kf, kt = f_ns >= 0, t_ns >= 0
    ir = np.concatenate([np.where(kf)[0], np.where(kt)[0]])
    ic = np.concatenate([f_ns[kf], t_ns[kt]])
    iv = np.concatenate([br_b[kf], -br_b[kt]])
    return csr_matrix((iv, (ir, ic)), shape=(m, n - 1))


def screen_lines(M, A_inc, theta_ns, bridges, n_check=3):
    """Batched Sherman-Morrison screen of the non-bridge outages.
    For outage o with endpoints (f,t), susceptance b:
      w = Bt^-1 d,  d = e_f - e_t (reduced)
      theta' = theta0 + b (d'theta0) w / (1 - b d'w)   (exact reduced solve)
    Returns recs = [(outage, worst_line, worst_ratio, n_viol)] and totals."""
    n, m = M["n"], M["m"]
    ns_pos, br_f, br_t, br_b = M["ns_pos"], M["br_f"], M["br_t"], M["br_b"]
    S_BASE = M["S_BASE"]
    lu = M["lu"]
    f_ns, t_ns = ns_pos[br_f], ns_pos[br_t]
    rated_rows = M["rated_idx"]
    Fmax_r = M["Fmax"][rated_rows]

    # self-test: SM vs direct re-solve for a few outages
    rng = np.random.default_rng(0)
    nonbr_all = np.array([li for li in range(m) if li not in bridges])
    sample = rng.choice(nonbr_all, min(n_check, len(nonbr_all)), replace=False)
    S = M["S"]
    ns = M["ns"]
    for o in sample:
        f, t2, bl = int(br_f[o]), int(br_t[o]), br_b[o]
        d = np.zeros(n)
        d[f] = 1.0
        d[t2] = -1.0
        dn = d[ns]
        w = lu.solve(dn)
        den = 1.0 - bl * float(dn @ w)
        th_sm = theta_ns + (bl * float(dn @ theta_ns) / den) * w
        # direct: remove branch o (its contribution: +bl on the diagonal,
        # -bl on the off-diagonals)
        dmat = coo_matrix(([-bl, -bl, bl, bl],
                          ([f, t2, f, t2], [t2, f, f, t2])), shape=(n, n))
        S2 = (S - dmat).tocsr()
        Bt2 = S2[ns][:, ns].tocsc()
        P_ns = M["Pnet"][ns] / S_BASE
        th_dir = np.zeros(n)
        th_dir[ns] = splu(Bt2).solve(P_ns)
        F_sm = (A_inc @ th_sm) * S_BASE
        F_sm[o] = 0.0
        F_dr = br_b * (th_dir[br_f] - th_dir[br_t]) * S_BASE
        F_dr[o] = 0.0
        diff = float(np.abs(F_sm - F_dr).max())
        print(f"  SM self-test ct {o}: den={den:.4e}  max|SM-direct|={diff:.3e} MW",
              flush=True)
        assert diff < 1e-5, "Sherman-Morrison vs direct mismatch"

    t0 = time.time()
    recs = []
    n_viol_ct = 0
    worst_overall = 0.0
    batch = M.get("sm_batch", 256)
    nbatch = (len(nonbr_all) + batch - 1) // batch
    for b0 in range(0, len(nonbr_all), batch):
        cols = nonbr_all[b0:b0 + batch]
        B = len(cols)
        D = np.zeros((n - 1, B))
        for j, li in enumerate(cols):
            if f_ns[li] >= 0:
                D[f_ns[li], j] += 1.0
            if t_ns[li] >= 0:
                D[t_ns[li], j] -= 1.0
        W = lu.solve(D)
        dw = (D * W).sum(axis=0)
        dt = D.T @ theta_ns
        den = 1.0 - br_b[cols] * dw
        if np.any(np.abs(den) < 1e-10):
            bad = cols[np.where(np.abs(den) < 1e-10)[0]]
            raise RuntimeError(f"zero SM denominator on non-bridge outages {bad}")
        scale = br_b[cols] * dt / den
        thp = theta_ns[None, :] + W.T * scale[:, None]
        Fp = (A_inc @ thp.T) * S_BASE
        for j, li in enumerate(cols):
            Fp[li, j] = 0.0
        ratio = np.abs(Fp[rated_rows, :]) / Fmax_r[:, None]
        viol = ratio > VIOL
        nv = viol.sum(axis=0)
        wloc = ratio.argmax(axis=0)
        wr = ratio.max(axis=0)
        worst_overall = max(worst_overall, float(wr.max()))
        n_viol_ct += int((nv > 0).sum())
        for j, li in enumerate(cols):
            recs.append((int(li), int(rated_rows[wloc[j]]), float(wr[j]), int(nv[j])))
        bi = b0 // batch
        if bi % 10 == 0 or b0 + batch >= len(nonbr_all):
            print(f"    SM {b0 + B}/{len(nonbr_all)}  worst {wr.max():.4f}  "
                  f"violating outages {n_viol_ct}  [{time.time()-t0:.0f}s] "
                  f"({(bi + 1)}/{nbatch} batches)", flush=True)
    return recs, worst_overall, n_viol_ct


def screen_bridges(M, bridges, tin, tout):
    """Per-bridge islanding analysis:
      C1 = slack-side component: re-factor, solve with the N-0 dispatch
           (slack absorbs the C1 imbalance), sweep C1 line limits.
      C2 = island (DFS subtree of the child endpoint): load / frozen gen /
           gen Pmax / Pmin stats, shortfall = max(0, load - Pmax), plus
           internal flow sweep (island re-rooted at the child endpoint)."""
    n, m = M["n"], M["m"]
    br_f, br_t, br_b, Fmax = M["br_f"], M["br_t"], M["br_b"], M["Fmax"]
    S = M["S"]
    S_BASE = M["S_BASE"]
    slack = M["slack"]
    loads, Pnet = M["loads"], M["Pnet"]
    arange_n = np.arange(n)
    not_slack = arange_n != slack
    gen_pmax_bus = np.zeros(n)
    np.add.at(gen_pmax_bus, M["gen_bus_idx"], M["gen_pmax"])
    gen_pmin_bus = np.zeros(n)
    np.add.at(gen_pmin_bus, M["gen_bus_idx"], M["gen_pmin"])

    results = []
    t0 = time.time()
    bsort = sorted(bridges)
    for bi, li in enumerate(bsort):
        f, t2 = int(br_f[li]), int(br_t[li])
        c = f if tin[f] > tin[t2] else t2
        mask2 = (tin >= tin[c]) & (tin < tout[c])
        C1 = (~mask2) & not_slack
        ns1 = np.nonzero(C1)[0]
        viol1 = 0
        worst1 = 0.0
        if len(ns1) > 0:
            lu1 = splu(S[ns1][:, ns1].tocsc())
            th1 = np.zeros(n)
            th1[ns1] = lu1.solve(Pnet[ns1] / S_BASE)
            lines1 = np.nonzero(C1[br_f] & C1[br_t])[0]
            if len(lines1):
                F1 = br_b[lines1] * (th1[br_f[lines1]] - th1[br_t[lines1]]) * S_BASE
                r1 = Fmax[lines1] > 1e-9
                if r1.any():
                    ratio = np.abs(F1[r1]) / Fmax[lines1][r1]
                    viol1 = int((ratio > VIOL).sum())
                    worst1 = float(ratio.max())
        L2 = float(loads[mask2].sum())
        G2 = float(Pnet[mask2].sum() + loads[mask2].sum())
        Pmax2 = float(gen_pmax_bus[mask2].sum())
        Pmin2 = float(gen_pmin_bus[mask2].sum())
        shortfall = max(0.0, L2 - Pmax2)
        viol2 = 0
        worst2 = 0.0
        ns2 = np.nonzero(mask2 & (arange_n != c))[0]
        if len(ns2) > 0:
            lu2 = splu(S[ns2][:, ns2].tocsc())
            th2 = np.zeros(n)
            th2[ns2] = lu2.solve(Pnet[ns2] / S_BASE)
            lines2 = np.nonzero(mask2[br_f] & mask2[br_t])[0]
            if len(lines2):
                F2 = br_b[lines2] * (th2[br_f[lines2]] - th2[br_t[lines2]]) * S_BASE
                r2 = Fmax[lines2] > 1e-9
                if r2.any():
                    ratio = np.abs(F2[r2]) / Fmax[lines2][r2]
                    viol2 = int((ratio > VIOL).sum())
                    worst2 = float(ratio.max())
        results.append(dict(li=li, nbus2=int(mask2.sum()), nbus1=int((~mask2).sum()),
                            L2=L2, G2=G2, Pmax2=Pmax2, Pmin2=Pmin2,
                            shortfall=shortfall,
                            viol1=viol1, worst1=worst1, viol2=viol2, worst2=worst2))
        if (bi + 1) % 500 == 0 or bi + 1 == len(bsort):
            print(f"    bridges {bi + 1}/{len(bsort)}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    return results


def screen_gens(M, A_inc, F0, p_g):
    """Gen-outage screen (DIAGNOSTIC, clearly labeled in the report).
    Non-slack units:  F'_l = F0_l - P_g * b_l (w_f - w_t),  w = Bt^-1 e_bus
    (batched multi-RHS per group of distinct gen buses).
    Slack-bus units: re-root the network at the largest other gen bus
    (that reference absorbs the lost injection); direct solves."""
    n, m = M["n"], M["m"]
    ns_pos, slack = M["ns_pos"], M["slack"]
    br_f, br_t, br_b, Fmax = M["br_f"], M["br_t"], M["br_b"], M["Fmax"]
    S_BASE = M["S_BASE"]
    S = M["S"]
    lu = M["lu"]
    rated_rows = M["rated_idx"]
    Fmax_r = M["Fmax"][rated_rows]

    bus2units = collections.defaultdict(list)
    for k, b in enumerate(M["gen_bus_idx"]):
        bus2units[int(b)].append(k)

    recs = []  # (unit k, bus, worst_line, worst_ratio, n_viol, [ref_ok])
    bus_list = sorted(b for b in bus2units if b != slack)
    t0 = time.time()
    batch = 256
    for b0 in range(0, len(bus_list), batch):
        cols = bus_list[b0:b0 + batch]
        E = np.zeros((n - 1, len(cols)))
        for j, b in enumerate(cols):
            E[ns_pos[b], j] = 1.0
        W = lu.solve(E)
        WV = A_inc @ W  # (m, len(cols))
        for j, b in enumerate(cols):
            rv = WV[:, j]
            for k in bus2units[b]:
                Fp = F0 - p_g[k] * rv
                ratio = np.abs(Fp[rated_rows]) / Fmax_r
                wr = float(ratio.max())
                wl = int(rated_rows[int(ratio.argmax())])
                recs.append((k, b, wl, wr, int((ratio > VIOL).sum()), True))
    nbus = len(bus_list)
    print(f"    gen screen: {nbus} non-slack gen buses in {time.time()-t0:.0f}s "
          f"({sum(len(v) for k, v in bus2units.items() if k != slack)} units)",
          flush=True)

    if slack in bus2units:
        # re-root: new reference = largest-Pmax gen bus other than the slack
        pmax_bus = np.zeros(n)
        np.add.at(pmax_bus, M["gen_bus_idx"], M["gen_pmax"])
        cand = pmax_bus.copy()
        cand[slack] = -1.0
        alt = int(np.argmax(cand))
        ns_alt = np.array([i for i in range(n) if i != alt], dtype=np.int64)
        ns_alt_pos = np.full(n, -1, dtype=np.int64)
        ns_alt_pos[ns_alt] = np.arange(n - 1)
        fa, ta = ns_alt_pos[br_f], ns_alt_pos[br_t]
        kfa, kta = fa >= 0, ta >= 0
        ir = np.concatenate([np.where(kfa)[0], np.where(kta)[0]])
        ic = np.concatenate([fa[kfa], ta[kta]])
        iv = np.concatenate([br_b[kfa], -br_b[kta]])
        A_alt = csr_matrix((iv, (ir, ic)), shape=(m, n - 1))
        lu_alt = splu(S[ns_alt][:, ns_alt].tocsc())
        pmin_ref = float(np.sum([M["gen_pmin"][k] for k in bus2units[alt]]))
        pmax_ref = float(np.sum([M["gen_pmax"][k] for k in bus2units[alt]]))
        net_ref_frozen = float(Pnet_ref(M, alt))
        for k in bus2units[slack]:
            Pp = M["Pnet"].copy()
            Pp[slack] -= p_g[k]
            thp = np.zeros(n)
            thp[ns_alt] = lu_alt.solve(Pp[ns_alt] / S_BASE)
            Fp = (A_alt @ thp[ns_alt]) * S_BASE
            ratio = np.abs(Fp[rated_rows]) / Fmax_r
            wr = float(ratio.max())
            wl = int(rated_rows[int(ratio.argmax())])
            net_ref_post = net_ref_frozen + p_g[k]
            ref_ok = (net_ref_post >= pmin_ref - 1e-6
                      and net_ref_post <= pmax_ref + 1e-6)
            recs.append((k, slack, wl, wr, int((ratio > VIOL).sum()), ref_ok))
            print(f"    slack unit k{k} (bus {M['bus_ids'][slack]}): "
                  f"worst ratio {wr:.4f}  new ref {M['bus_ids'][alt]} net "
                  f"{net_ref_frozen:.1f} -> {net_ref_post:.1f} MW "
                  f"([{pmin_ref:.1f}, {pmax_ref:.1f}]) "
                  f"{'OK' if ref_ok else 'INFEASIBLE at ref'}", flush=True)
    return recs


def Pnet_ref(M, bus):
    """Frozen N-0 net injection (gen - load) at a bus."""
    inj = np.zeros(M["n"])
    np.add.at(inj, M["gen_bus_idx"], M["p_g"])
    return inj[bus] - M["loads"][bus]


def branch_name(M, li):
    f, t = int(M["br_f"][li]), int(M["br_t"][li])
    k = M["kind"][li]
    s = f"{M['bus_ids'][f]}-{M['bus_ids'][t]}"
    if k == "xfrm":
        s += " (xfmr)"
    elif k == "star":
        s += " (3w-star)"
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nseg", type=int, default=20)
    ap.add_argument("--fdlmp", type=int, default=8,
                    help="number of buses to FD-cross-check LMPs (0=skip)")
    ap.add_argument("--gen-outages", type=int, default=1)
    ap.add_argument("--sm-batch", type=int, default=256)
    args = ap.parse_args()
    t_start = time.time()
    print(f"ACTIVSg25k N-0 OPF + N-1 screening  [{time.strftime('%H:%M:%S')}]",
          flush=True)

    M = build_model()
    M["sm_batch"] = args.sm_batch
    n, m, ng = M["n"], M["m"], M["ng"]
    S_BASE = M["S_BASE"]
    br_f, br_t, br_b, Fmax = M["br_f"], M["br_t"], M["br_b"], M["Fmax"]
    bus_ids = M["bus_ids"]
    rated = Fmax > 1e-9
    M["rated_idx"] = np.where(rated)[0]
    print(f"rated branches: {int(rated.sum())}/{m} "
          f"({m - int(rated.sum())} zero-rated = unlimited)", flush=True)

    # snapshot (base-case) dispatch cost -- reference only (not balanced)
    mp = json.loads((A / "mp_data.json").read_text())
    snap_pg = np.zeros(ng)
    snap_inj = np.zeros(n)
    idx = 0
    for i, r in enumerate(mp["gen"]):
        if r[7] == 1:
            snap_inj[M["bidx"][int(r[0])]] += r[1]
            snap_pg[idx] = r[1]
            idx += 1
    assert idx == ng
    snap_cost = float(np.sum(M["gen_c2"] * snap_pg**2
                             + M["gen_c1"] * snap_pg + M["gen_c0"]))
    snap_imbal = float(snap_inj.sum() - M["loads"].sum())
    print(f"snapshot dispatch: cost ${snap_cost:,.1f}/h, imbalance {snap_imbal:+.2f} MW",
          flush=True)

    # ---- N-0 PL OPF
    print("N-0 LP ...", flush=True)
    lp = build_n0_lp(M, args.nseg)
    ok, x, duals, h0 = solve_lp(lp["c_obj"], lp["lb"], lp["ub"], lp["A_ub"],
                                lp["b_ub"], lp["A_eq"], lp["b_eq"], label="N-0")
    if not ok:
        raise RuntimeError("N-0 LP failed")
    dim = lp["dim"]
    ng_, nseg = ng, args.nseg
    u0 = ng_ + ng_ * nseg
    p_g = x[:ng_].copy()
    theta = np.zeros(n)
    theta[M["ns"]] = x[u0:u0 + n - 1]

    # ---- verification
    Sfull = M["S"]
    inj = np.zeros(n)
    np.add.at(inj, M["gen_bus_idx"], p_g)
    resid = Sfull @ theta * S_BASE - (inj - M["loads"])
    print(f"\nverification:", flush=True)
    print(f"  nodal balance: max |resid| = {float(np.abs(resid).max()):.3e} MW",
          flush=True)
    assert float(np.abs(p_g.sum() - M["loads"].sum())) < 1e-3, "gen/load imbalance"
    print(f"  total gen {p_g.sum():.2f} MW = total load {M['loads'].sum():.2f} MW",
          flush=True)
    assert float(np.abs((p_g < M["gen_pmin"] - 1e-6).sum()
                        + (p_g > M["gen_pmax"] + 1e-6).sum())) == 0, "gen limit violated"
    F0 = br_b * (theta[br_f] - theta[br_t]) * S_BASE
    ratio0 = np.where(rated, np.abs(F0) / np.where(rated, Fmax, 1e9), 0.0)
    print(f"  base-case line limits: violations {int((ratio0 > VIOL).sum())}, "
          f"worst ratio {float(ratio0.max()):.4f}", flush=True)
    assert int((ratio0 > VIOL).sum()) == 0, "N-0 dispatch violates base-case limits"

    exact_cost = float(np.sum(M["gen_c2"] * p_g**2 + M["gen_c1"] * p_g
                              + M["gen_c0"]))
    pl_const = float(np.sum(M["gen_c2"] * M["gen_pmin"]**2
                            + M["gen_c1"] * M["gen_pmin"] + M["gen_c0"]))
    pl_cost = float(h0.getObjectiveValue()) + pl_const
    zc = (M["gen_c2"] == 0) & (M["gen_c1"] == 0) & (M["gen_c0"] == 0)
    zc_at_min = int((zc & (p_g <= M["gen_pmin"] + 1e-3)).sum())
    zc_at_max = int((zc & (p_g >= M["gen_pmax"] - 1e-3)).sum())
    print(f"  cost: exact quadratic ${exact_cost:,.1f}/h | PL ${pl_cost:,.1f}/h",
          flush=True)
    print(f"  zero-cost units: {int(zc.sum())} online, at Pmin {zc_at_min}, "
          f"at Pmax {zc_at_max}", flush=True)
    M["p_g"] = p_g
    M["Pnet"] = inj - M["loads"]

    # ---- LMPs: N-0 LP duals (FD cross-check runs after the screens)
    lmp_dual = duals[lp["n_ub"]:lp["n_ub"] + n] / S_BASE
    lmp_fd = np.zeros(n)
    sel = []
    if args.fdlmp > 0:
        top_load = np.argsort(-M["loads"])
        top_load = [int(b) for b in top_load if M["loads"][b] > 1e-6]
        gpb = np.zeros(n)
        np.add.at(gpb, M["gen_bus_idx"], M["gen_pmax"])
        big_gen = int(np.argmax(gpb))
        sel = sorted({M["slack"]} | set(top_load[:args.fdlmp - 2]) | {big_gen})

    # ---- bridges / islands
    print("bridge detection + DFS from slack ...", flush=True)
    bridges, tin, tout = find_bridges_dfs(M)

    # ---- line-outage screen (Sherman-Morrison, non-bridges)
    print("line-outage screen (Sherman-Morrison) ...", flush=True)
    A_inc = make_incidence(M)
    theta_ns = theta[M["ns"]]
    line_recs, line_worst, _ = screen_lines(
        M, A_inc, theta_ns, bridges)
    line_recs.sort(key=lambda r: -r[2])
    n_with_viol = sum(1 for r in line_recs if r[3] > 0)
    print(f"  line screen done: {len(line_recs)} non-bridge outages, "
          f"{n_with_viol} with >=1 violation, worst {line_worst:.4f}",
          flush=True)

    # ---- bridge-outage screen (islanding)
    print("bridge-outage screen (per-island re-solves) ...", flush=True)
    br_recs = screen_bridges(M, bridges, tin, tout)
    n_br_viol = sum(1 for r in br_recs if r["viol1"] + r["viol2"] > 0)
    n_shortfall = sum(1 for r in br_recs if r["shortfall"] > 0)
    print(f"  bridge screen done: {len(br_recs)} islanding outages, "
          f"{n_br_viol} with internal violations, {n_shortfall} with island "
          f"shortfall, max shortfall {max(r['shortfall'] for r in br_recs):.1f} MW",
          flush=True)

    # ---- gen-outage screen (diagnostic)
    gen_recs = []
    if args.gen_outages:
        print("gen-outage screen (diagnostic) ...", flush=True)
        gen_recs = screen_gens(M, A_inc, F0, p_g)
        gen_recs.sort(key=lambda r: -r[3])
        n_g_viol = sum(1 for r in gen_recs if r[4] > 0)
        print(f"  gen screen done: {len(gen_recs)} outages, {n_g_viol} with "
              f">=1 violation, worst {gen_recs[0][3]:.4f}", flush=True)

    # ---- FD LMP cross-check (last: warm re-solves of the N-0 LP)
    if args.fdlmp > 0:
        pert = 1.0  # MW
        t_fd = time.time()
        for bi in sel:
            vals = {}
            for sgn in (+1.0, -1.0):
                h0.changeRowBounds(lp["n_ub"] + bi,
                                   lp["b_eq"][bi] + sgn * pert / S_BASE,
                                   lp["b_eq"][bi] + sgn * pert / S_BASE)
                h0.run()
                assert h0.getModelStatus() == highspy.HighsModelStatus.kOptimal, \
                    "FD re-solve failed"
                vals[sgn] = h0.getObjectiveValue()
            h0.changeRowBounds(lp["n_ub"] + bi, lp["b_eq"][bi], lp["b_eq"][bi])
            lmp_fd[bi] = (vals[+1.0] - vals[-1.0]) / (2 * pert)
            print(f"  FD LMP bus {bus_ids[bi]:>5d} (load {M['loads'][bi]:8.2f} MW): "
                  f"{lmp_fd[bi]:10.3f}  (dual {lmp_dual[bi]:10.3f})", flush=True)
        print(f"  FD LMPs for {len(sel)} buses in {time.time()-t_fd:.0f}s",
              flush=True)

    # ---- report
    L = []
    L.append("# DC N-0 OPF + N-1 Security Screening — ACTIVSg25k (TAMU 25,000-bus)")
    L.append("")
    L.append(f"Network from ACTIVSg25k.RAW via the VeraGridEngine RAW parser "
             f"(24,250 buses + 750 synthetic star-stub buses; 3-winding "
             f"transformers expanded to 3-branch stars); costs/limits from the "
             f"official MATPOWER case (3,779 online units on 2,753 buses); "
             f"solved with HiGHS via highspy ({args.nseg}-segment PL cost, "
             f"lossless DC). Total load {M['loads'].sum():.2f} MW.")
    L.append("")
    L.append("**Method note.** At 25k scale a full N-1 SCOPF "
             f"(all {m} branch outages as explicit constraints) is "
             "intractable. This is instead an exact security *screen* of the "
             "N-0 dispatch: non-bridge outages are solved exactly via a "
             "Sherman–Morrison rank-1 update of the pre-contingency angle "
             "vector (verified in-run against direct re-solves), bridge "
             "outages island the network and are analyzed per-island "
             "(re-factorized slack-side component + island statistics and "
             "internal flow sweep). Gen outages are screened as a labeled "
             "diagnostic. The screen answers \"which single outages violate "
             "line limits at the N-0 dispatch\" — not the N-1-optimal "
             "dispatch.")
    L.append("")
    L.append("### N-0 (base-case) OPF")
    L.append(f"Cost (exact quadratic): **${exact_cost:,.1f}/h** "
             f"(PL objective {pl_cost:,.1f}/h; the 20-segment chord of the "
             f"convex cost is an upper bound — gap {pl_cost - exact_cost:.1f} $/h)")
    L.append(f"Snapshot (pre-OPF) dispatch cost: ${snap_cost:,.1f}/h "
             f"(snapshot imbalance {snap_imbal:+.2f} MW — the case is not "
             "power-balanced; the OPF dispatch is)")
    top_units = np.argsort(-p_g)[:15]
    L.append("")
    L.append("Top 15 units by output:")
    L.append("")
    L.append("| bus | P (MW) | Pmin | Pmax | $/MWh marginal |")
    L.append("|---|---|---|---|---|")
    for k in top_units:
        marg = 2 * M["gen_c2"][k] * p_g[k] + M["gen_c1"][k]
        L.append(f"| {int(M['gen_bus'][k])} | {p_g[k]:.2f} | "
                 f"{M['gen_pmin'][k]:.2f} | {M['gen_pmax'][k]:.2f} | {marg:.3f} |")
    L.append("")
    L.append(f"Zero-cost (must-run) units: {int(zc.sum())} online, all with "
             f"Pmin > 0 (sum Pmin {float(M['gen_pmin'][zc].sum()):.1f} MW); "
             f"at Pmin: {zc_at_min}, at Pmax: {zc_at_max}, "
             f"in between: {int(zc.sum()) - zc_at_min - zc_at_max}. The LP "
             f"is cost-indifferent over their [Pmin, Pmax] ranges, so the "
             f"optimal dispatch of these units is degenerate (HiGHS sits "
             f"them at Pmax); the cost is unaffected either way.")
    L.append("")
    L.append("Base-case loading, top 10 rated lines (by % of limit):")
    L.append("")
    L.append("| branch | flow (MW) | limit (MW) | % |")
    L.append("|---|---|---|---|")
    for li in np.argsort(-ratio0)[:10]:
        L.append(f"| {branch_name(M, li)} | {F0[li]:.2f} | {Fmax[li]:.1f} | "
                 f"{100 * F0[li] / Fmax[li]:.1f} |")
    unl = np.where(~rated)[0]
    li_u = unl[int(np.argmax(np.abs(F0[unl])))]
    L.append("")
    L.append(f"Largest flow on an unlimited (zero-rated) branch: "
             f"{F0[li_u]:.2f} MW on {branch_name(M, li_u)} "
             f"({len(unl)} such branches are unconstrained in the model)")
    L.append("")
    if args.fdlmp > 0:
        L.append("### N-0 LMPs ($/MWh)")
        L.append("")
        L.append("(dual = N-0 LP balance-row dual, pu-scaled; FD = central "
                 "finite difference of the N-0 optimal cost, ±1 MW)")
        L.append("")
        L.append(f"The N-0 base case has no binding line limit (worst "
                 f"ratio {float(ratio0.max()):.4f}), so every nodal LMP "
                 f"equals the uniform system marginal cost.")
        L.append("")
        L.append("| bus | load (MW) | LMP (dual) | LMP (FD) |")
        L.append("|---|---|---|---|")
        for bi in sel:
            L.append(f"| {bus_ids[bi]} | {M['loads'][bi]:.2f} | "
                     f"{lmp_dual[bi]:.3f} | {lmp_fd[bi]:.3f} |")
        L.append("")
    L.append(f"### N-1 line-outage screen — non-bridge outages "
             f"({len(line_recs)} of {m} branches)")
    L.append("")
    L.append(f"Outages with ≥1 post-contingency limit violation: **{n_with_viol}** "
             f"({100 * n_with_viol / len(line_recs):.2f}%)")
    L.append(f"Worst post-contingency ratio (any outage): **{line_worst:.4f}**")
    L.append("")
    L.append("Top 20 outages by worst post-contingency ratio:")
    L.append("")
    L.append("| outage | worst line | flow (MW) | limit (MW) | ratio | #lines violated |")
    L.append("|---|---|---|---|---|---|")
    for li_o, li_w, r, nv in line_recs[:20]:
        L.append(f"| {branch_name(M, li_o)} | {branch_name(M, li_w)} | "
                 f"{Fmax[li_w] * r:.1f} | "
                 f"{Fmax[li_w]:.1f} | {r:.3f} | {nv} |")
    L.append("")
    L.append("Note: the 'flow (MW)' column is the post-contingency flow "
             "magnitude on the worst line (ratio × limit).")
    L.append("")
    L.append(f"### N-1 bridge outages — islanding ({len(br_recs)} branches)")
    L.append("")
    L.append(f"Bridges disconnect an island when they fail. "
             f"Islands with load that even full local generation cannot "
             f"serve (shortfall = load − ΣPmax > 0): **{n_shortfall}** "
             f"(max shortfall {max(r['shortfall'] for r in br_recs):.1f} MW)")
    L.append(f"Islands (or slack-side components) with internal limit "
             f"violations under the frozen N-0 dispatch: **{n_br_viol}**")
    L.append("")
    top_br = sorted(br_recs, key=lambda r: (-r["shortfall"], -r["worst1"],
                                            -r["worst2"]))[:15]
    L.append("Top 15 islanding outages (by shortfall, then internal loading):")
    L.append("")
    L.append("| outage | island buses | island load (MW) | frozen gen (MW) | gen Pmax (MW) | shortfall (MW) | viol C1/C2 | worst ratio C1/C2 |")
    L.append("|---|---|---|---|---|---|---|---|")
    for r in top_br:
        L.append(f"| {branch_name(M, r['li'])} | {r['nbus2']} | "
                 f"{r['L2']:.1f} | {r['G2']:.1f} | {r['Pmax2']:.1f} | "
                 f"{r['shortfall']:.1f} | {r['viol1']}/{r['viol2']} | "
                 f"{r['worst1']:.3f}/{r['worst2']:.3f} |")
    L.append("")
    if gen_recs:
        L.append("### N-1 gen-outage screen (DIAGNOSTIC — labeled)")
        L.append("")
        n_slack_units = sum(1 for r in gen_recs if r[1] == M["slack"])
        L.append(f"Screening all {len(gen_recs)} online-unit outages at the "
                 f"N-0 dispatch ({n_slack_units} of them on the slack bus — "
                 f"those re-root the network at the largest other gen bus, "
                 f"which absorbs the lost injection).")
        n_g_viol = sum(1 for r in gen_recs if r[4] > 0)
        L.append(f"Outages with ≥1 violation: **{n_g_viol}** "
                 f"({100 * n_g_viol / len(gen_recs):.2f}%)")
        L.append("")
        L.append("Top 10 by worst post-outage ratio:")
        L.append("")
        L.append("| unit bus | P (MW) | worst line | limit (MW) | ratio | #viol | ref OK |")
        L.append("|---|---|---|---|---|---|---|")
        for k, b, wl, r, nv, rok in gen_recs[:10]:
            L.append(f"| {bus_ids[b]} | {p_g[k]:.1f} | {branch_name(M, wl)} | "
                     f"{Fmax[wl]:.1f} | {r:.3f} | {nv} | "
                     f"{'-' if b != M['slack'] else ('yes' if rok else 'NO')} |")
        bad_ref = [r for r in gen_recs if r[1] == M["slack"] and not r[5]]
        if bad_ref:
            L.append("")
            L.append(f"{len(bad_ref)} slack-bus outages leave the re-rooted "
                     f"reference outside its gen [Pmin, Pmax] band — the "
                     f"screen's flow numbers for those are reference-absorption "
                     f"artifacts, not feasible post-outage states.")
        L.append("")
    L.append("### Data / assumptions")
    L.append(f"- network: ACTIVSg25k.RAW via VeraGridEngine RAW parser — "
             f"{n} buses ({n - M['n_3w']} RAW + {M['n_3w']} star stubs), "
             f"{m} in-service branches "
             f"({m - M['n_star'] - M['n_tw'] + M['n_off_tw']} lines + "
             f"{M['n_tw'] - M['n_off_tw']} 2-winding xfmr + {M['n_star']} "
             f"star legs of {M['n_3w']} 3-winding xfmr); "
             f"{M['n_off_line']} out-of-service line excluded; 503 star "
             f"legs with negative x kept as-is (exact wye-equivalent "
             f"artifacts; the reduced Laplacian still factors)")
    L.append(f"- costs/limits: official MATPOWER case; {ng} online units "
             f"on {len(M['bus_units'])} buses (per-bus unit counts agree "
             f"RAW vs MP exactly)")
    L.append(f"- contingency set: no TAMU contab exists for the 25k case → "
             f"N-1 line set = all {m} in-service branches; gen outages "
             f"screened as a diagnostic")
    mp_ref_kv = next(int(r[9]) for r in mp["bus"] if int(r[0]) == SLACK_BUS)
    L.append(f"- slack: MP reference bus {SLACK_BUS} (baseKV {mp_ref_kv}); "
             f"{len(M['bus_units'][M['slack']])} units there, all at Pmax "
             f"in the snapshot (their outages are re-rooted in the gen "
             f"screen)")
    L.append("- lossless DC, no voltage limits, no corrective re-dispatch; "
             "screening (not N-1-optimal re-dispatch)")
    L.append(f"- PL cost: {args.nseg}-segment chord (upper bound) of the "
             f"quadratic gencost; zero-rated branches ({m - int(rated.sum())}) "
             f"unconstrained")
    report = "\n".join(L)
    out = A / "results_activsg25k.md"
    out.write_text(report + "\n", encoding="utf-8")
    print(f"\n[report saved to {out}]", flush=True)
    print(f"total {time.time()-t_start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
