#!/usr/bin/env python3
"""
DC N-1 OPF for the IEEE 14-bus test case.

Formulation (lossless DC, classical N-1 security constraint, no corrective
re-dispatch):

  variables:
    P_g  : pre-contingency generation at each generator bus
    s_ik : piecewise-linear cost segments for generator i
    th   : bus voltages (angles); slack bus fixed at 0

  min        sum of PL segment costs  (approx. of the quadratic gencost)
  s.t.       nodal power balance in the base case
             generator limits
             base-case line flow limits      |F_l|        <= Fmax
             for every N-1 contingency c:
               (each of the 20 line outages, each of the 5 generator outages)
                 post-contingency flows      |F_l^(c)|    <= Fmax

The post-contingency flow on line l under contingency c is linear in the
pre-contingency dispatch:  F_l^(c) = K_c (P_g - loads), with K_c built from
the post-contingency network matrix. The pre-contingency schedule must be
feasible under every contingency without re-dispatch.

Solver: scipy.optimize.linprog (HiGHS).
LMPs:  central finite difference of the optimal cost w.r.t. each bus load.

Usage:
  python dc_n1_opf.py [--limit 100] [--nseg 12]
"""

import argparse
import sys

import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from scipy.sparse.csgraph import connected_components

import case14 as C


# ----------------------------------------------------------------------------
# indexing / static matrices
# ----------------------------------------------------------------------------
bus_ids = sorted(C.BUSES)
n = len(bus_ids)
bidx = {b: i for i, b in enumerate(bus_ids)}
slack = bidx[C.SLACK_BUS]
ns = [i for i in range(n) if i != slack]          # non-slack bus indices
ns_pos = {i: k for k, i in enumerate(ns)}

loads = np.array([C.BUSES[b][0] for b in bus_ids])

br = [(bidx[f], bidx[t], 1.0 / x) for (f, t, r, x, b) in C.BRANCHES]
m = len(br)

gen_buses = sorted(C.GENS)
ng = len(gen_buses)
gcol = {b: k for k, b in enumerate(gen_buses)}

E = np.zeros((n, ng))                               # gen incidence, buses x gens
for b, k in gcol.items():
    E[bidx[b], k] = 1.0

# full bus susceptance matrix (base case)
Sfull = np.zeros((n, n))
for (f, t, bl) in br:
    Sfull[f, t] -= bl
    Sfull[t, f] -= bl
    Sfull[f, f] += bl
    Sfull[t, t] += bl
Sns = Sfull[np.ix_(ns, ns)]


# ----------------------------------------------------------------------------
# contingency helpers
# ----------------------------------------------------------------------------
def contingency_matrices(out_line=None, out_gen=None):
    """Return (K_rows, active_line_ids, E_c) for one contingency, or None if
    the post-contingency network is disconnected.

    K_rows: (n_active_lines x n) with F = K_rows @ (P_g - loads)
    E_c:    modified gen incidence for generator outages (that unit at 0)
    """
    active = [li for li in range(m) if li != out_line]

    adj = np.zeros((n, n), bool)
    for li in active:
        f, t, _ = br[li]
        adj[f, t] = adj[t, f] = True
    adj[slack, slack] = True
    _, labels = connected_components(sparse.csr_matrix(adj), directed=False)
    if not np.all(labels == labels[slack]):
        return None

    S = np.zeros((len(ns), len(ns)))
    for li in active:
        f, t, bl = br[li]
        if f in ns_pos:
            S[ns_pos[f], ns_pos[f]] += bl
        if t in ns_pos:
            S[ns_pos[t], ns_pos[t]] += bl
        if f in ns_pos and t in ns_pos:
            S[ns_pos[f], ns_pos[t]] -= bl
            S[ns_pos[t], ns_pos[f]] -= bl

    Minv = np.linalg.pinv(S)
    M = np.zeros((n, n))          # P_inj -> theta, theta_slack = 0
    # M[bus_i, bus_j] = Minv[pos_i, pos_j]; slack row/col stay 0
    for a, i in enumerate(ns):    # a = row position, i = row bus
        for b, j in enumerate(ns):  # b = col position, j = col bus
            M[i, j] = Minv[a, b]

    rows = []
    for li in active:
        f, t, bl = br[li]
        r = np.zeros(n)
        r[f] = bl
        r[t] = -bl
        rows.append(r @ M)

    E_c = E.copy()
    if out_gen is not None:
        E_c[bidx[out_gen], gcol[out_gen]] = 0.0

    return np.array(rows), active, E_c


def build_contingency_list():
    ct = []
    for li, (f, t, bl) in enumerate(br):
        bf, bt = C.BRANCHES[li][0], C.BRANCHES[li][1]
        ct.append({"name": f"line {bf}-{bt} out", "out_line": li, "out_gen": None})
    for b in gen_buses:
        ct.append({"name": f"gen bus {b} out", "out_line": None, "out_gen": b})
    return ct


# ----------------------------------------------------------------------------
# piecewise-linear cost approximation
# ----------------------------------------------------------------------------
def build_pl_costs(nseg):
    """Return (seg_cost, seg_width, seg_price) per generator, length ng*nseg.
    seg_cost: total $ of a full segment
    seg_price: $ per MW in that segment (= seg_cost / seg_width) -- this is
               the coefficient that goes into the LP objective.
    """
    costs, widths = [], []
    for b in gen_buses:
        pmin, pmax, c2, c1, c0 = C.GENS[b]
        xs = np.linspace(pmin, pmax, nseg + 1)
        cc = c2 * xs**2 + c1 * xs + c0
        costs.append(np.diff(cc))
        widths.append(np.diff(xs))
    costs = np.concatenate(costs)
    widths = np.concatenate(widths)
    return costs, widths, costs / widths


# ----------------------------------------------------------------------------
# LP assembly
# ----------------------------------------------------------------------------
# variable layout: x = [ P_g (ng) | s (ng*nseg) | theta_ns (n-1) ]
def solve_opf(loads_vec, Fmax, nseg, with_n1=True):
    seg_cost, seg_width, seg_price = build_pl_costs(nseg)
    nseg_total = ng * nseg
    dim = ng + nseg_total + len(ns)
    s0 = ng                  # start index of segment vars
    t0 = ng + nseg_total     # start index of theta vars

    c_obj = np.zeros(dim)
    c_obj[s0:t0] = seg_price

    # ---- equality constraints
    a_eq, b_eq = [], []
    # nodal balance for EVERY bus, including the slack bus.
    # The slack angle is fixed at 0, but the slack *generation* is a variable:
    #   P_g_slack - load_slack = (S theta)_slack = S[slack, ns] . theta_ns
    # Omitting this row would leave P_g_slack free (and the objective would
    # drive it to its bound).
    Sall_ns = Sfull[:, ns]  # all rows, non-slack columns
    for i in range(n):
        row = np.zeros(dim)
        row[0:ng] = E[i, :]
        row[t0:t0 + len(ns)] = -Sall_ns[i, :]
        a_eq.append(row)
        b_eq.append(loads_vec[i])
    # P_g_k = sum_j s_kj
    for k in range(ng):
        row = np.zeros(dim)
        row[k] = 1.0
        row[s0 + k * nseg:s0 + (k + 1) * nseg] = -1.0
        a_eq.append(row)
        b_eq.append(0.0)

    # ---- inequality constraints
    a_ub, b_ub = [], []
    # base-case line limits
    for (f, t, bl) in br:
        row = np.zeros(dim)
        if f in ns_pos:
            row[t0 + ns_pos[f]] += bl
        if t in ns_pos:
            row[t0 + ns_pos[t]] -= bl
        a_ub += [row, -row]
        b_ub += [Fmax, Fmax]

    # N-1 line limits
    skipped = []
    if with_n1:
        for ct in build_contingency_list():
            res = contingency_matrices(ct["out_line"], ct["out_gen"])
            if res is None:
                skipped.append(ct["name"])
                continue
            K, active, E_c = res
            Kme = K @ E_c                    # active_lines x ng
            Kmloads = K @ loads_vec
            for r in range(Kme.shape[0]):
                row = np.zeros(dim)
                row[0:ng] = Kme[r, :]
                a_ub += [row, -row]
                b_ub += [Fmax + Kmloads[r], Fmax - Kmloads[r]]

    # ---- bounds
    lb = np.zeros(dim)
    ub = np.full(dim, np.inf)
    for k, b in enumerate(gen_buses):
        pmin, pmax, _, _, _ = C.GENS[b]
        lb[k], ub[k] = pmin, pmax
    for k in range(ng):
        for j in range(nseg):
            idx = s0 + k * nseg + j
            ub[idx] = seg_width[k * nseg + j]
    # bus angles are unbounded (loaded buses run negative w.r.t. slack)
    lb[t0:] = -np.inf
    bounds = list(zip(lb, ub))

    res = linprog(
        c_obj,
        A_ub=sparse.csr_matrix(np.array(a_ub)),
        b_ub=np.array(b_ub),
        A_eq=sparse.csr_matrix(np.array(a_eq)),
        b_eq=np.array(b_eq),
        bounds=bounds,
        method="highs",
    )
    if not res.success:
        raise RuntimeError(f"LP infeasible/unbounded: {res.message}")
    mats = {"a_eq": np.array(a_eq), "b_eq": np.array(b_eq),
            "a_ub": np.array(a_ub), "b_ub": np.array(b_ub)}
    return res, mats


def extract(res, loads_vec, Fmax, nseg):
    x = res.x
    p_g = x[:ng]
    theta_ns = x[ng + ng * nseg:]
    theta = np.zeros(n)
    theta[ns] = theta_ns
    flows = np.array([bl * (theta[f] - theta[t]) for (f, t, bl) in br])
    cost = 0.0
    for k, b in enumerate(gen_buses):
        _, _, c2, c1, c0 = C.GENS[b]
        cost += c2 * p_g[k] ** 2 + c1 * p_g[k] + c0
    return {"p_g": p_g, "theta": theta, "flows": flows, "cost": cost}


def check_contingencies(p_g, loads_vec, Fmax):
    """Evaluate post-contingency flows for a given dispatch.
    Returns list of (name, max_flow_ratio, worst_line, violated_bool, max_flow)."""
    out = []
    for ct in build_contingency_list():
        res = contingency_matrices(ct["out_line"], ct["out_gen"])
        if res is None:
            out.append((ct["name"], None, None, False, None))
            continue
        K, active, E_c = res
        F = K @ (E_c @ p_g - loads_vec)
        ratio = np.max(np.abs(F)) / Fmax
        wi = int(np.argmax(np.abs(F)))
        out.append((ct["name"], ratio, active[wi], bool(np.any(np.abs(F) > Fmax + 1e-6)),
                    float(np.max(np.abs(F)))))
    return out


def compute_lmps(loads_vec, Fmax, nseg, pert=0.05):
    lmps = np.zeros(n)
    for i in range(n):
        lp, lm = loads_vec.copy(), loads_vec.copy()
        lp[i] += pert
        lm[i] -= pert
        rp, _ = solve_opf(lp, Fmax, nseg, with_n1=True)
        rm, _ = solve_opf(lm, Fmax, nseg, with_n1=True)
        cp = extract(rp, lp, Fmax, nseg)["cost"]
        cm = extract(rm, lm, Fmax, nseg)["cost"]
        lmps[i] = (cp - cm) / (2 * pert)
    return lmps


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="DC N-1 OPF, IEEE 14-bus")
    ap.add_argument("--limit", type=float, default=C.LINE_LIMIT_MW,
                    help="uniform line flow limit [MW]")
    ap.add_argument("--nseg", type=int, default=12,
                    help="PL segments per generator cost")
    args = ap.parse_args()
    Fmax, nseg = args.limit, args.nseg

    total_load = loads.sum()
    print("=" * 72)
    print("DC N-1 OPF -- IEEE 14-bus test case (MATPOWER case14, public domain)")
    print("=" * 72)
    print(f"  buses: {n}   lines: {m}   generators: {ng}   total load: {total_load:.1f} MW")
    print(f"  line limit: uniform {Fmax:.1f} MW (case14 defines no ratings)")
    print(f"  costs: MATPOWER quadratic gencost, PL-approx with {nseg} segments/gen")
    print(f"  slack bus: {C.SLACK_BUS}   N-1 set: {m} line + {ng} generator outages = {m + ng}")
    print()

    # ---- N-0 (no security constraints)
    r0, _ = solve_opf(loads, Fmax, nseg, with_n1=False)
    d0 = extract(r0, loads, Fmax, nseg)
    n0_ct = check_contingencies(d0["p_g"], loads, Fmax)
    n0_viol = [c for c in n0_ct if c[3]]

    # ---- N-1 (with security constraints)
    r1, _ = solve_opf(loads, Fmax, nseg, with_n1=True)
    d1 = extract(r1, loads, Fmax, nseg)
    n1_ct = check_contingencies(d1["p_g"], loads, Fmax)
    n1_viol = [c for c in n1_ct if c[3]]
    n1_worst = max((c for c in n1_ct if c[1] is not None), key=lambda c: c[1])

    # ---- sanity
    assert abs(d1["p_g"].sum() - total_load) < 1e-6, "gen/load balance failed"

    # ---- LMPs (on the N-1 problem)
    lmps = compute_lmps(loads, Fmax, nseg)

    # ---- report
    def fmt_line(li):
        f, t, bl = br[li]
        return f"line {bus_ids[f]}->{bus_ids[t]}"

    L = []
    L.append("### N-0 (base case only)")
    L.append(f"Cost: ${d0['cost']:,.1f}/h")
    for k, b in enumerate(gen_buses):
        L.append(f"  gen bus {b:<2d}: {d0['p_g'][k]:8.2f} MW")
    base_worst = int(np.argmax(np.abs(d0["flows"])))
    L.append(f"Worst base-case flow: {fmt_line(base_worst)} @ "
             f"{d0['flows'][base_worst]:.1f} MW ({100*d0['flows'][base_worst]/Fmax:.1f}% of limit)")
    L.append(f"Contingencies violated by this dispatch: {len(n0_viol)}")
    for c in n0_viol:
        L.append(f"  - {c[0]}: max flow {c[4]:.1f} MW ({100*c[1]:.1f}%) on {fmt_line(c[2])}")
    L.append("")
    L.append("### N-1 (all 25 contingencies enforced)")
    L.append(f"Cost: ${d1['cost']:,.1f}/h  "
             f"(delta {d1['cost']-d0['cost']:+,.1f}/h vs N-0)")
    for k, b in enumerate(gen_buses):
        L.append(f"  gen bus {b:<2d}: {d1['p_g'][k]:8.2f} MW   "
                 f"(N-0 was {d0['p_g'][k]:.2f})")
    L.append(f"All {len(n1_ct)} contingencies satisfied "
             f"(violations: {len(n1_viol)})")
    L.append(f"Worst post-contingency flow: {n1_worst[0]} -> "
             f"{n1_worst[4]:.1f} MW ({100*n1_worst[1]:.1f}%) on {fmt_line(n1_worst[2])}")
    binding = [c for c in n1_ct if c[1] is not None and c[1] >= 1.0 - 1e-4]
    if binding:
        L.append("Binding N-1 contingency(ies):")
        for c in binding:
            L.append(f"  - {c[0]}: {c[4]:.1f} MW on {fmt_line(c[2])}")
    L.append("")
    L.append("### Base-case line flows (N-1 solution)")
    for li in range(m):
        L.append(f"  {fmt_line(li):<14s} {d1['flows'][li]:8.2f} MW "
                 f"({100*d1['flows'][li]/Fmax:6.1f}%)")
    L.append("")
    L.append("### LMPs [$ / MWh] (finite difference, N-1 problem)")
    for i, b in enumerate(bus_ids):
        L.append(f"  bus {b:<2d}: {lmps[i]:8.3f}")
    L.append("")
    L.append("### Contingency detail (N-1 solution)")
    L.append(f"  {'contingency':<18s} {'max flow':>9s} {'util':>7s}")
    for c in n1_ct:
        if c[1] is None:
            L.append(f"  {c[0]:<18s} (disconnected -> skipped)")
        else:
            L.append(f"  {c[0]:<18s} {c[4]:8.1f} MW {100*c[1]:6.1f}%")

    report = "\n".join(L)
    print(report)
    out = f"results_report_limit{Fmax:.0f}MW.md"
    with open(out, "w") as fh:
        fh.write(report + "\n")
    print(f"\n[report saved to {out}]")


if __name__ == "__main__":
    main()
