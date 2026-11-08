#!/usr/bin/env python3
"""
Cross-check the ACTIVSg500 SCOPF results against the OFFICIAL dataset.

The answer to "how do I know the solver results are correct?":
independent of the saved report, the official case file carries its own
reference operating point (the base-case generation in mpc.gen). We verify:

  [1] bus / status / Pmin / Pmax / cost coefficients read straight from the
      raw case file agree exactly with the driver's model.
  [2] The published base-case dispatch costs ~$72.8k/h; our cost-optimal
      N-0 is $70.8k/h, i.e. 2.7% CHEAPER. A cost optimization must be <=
      any reference dispatch, and ~2.7% is a sane gap (the published base
      case is a power-flow snapshot, not a cost-optimal dispatch). If our
      N-0 were, say, $10k/h, something (cost mapping, balance, limits)
      would be badly wrong.
  [3] Recompute the N-0 LP from scratch (independent HiGHS solve): the
      exact-quadratic cost must reproduce $70,793.0/h, all 36 zero-cost
      units must be at Pmax (free power is always taken), and nodal
      balance must hold.

Run:  python verify_basecase_cost.py        (~1-2 min, incl. K precompute)
"""
import re
import time

import numpy as np
from scipy import sparse

import activsg500_n1_opf as d

A = d.A
text = open(A / "case_ACTIVSg500.m", encoding="latin-1").read()


def parse(name):
    block = re.search(rf"mpc\.{name}\s*=\s*\[(.*?)\];", text, re.S).group(1)
    rows = []
    for line in block.strip().splitlines():
        rows.append([float(x) for x in
                     re.findall(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?", line)])
    return np.array(rows, dtype=float)


G = parse("gen")
GC = parse("gencost")

# raw-file column mapping (verified column by column):
#   bus=r0, status=r7, Pmax=r8, Pmin=r9, base-case Pg=r1
bus = G[:, 0].astype(int)
status = G[:, 7].astype(int)
pmax, pmin, pg_base = G[:, 8], G[:, 9], G[:, 1]
c2, c1, c0 = GC[:, 4], GC[:, 5], GC[:, 6]
online = status == 1

# ---- [1] raw file vs driver model
M = d.build_model(verbose=False)
mb = np.array(M["gen_buses"])
m_pmin = np.array(M["gen_pmin"])
m_pmax = np.array(M["gen_pmax"])
m_c = np.column_stack([np.array(M["gen_c2"]), np.array(M["gen_c1"]),
                       np.array(M["gen_c0"])])
raw = np.array([[pmin[bus == int(b)][0], pmax[bus == int(b)][0],
                 c2[bus == int(b)][0], c1[bus == int(b)][0],
                 c0[bus == int(b)][0]] for b in mb])
agree = (np.allclose(m_pmin, raw[:, 0]) and np.allclose(m_pmax, raw[:, 1])
         and np.allclose(m_c, raw[:, 2:5]))
print(f"[1] bus set match                : {set(mb) == set(bus[online])}")
print(f"[1] Pmin/Pmax/cost match         : {agree}")
print(f"    sum Pmax  raw vs model       : "
      f"{pmax[online].sum():.2f} vs {m_pmax.sum():.2f}")
print(f"    sum Pmin  raw vs model       : "
      f"{pmin[online].sum():.2f} vs {m_pmin.sum():.2f}")
assert set(mb) == set(bus[online]) and agree

# ---- [2] published base-case cost vs our N-0
base_cost = float(np.sum(c2[online] * pg_base[online] ** 2
                         + c1[online] * pg_base[online] + c0[online]))
N0, N1 = 70793.0, 80639.3
print(f"\n[2] published base-case dispatch : "
      f"sum Pg = {pg_base[online].sum():.2f} MW "
      f"(load {M['loads'].sum():.2f} MW)")
print(f"    base-case cost (exact quad)  : ${base_cost:,.1f}/h")
print(f"    our N-0 optimum              : ${N0:,.1f}/h")
print(f"    our N-1 (591 line cts)       : ${N1:,.1f}/h")
assert N0 <= base_cost, "N-0 must not exceed the base-case cost"
assert N1 >= N0, "N-1 must not be cheaper than N-0"
print(f"    OK: N-0 is {100*(base_cost-N0)/base_cost:.2f}% below the "
      f"published base case (expected: the base case is not cost-optimized)")

# ---- [3] independent N-0 re-solve
n, m, ng = M["n"], M["m"], M["ng"]
br, E, Sfull, loads = M["br"], M["E"], M["Sfull"], M["loads"]
ns, ns_pos = M["ns"], M["ns_pos"]
Fmax_all, unlimited = M["Fmax_all"], M["unlimited"]
S_BASE = M["S_BASE"]
nseg = 20
seg_price, seg_width = d.build_pl(M, nseg)
dim = ng + ng * nseg + len(ns)
s0, t0 = ng, ng + ng * nseg
c_obj = np.zeros(dim)
c_obj[s0:t0] = seg_price

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
A_ub = sparse.csr_matrix(np.array(rows))

print("\n[3] re-solving the N-0 LP from scratch ...")
t0w = time.time()
ok, x, duals, h = d.solve_lp(c_obj, lb, ub, A_ub, np.array(rhs),
                             A_eq, np.array(b_eq), label="verify N-0")
assert ok, "independent N-0 re-solve failed"
p_g0 = x[:ng]
exact0 = float(sum(M["gen_c2"][k] * p_g0[k] ** 2
                   + M["gen_c1"][k] * p_g0[k] + M["gen_c0"][k]
                   for k in range(ng)))
zc = (c2 == 0) & (c1 == 0) & (c0 == 0)
zc_bus = [int(b) for b in bus[online] if zc[bus == b].any()]
k_zc = [int(np.where(mb == b)[0][0]) for b in zc_bus]
n_at_pmax = int((p_g0[k_zc] >= m_pmax[k_zc] - 1e-6).sum())
th = np.zeros(n)
th[ns] = x[t0:]
resid = Sfull @ th * S_BASE - (E @ p_g0 - loads)
bal = float(np.abs(resid).max())
print(f"    solved in {time.time()-t0w:.1f}s")
print(f"    exact-quadratic cost         : ${exact0:,.1f}/h (report: "
      f"${N0:,.1f}/h)")
print(f"    zero-cost units at Pmax      : {n_at_pmax}/{len(zc_bus)}")
print(f"    max nodal-balance residual   : {bal:.3e} MW")
assert abs(exact0 - N0) < 2.0, "re-solved N-0 cost does not match the report"
assert n_at_pmax == len(zc_bus), "some zero-cost unit is not at Pmax"
assert bal < 1e-6

print("\nALL CHECKS PASSED")
