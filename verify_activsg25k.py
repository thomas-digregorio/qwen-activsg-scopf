"""Independent verification of the 25k driver results:
 1. N-0: re-solve the PL OPF, recompute exact cost, balance, limits
    (incl. check that the PL chord cost is an upper bound on the exact)
 2. Line screen: worst outages checked by DIRECT re-solve (no SM)
 3. Island stats: recompute island load / gen Pmax / frozen gen from the
    RAW load blocks, the MP case, and the re-solved dispatch (validates
    DFS masks + bus mapping + aggregation)
 4. Gen screen: direct re-solve for the worst gen outages
Runs from its own directory (repo root or OPF-scratch) next to
activsg25k_n1_screen.py.
"""
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "drv", str(HERE / "activsg25k_n1_screen.py"))
drv = importlib.util.module_from_spec(spec)
sys.modules["drv"] = drv
spec.loader.exec_module(drv)

M = drv.build_model()
n, m, ng = M["n"], M["m"], M["ng"]
S_BASE = M["S_BASE"]
br_f, br_t, br_b, Fmax = M["br_f"], M["br_t"], M["br_b"], M["Fmax"]
rated = Fmax > 1e-9
M["rated_idx"] = np.where(rated)[0]
mp = json.loads((drv.A / "mp_data.json").read_text())

# ---- 1. N-0 re-solve
lp = drv.build_n0_lp(M, 20)
t0 = time.time()
ok, x, duals, h = drv.solve_lp(lp["c_obj"], lp["lb"], lp["ub"], lp["A_ub"],
                               lp["b_ub"], lp["A_eq"], lp["b_eq"], label="V-N0")
assert ok, "N-0 re-solve failed"
p_g = x[:ng]
M["p_g"] = p_g
theta = np.zeros(n)
theta[M["ns"]] = x[ng + ng * 20:]
inj = np.zeros(n)
np.add.at(inj, M["gen_bus_idx"], p_g)
res = M["S"] @ theta * S_BASE - (inj - M["loads"])
exact = float(np.sum(M["gen_c2"] * p_g**2 + M["gen_c1"] * p_g + M["gen_c0"]))
pl_const = float(np.sum(M["gen_c2"] * M["gen_pmin"]**2
                        + M["gen_c1"] * M["gen_pmin"] + M["gen_c0"]))
pl = float(h.getObjectiveValue()) + pl_const
F0 = br_b * (theta[br_f] - theta[br_t]) * S_BASE
viol0 = int((np.abs(F0[rated]) / Fmax[rated] > 1 + 1e-9).sum())
print(f"\n1. N-0 re-solve: max|balance res|={float(np.abs(res).max()):.3e} MW, "
      f"gen-load gap={float(p_g.sum() - M['loads'].sum()):.3e} MW, "
      f"limit viol={viol0}  [{time.time() - t0:.0f}s]")
print(f"   exact quadratic ${exact:,.4f}/h | PL ${pl:,.4f}/h | "
      f"PL - exact = {pl - exact:+.6f}  (theory: chord of convex cost is "
      f"an upper bound, so PL >= exact)")
assert pl - exact > -1e-6, "PL below exact: chord property violated"
assert pl - exact < 0.01 * exact, "PL gap implausibly large"
assert abs(p_g.sum() - M["loads"].sum()) < 1e-3
assert viol0 == 0
assert np.abs(res).max() < 1e-3
assert (p_g < M["gen_pmin"] - 1e-6).sum() == 0
assert (p_g > M["gen_pmax"] + 1e-6).sum() == 0

# ---- 2. line screen: worst outages by direct re-solve
bridges, tin, tout = drv.find_bridges_dfs(M)
A_inc = drv.make_incidence(M)
theta_ns = theta[M["ns"]]
recs, worst, nviol = drv.screen_lines(M, A_inc, theta_ns, bridges, n_check=0)
recs.sort(key=lambda r: -r[2])
print(f"\n2. line screen re-run: {len(recs)} outages, "
      f"{nviol} with violations, worst {worst:.4f}")
Pnet = inj - M["loads"]
M["Pnet"] = Pnet
S = M["S"]
ns = M["ns"]
for li_o, li_w, r_sm, nv in recs[:5]:
    f, t2, bl = int(br_f[li_o]), int(br_t[li_o]), br_b[li_o]
    dmat = coo_matrix(([-bl, -bl, bl, bl],
                       ([f, t2, f, t2], [t2, f, f, t2])), shape=(n, n))
    S2 = (S - dmat).tocsr()
    th = np.zeros(n)
    th[ns] = splu(S2[ns][:, ns].tocsc()).solve(Pnet[ns] / S_BASE)
    F = br_b * (th[br_f] - th[br_t]) * S_BASE
    F[li_o] = 0.0
    ratio = np.where(rated, np.abs(F) / np.where(rated, Fmax, 1e9), 0.0)
    w = float(ratio.max())
    wi = int(np.argmax(ratio))
    status = "OK" if abs(w - r_sm) < 1e-6 else "MISMATCH"
    print(f"   ct {li_o} ({drv.branch_name(M, li_o)}): screen worst "
          f"{r_sm:.6f} (line {li_w}) vs direct {w:.6f} (line {wi}) "
          f"-> {status}")
    assert status == "OK" and wi == li_w
    # violated-line count check
    nv_direct = int((ratio[rated] > drv.VIOL).sum())
    assert nv_direct == nv, (nv_direct, nv)

# ---- 3. island stats: validate mask + bus mapping + aggregation
# NOTE: RAW and MP are different-time snapshots, so island load is NOT
# compared against MP loads. Instead: load re-summed from RAW load blocks
# (vg_data.json) independently of the model array; gen Pmax re-summed
# from the MP case (official source); frozen gen re-summed from the
# re-solved N-0 dispatch via the per-bus unit map.
gen_pmax_bus = np.zeros(n)
np.add.at(gen_pmax_bus, M["gen_bus_idx"], M["gen_pmax"])
rows = []
for li in sorted(bridges):
    f, t2 = int(br_f[li]), int(br_t[li])
    c = f if tin[f] > tin[t2] else t2
    mask2 = (tin >= tin[c]) & (tin < tout[c])
    L2 = float(M["loads"][mask2].sum())
    Pmax2 = float(gen_pmax_bus[mask2].sum())
    rows.append((max(0.0, L2 - Pmax2), Pmax2, li, c, mask2))
rows.sort(key=lambda r: (-r[0], -r[1]))
sub = [r[2] for r in rows[:8]]
genrich = [r[2] for r in rows if r[1] > 100.0][:4]
sub = list(set(sub + genrich))
print(f"\n3. island stats check ({len(sub)} bridges, "
      f"top by shortfall + {len(genrich)} gen-rich):")
vg = json.loads((drv.A / "vg_data.json").read_text())
raw_load_by_bus = {}
for l in vg["loads"]:
    if l["STATUS"] in (1, None):
        b = int(l["I"])
        raw_load_by_bus[b] = raw_load_by_bus.get(b, 0.0) + (l["PL"] or 0.0)
for b in vg["buses"]:
    if b.get("BL"):
        raw_load_by_bus[int(b["I"])] = \
            raw_load_by_bus.get(int(b["I"]), 0.0) + b["BL"]
br_recs = {r["li"]: r for r in drv.screen_bridges(M, sub, tin, tout)}
mp_bus = {int(r[0]): r for r in mp["bus"]}
mp_pmax = {}
for r in mp["gen"]:
    b = int(r[0])
    if r[7] == 1:
        mp_pmax[b] = mp_pmax.get(b, 0.0) + r[8]
inj2 = np.zeros(n)
np.add.at(inj2, M["gen_bus_idx"], p_g)
for li in sub:
    r = br_recs[li]
    f, t2 = int(br_f[li]), int(br_t[li])
    c = f if tin[f] > tin[t2] else t2
    mask2 = (tin >= tin[c]) & (tin < tout[c])
    raw_ids = [M["bus_ids"][i] for i in np.nonzero(mask2)[0]
               if M["bus_ids"][i] < 10**6]
    load_raw = sum(raw_load_by_bus.get(b, 0.0) for b in raw_ids)
    pmax_mp = sum(mp_pmax.get(b, 0.0) for b in raw_ids)
    g_model = float(inj2[mask2].sum())
    ok_ld = abs(load_raw - r["L2"]) < 1e-6
    ok_pm = abs(pmax_mp - r["Pmax2"]) < 1e-6
    ok_g = abs(g_model - r["G2"]) < 1e-3
    print(f"   ct {li} ({drv.branch_name(M, li)}): island {r['nbus2']} buses | "
          f"load model {r['L2']:.3f} vs RAW-resum {load_raw:.3f} "
          f"{'OK' if ok_ld else 'MISMATCH'} | Pmax model {r['Pmax2']:.3f} "
          f"vs MP-resum {pmax_mp:.3f} {'OK' if ok_pm else 'MISMATCH'} | "
          f"frozen gen model {g_model:.3f} vs N0-resum {r['G2']:.3f} "
          f"{'OK' if ok_g else 'MISMATCH'} | shortfall {r['shortfall']:.3f}")
    assert ok_ld and ok_pm and ok_g

# ---- 4. gen screen: direct re-solve for the 3 worst
gen_recs = drv.screen_gens(M, A_inc, F0, p_g)
gen_recs.sort(key=lambda r: -r[3])
ns_pos = M["ns_pos"]
print(f"\n4. gen screen re-run: {len(gen_recs)} outages, worst {gen_recs[0][3]:.4f}")
slack = M["slack"]
checked = 0
for k, bus, wl, r_sm, nv, rok in gen_recs[:10]:
    if bus == slack:
        continue  # slack outages are re-rooted; direct e-bus check invalid
    d = np.zeros(n)
    d[bus] = 1.0
    wvec = M["lu"].solve(d[ns])  # reduced (n-1) space
    wfull = np.zeros(n)
    wfull[ns] = wvec
    th = theta + (-p_g[k] / S_BASE) * wfull
    F = br_b * (th[br_f] - th[br_t]) * S_BASE
    ratio = np.where(rated, np.abs(F) / np.where(rated, Fmax, 1e9), 0.0)
    wd = float(ratio.max())
    status = "OK" if abs(wd - r_sm) < 1e-6 else "MISMATCH"
    print(f"   unit k{k} bus {M['bus_ids'][bus]} (P {p_g[k]:.2f} MW): "
          f"screen {r_sm:.6f} vs direct {wd:.6f} -> {status}")
    assert status == "OK"
    checked += 1
    if checked >= 3:
        break

print("\nALL VERIFICATIONS PASSED")
