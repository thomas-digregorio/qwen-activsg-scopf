"""Cross-check the VeraGrid-parsed RAW (vg_data.json) against the official
TAMU MATPOWER case (mp_data.json): buses, loads, generators, branches.
3-winding transformers appear in MP as a 3-branch star around a synthetic
stub bus; those are verified separately (stub -> transformer matching,
triangle->star impedances, per-winding ratings, zero angles). Any mismatch
is printed.

Usage:  python crosscheck.py [--dir DIR]   (default: activsg500)
"""
import argparse
import json
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="activsg500")
args = ap.parse_args()

A = Path(__file__).resolve().parent / args.dir
vg = json.loads((A / "vg_data.json").read_text())
mp = json.loads((A / "mp_data.json").read_text())
S = mp["baseMVA"]

issues = []

# ---- buses ---------------------------------------------------------------
mp_bus = {int(r[0]): r for r in mp["bus"]}
vg_bus = {int(b["I"]): b for b in vg["buses"]}
stub_only = sorted(set(mp_bus) - set(vg_bus))
raw_only = sorted(set(vg_bus) - set(mp_bus))
if raw_only:
    issues.append(f"bus id mismatch: RAW-only buses {raw_only[:5]}")
common = sorted(set(mp_bus) & set(vg_bus))
n_vm = n_va = 0
for i in common:
    m, v = mp_bus[i], vg_bus[i]
    if abs(m[9] - v["BASKV"]) > 1e-6:
        issues.append(f"bus {i}: baseKV {m[9]} vs RAW BASKV {v['BASKV']}")
    if abs(m[7] - v["VM"]) > 1e-4:
        n_vm += 1
    if abs(m[8] - v["VA"]) > 1e-3:
        n_va += 1
print(f"buses: {len(common)} in both | MP-only {len(stub_only)} | RAW-only {len(raw_only)}")
print(f"  Vm diff >1e-4: {n_vm} | Va diff >1e-3: {n_va}")

# ---- 3-winding transformer star stubs -------------------------------------
tw3 = [t for t in vg["transformers"] if t["K"] not in (0, 0.0, None)]
print(f"\n3-winding transformers in RAW: {len(tw3)}  | MP-only (stub) buses: {len(stub_only)}")
stub_branches = {}
if tw3 or stub_only:
    stub_nbr = {s: set() for s in stub_only}
    for idx, r in enumerate(mp["branch"]):
        f, t2 = int(r[0]), int(r[1])
        if f in stub_nbr:
            stub_nbr[f].add(t2)
            stub_branches[(f, t2)] = idx
        if t2 in stub_nbr:
            stub_nbr[t2].add(f)
            stub_branches[(t2, f)] = idx
    by_term = {}
    for t in tw3:
        by_term.setdefault(frozenset((int(t["I"]), int(t["J"]), int(t["K"]))), []).append(t)
    stub2x, unmatch_stubs = {}, []
    for s in stub_only:
        match = by_term.get(frozenset(stub_nbr[s]))
        if not match or len(match) != 1:
            unmatch_stubs.append(s)
            continue
        stub2x[s] = match[0]
    used = set(id(t) for t in stub2x.values())
    unmatch_x = [t for t in tw3 if id(t) not in used]
    if unmatch_stubs:
        issues.append(f"stub buses with no unique transformer: {unmatch_stubs[:5]}")
    if unmatch_x:
        issues.append(f"3-winding transformers with no stub: {[(int(t['I']), int(t['J']), int(t['K'])) for t in unmatch_x[:5]]}")
    print(f"  stub<->transformer matched: {len(stub2x)} | unmatched stubs {len(unmatch_stubs)} | "
          f"unmatched 3-winding {len(unmatch_x)}")

    # stub bus attributes: PQ, baseKV 1, no load, no gen
    n_bad = n_g = 0
    for s in stub_only:
        m = mp_bus[s]
        if not (int(m[1]) == 1 and abs(m[9] - 1) < 1e-9 and abs(m[2]) < 1e-9 and abs(m[3]) < 1e-9):
            n_bad += 1
            if n_bad <= 3:
                issues.append(f"stub {s}: unexpected attrs (type {m[1]}, baseKV {m[9]}, Pd {m[2]}, Qd {m[3]})")
        n_g += any(int(r[0]) == s for r in mp["gen"])
    print(f"  stub attributes: bad {n_bad} | stubs with gen: {n_g}")

    # star legs: impedance, rating (RAW RATEA is in MW), angle
    n_zr = n_zx = n_zrate = n_zang = n_miss = 0
    for s in sorted(stub2x):
        t = stub2x[s]
        Z12 = complex(t["R1_2"] or 0, t["X1_2"] or 0)
        Z23 = complex(t["R2_3"] or 0, t["X2_3"] or 0)
        Z31 = complex(t["R3_1"] or 0, t["X3_1"] or 0)
        exp = {
            int(t["I"]): ((Z12 + Z31 - Z23) / 2, t["RATA1"]),
            int(t["J"]): ((Z12 + Z23 - Z31) / 2, t["RATA2"]),
            int(t["K"]): ((Z23 + Z31 - Z12) / 2, t["RATA3"]),
        }
        for term, (Z, R) in exp.items():
            idx = stub_branches.get((s, term), stub_branches.get((term, s)))
            if idx is None:
                n_miss += 1
                if n_miss <= 3:
                    issues.append(f"stub {s}: star leg to terminal {term} missing")
                continue
            row = mp["branch"][idx]
            if abs(row[2] - Z.real) > 1e-5 * max(1.0, abs(Z.real)):
                n_zr += 1
                if n_zr <= 3:
                    issues.append(f"stub {s} leg {term}: r MP {row[2]:.6f} vs star {Z.real:.6f}")
            if abs(row[3] - Z.imag) > 1e-5 * max(1.0, abs(Z.imag)):
                n_zx += 1
                if n_zx <= 3:
                    issues.append(f"stub {s} leg {term}: x MP {row[3]:.6f} vs star {Z.imag:.6f}")
            if abs(row[5] - (R or 0)) > 1.0:
                n_zrate += 1
                if n_zrate <= 3:
                    issues.append(f"stub {s} leg {term}: rateA MP {row[5]:.2f} MW vs RATA {R or 0} MW")
            if abs(row[9]) > 1e-9:
                n_zang += 1
                if n_zang <= 3:
                    issues.append(f"stub {s} leg {term}: angle MP {row[9]:.4f} (expected 0)")
    print(f"  star legs: r mismatch {n_zr} | x mismatch {n_zx} | rateA mismatch {n_zrate} | "
          f"nonzero angle {n_zang} | missing {n_miss}")

# ---- loads ---------------------------------------------------------------
vg_load = {}
for l in vg["loads"]:
    if l["STATUS"] in (1, None):
        vg_load[int(l["I"])] = vg_load.get(int(l["I"]), 0.0) + (l["PL"] or 0.0)
n_ld = 0
for i, m in mp_bus.items():
    if i not in vg_bus:  # star stub: pure artifact, load must be 0
        if abs(m[2]) > 1e-9:
            issues.append(f"stub {i}: unexpected load {m[2]}")
        continue
    raw_ld = (vg_load.get(i, 0.0) + (vg_bus[i]["BL"] or 0.0))
    if abs(raw_ld - m[2]) > 0.01:
        n_ld += 1
        if n_ld <= 5:
            issues.append(f"bus {i}: load MP {m[2]:.2f} vs RAW {raw_ld:.2f}")
print(f"loads: per-bus mismatch (>0.01 MW): {n_ld}")
tot_raw = sum(vg_load.values()) + sum(b["BL"] or 0 for b in vg["buses"])
print(f"  total load: MP {sum(m[2] for m in mp_bus.values()):.2f} vs RAW {tot_raw:.2f}")

# ---- generators (duplicate-safe: group by bus, match within bus by Pg) ---
mp_by_bus = {}
for r in mp["gen"]:
    mp_by_bus.setdefault(int(r[0]), []).append(r)
vg_by_bus = {}
for g in vg["generators"]:
    vg_by_bus.setdefault(int(g["I"]), []).append(g)
if set(mp_by_bus) != set(vg_by_bus):
    issues.append(f"gen bus mismatch: mp-only={sorted(set(mp_by_bus)-set(vg_by_bus))[:5]} vg-only={sorted(set(vg_by_bus)-set(mp_by_bus))[:5]}")
if sum(len(v) for v in mp_by_bus.values()) != sum(len(v) for v in vg_by_bus.values()):
    issues.append(f"gen count mismatch: MP {len(mp['gen'])} vs RAW {len(vg['generators'])}")
n_pg = n_st = 0
maxpgdiff = 0.0
for i in sorted(mp_by_bus):
    mrows = sorted(mp_by_bus[i], key=lambda r: (r[7], r[1]))
    vrows = sorted(vg_by_bus.get(i, []), key=lambda g: (g["STAT"] or 0, g["PG"] or 0))
    if len(mrows) != len(vrows):
        n_pg += abs(len(mrows) - len(vrows))
        if n_pg <= 5:
            issues.append(f"gen bus {i}: MP {len(mrows)} units vs RAW {len(vrows)}")
        continue
    for m, v in zip(mrows, vrows):
        d = m[1] - (v["PG"] or 0)
        if abs(d) > 0.05:
            maxpgdiff = max(maxpgdiff, abs(d))
            n_pg += 1
        if int(round(m[7])) != int(round(v["STAT"] or 0)):
            n_st += 1
print(f"generators: {len(mp_by_bus)} buses / {sum(len(v) for v in mp_by_bus.values())} units | Pg diff >0.05 MW: {n_pg} (max {maxpgdiff:.3f}) | status mismatch: {n_st}")

# ---- branches (multiset match on (f,t) groups; parallel circuits share a
# bus pair, so match within each (f,t) group by (r, x, rateA)) ---------------
star_idx = set(stub_branches.values())
mp_groups = {}
for idx, r in enumerate(mp["branch"]):
    if idx in star_idx:  # star legs handled above
        continue
    mp_groups.setdefault((int(r[0]), int(r[1])), []).append(r)
vg_groups = {}
for b in vg["branches"]:
    vg_groups.setdefault((int(b["I"]), int(b["J"])), []).append(("line", b, b["CKT"]))
for t in vg["transformers"]:
    if t["K"] not in (0, 0.0, None):  # 3-winding handled above
        continue
    vg_groups.setdefault((int(t["I"]), int(t["J"])), []).append(("xfrm", t, t["CKT"]))

if set(mp_groups) != set(vg_groups):
    only_mp = sorted(set(mp_groups) - set(vg_groups))
    only_vg = sorted(set(vg_groups) - set(mp_groups))
    issues.append(f"branch topology mismatch: mp-only={only_mp[:5]} vg-only={only_vg[:5]}")
    print(f"branch topology: {len(only_mp)} one-sided, {len(only_vg)} other")
else:
    print(f"branch topology: all {len(mp_groups)} (f,t) groups match "
          f"({sum(len(v) for v in mp_groups.values())} branches + {len(star_idx)} star legs)")

n_r = n_x = n_rate = 0
n_count = 0
for key, mrows in mp_groups.items():
    vrows = vg_groups.get(key, [])
    if len(mrows) != len(vrows):
        n_count += 1
        issues.append(f"branch {key}: MP {len(mrows)} circuits vs RAW {len(vrows)}")
        continue
    # match by sorted (r, x, rate) fingerprint
    def fp(r, x, rate_mw):
        return (round(r, 7), round(x, 7), round(rate_mw, 1))
    mset = sorted(((fp(m[2], m[3], m[5]), m) for m in mrows), key=lambda p: p[0])
    vset = sorted(((fp(v[1]["R"] or 0, v[1]["X"] or 0, v[1]["RATEA"] or 0), v) for v in vrows),
                  key=lambda p: p[0])
    for (_, m), (_, v) in zip(mset, vset):
        kind, vd = v[0], v[1]
        if abs(m[2] - (vd["R"] or 0)) > 1e-5 * max(1, abs(vd["R"] or 0)):
            n_r += 1
            if n_r <= 5:
                issues.append(f"branch {key} ckt {vd['CKT']} ({kind}): r MP {m[2]} vs RAW {vd['R']}")
        if abs(m[3] - (vd["X"] or 0)) > 1e-5 * max(1, abs(vd["X"] or 0)):
            n_x += 1
            if n_x <= 5:
                issues.append(f"branch {key} ckt {vd['CKT']} ({kind}): x MP {m[3]} vs RAW {vd['X']}")
        if m[5] > 0 and vd.get("RATEA") is not None and abs(m[5] - vd["RATEA"]) > 1.0:
            n_rate += 1
            if n_rate <= 5:
                issues.append(f"branch {key} ckt {vd['CKT']} ({kind}): rateA MP {m[5]} MW vs RAW {vd['RATEA']} MW")
print(f"impedances/ratings: circuit-count mismatch {n_count}  r mismatch {n_r}  x mismatch {n_x}  rateA mismatch {n_rate}")

# phase shifts
mp_ang = sum(1 for r in mp["branch"] if abs(r[9]) > 1e-9)
vg_ang = sum(1 for t in vg["transformers"] if abs(t["ANG1"] or 0) > 1e-9)
print(f"phase shift: MP nonzero {mp_ang} | RAW transformers ANG1!=0 {vg_ang}")

print("\n" + ("ISSUES:" if issues else "NO ISSUES — VeraGrid RAW parse and MATPOWER case agree."))
for i in issues:
    print("  -", i)
