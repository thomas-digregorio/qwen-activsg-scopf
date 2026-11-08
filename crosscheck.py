"""Cross-check the VeraGrid-parsed RAW (vg_data.json) against the official
TAMU MATPOWER case (mp_data.json): topology, impedances, ratings, loads,
generators. Any mismatch is printed."""
import json
from pathlib import Path

A = Path(__file__).resolve().parent / "activsg500"
vg = json.loads((A / "vg_data.json").read_text())
mp = json.loads((A / "mp_data.json").read_text())
S = mp["baseMVA"]

issues = []

# ---- buses ---------------------------------------------------------------
mp_bus = {int(r[0]): r for r in mp["bus"]}
vg_bus = {int(b["I"]): b for b in vg["buses"]}
if set(mp_bus) != set(vg_bus):
    issues.append(f"bus id mismatch: mp-only={sorted(set(mp_bus)-set(vg_bus))[:5]} vg-only={sorted(set(vg_bus)-set(mp_bus))[:5]}")
n_vm = n_va = 0
for i in sorted(mp_bus):
    m, v = mp_bus[i], vg_bus[i]
    if abs(m[9] - v["BASKV"]) > 1e-6:
        issues.append(f"bus {i}: baseKV {m[9]} vs RAW BASKV {v['BASKV']}")
    if abs(m[7] - v["VM"]) > 1e-4:
        n_vm += 1
    if abs(m[8] - v["VA"]) > 1e-3:
        n_va += 1
print(f"buses: {len(mp_bus)} match | Vm max diff count>{'n/a' if n_vm==0 else n_vm}: {n_vm}  Va: {n_va}")
if n_vm or n_va:
    worst = max((i for i in mp_bus if abs(mp_bus[i][7] - vg_bus[i]['VM']) > 1e-4),
                key=lambda i: abs(mp_bus[i][7] - vg_bus[i]['VM']), default=None)
    print(f"  worst Vm bus {worst}: mp {mp_bus[worst][7]} vs raw {vg_bus[worst]['VM']}")

# ---- loads ---------------------------------------------------------------
vg_load = {}
for l in vg["loads"]:
    if l["STATUS"] in (1, None):
        vg_load[int(l["I"])] = vg_load.get(int(l["I"]), 0.0) + (l["PL"] or 0.0)
n_ld = 0
for i, m in mp_bus.items():
    raw_ld = (vg_load.get(i, 0.0) + (vg_bus[i]["BL"] or 0.0))
    if abs(raw_ld - m[2]) > 0.01:
        n_ld += 1
        if n_ld <= 5:
            issues.append(f"bus {i}: load MP {m[2]:.2f} vs RAW {raw_ld:.2f}")
print(f"loads: per-bus mismatch (>0.01 MW): {n_ld}")
tot_raw = sum(vg_load.values()) + sum(b["BL"] or 0 for b in vg["buses"])
print(f"  total load: MP {sum(m[2] for m in mp_bus.values()):.2f} vs RAW {tot_raw:.2f}")

# ---- generators ----------------------------------------------------------
mp_gen = {int(r[0]): r for r in mp["gen"]}
vg_gen = {int(g["I"]): g for g in vg["generators"]}
if set(mp_gen) != set(vg_gen):
    issues.append(f"gen bus mismatch: mp-only={sorted(set(mp_gen)-set(vg_gen))} vg-only={sorted(set(vg_gen)-set(mp_gen))}")
diffs = []
for i in sorted(mp_gen):
    d = mp_gen[i][1] - (vg_gen[i]["PG"] or 0)
    if abs(d) > 0.05:
        diffs.append((i, d, int(round(mp_gen[i][7]))))
print(f"generators: {len(mp_gen)} buses match | Pg diff >0.05 MW: {len(diffs)} "
      f"(max {max((abs(d) for _, d, _ in diffs), default=0):.3f})")
if diffs:
    off = sum(1 for _, _, st in diffs if st == 0)
    print(f"  of those, {off} are status=0 (offline units keep a stale PG in the RAW; "
          f"MP zeros them — expected)")
    print(f"  online (status=1) Pg diffs: {[ (i, round(d,2)) for i, d, st in diffs if st == 1 ][:10]}")
n_st = 0
for i in sorted(mp_gen):
    if int(round(mp_gen[i][7])) != int(round(vg_gen[i]["STAT"])):
        n_st += 1
        if n_st <= 5:
            issues.append(f"gen bus {i}: status MP {mp_gen[i][7]} vs RAW {vg_gen[i]['STAT']}")
print(f"  status mismatch: {n_st}")

# ---- branches (multiset match on (f,t) groups; parallel circuits share a
# bus pair, so match within each (f,t) group by (r, x, rateA)) ---------------
mp_groups = {}
for r in mp["branch"]:
    mp_groups.setdefault((int(r[0]), int(r[1])), []).append(r)
vg_groups = {}
for b in vg["branches"]:
    vg_groups.setdefault((int(b["I"]), int(b["J"])), []).append(("line", b, b["CKT"]))
for t in vg["transformers"]:
    vg_groups.setdefault((int(t["I"]), int(t["J"])), []).append(("xfrm", t, t["CKT"]))

if set(mp_groups) != set(vg_groups):
    only_mp = sorted(set(mp_groups) - set(vg_groups))
    only_vg = sorted(set(vg_groups) - set(mp_groups))
    issues.append(f"branch topology mismatch: mp-only={only_mp[:5]} vg-only={only_vg[:5]}")
    print(f"branch topology: {len(only_mp)} one-sided, {len(only_vg)} other")
else:
    print(f"branch topology: all {len(mp_groups)} (f,t) groups match "
          f"({sum(len(v) for v in mp_groups.values())} total branches)")

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
    vset = sorted(((fp(v[1]["R"] or 0, v[1]["X"] or 0, (v[1]["RATEA"] or 0) * S), v) for v in vrows),
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
        if m[5] > 0 and vd.get("RATEA") is not None and abs(m[5] - vd["RATEA"] * S) > 1.0:
            n_rate += 1
            if n_rate <= 5:
                issues.append(f"branch {key} ckt {vd['CKT']} ({kind}): rateA MP {m[5]} MW vs RAW {vd['RATEA']} pu*{S}")
print(f"impedances/ratings: circuit-count mismatch {n_count}  r mismatch {n_r}  x mismatch {n_x}  rateA mismatch {n_rate}")

# phase shifts
mp_ang = sum(1 for r in mp["branch"] if abs(r[9]) > 1e-9)
vg_ang = sum(1 for t in vg["transformers"] if abs(t["ANG1"] or 0) > 1e-9)
print(f"phase shift: MP nonzero {mp_ang} | RAW transformers ANG1!=0 {vg_ang}")

print("\n" + ("ISSUES:" if issues else "NO ISSUES — VeraGrid RAW parse and MATPOWER case agree."))
for i in issues:
    print("  -", i)
