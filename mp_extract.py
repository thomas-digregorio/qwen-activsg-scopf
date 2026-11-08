"""Extract all arrays from case_ACTIVSg500.m (TAMU's official MATPOWER
representation) into JSON + print OPF-relevant summaries."""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
MP = HERE / "activsg500" / "case_ACTIVSg500.m"
CT = HERE / "activsg500" / "contab_ACTIVSg500.m"
OUT = HERE / "activsg500" / "mp_data.json"

text = MP.read_text(errors="ignore")


def matp_array(name, src):
    m = re.search(name + r"\s*=\s*\[(.*?)\];", src, re.S)
    if not m:
        return None
    rows = []
    for ln in m.group(1).splitlines():
        ln = ln.strip().rstrip(";")
        if not ln or ln.startswith("%"):
            continue
        rows.append([float(x) for x in ln.split()])
    return rows


bus = matp_array("mpc.bus", text)
gen = matp_array("mpc.gen", text)
branch = matp_array("mpc.branch", text)
gencost = matp_array("mpc.gencost", text)
baseMVA = float(re.search(r"mpc\.baseMVA\s*=\s*([\d.]+)", text).group(1))

ct_text = CT.read_text(errors="ignore")
contab = matp_array("mpc.contab", ct_text)

data = {
    "baseMVA": baseMVA,
    "bus": bus,        # [i, type, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, zone, Vmax, Vmin, ...]
    "gen": gen,        # [bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2, ...]
    "branch": branch,  # [fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, ...]
    "gencost": gencost,  # [2, startup, shutdown, n, c2, c1, c0]
    "contab": contab,
}
OUT.write_text(json.dumps(data))
print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")

# ---- summaries -------------------------------------------------------------
total_load = sum(r[2] for r in bus)
total_pg = sum(r[1] for r in gen)
total_pmax = sum(r[8] for r in gen)
print(f"\nbus: {len(bus)}  gen: {len(gen)}  branch: {len(branch)}  gencost: {len(gencost)}")
print(f"total load Pd: {total_load:.2f} MW | total initial Pg: {total_pg:.2f} | total Pmax: {total_pmax:.2f}")

# status / online units
on = [r for r in gen if r[7] == 1]
print(f"gen status=1: {len(on)}  (status values: {sorted(set(r[7] for r in gen))})")

# units with zero gencost: what are their Pmax/Pmin/status?
print("\nzero-gencost units (bus, status, Pmin, Pmax):")
nzc = 0
for i, gc in enumerate(gencost):
    coeffs = gc[4:4 + int(gc[3])]
    if all(abs(c) < 1e-12 for c in coeffs):
        g = gen[i]
        if not (g[7] == 0 and g[8] == 0.0 and g[9] == 0.0):
            print(f"  unit {i+1}: bus {int(g[0])}, status {g[7]}, Pmin {g[9]:.1f}, Pmax {g[8]:.1f}")
            nzc += 1
print(f"  -> {nzc} zero-cost units that are NOT (status=0, Pmax=0, Pmin=0)")

# branch angle / ratio / status
angles = [b[9] for b in branch]
ratios = [b[8] for b in branch]
status = [b[10] for b in branch]
nz_ang = [(int(b[0]), int(b[1]), a) for b, a in zip(branch, angles) if abs(a) > 1e-9]
nz_rat = [(int(b[0]), int(b[1]), r) for b, r in zip(branch, ratios) if abs(r - 1.0) > 1e-9]
print(f"\nbranch status values: {sorted(set(s for s in status))}")
print(f"nonzero angle: {len(nz_ang)}  max|angle|={max(abs(a) for a in angles):.4f} deg")
print(f"non-unity ratio: {len(nz_rat)}  max|r-1|={max(abs(r-1) for r in ratios):.4f}")
print(f"rateA: min={min(b[5] for b in branch):.2f} max={max(b[5] for b in branch):.2f} (zero-rated: {sum(1 for b in branch if b[5] == 0)})")

# gencost coefficient spread
import numpy as np
c2 = [gc[4] for gc in gencost]
c1 = [gc[5] for gc in gencost]
c0 = [gc[6] for gc in gencost]
print(f"\ngencost: c2 in [{min(c2):.5f}, {max(c2):.5f}]  c1 in [{min(c1):.2f}, {max(c1):.2f}]  c0 in [{min(c0):.1f}, {max(c0):.1f}]")

# contingency table
if contab is not None:
    nct = len(contab)
    # each row: [br, bus, br, bus, ..., -1]
    n_br = sum(1 for row in contab for v in row if v > 0 and v <= 597 and row[row.index(v, row.index(-1) - 1) if -1 in row else 0] >= 0) if False else None
    # simpler count: split each row at -1
    n_lines = n_gens = 0
    for row in contab:
        vals = row[row.index(-1) - 1:] if -1 in row else row
        # actually: entries before -1 are alternating? MATPOWER contab row = [br bus ...]
        pre = row[:row.index(-1)] if -1 in row else row
        n_lines += sum(1 for v in pre if v > 0 and v <= len(branch))
        n_gens += sum(1 for v in pre if -(len(gen)) <= v < 0)
    print(f"\ncontab: {nct} rows -> counted {n_lines} line outages, {n_gens} gen outages (heuristic)")
    print("  first row sample:", contab[0][:8], "..." if len(contab[0]) > 8 else "")
