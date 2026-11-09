"""Parse a TAMU ACTIVSg .RAW through VeraGridEngine's RAW parser and dump
the network data to JSON.

Usage:  python vg_extract.py [--dir DIR] [--stem STEM]
  --dir   data directory (default: activsg500)
  --stem  RAW filename without extension (default: ACTIVSg500)

v33 object attribute map (verified by introspection):
  branch:      I, J, CKT, R, X, B, RATEA, ST, MET
  transformer: I, J, K, CKT, R1_2, R2_3, X1_2, X2_3, B, MAG1/2, ANG1, ANSTAR,
               STAT, RATA1, NOMV1/2
  generator:   I, ID, MBASE, PB, PT, PG, QG, STAT, VS
  load:        I, ID, PL, QL, STATUS
  bus:         I, NAME, BASKV, AREA, ZONE, OWNER, VM, VA, BL, GL, EVHI, EVLO
"""
import argparse
import json
from pathlib import Path

from VeraGridEngine.IO.raw.raw_parser_writer import read_raw

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="activsg500")
ap.add_argument("--stem", default="ACTIVSg500")
args = ap.parse_args()

HERE = Path(__file__).resolve().parent
RAW = HERE / args.dir / (args.stem + ".RAW")
OUT = HERE / args.dir / "vg_data.json"

log = read_raw(str(RAW))
c = log.result if hasattr(log, "result") else log
print("parsed:", type(c).__name__, "version", c.version, "SBASE", c.SBASE)

buses = [{k: getattr(b, k) for k in
          ("I", "NAME", "BASKV", "AREA", "ZONE", "OWNER", "VM", "VA", "BL", "GL", "EVHI", "EVLO")}
         for b in c.buses]

branches = [{k: getattr(b, k) for k in
             ("I", "J", "CKT", "R", "X", "B", "RATEA", "RATEB", "RATEC", "ST", "MET")}
            for b in c.branches]

transformers = []
for t in c.transformers:
    transformers.append({
        "I": t.I, "J": t.J, "K": t.K, "CKT": t.CKT,
        "R": t.R1_2 + t.R2_3,          # 2-winding equivalent series R
        "X": t.X1_2 + t.X2_3,          # 2-winding equivalent series X
        "R1_2": t.R1_2 or 0.0, "X1_2": t.X1_2 or 0.0,
        "R2_3": t.R2_3 or 0.0, "X2_3": t.X2_3 or 0.0,
        "R3_1": t.R3_1 or 0.0, "X3_1": t.X3_1 or 0.0,
        "ANG1": t.ANG1 or 0.0, "ANG2": t.ANG2 or 0.0, "ANG3": t.ANG3 or 0.0,
        "ANSTAR": t.ANSTAR or 0.0,
        "STAT": t.STAT, "RATEA": t.RATA1 or 0.0,
        "RATA1": t.RATA1 or 0.0, "RATA2": t.RATA2 or 0.0, "RATA3": t.RATA3 or 0.0,
        "NOMV1": t.NOMV1, "NOMV2": t.NOMV2, "NOMV3": t.NOMV3,
    })

generators = [{k: getattr(g, k) for k in
               ("I", "ID", "MBASE", "PB", "PT", "PG", "QG", "STAT", "VS")}
              for g in c.generators]

loads = [{k: getattr(l, k) for k in ("I", "ID", "PL", "QL", "STATUS", "LOADTYPE", "SCALE")}
         for l in c.loads]

data = {"version": c.version, "SBASE": c.SBASE,
        "buses": buses, "branches": branches, "transformers": transformers,
        "generators": generators, "loads": loads}
OUT.write_text(json.dumps(data))
print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")

print(f"\ncounts: buses {len(buses)}  branches {len(branches)}  transformers {len(transformers)}  "
      f"generators {len(generators)}  loads {len(loads)}")
n3w = sum(1 for t in transformers if t["K"] not in (0, 0.0, None))
print(f"transformers: 2-winding {len(transformers)-n3w}  3-winding {n3w}  | STAT!=1: {sum(1 for t in transformers if t['STAT'] != 1)}")
bad_x = [t for t in transformers if abs(t["X"]) < 1e-12]
print(f"transformers with X==0: {len(bad_x)}  | STAT!=1: {sum(1 for t in transformers if t['STAT'] != 1)}")
print(f"branch ST values: {sorted(set(b['ST'] for b in branches))}")
print(f"total PL: {sum(l['PL'] or 0 for l in loads):.2f} MW | total BL: {sum(b['BL'] or 0 for b in buses):.2f} MW")
print(f"total PG: {sum(g['PG'] or 0 for g in generators):.2f} MW | gen STAT values: {sorted(set(g['STAT'] for g in generators))}")
