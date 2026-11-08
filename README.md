# OPF-scratch — DC N-1 OPF: IEEE 14-bus and ACTIVSg500

Scratch workspace for power-system OPF work. All data public-domain,
all computation local (Python 3.13 + numpy/scipy, HiGHS 1.15.1 via
`highspy`; the 14-bus run also works through `scipy.optimize.linprog`).

## Files

| file | purpose |
|---|---|
| `case14.py` | IEEE 14-bus network data (buses, branches, generators, costs) |
| `dc_n1_opf.py` | DC N-1 OPF solver + report generator (14-bus) |
| `test_contingency_k.py` | regression test: contingency K-matrices vs independent DC power flow |
| `test_highspy_lmps.py` | LMP cross-check (duals vs finite differences) on the 14-bus case |
| `results_report_limit100MW.md` | latest 14-bus run report (uniform 100 MW line limit) |
| `vg_extract.py` | ACTIVSg500: RAW → `activsg500/vg_data.json` via VeraGridEngine parser |
| `mp_extract.py` | ACTIVSg500: official MATPOWER case → `activsg500/mp_data.json` |
| `crosscheck.py` | ACTIVSg500: cross-validates the two parses (impedances, ratings, costs) |
| `activsg500_n1_opf.py` | ACTIVSg500: DC N-1 SCOPF driver (HiGHS, constraint generation) |
| `test_contingency_k_500.py` | ACTIVSg500: K-matrix regression vs independent DC power flow |
| `verify_basecase_cost.py` | ACTIVSg500: cross-check vs official case file + independent N-0 re-solve |
| `activsg500/` | dataset + generated JSON + `results_activsg500.md` report |
| `activsg500_dataset.zip` | TAMU ACTIVSg500 dataset (RAW, MATPOWER case, contab) |
| `download_activsg500.ps1` | one-shot download script for the dataset |

Run:

```
python dc_n1_opf.py            # default: 100 MW uniform line limit
python dc_n1_opf.py --limit 70 # tighter limits (N-1 binds harder)
python test_contingency_k.py   # verification
```

## Data provenance

`case14.py` was transcribed on 2026-11-07 from MATPOWER `case14`
(retrieved from `raw.githubusercontent.com/MATPOWER/matpower/master/data/case14.m`),
itself converted from the IEEE Common Data Format (UW archive).
14 buses, 20 lines, 5 generators, total load 259.0 MW, 100 MVA base.

## Formulation

Lossless DC, classical N-1 security constraint (no corrective re-dispatch):

- variables: `P_g` per generator bus, bus angles (slack bus 1 fixed at 0),
  piecewise-linear (PL) cost-segment variables
- objective: PL approximation (12 segments) of the MATPOWER quadratic
  generator costs `c2*P^2 + c1*P`
- constraints:
  - nodal power balance at **every** bus (incl. slack — the slack row pins
    the slack generation to the network residual; omitting it leaves
    `P_g_slack` free and the objective drives it to its bound)
  - generator limits, base-case line flow limits `|F_l| <= Fmax`
  - for each of the 25 N-1 contingencies (20 line outages + 5 generator
    outages): post-contingency flows `|F_l^(c)| <= Fmax`, with
    `F_l^(c) = K_c (P_g - loads)` linear in the pre-contingency dispatch;
    disconnected contingencies are skipped
- solver: HiGHS (sparse). LMPs by central finite difference of the optimal
  cost w.r.t. each bus load (±0.05 MW).

## Assumptions / caveats

- **Line limits:** the raw case14 data defines no thermal ratings
  (all `rateA = 0` in MATPOWER). A **uniform 100 MW** (1.0 pu) limit is
  assumed — change with `--limit`.
- Lossless DC: no losses, voltages at 1.0 pu, off-nominal transformer taps
  (lines 4-7, 4-9, 5-6) ignored.
- No corrective re-dispatch under contingencies (pre-contingency schedule
  must be post-contingency feasible — the standard SCOPF N-1 constraint).
- PL cost is the 12-segment chord approximation of the quadratic cost.
  For a convex cost curve the chord lies ABOVE the curve, so it is an
  upper bound: the LP objective is >= the exact quadratic cost at the
  same dispatch (max gap ~0.2% on this case with 12 segments).

## Verification

- nodal balance asserted on the solution (gen total = load total)
- `test_contingency_k.py` checks the K-matrix of **every** contingency
  against an independently built post-contingency DC power flow
  (24 connected contingencies verified to 1e-7; the 7-8 outage correctly
  identified as disconnecting bus 8)
- KKT sanity: at the N-1 optimum, gen1 marginal ($28.6/MWh at 100 MW)
  equals the bus-1 LMP; N-1 constraints bind at 100% utilization

## Results (uniform 100 MW limit, 2026-11-07)

- N-0: $7,931.9/h — gen1 154.8 (export capped by base line limits),
  gen2 46.7, gen3 50.0, gen8 7.6; violates **14 of 25** contingencies
  (e.g. line 1-2 outage forces 154.8 MW through line 1-5).
- N-1: $8,438.4/h (+$506.4/h = price of N-1 security) — gen1 capped at
  exactly 100 MW (each of lines 1-2 / 1-5 must alone carry bus-1 export
  under the other's outage), gen6 36.4 and gen8 41.7 brought in;
  all 25 contingencies satisfied, worst post-contingency flow 100%.
- Binding contingencies: line 1-2 out, line 1-5 out, gen bus 2 out.
- LMPs: bus 1 $28.61/MWh, buses 2-14 ≈ $40.7-41.1/MWh.

## ACTIVSg500 (TAMU 500-bus test case)

DC N-1 SCOPF on the public synthetic ACTIVSg500 (electricgrids.engr.tamu.edu;
free for commercial or non-commercial use — cite A. Birchfield, T. Xu,
K. Gegner, K. Shetye, T. Overbye, "Grid Structural Characteristics as
Validation Criteria for Synthetic Networks," IEEE Trans. Power Syst. 32(4),
2017, and T. Xu et al., "Application of Large-Scale Synthetic Power System
Models for Energy Economic Studies," HICSS 2017, for the grid + cost data).
Network from `ACTIVSg500.RAW` via the **VeraGridEngine** RAW parser; costs and
generator limits from the official MATPOWER case `case_ACTIVSg500.m` (the RAW
carries no P records); contingency set from TAMU's `contab_ACTIVSg500.m`
(591 line outages + 90 generator outages, of which 34 are offline units and
therefore redundant → 647 effective). Solved with **HiGHS** via `highspy`.

### Data pipeline

```
python vg_extract.py    # RAW -> activsg500/vg_data.json (needs VeraGridEngine)
python mp_extract.py    # case_ACTIVSg500.m -> activsg500/mp_data.json
python crosscheck.py    # 0 impedance/rating mismatches
```

500 buses, 466 lines, 131 transformers, 90 generators (56 online, 34 offline),
206 loads, total load **7750.72 MW**, 100 MVA base, 36 zero-cost must-run units.

### Run

```
python activsg500_n1_opf.py             # full run, ~12 min
python activsg500_n1_opf.py --gentest 0 # skip the ~8-min gen diagnostic
python test_contingency_k_500.py        # K regression (passes)
python verify_basecase_cost.py          # cross-checks vs official case file
```

Formulation is identical to the 14-bus case (lossless DC, classical N-1,
no corrective re-dispatch, 20-segment PL chord cost). The full N-1 LP
(705,666 flow rows + 556 equality rows × 1675 vars, 38.8M nnz) is built
directly in CSR; solved by constraint generation (6 iterations / 6 solves)
and verified against the full LP (same objective, 0 violations).

### Headline results (2026-11-08, `activsg500/results_activsg500.md`)

- **N-0: $70,793.0/h** — dispatches 318 of 591 line outages into violation
  (worst post-contingency ratio 1.79).
- **N-1 (all 591 line outages): $80,639.3/h** (+$9,846.3/h = price of
  N-1 security); worst post-contingency ratio exactly 1.0000; binding
  contingencies: line 423-496 out, the two 9-7 transformers, the two 17-14
  transformers, the two 144-143 transformers.
- Dispatch: security curtails the two largest units — bus 9
  771.8→680.0 MW, bus 17 888.9→674.0 MW — all 36 zero-cost units
  sit at Pmax, and gen 144/145/456 settle interior.
- **LMPs: uniform $27.280/MWh** at all 12 sampled buses (duals of the
  balance rows — pu-scaled, divided by S_BASE — agree with ±0.05 MW central
  finite differences to 4 digits).
- **Generator outages (diagnostic):** the joint line+gen contingency set is
  **infeasible** (proven by HiGHS). Per-outage testing isolates the cause to
  exactly 2 of 56 units — bus 9 (Pmax 771.8 MW) and bus 17 (Pmax 888.9 MW),
  the two largest generators in the case: tripping either forces line 44-1
  (limit 320.3 MW) to carry ~443-453 MW. Excluding those two, the
  645-contingency SCOPF is feasible at **$81,955.3/h**.
- 253 of 647 contingencies disconnect a radial island (assumed shed):
  max island load 161.3 MW, max island gen Pmax 888.9 MW, 36 islands
  carry generation.

### How the results are verified

There is no published DC N-1 SCOPF cost benchmark for ACTIVSg500 to match
(the TAMU references — Birchfield et al., Xu et al. — define the synthetic
grid and its cost data, not a reference SCOPF solution). Instead the results
rest on independent checks:

- **Data correctness** — `crosscheck.py`: the VeraGrid RAW parse and the
  official MATPOWER case agree on every branch impedance, rating, and cost
  (0 mismatches). `verify_basecase_cost.py` [1] re-reads bus/status/Pmin/Pmax/
  cost straight from the raw case file and matches the model exactly
  (sum Pmax 8863.65, sum Pmin 2659.06).
- **Published base-case anchor** — the case file carries its own reference
  operating point (base-case `Pg`). Its exact-quadratic cost is
  **$72,792.8/h**; our cost-optimal N-0 is **$70,793.0/h**, i.e. 2.75%
  *cheaper*. A cost optimization must be ≤ any reference dispatch, and a
  ~2.7% gap is exactly what you expect from a power-flow snapshot that was
  never cost-optimized. (If N-0 came out at, say, $10k/h, the cost mapping
  or balance would be badly wrong.)
- **Independent re-solve** — `verify_basecase_cost.py` [3] rebuilds and solves
  the N-0 LP from scratch (separate HiGHS call): it reproduces $70,793.0/h
  to the dollar, all 36 zero-cost units sit at Pmax (free power is always
  taken), and the nodal-balance residual is 4e-12 MW.
- **Contingency math** — `test_contingency_k_500.py` checks the K-matrix of
  35 contingencies against an independently built post-contingency DC power
  flow (agrees to 1e-7).
- **CG vs full LP** — the constraint-generation solution equals the full
  705,666-row N-1 LP to the dollar (same objective), with 0 violations across
  all 591 line outages.
- **LMPs** — HiGHS equality duals (÷ S_BASE) agree with ±0.05 MW central
  finite differences to 4 digits at all 12 sampled buses.
- **Infeasibility proof** — the joint line+gen set is proven infeasible by
  HiGHS and confirmed by a second, independent HiGHS call; per-outage testing
  isolates the cause to exactly the two largest generators (bus 9, bus 17).
