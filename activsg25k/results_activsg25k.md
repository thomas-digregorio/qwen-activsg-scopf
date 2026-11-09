# DC N-0 OPF + N-1 Security Screening — ACTIVSg25k (TAMU 25,000-bus)

Network from ACTIVSg25k.RAW via the VeraGridEngine RAW parser (24,250 buses + 750 synthetic star-stub buses; 3-winding transformers expanded to 3-branch stars); costs/limits from the official MATPOWER case (3,779 online units on 2,753 buses); solved with HiGHS via highspy (20-segment PL cost, lossless DC). Total load 234527.57 MW.

**Method note.** At 25k scale a full N-1 SCOPF (all 32229 branch outages as explicit constraints) is intractable. This is instead an exact security *screen* of the N-0 dispatch: non-bridge outages are solved exactly via a Sherman–Morrison rank-1 update of the pre-contingency angle vector (verified in-run against direct re-solves), bridge outages island the network and are analyzed per-island (re-factorized slack-side component + island statistics and internal flow sweep). Gen outages are screened as a labeled diagnostic. The screen answers "which single outages violate line limits at the N-0 dispatch" — not the N-1-optimal dispatch.

### N-0 (base-case) OPF
Cost (exact quadratic): **$5,856,237.3/h** (PL objective 5,856,238.7/h; the 20-segment chord of the convex cost is an upper bound — gap 1.4 $/h)
Snapshot (pre-OPF) dispatch cost: $6,027,935.0/h (snapshot imbalance +5160.93 MW — the case is not power-balanced; the OPF dispatch is)

Top 15 units by output:

| bus | P (MW) | Pmin | Pmax | $/MWh marginal |
|---|---|---|---|---|
| 42018 | 1299.00 | 441.72 | 1299.00 | 8.000 |
| 62256 | 1298.00 | 169.75 | 1298.00 | 8.000 |
| 62255 | 1298.00 | 272.80 | 1298.00 | 8.000 |
| 14645 | 1290.70 | 377.19 | 1290.70 | 8.000 |
| 55528 | 1269.90 | 573.73 | 1269.90 | 8.000 |
| 35457 | 1242.00 | 403.03 | 1242.00 | 8.000 |
| 55587 | 1220.50 | 389.89 | 1220.50 | 8.000 |
| 55588 | 1220.50 | 289.03 | 1220.50 | 8.000 |
| 25692 | 1220.30 | 162.61 | 1220.30 | 8.000 |
| 25691 | 1220.30 | 457.33 | 1220.30 | 8.000 |
| 11434 | 1205.10 | 391.05 | 1205.10 | 8.000 |
| 11433 | 1205.10 | 277.37 | 1205.10 | 8.000 |
| 14472 | 1170.00 | 314.52 | 1170.00 | 8.000 |
| 14471 | 1170.00 | 504.95 | 1170.00 | 8.000 |
| 62080 | 1138.50 | 160.11 | 1138.50 | 8.000 |

Zero-cost (must-run) units: 1483 online, all with Pmin > 0 (sum Pmin 7982.5 MW); at Pmin: 0, at Pmax: 1483, in between: 0. The LP is cost-indifferent over their [Pmin, Pmax] ranges, so the optimal dispatch of these units is degenerate (HiGHS sits them at Pmax); the cost is unaffected either way.

Base-case loading, top 10 rated lines (by % of limit):

| branch | flow (MW) | limit (MW) | % |
|---|---|---|---|
| 63405-1000689 (3w-star) | -131.10 | 141.1 | -92.9 |
| 55584-55529 | -294.48 | 320.0 | -92.0 |
| 55584-55529 | -294.48 | 320.0 | -92.0 |
| 55584-55529 | -294.48 | 320.0 | -92.0 |
| 63115-63422 | -294.15 | 320.0 | -91.9 |
| 63115-63422 | -294.15 | 320.0 | -91.9 |
| 63115-63422 | -294.15 | 320.0 | -91.9 |
| 55288-55584 | -262.83 | 288.7 | -91.0 |
| 55288-55584 | -262.83 | 288.7 | -91.0 |
| 55288-55584 | -262.83 | 288.7 | -91.0 |

Largest flow on an unlimited (zero-rated) branch: 2412.95 MW on 11378-11383 (8899 such branches are unconstrained in the model)

### N-0 LMPs ($/MWh)

(dual = N-0 LP balance-row dual, pu-scaled; FD = central finite difference of the N-0 optimal cost, ±1 MW)

The N-0 base case has no binding line limit (worst ratio 0.9294), so every nodal LMP equals the uniform system marginal cost.

| bus | load (MW) | LMP (dual) | LMP (FD) |
|---|---|---|---|
| 13189 | 99.59 | 30.028 | 30.028 |
| 13262 | 106.24 | 30.028 | 30.028 |
| 31072 | 99.31 | 30.028 | 30.028 |
| 46073 | 0.00 | 30.028 | 30.028 |
| 57337 | 96.28 | 30.028 | 30.028 |
| 62120 | 0.00 | 30.028 | 30.028 |
| 63232 | 94.96 | 30.028 | 30.028 |
| 69393 | 113.36 | 30.028 | 30.028 |

### N-1 line-outage screen — non-bridge outages (21290 of 32229 branches)

Outages with ≥1 post-contingency limit violation: **477** (2.24%)
Worst post-contingency ratio (any outage): **12.3311**

Top 20 outages by worst post-contingency ratio:

| outage | worst line | flow (MW) | limit (MW) | ratio | #lines violated |
|---|---|---|---|---|---|
| 53794-53801 | 53796-53795 (xfmr) | 925.2 | 75.0 | 12.331 | 20 |
| 25884-25897 | 25885-25811 | 2281.8 | 283.3 | 8.054 | 15 |
| 46153-46159 | 46228-46154 | 901.9 | 130.3 | 6.920 | 4 |
| 11378-11383 | 11379-11463 | 1433.0 | 225.2 | 6.362 | 17 |
| 55585-55589 | 55105-55123 | 1031.4 | 185.3 | 5.566 | 18 |
| 16262-16271 | 45069-45058 | 911.7 | 169.0 | 5.395 | 78 |
| 25922-25929 | 40057-25924 | 352.0 | 72.3 | 4.870 | 25 |
| 16288-16299 | 16297-15069 | 835.0 | 175.8 | 4.749 | 38 |
| 16250-16259 | 15358-16251 | 939.0 | 198.4 | 4.732 | 13 |
| 17453-17450 | 17594-17481 | 218.1 | 48.4 | 4.507 | 6 |
| 39841-39845 | 39842-39841 (xfmr) | 2003.2 | 449.1 | 4.461 | 17 |
| 16014-16019 | 15296-16015 | 251.8 | 59.9 | 4.204 | 15 |
| 14468-14475 | 14643-14469 | 2240.6 | 537.9 | 4.166 | 9 |
| 25912-25920 | 33016-25913 | 717.3 | 196.3 | 3.655 | 11 |
| 35340-17357 | 35254-35224 | 352.7 | 100.6 | 3.506 | 44 |
| 25801-25811 | 25802-1000239 (3w-star) | 343.8 | 98.7 | 3.483 | 12 |
| 46064-46070 | 46065-46064 (xfmr) | 1684.8 | 500.5 | 3.366 | 21 |
| 45992-46001 | 45993-1000477 (3w-star) | 981.4 | 295.6 | 3.321 | 15 |
| 43751-43757 | 43752-43751 (xfmr) | 1846.8 | 602.4 | 3.066 | 3 |
| 16227-16232 | 16019-16138 | 734.3 | 239.6 | 3.065 | 4 |

Note: the 'flow (MW)' column is the post-contingency flow magnitude on the worst line (ratio × limit).

### N-1 bridge outages — islanding (10939 branches)

Bridges disconnect an island when they fail. Islands with load that even full local generation cannot serve (shortfall = load − ΣPmax > 0): **357** (max shortfall 77.2 MW)
Islands (or slack-side components) with internal limit violations under the frozen N-0 dispatch: **1315**

Top 15 islanding outages (by shortfall, then internal loading):

| outage | island buses | island load (MW) | frozen gen (MW) | gen Pmax (MW) | shortfall (MW) | viol C1/C2 | worst ratio C1/C2 |
|---|---|---|---|---|---|---|---|
| 37187-37189 | 2 | 77.2 | 0.0 | 0.0 | 77.2 | 6/0 | 3.022/0.000 |
| 57432-1000585 (3w-star) | 2 | 73.0 | 0.0 | 0.0 | 73.0 | 0/0 | 0.929/0.398 |
| 57432-57436 | 3 | 73.0 | 0.0 | 0.0 | 73.0 | 0/0 | 0.929/0.398 |
| 41312-41314 | 2 | 72.9 | 0.0 | 0.0 | 72.9 | 12/0 | 3.199/0.000 |
| 13248-1000021 (3w-star) | 2 | 69.1 | 0.0 | 0.0 | 69.1 | 15/0 | 3.717/0.146 |
| 69289-69291 | 2 | 69.1 | 0.0 | 0.0 | 69.1 | 8/0 | 3.167/0.000 |
| 57341-57342 | 1 | 68.4 | 0.0 | 0.0 | 68.4 | 5/0 | 1.945/0.000 |
| 41080-41081 | 1 | 66.4 | 0.0 | 0.0 | 66.4 | 13/0 | 3.745/0.000 |
| 41083-41088 | 4 | 66.4 | 0.0 | 0.0 | 66.4 | 16/0 | 2.573/0.000 |
| 33097-33098 | 1 | 65.0 | 0.0 | 0.0 | 65.0 | 36/0 | 10.769/0.000 |
| 57357-57361 | 4 | 63.3 | 0.0 | 0.0 | 63.3 | 2/0 | 1.379/0.000 |
| 39105-39110 | 4 | 62.7 | 0.0 | 0.0 | 62.7 | 14/0 | 4.123/0.000 |
| 69050-69055 | 4 | 62.3 | 0.0 | 0.0 | 62.3 | 16/0 | 5.801/0.000 |
| 39764-39767 | 2 | 61.5 | 0.0 | 0.0 | 61.5 | 19/0 | 4.544/0.000 |
| 39170-39171 | 1 | 61.1 | 0.0 | 0.0 | 61.1 | 25/0 | 17.965/0.000 |

### N-1 gen-outage screen (DIAGNOSTIC — labeled)

Screening all 3779 online-unit outages at the N-0 dispatch (6 of them on the slack bus — those re-root the network at the largest other gen bus, which absorbs the lost injection).
Outages with ≥1 violation: **1** (0.03%)

Top 10 by worst post-outage ratio:

| unit bus | P (MW) | worst line | limit (MW) | ratio | #viol | ref OK |
|---|---|---|---|---|---|---|
| 53799 | 925.2 | 53333-53344 | 244.1 | 1.024 | 1 | - |
| 22117 | 67.8 | 21973-22117 | 193.1 | 0.969 | 0 | - |
| 53914 | 655.8 | 59684-60014 | 310.1 | 0.964 | 0 | - |
| 49520 | 12.0 | 49516-49518 | 54.8 | 0.947 | 0 | - |
| 53812 | 272.4 | 59684-60014 | 310.1 | 0.935 | 0 | - |
| 63395 | 14.5 | 63405-1000689 (3w-star) | 141.1 | 0.931 | 0 | - |
| 63395 | 14.5 | 63405-1000689 (3w-star) | 141.1 | 0.931 | 0 | - |
| 53860 | 246.5 | 59684-60014 | 310.1 | 0.931 | 0 | - |
| 53861 | 246.5 | 59684-60014 | 310.1 | 0.931 | 0 | - |
| 53862 | 246.5 | 59684-60014 | 310.1 | 0.931 | 0 | - |

### Data / assumptions
- network: ACTIVSg25k.RAW via VeraGridEngine RAW parser — 25000 buses (24250 RAW + 750 star stubs), 32229 in-service branches (24699 lines + 5280 2-winding xfmr + 2250 star legs of 750 3-winding xfmr); 1 out-of-service line excluded; 503 star legs with negative x kept as-is (exact wye-equivalent artifacts; the reduced Laplacian still factors)
- costs/limits: official MATPOWER case; 3779 online units on 2753 buses (per-bus unit counts agree RAW vs MP exactly)
- contingency set: no TAMU contab exists for the 25k case → N-1 line set = all 32229 in-service branches; gen outages screened as a diagnostic
- slack: MP reference bus 62120 (baseKV 500); 6 units there, all at Pmax in the snapshot (their outages are re-rooted in the gen screen)
- lossless DC, no voltage limits, no corrective re-dispatch; screening (not N-1-optimal re-dispatch)
- PL cost: 20-segment chord (upper bound) of the quadratic gencost; zero-rated branches (8899) unconstrained
