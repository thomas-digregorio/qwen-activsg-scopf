# DC N-1 SCOPF — ACTIVSg500 (TAMU 500-bus test case)

Network from ACTIVSg500.RAW via the VeraGridEngine RAW parser; costs/limits from the official MATPOWER case; contingencies from contab_ACTIVSg500.m; solved with HiGHS via highspy (constraint generation; 20-segment PL cost, lossless DC, no corrective re-dispatch). Total load 7750.72 MW.

### N-0 (base case only)
Cost (exact quadratic): $70,793.0/h
  gen bus 9   :    771.80 MW
  gen bus 16  :    444.45 MW
  gen bus 17  :    888.90 MW
  gen bus 18  :    157.60 MW
  gen bus 49  :      1.30 MW
  gen bus 50  :      2.60 MW
  gen bus 71  :      4.00 MW
  gen bus 72  :      2.00 MW
  gen bus 73  :      1.40 MW
  gen bus 82  :    361.90 MW
  gen bus 127 :    162.00 MW
  gen bus 128 :    162.00 MW
  gen bus 144 :    602.55 MW
  gen bus 145 :    285.43 MW
  gen bus 167 :    193.50 MW
  gen bus 168 :    193.50 MW
  gen bus 169 :    387.00 MW
  gen bus 197 :    293.60 MW
  gen bus 198 :    207.30 MW
  gen bus 222 :      5.20 MW
  gen bus 223 :     18.00 MW
  gen bus 224 :    511.20 MW
  gen bus 225 :    514.80 MW
  gen bus 231 :      2.90 MW
  gen bus 258 :      1.10 MW
  gen bus 301 :      2.90 MW
  gen bus 302 :      2.90 MW
  gen bus 305 :      1.00 MW
  gen bus 306 :      2.40 MW
  gen bus 319 :     42.30 MW
  gen bus 350 :     15.00 MW
  gen bus 351 :      6.20 MW
  gen bus 352 :     14.40 MW
  gen bus 353 :     74.00 MW
  gen bus 410 :     53.56 MW
  gen bus 411 :     55.46 MW
  gen bus 412 :     55.46 MW
  gen bus 413 :     55.46 MW
  gen bus 430 :      1.00 MW
  gen bus 431 :      3.30 MW
  gen bus 432 :      1.20 MW
  gen bus 433 :      6.70 MW
  gen bus 434 :      4.30 MW
  gen bus 437 :    197.00 MW
  gen bus 438 :    197.00 MW
  gen bus 439 :    197.00 MW
  gen bus 455 :     66.85 MW
  gen bus 456 :    133.70 MW
  gen bus 458 :    109.60 MW
  gen bus 480 :     56.00 MW
  gen bus 481 :     45.00 MW
  gen bus 482 :     45.00 MW
  gen bus 483 :     24.00 MW
  gen bus 484 :     28.00 MW
  gen bus 497 :     60.00 MW
  gen bus 498 :     18.00 MW
N-1 violations of this dispatch: 318 of 591 line outages (worst ratio 1.7933)

### N-1 (SCOPF, classical: all 591 line outages enforced)
Cost (exact quadratic): $80,639.3/h (PL obj 80,639.8; delta vs N-0 +9,846.3/h)
  gen bus 9   :    680.00 MW 
  gen bus 16  :    444.45 MW  **
  gen bus 17  :    674.00 MW 
  gen bus 18  :    157.60 MW  **
  gen bus 49  :      1.30 MW  **
  gen bus 50  :      2.60 MW  **
  gen bus 71  :      4.00 MW  **
  gen bus 72  :      2.00 MW  **
  gen bus 73  :      1.40 MW  **
  gen bus 82  :    361.90 MW  **
  gen bus 127 :    162.00 MW  **
  gen bus 128 :    162.00 MW  **
  gen bus 144 :    336.00 MW 
  gen bus 145 :    225.17 MW 
  gen bus 167 :    193.50 MW  **
  gen bus 168 :    193.50 MW  **
  gen bus 169 :    387.00 MW  **
  gen bus 197 :    293.60 MW  **
  gen bus 198 :    207.30 MW  **
  gen bus 222 :      5.20 MW  **
  gen bus 223 :     18.00 MW  **
  gen bus 224 :    511.20 MW  **
  gen bus 225 :    514.80 MW  **
  gen bus 231 :      2.90 MW  **
  gen bus 258 :      1.10 MW  **
  gen bus 301 :      2.90 MW  **
  gen bus 302 :      2.90 MW  **
  gen bus 305 :      1.00 MW  **
  gen bus 306 :      2.40 MW  **
  gen bus 319 :     42.30 MW  **
  gen bus 350 :     15.00 MW  **
  gen bus 351 :      6.20 MW  **
  gen bus 352 :     14.40 MW  **
  gen bus 353 :     74.00 MW  **
  gen bus 410 :    136.95 MW  **
  gen bus 411 :    136.95 MW  **
  gen bus 412 :    136.95 MW  **
  gen bus 413 :    136.95 MW  **
  gen bus 430 :      1.00 MW  **
  gen bus 431 :      3.30 MW  **
  gen bus 432 :      1.20 MW  **
  gen bus 433 :      6.70 MW  **
  gen bus 434 :      4.30 MW  **
  gen bus 437 :    197.00 MW  **
  gen bus 438 :    197.00 MW  **
  gen bus 439 :    197.00 MW  **
  gen bus 455 :    222.83 MW  **
  gen bus 456 :    283.37 MW 
  gen bus 458 :    109.60 MW  **
  gen bus 480 :     56.00 MW  **
  gen bus 481 :     45.00 MW  **
  gen bus 482 :     45.00 MW  **
  gen bus 483 :     24.00 MW  **
  gen bus 484 :     28.00 MW  **
  gen bus 497 :     60.00 MW  **
  gen bus 498 :     18.00 MW  **
All 591 line-outage contingencies satisfied (violations: 0); worst post-contingency ratio 1.0000
Binding N-1 contingencies (>=99.9% of limit), top 15:
  - line 423-496 out: ratio 1.0000
  - line 9-7 (xfmr 1) out: ratio 1.0000
  - line 9-7 (xfmr 2) out: ratio 1.0000
  - line 17-14 (xfmr 1) out: ratio 1.0000
  - line 17-14 (xfmr 2) out: ratio 1.0000
  - line 144-143 (xfmr 1) out: ratio 1.0000
  - line 144-143 (xfmr 2) out: ratio 1.0000

### Base-case flows (N-1 dispatch), top 15 by |flow|
  220-471               529.17 MW (  35.4% of 1494 MW)
  440-386              -511.24 MW ( -34.2% of 1494 MW)
  16-14 (xfmr 1)        444.45 MW (  74.6% of 596 MW)
  386-14               -435.45 MW ( -29.1% of 1494 MW)
  386-14 (2)           -435.12 MW ( -29.1% of 1494 MW)
  263-262 (xfmr 1)     -411.02 MW ( -62.7% of 656 MW)
  57-123               -406.66 MW ( -27.2% of 1494 MW)
  15-14 (xfmr 1)       -405.49 MW ( -59.5% of 682 MW)
  429-428 (xfmr 1)     -388.03 MW ( -60.1% of 646 MW)
  169-165 (xfmr 1)      387.00 MW (  71.3% of 543 MW)
  82-80 (xfmr 1)        361.90 MW (  57.6% of 628 MW)
  9-7 (xfmr 1)          340.00 MW (  50.0% of 680 MW)
  9-7 (xfmr 2)          340.00 MW (  50.0% of 680 MW)
  17-14 (xfmr 2)        337.00 MW (  50.0% of 674 MW)
  17-14 (xfmr 1)        337.00 MW (  50.0% of 674 MW)

### LMPs [$ / MWh]
(dual = optimal dual of the N-1 LP balance rows -- pu-scaled, so divided by S_BASE --; FD = central finite difference of the N-1 optimal cost, +/-0.05 MW)
    bus     load    LMP(FD)    LMP(dual)
      1     0.00     27.280       27.280
      4   132.69     27.280       27.280
     22   148.70     27.280       27.280
     59   140.94     27.280       27.280
    142   161.34     27.280       27.280
    321   150.42     27.280       27.280
    327   121.07     27.280       27.280
    424   157.09     27.280       27.280
    469   128.53     27.280       27.280
    474   167.22     27.280       27.280
    491   119.53     27.280       27.280
    499   122.33     27.280       27.280

### N-1 including generator outages (TAMU contab CT_TGEN)
Infeasible: no dispatch keeps post-contingency flows within limits for every line AND every generator outage (proven infeasible by HiGHS, confirmed by an independent HiGHS call through scipy).
Per-outage test (all 591 line outages + each single gen outage) isolates the cause to 2 of 56:
  - gen bus 9 (Pmin 231.5 / Pmax 771.8 MW, must-run): at the N-1 dispatch its outage forces line 44-1 to carry 443 MW vs limit 320.3 MW (ratio 1.384)
  - gen bus 17 (Pmin 266.7 / Pmax 888.9 MW, must-run): at the N-1 dispatch its outage forces line 44-1 to carry 453 MW vs limit 320.3 MW (ratio 1.413)
They are the two largest generators in the case (both with Pmin > 0); their post-outage imports at their buses exceed what the connector paths can carry simultaneously with the line-outage security requirements.
Excluding those 2 gen outages, the 645-contingency SCOPF is feasible: $81,955.3/h (+$11,162.3/h vs N-0).

### Data / assumptions
- network: ACTIVSg500.RAW via VeraGridEngine RAW parser (500 buses, 466 lines, 131 transformers, 90 generators, 206 loads)
- costs/limits: official MATPOWER case (the RAW has no P records); 56 online units, 34 offline
- primary contingency set: 591 TAMU line outages (classical N-1 security constraint; gen outages -- see diagnostic above)
- 253 of 647 contingencies disconnect an island (radial cut-edge outages); post-contingency flow limits are checked on the network that remains connected to the slack, i.e. the island is assumed shed/handled locally (standard contingency practice)
  island sizes: max load 161.3 MW, max gen Pmax 888.9 MW, islands with generation: 36
- lossless DC, no voltage limits, no corrective re-dispatch
- PL cost: 20-segment chord (upper bound) of quadratic gencost
