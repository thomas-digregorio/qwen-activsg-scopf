"""
IEEE 14-bus test case (public domain).

Provenance: MATPOWER `case14`, retrieved 2026-11-07 from
https://raw.githubusercontent.com/MATPOWER/matpower/master/data/case14.m
It was converted from the IEEE Common Data Format (ieee14cdf.txt, UW archive)
by cdf2matp rev. 2393.

Modeling notes:
- baseMVA = 100.
- case14 defines NO branch thermal ratings (all rateA = 0, i.e. unlimited in
  MATPOWER). A uniform limit of 100 MW (1.0 pu) per branch is ASSUMED for the
  OPF / N-1 study; see LINE_LIMIT_MW.
- gencost rows follow MATPOWER model 2:
  [2, startup, shutdown, n, c2, c1, c0] with cost = c2*P**2 + c1*P + c0 ($/h).
"""

BASE_MVA = 100.0

# bus_id: (Pd [MW], Qd [Mvar], type)   type: 2=PQ load bus, 1=PV, 3=slack
BUSES = {
    1:  (0.0,  0.0,  3),
    2:  (21.7, 12.7, 2),
    3:  (94.2, 19.0, 2),
    4:  (47.8, -3.9, 1),
    5:  (7.6,  1.6,  1),
    6:  (11.2, 7.5,  2),
    7:  (0.0,  0.0,  1),
    8:  (0.0,  0.0,  2),
    9:  (29.5, 16.6, 1),
    10: (9.0,  5.8,  1),
    11: (3.5,  1.8,  1),
    12: (6.1,  1.6,  1),
    13: (13.5, 5.8,  1),
    14: (14.9, 5.0,  1),
}

# (fbus, tbus, r, x, b)  -- b = line charging susceptance [Mvar]
BRANCHES = [
    (1,  2,  0.01938, 0.05917, 0.0528),
    (1,  5,  0.05403, 0.22304, 0.0492),
    (2,  3,  0.04699, 0.19797, 0.0438),
    (2,  4,  0.05811, 0.17632, 0.0340),
    (2,  5,  0.05695, 0.17388, 0.0346),
    (3,  4,  0.06701, 0.17103, 0.0128),
    (4,  5,  0.01335, 0.04211, 0.0000),
    (4,  7,  0.00000, 0.20912, 0.0000),
    (4,  9,  0.00000, 0.55618, 0.0000),
    (5,  6,  0.00000, 0.25202, 0.0000),
    (6,  11, 0.09498, 0.19890, 0.0000),
    (6,  12, 0.12291, 0.25581, 0.0000),
    (6,  13, 0.06615, 0.13027, 0.0000),
    (7,  8,  0.00000, 0.17615, 0.0000),
    (7,  9,  0.00000, 0.11001, 0.0000),
    (9,  10, 0.03181, 0.08450, 0.0000),
    (9,  14, 0.12711, 0.27038, 0.0000),
    (10, 11, 0.08205, 0.19207, 0.0000),
    (12, 13, 0.22092, 0.19988, 0.0000),
    (13, 14, 0.17093, 0.34802, 0.0000),
]

# bus: (Pmin [MW], Pmax [MW], c2 [$/(MWh^2)], c1 [$/(MWh)], c0 [$])
GENS = {
    1: (0.0, 332.4, 0.0430292599, 20.0, 0.0),
    2: (0.0, 140.0, 0.25,         20.0, 0.0),
    3: (0.0, 100.0, 0.01,         40.0, 0.0),
    6: (0.0, 100.0, 0.01,         40.0, 0.0),
    8: (0.0, 100.0, 0.01,         40.0, 0.0),
}

LINE_LIMIT_MW = 100.0  # assumed uniform rating (case14 defines none)
SLACK_BUS = 1
