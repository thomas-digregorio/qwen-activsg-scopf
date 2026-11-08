### N-0 (base case only)
Cost: $7,931.9/h
  gen bus 1 :   154.78 MW
  gen bus 2 :    46.67 MW
  gen bus 3 :    50.00 MW
  gen bus 6 :    -0.00 MW
  gen bus 8 :     7.56 MW
Worst base-case flow: line 1->2 @ 100.0 MW (100.0% of limit)
Contingencies violated by this dispatch: 14
  - line 1-2 out: max flow 154.8 MW (154.8%) on line 1->5
  - line 1-5 out: max flow 154.8 MW (154.8%) on line 1->2
  - line 3-4 out: max flow 100.3 MW (100.3%) on line 1->2
  - line 4-5 out: max flow 112.8 MW (112.8%) on line 1->2
  - line 5-6 out: max flow 102.4 MW (102.4%) on line 1->2
  - line 6-11 out: max flow 100.2 MW (100.2%) on line 1->2
  - line 6-12 out: max flow 100.0 MW (100.0%) on line 1->2
  - line 6-13 out: max flow 100.2 MW (100.2%) on line 1->2
  - line 10-11 out: max flow 100.0 MW (100.0%) on line 1->2
  - line 12-13 out: max flow 100.0 MW (100.0%) on line 1->2
  - line 13-14 out: max flow 100.1 MW (100.1%) on line 1->2
  - gen bus 2 out: max flow 139.1 MW (139.1%) on line 1->2
  - gen bus 3 out: max flow 137.3 MW (137.3%) on line 1->2
  - gen bus 8 out: max flow 105.0 MW (105.0%) on line 1->2

### N-1 (all 25 contingencies enforced)
Cost: $8,438.4/h  (delta +506.4/h vs N-0)
  gen bus 1 :   100.00 MW   (N-0 was 154.78)
  gen bus 2 :    39.30 MW   (N-0 was 46.67)
  gen bus 3 :    41.67 MW   (N-0 was 50.00)
  gen bus 6 :    36.37 MW   (N-0 was -0.00)
  gen bus 8 :    41.67 MW   (N-0 was 7.56)
All 25 contingencies satisfied (violations: 0)
Worst post-contingency flow: gen bus 2 out -> 100.0 MW (100.0%) on line 1->2
Binding N-1 contingency(ies):
  - line 1-2 out: 100.0 MW on line 1->5
  - line 1-5 out: 100.0 MW on line 1->2
  - gen bus 2 out: 100.0 MW on line 1->2

### Base-case line flows (N-1 solution)
  line 1->2         67.07 MW (  67.1%)
  line 1->5         32.93 MW (  32.9%)
  line 2->3         37.57 MW (  37.6%)
  line 2->4         27.67 MW (  27.7%)
  line 2->5         19.42 MW (  19.4%)
  line 3->4        -14.96 MW ( -15.0%)
  line 4->5        -35.66 MW ( -35.7%)
  line 4->7         -4.87 MW (  -4.9%)
  line 4->9          5.44 MW (   5.4%)
  line 5->6          9.10 MW (   9.1%)
  line 6->11         8.34 MW (   8.3%)
  line 6->12         7.84 MW (   7.8%)
  line 6->13        18.08 MW (  18.1%)
  line 7->8        -41.67 MW ( -41.7%)
  line 7->9         36.79 MW (  36.8%)
  line 9->10         4.16 MW (   4.2%)
  line 9->14         8.58 MW (   8.6%)
  line 10->11       -4.84 MW (  -4.8%)
  line 12->13        1.74 MW (   1.7%)
  line 13->14        6.32 MW (   6.3%)

### LMPs [$ / MWh] (finite difference, N-1 problem)
  bus 1 :   28.606
  bus 2 :   41.083
  bus 3 :   40.927
  bus 4 :   40.792
  bus 5 :   40.694
  bus 6 :   40.727
  bus 7 :   40.774
  bus 8 :   40.774
  bus 9 :   40.765
  bus 10:   40.758
  bus 11:   40.743
  bus 12:   40.730
  bus 13:   40.733
  bus 14:   40.751

### Contingency detail (N-1 solution)
  contingency         max flow    util
  line 1-2 out          100.0 MW  100.0%
  line 1-5 out          100.0 MW  100.0%
  line 2-3 out           59.3 MW   59.3%
  line 2-4 out           59.5 MW   59.5%
  line 2-5 out           60.1 MW   60.1%
  line 3-4 out           70.2 MW   70.2%
  line 4-5 out           77.5 MW   77.5%
  line 4-7 out           67.2 MW   67.2%
  line 4-9 out           66.9 MW   66.9%
  line 5-6 out           67.6 MW   67.6%
  line 6-11 out          67.4 MW   67.4%
  line 6-12 out          67.1 MW   67.1%
  line 6-13 out          67.2 MW   67.2%
  line 7-8 out       (disconnected -> skipped)
  line 7-9 out           66.0 MW   66.0%
  line 9-10 out          66.9 MW   66.9%
  line 9-14 out          66.8 MW   66.8%
  line 10-11 out         67.2 MW   67.2%
  line 12-13 out         67.1 MW   67.1%
  line 13-14 out         67.3 MW   67.3%
  gen bus 1 out          67.1 MW   67.1%
  gen bus 2 out         100.0 MW  100.0%
  gen bus 3 out          98.2 MW   98.2%
  gen bus 6 out          90.0 MW   90.0%
  gen bus 8 out          94.5 MW   94.5%
