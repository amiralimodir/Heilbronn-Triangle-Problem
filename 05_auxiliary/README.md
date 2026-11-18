# Auxiliary computations

The available main file runs the nine-point complementary Case B for
Proposition 5:

```console
python 05_auxiliary/proposition5_case_b.py
```

| Input or parameter | Value |
|---|---:|
| n | 9 |
| LB | 0 |
| UB | 0.0723767 |
| epsilon | 0.01 |
| MIPGap | 1e-6 |
| MIPGapAbs | 0 |
| Symmetry | 2 |

LB=0 remains zero. No explicit time limit or thread count is set; remaining
parameters use the installed Gurobi version's defaults. The restriction
`x[1] >= x[0]` is inactive; `x[0] <= 0.5` remains active.

Each run exports its model, parameters, solver log, versions, status, objective,
dual bound, gap, solution, and available plot into `05_auxiliary/runs/`.
Use an ASCII-only path for these outputs. `--help` and importing the file do
not start optimization.

The reported reference run has objective/bound approximately
0.05463409383500306 and elapsed time 6576.674431800842 seconds. Its complete
log/model record is not yet included. Numerical solver bounds retain the
qualifications stated in the paper.

The final max-random lower-bound and y-bound computation files will be added
when supplied by the authors.

