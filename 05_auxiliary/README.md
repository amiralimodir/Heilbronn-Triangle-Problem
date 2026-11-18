# Auxiliary computations

Three standalone procedural files provide the auxiliary computations.

| Computation | File |
|---|---|
| Random feasible input lower bounds | [max_random.py](max_random.py) |
| Rectangular subproblems and strip-derived y-bounds | [heilbronn_in_rectangle.py](heilbronn_in_rectangle.py) |
| Proposition 5, complementary Case B | [proposition5_case_b.py](proposition5_case_b.py) |

## Random lower bounds

```console
python 05_auxiliary/max_random.py --list
python 05_auxiliary/max_random.py --n 9 --trials 500000 --seed 12345
python 05_auxiliary/max_random.py --n 3 4 5 6 7 8 9 10 11 12 --seed 12345
```

With no filters, every n=3,...,12 is processed sequentially using 500,000
attempted configurations per n. The sampling order is left, right, left,
right, bottom, bottom, top, top, then exactly n-8 interior points when n>8.
The distance rejection filter remains 0.167718; `--min-distance` changes it.
An attempt limit prevents a rejection loop from running indefinitely.

The supplied loop always appended two interior points. That made requested
n=9,11,12 configurations contain ten points. This script keeps n fixed and
generates exactly the requested number. It does not retroactively reproduce
those incorrect point-count runs.

Each instance uses seed `base_seed + n`. If `--seed` is omitted, the generated
base seed is saved. Results go under `05_auxiliary/results/max_random/` unless
`--output` is supplied. Every run saves its settings, best coordinates, exact
minimum triangle area, and a feasible lower bound rounded down to seven
decimal places (`feasible_lower_bound_7dp`). The exact calculation uses the
saved decimal coordinates. Batch summaries are updated after each instance.

No solver or third-party package is needed. The boundary sampling and distance
filter are heuristic choices; the saved feasible configuration supplies the
lower bound. They do not prove optimality. To use a result, copy the selected
lower-bound value into the appropriate experiment's input file. No bound file
is overwritten automatically.

## Rectangular subproblems

```console
python 05_auxiliary/heilbronn_in_rectangle.py --list
python 05_auxiliary/heilbronn_in_rectangle.py --n 5 --m 4
python 05_auxiliary/heilbronn_in_rectangle.py --n 3 4 5 6 7 8 9 --m 2 3 4 5 6 7 8 9 10 --threshold 0.0548756 --target-n 9
```

Here `--n` is the number of points in the auxiliary rectangle and `--m` is M
in its dimensions 1 by 1/M. The default selection is n=3,...,9 and M=1,...,7.
The objective bounds are
zero and 1/(2M). No time limit is set unless `--time-limit` is supplied; solver
settings otherwise use the installed Gurobi version's defaults.

The file contains one model: the unrestricted rectangular problem. It orders
the y-coordinates to fix labels and defines every product w[i,j]=x[i]*y[j].
It does not fix any point to a boundary or force two points onto an edge.

Provide `--threshold` only when it is a feasible lower bound for the intended
target instance. For each M, the sweep stops when a finite numerical **ObjBound**
is strictly below this threshold. Incumbent objectives are saved separately and
never used to establish exclusion. A time-limit run can establish a numerical
exclusion if its finite dual bound satisfies the comparison. Missing or
unusable bounds do not yield a capacity.

With `--target-n`, successful comparisons can generate
strip capacities and bounds on the ordered y-coordinates. These remain
numerical solver results under the paper's stated qualifications. Endpoint and
near-boundary equalities belong to Group 2 and are not added by this step.
Outputs include per-run settings, model, log, status, primal/dual values and
coordinates, plus a batch summary. Generated bounds do not modify the fixed
G3 tables or the comparison/final input files.

Both scripts support `--help` and `--list` without starting a computation.

## Proposition 5, Case B

The third file runs the nine-point complementary Case B:

```console
python 05_auxiliary/proposition5_case_b.py --check
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

Case A is the n=9 instance of `04_final_results/final_results.py`. Its saved
coordinates supply the independently evaluated configuration value; this
file supplies the complementary Case-B solver bound.

The reported reference run has objective/bound approximately
0.05463409383500306 and elapsed time 6576.674431800842 seconds. Its complete
log/model record is not yet included. Numerical solver bounds retain the
qualifications stated in the paper.


