# Solver comparison

`solver_comparison.py` runs every row of the paper's solver-comparison table:
Gurobi, BARON, SCIP, and Couenne for each of n=3,4,5,6,7 (20 runs).

```console
python 02_solver_comparison/solver_comparison.py --list
python 02_solver_comparison/solver_comparison.py
```

Runs are sequential. Select a smaller set when needed:

```console
python 02_solver_comparison/solver_comparison.py --solvers gurobi baron --n 3 4 5 6 7
python 02_solver_comparison/solver_comparison.py --solvers scip couenne --n 7
python 02_solver_comparison/solver_comparison.py --missing-duals
```

`--missing-duals` selects only Couenne n=3…7 and SCIP n=7; it is not the default.
`--all` is an optional alias for the default full schedule. `--help` and `--list`
do not load solver packages or start optimization.

## Model and settings

All cases use LB=0, UB=0.5, no enhancement groups, and a requested time limit
of 86400 seconds per case. Other solver settings retain their defaults.
The comparison environment uses Gurobi 11.0.1, Pyomo 6.8.2, SCIP 10.1.0,
Couenne 0.5.8, and BARON 22.4.20 with Python 3.8.12.

Gurobi uses the native Approach-1 implementation, including the two
sign-strengthening inequalities in the paper. BARON, SCIP, and Couenne use
the supplied Pyomo implementation, which contains the two big-M disjunction
inequalities but omits those additional strengthening inequalities. Both model
variants express the same integer feasible configurations, but their continuous
relaxations and input representations differ. This distinction must be retained
when interpreting the solver comparison. This script produces new runs; it does
not promise identical search paths or timings to previously reported runs.

Edit `SCIP_EXECUTABLE`, `COUENNE_EXECUTABLE`, and `BARON_EXECUTABLE` near the top
of the file to match the installed executables. Gurobi uses `gurobipy` and its
installed license. The script checks the selected dependencies before the batch.

## Outputs

The default output location is `~/heilbronn_solver_comparison`. Override it with
`--output C:/heilbronn_results/solvers`; use a writable path containing ASCII
characters only. Each invocation creates a new folder, and every solver/n case
gets its own files:

- `run_config.json`: input bounds, settings, model variant and timing convention;
- `solver.log` and exported model/interface files;
- `result.json`: objective, dual bound, coordinates, elapsed time, termination
  status, versions and any recovery errors;
- Gurobi additionally saves `model.mps`, `parameters.prm` and, if available,
  `solution.sol`; Pyomo retains its solver-interface files and result record.

After every case the script updates `solver_results.json` and `solver_table.csv`.
The CSV contains n, solver, dual bound, recorded seconds, and termination status;
it never substitutes the incumbent objective for a missing dual bound. Missing
bounds remain JSON `null` and blank CSV cells. Errors are saved and the batch
continues; use `--stop-on-error` to stop after the first failure or missing bound.

Gurobi time measures the `model.optimize()` call, while the external-solver time
measures `opt.solve()`, including interface export and solution loading. Requested
time limits do not imply exactly 86400 elapsed seconds. Gurobi's own `Runtime`
is also stored separately.

Couenne internally minimizes `-z`; its final global lower bound is negated to
recover the maximization upper bound. The reader includes an allowance for the
printed bound's precision. This accounts for output rounding, not all numerical
error in the optimization algorithm.

