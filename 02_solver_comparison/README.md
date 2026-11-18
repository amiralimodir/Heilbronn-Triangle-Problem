# Solver comparison

`solver_comparison.py` contains the complete Pyomo Approach-1 model and the
corrected dual-bound reader.

```console
python 02_solver_comparison/solver_comparison.py
python 02_solver_comparison/solver_comparison.py --all
```

Edit the solver executable paths and `RUNS_TO_REPEAT` near the beginning of the
file. The default selection repeats only the cases whose dual bounds remain
outstanding: Couenne n=3…7 and SCIP n=7. Each case has a requested time limit
of 86400 seconds. The model uses LB=0 and UB=0.5 without enhancement groups.
Other solver parameters use their respective defaults.

Use `--all` for the full set: SCIP n=3…7, Couenne n=3…7, and BARON n=7.
The executable paths for all three solvers are editable constants near the top
of the file.

The comparison environment is Python 3.8.12, Pyomo 6.8.2, SCIP 10.1.0,
Couenne 0.5.8, and BARON 22.4.20. The output directory is configured by
`OUTPUT_ROOT`; use an ASCII-only writable path.

Each run saves its log, model/interface files, parameters, versions, termination
status, objective, available dual bound, coordinates, and elapsed time. Results
are checkpointed after each case. Bound-recovery failures stop the batch and
remain missing values; they are not replaced by incumbent objectives.

Couenne solves an internal minimization of `-z`, so its final global lower bound
is negated to obtain the maximization upper bound. The reader accounts for the
printed bound's precision. This does not independently certify all numerical
error in the optimization algorithm.

`--help` prints usage without starting a solver.

