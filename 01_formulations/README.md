# Formulation comparison

`formulations.py` runs both formulation-comparison tables from one file:

| Suite | Cases | Model settings |
|---|---|---|
| `baseline` | Approaches 1, 2, 3 for n=3,...,7: 15 runs | No geometric enhancements; objective bounds 0 and 0.5; default Gurobi stopping gaps |
| `shared` | Approaches 1 and 3 for n=7,8: 4 runs | The shared geometric subset described below |

The default is all 19 runs, each with TimeLimit=86400 seconds. Approach 3 uses
the MILP relaxation with P=10. Install `requirements-comparison.txt` and a
Gurobi license before solving.
Approach 2 retains its original `z >= 1e-10` non-collinearity constraint;
this effective lower limit is recorded separately from the zero input bound.

```console
python 01_formulations/formulations.py --list
python 01_formulations/formulations.py
python 01_formulations/formulations.py --suite baseline
python 01_formulations/formulations.py --suite shared
python 01_formulations/formulations.py --suite baseline --n 7 --approaches 1 3
```

`--list` and `--help` do not require Gurobi and do not solve anything.
Use `--output PATH` to choose the results directory.

## Shared comparison

Both approaches impose ordered y coordinates, y1=y2=0, yn=1, the n-specific
y bounds, and monotonicity of the product variables, with wi1=wi2=0.
The objective bounds are LB=0 and UB=0.5 in both models. Zero remains zero;
no automatic lower-bound formula is applied.

Both approaches use MIPGap=1e-6 and MIPGapAbs=0, with a requested
TimeLimit=86400 seconds and otherwise default Gurobi settings. These are
the settings of the shared comparison, not the final-run profile.

This selected subset does not include the rectangle or square packing
constraints of Group 4. Group 4 is implemented in
[the enhancement comparison](../03_enhancement_comparison/README.md).

## Outputs and optional runs

Each run gets its own directory containing `run_config.json`, `solver.log`,
`model.mps`, `parameters.prm`, `result.json`, and, when available, `solution.sol`.
The batch's `table_results.csv` collects dual bounds, times, statuses, and
configuration values. Dual bounds and termination statuses are recorded even
if no incumbent is found. Coordinates are stored as x/y arrays for independent
verification. Missing or nonfinite solver bounds are written as JSON null.

For Approach 3, the incumbent of the MILP relaxation need not equal the minimum
triangle area of its coordinates. `recomputed_minimum_area` is a separate
floating-point evaluation; use the area-verification script for exact decimal
coordinate evaluation.

`--approach3-mode miqcp` selects the exact version of Approach 3, and `--p`
changes the expansion depth. These optional variants are marked in the run
configuration and are not the MILP P=10 entries in the paper tables.
`--time-limit` can be used for a shorter exploratory run. All three formulations
receive the selected time limit, including Approach 2 and exact Approach 3.

