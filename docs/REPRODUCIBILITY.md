# Reproducing the experiments

Each experiment has one main model or validation file in a numbered directory.
Their READMEs identify inputs, solver settings, execution commands, and outputs.
The comparison and final-run parameter files are intentionally separate.

| Paper experiment | Cases supplied by the code |
|---|---|
| Unenhanced formulation comparison | Approaches 1, 2, 3 (MILP, P=10), n=3…7 |
| Shared-enhancement comparison | Approaches 1 and 3 (MILP, P=10), n=7,8 |
| Solver comparison | Gurobi, BARON, SCIP, Couenne, n=3…7 |
| Enhancement comparison | All 16 G1–G4 combinations, n=7,8,9 |
| Final results | G1–G3, n=3…12 |
| Proposition 5 | Case B; Case A uses the final n=9 script |
| Random feasible lower bounds | Boundary-first random configurations, n=3…12 |
| Auxiliary rectangles | Unrestricted rectangle model for selected point counts and heights 1/M |

Use `--list` to inspect the first four experiment selections. Their default
commands cover the complete tables, and filters permit individual reruns.
Their READMEs document model variants, settings, and generated summaries.
The exact coordinate-area verifier accepts user-supplied files and saved run
outputs. The repository does not require bundled numerical data.

The auxiliary scripts now generate random feasible lower bounds and solve
rectangular subproblems. The random script saves its seed and the actual best
configuration; its point count remains exactly the requested n. Rectangle
exclusions use the numerical dual bound, not the incumbent objective.

The rectangle script solves the unrestricted problem, with label ordering
and no boundary fixing. Its numerical dual comparisons can generate the
strip-derived y-bounds under the paper's numerical qualifications. Generated
outputs never overwrite the fixed G3 tables or existing LB/UB settings.

Executable coverage does not imply that every numerical run has already
finished or that every paper entry has been regenerated.

Each reported result must retain its own settings, inputs, status, primal/dual
bounds, gap, runtime, and solver version. Do not assign a different run's timing
or settings to it. Missing values remain unfilled until the corresponding data
is available.

The selected final script uses G1–G3. Unproved restrictions for n=10…12 mean
that its numerical dual bounds concern the restricted models. Explicit feasible
coordinates can independently establish lower bounds for the full problem.
