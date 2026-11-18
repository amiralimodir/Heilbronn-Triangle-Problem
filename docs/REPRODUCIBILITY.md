# Reproducing the experiments

The numbered directories each provide one main model or validation file.
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

Use `--list` to inspect the first four experiment selections. Their default
commands cover the complete tables, and filters permit individual reruns.
Their READMEs document model variants, settings, and generated summaries.
The exact coordinate-area verifier accepts user-supplied files and saved run
outputs. The repository does not require bundled numerical data.

The final max-random lower-bound and auxiliary y-bound computation sources
are still to be supplied. Until then, the main models consume the specified
bounds; they do not generate those inputs. Executable coverage does not imply
that every numerical run has already finished or that every paper entry has
been regenerated.

Each reported result must retain its own settings, inputs, status, primal/dual
bounds, gap, runtime, and solver version. Do not assign a different run's timing
or settings to it. Missing values remain unfilled until the corresponding data
is available.

The selected final script uses G1–G3. Unproved restrictions for n=10…12 mean
that its numerical dual bounds concern the restricted models. Explicit feasible
coordinates can independently establish lower bounds for the full problem.
