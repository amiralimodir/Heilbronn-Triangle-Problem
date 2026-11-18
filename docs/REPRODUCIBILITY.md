# Reproducibility checklist

The numbered directories each provide one main model or validation file.
Their READMEs identify inputs, solver settings, execution commands, and outputs.
The comparison and final-run parameter files are intentionally separate.

Available:

- Complete procedural model scripts for the listed experiments.
- Recorded comparison settings and selected final-run settings.
- Available solver logs, models, solutions, and numerical results under `data/`.
- Printed n=6…12 coordinates and an independent exact area verifier.

Still to be added:

- The complete enhancement-comparison result collection.
- Corrected Couenne n=3…7 and SCIP n=7 dual-bound results.
- Final max-random and y-bound computation files.
- Missing exact run inputs or complete logs for results currently represented
  only by console summaries, including the main n=9 and Proposition-5 runs.

Each reported result must retain its own settings, inputs, status, primal/dual
bounds, gap, runtime, and solver version. Do not assign a different run's timing
or settings to it. Missing values remain unfilled until the corresponding data
is available.

The selected final script uses G1–G3. Unproved restrictions for n=10…12 mean
that its numerical dual bounds concern the restricted models. Explicit feasible
coordinates can independently establish lower bounds for the full problem.
