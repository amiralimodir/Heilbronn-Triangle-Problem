# Reproducibility coverage

This file maps the requested evidence to what is actually available. A runnable
new recipe, a reconstructed historical recipe, and a saved execution archive
are different objects and are identified separately.

| Requested material | Available evidence | Remaining material |
|---|---|---|
| Solver versions and parameters | Section READMEs, configuration JSONs, two complete n=9 solver logs | Exact versions for historical console-only rows where not recorded |
| Source/model files | Complete procedural scripts; MPS and executed source for the two archived n=9 runs | Historical main 633.331-second run's exact source/model archive |
| Main n=9 final status, objective, dual, gap | Two clearly identified additional n=9 archives; manuscript coordinates independently verifiable | Complete archive for the specific paper final n=9 timing; do not substitute either additional archive |
| Proposition 5 B experiment | Author-reported output summary and documented reconstruction with LB=0, UB=.0723767; new executions export logs/model/parameters/status | Full original log, model, and exact version for the reported 6576.674-second run |
| Proposition 8 / strip-derived y bounds | Historical auxiliary examples retained as legacy | Final executed y-bound source and associated bound/status records, to be supplied by the author |
| Input feasible lower-bound generation | Explicit editable final-run inputs and recorded comparison inputs | Final max-random source and frozen input choices, to be supplied by the author |
| Accurate coordinates | Appendix v15 decimal strings, figure coordinates, saved n=10…12 decimal coordinates | Other historical full-precision solution archives where absent |
| Independent area recomputation | Standard-library rational verifier and checked expected values for every appendix configuration | No solver required for these checks |
| Enhancement ablation | All 16 combinations per n supported; recorded comparison profile and one complete comparison archive | Complete 48-run result collection and revised paper table |
| Multi-solver comparison | Recorded finite bounds and times; targeted rerun code | Corrected SCIP n=7 and Couenne n=3…7 dual-bound extraction/results |

## Interpretation

* A solver incumbent objective can differ from the minimum area recomputed
  from the returned coordinates. Both are retained as distinct quantities.
* A numerical dual bound applies to the model actually constructed, subject
  to numerical tolerances. It is not automatically an exact mathematical
  certificate for a broader model with additional restrictions removed.
* The printed appendix coordinates form independently checkable feasible
  configurations. Every one of their triangle areas is evaluated exactly.
* The selected G1–G3 recipe disables G4. The archived n=10…12 source retains
  its historical packing implementation. The repository does not relabel
  those old results as executions of the newly selected recipe.
* Null final LB entries require a deliberate input. In the selected final
  script, zero means zero. In several preserved historical sources, entering
  zero triggers log(n)/n²; their READMEs document that behavior.
* A missing or infinite value returned by a solver interface is not evidence
  that the optimizer found no finite mathematical bound. Missing duals are
  recorded as unavailable until the corrected extraction is complete.

The paper's n=7 timing previously combined evidence from different runs, and
its current n=11/12 displayed time caps differ from the saved elapsed times.
Its n=10 figure also still uses the previous run's coordinates, while the
updated appendix uses the new run.
The repository retains the actual available elapsed times and does not infer
missing row-level provenance from a matching printed objective.

## Completing the archive

Add each finished run in a distinct directory with its exact executed source,
input/parameter snapshots, solver log, exported model and solution, termination
status, primal/dual bounds, gap, solver/runtime versions, and recomputed
configuration value. Keep comparison inputs/settings frozen across its design.
Update the relevant coverage row after the missing evidence is added.

