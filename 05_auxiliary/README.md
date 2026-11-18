# Auxiliary computations

These files relate to the **Enhancements** section, especially the nine-point
near-boundary exclusion experiment (Proposition 5) and the ordered-y-coordinate
strip arguments. Each runnable file retains its procedural source structure.

**Pending author files:** the final max-random lower-bound source and the final
strip/y-bound auxiliary source will be supplied separately. The old files below
are retained only under `legacy/` and `source_as_received/`; they are not active
entry points for reproducing the revised paper.

## Nine-point complementary Case B

```console
python 05_auxiliary/proposition5_case_b.py
```

`--help` prints usage without importing Gurobi or starting a solve. Unknown
arguments are rejected. Importing the file does not start a computation.

This file is a **reconstructed reproduction recipe** based on the author's
supplied source and subsequent explicit clarification. It is not presented as a
byte-for-byte snapshot of the historical server process. The received source is
preserved unchanged in `source_as_received/proposition5_case_b.txt`, with its
SHA-256 in `source_provenance.json`.

The fixed inputs and explicit solver choices are:

| Setting | Value |
| --- | --- |
| n | 9 |
| Input and effective objective LB | 0 |
| Objective UB | 0.0723767 |
| Boundary distance epsilon | 0.01 |
| MIPGap | 1e-6 |
| MIPGapAbs | 0 |
| Symmetry | 2 |
| Explicit time limit | None |
| Explicit threads, seed, feasibility tolerances | Unset; Gurobi defaults |

The source attachment included an automatic `LB=0` to `log(n)/n**2`
substitution. The author clarified that the executed final run used **LB=0**;
that substitution is therefore removed in this recipe. This is a documented
author clarification, not something recovered from a saved historical model.

The source already comments out `x[1] >= x[0]`; it stays inactive. The active
reflection restriction is `x[0] <= 0.5`. The model keeps the supplied ordered-y
and product constraints, `y[0]=0`, `y[1]>=0.01`, `y[8]=1`, and hardcoded n=9
ordered-coordinate bounds. The old vertical indicator links with sums at most
two remain in the source: setting all indicators to zero makes them vacuous.
The square/rectangle occupancy calls are inactive. No additional boundary-count
or point-swapping restriction has been introduced during packaging.

Changes beyond the clarified LB are output/documentation changes: fixed driver
inputs for this experiment, unique run directories, raw logging, model and
parameter export, environment/status/primal/dual/gap output, and a noninteractive
plot by default, plus a main guard and help-only argument handling. The mathematical constraint-creation statements and their order
are preserved. Unused original Approach-2/3 functions are retained so that the
base code has not been split or rewritten into another architecture.

The author-reported reference computation terminated with numerical status
`OPTIMAL`, objective and best bound approximately **0.05463409383500306**, and
elapsed time **6576.674431800842 seconds**. These reference values came from the
provided console excerpt; the complete log/model snapshot for that exact
historical run has not been recovered here. Do not relabel another Case-B run's
log as this one. In particular, the separately available file named
`Proposition 5 - Case B new - 1.txt` still enables the point-swapping inequality
and the automatic logarithmic LB; it is not an interchangeable executed snapshot.

The new script records its actual Gurobi/Python versions when run. It exports
`solver.log`, `model.mps`, `parameters.prm` (all nondefault parameters),
`run_config.json` (selected effective parameters and environment), `result.json`,
`result.txt`, the executed source, and available solution/plot files into
`05_auxiliary/runs/case_b_<timestamp>/`. The Gurobi log reports hardware and
thread use. Remaining parameters use the recorded Gurobi version's defaults.
An ASCII-only output path is required for reliable file export.

The relevant comparison is a solver numerical dual bound against an independently
evaluated feasible configuration value. A numerical `OPTIMAL` status or a zero
displayed gap is not an exact-arithmetic global certificate.

## Strip/y-bound computation: evidence gap remains

`legacy/rectangle_legacy.py` preserves the available old `heilbronn in rectangle.py`
mathematical code, with an explanatory header. It is a **legacy reference**, not
a certified or confirmed reproduction of the revised auxiliary experiment.
Specifically it:

- forces two points onto the bottom edge for every n;
- uses `x >= c2/M` for a selected right-side point;
- returns the incumbent objective on numerical `OPTIMAL` status, rather than
  separately exporting the global dual bound;
- uses a hardcoded threshold 0.083859 in its batch loop;
- lacks saved per-run versions, parameter snapshots, and model/log exports.

Those historical choices cannot simply be described as the unrestricted
rectangle auxiliary problem from the current manuscript. The exact executed
source and solver logs used to derive the reported strip thresholds remain
needed. Do not generate purported historical bounds or fill this gap by silently
changing the model. Raw reported strip-result files, where supplied, should be
kept as reported results with that provenance limitation. This source is kept
whole to expose the issue, not hidden behind a refactored implementation.

## Random lower-bound exploration: legacy helper

`legacy/random_search_legacy.py` preserves the supplied `max random.py` algorithm.
This particular version generates **11 points**, uses a fixed pairwise-distance
threshold **0.167718**, takes **100000 samples**, and has no fixed random seed.
It starts with five boundary points and adds six random points. It does not save
the best coordinates. Thus it is not the documented eight-boundary-point recipe
and cannot reproduce or substantiate the frozen input LBs used in the paper's
comparisons.

Historical bounds should instead be accompanied by their actual feasible
coordinates and an independent area evaluation. A new seeded random-search
experiment would be a new experiment and must be labelled accordingly. The
legacy helper requires `matplotlib` and `tqdm`; it is not run automatically by
the repository's principal experiments.

## Validation performed while packaging

Python 3.8 syntax checks were applied without starting any optimization job.
The Case-B model construction and all unchanged model helper functions were
compared structurally with the supplied source. Execution on the original server
and a repeat of its solve time have not been claimed. Source files in
`source_as_received/` preserve their received bytes; no license files or solver
binaries are included.
