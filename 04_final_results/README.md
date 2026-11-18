# Final Results: selected Approach 1 configuration

`final_results.py` contains the complete procedural model and all its helpers in
one file. Select one point count from 3 through 12 with `--n`. It is derived from
the prepared `final_G123_all_constraints_2026-11-16` package; the model/helper
function bodies and their order are unchanged. The ten instance-specific headers
were consolidated into point-count-dependent flags. See `provenance.json` for
the SHA-256 hashes of all ten source files.

This is the selected configuration for new final runs. It is not claimed to be
the exact executed source/settings behind every existing table row, nor the
fastest configuration for every possible instance. In particular, the archived
10-, 11-, and 12-point results used the separate historical variant described
below. Enhancement-comparison experiments retain their own settings in the
corresponding repository section.

## Inputs and execution

Use Python 3.8 or newer and the reference Gurobi version 11.0.3, with a valid
Gurobi license. Matplotlib is optional for configuration images. Run from this
directory, or supply the full path to the script; JSON inputs are located beside
the script rather than relative to the terminal's working directory.

1. Edit `bounds_by_n.json`. For every instance you intend to run, replace
   `lower_bound: null` with a feasible lower bound you have independently
   established. Missing values deliberately cause an error; they are not guessed.
2. Check the corresponding `upper_bound` and its documented source in the same
   file. These defaults are prepared inputs, not claimed historical run inputs.
3. Validate and then launch the selected instance:

```console
python final_results.py --n 9 --check
python final_results.py --n 9
```

To inspect/export a model without optimizing, use:

```console
python final_results.py --n 12 --build-only
```

`--help` and `--check` do not import Gurobi or start optimization. `--check` still
requires a completed bounds entry. Importing the file does not launch a run.
`--build-only` requires Gurobi and constructs/exports the model, but does not
optimize it. Repeat the command with the desired `--n` for other instances.
Runs are sequential; a process-held lock prevents concurrent runs from this
folder. It does not control unrelated jobs on the computer.

Zero LB is accepted literally and remains zero. There is **no automatic
`log(n)/n**2` substitution** in this prepared profile. The fixed y-bound table is
not recomputed when the LB changes. For a controlled comparison, preserve the
same per-instance inputs and record any change explicitly.

## Solver profile

`solver_settings.json` retains the prepared reference profile exactly:

| Parameter | Value |
|---|---:|
| TimeLimit | 86400 seconds |
| MIPGap | 1e-6 |
| MIPGapAbs | 0 |
| Threads | 32 |
| Seed | 0 |
| Symmetry | 2 |
| ConcurrentMIP | 1 |
| NonConvex | -1 |
| FeasibilityTol | 1e-6 |
| IntFeasTol | 1e-5 |
| OptimalityTol | 1e-6 |

The script records the requested profile and actual parameter values. Its
one-day protocol requires `TimeLimit=86400`. Using another Gurobi version is
recorded as a deviation from the 11.0.3 reference. Runtime depends on the machine,
software version, model construction, parameters, and competing workloads.

## Active restrictions

Groups 1, 2, and 3 are selected; group 4 (square occupancy and both strip
orientations) is disabled. No group-4 variables or constraints are constructed.
The original point-count eligibility rules are preserved:

| n | Instance-dependent restrictions |
|---|---|
| 3 | Coordinate/product ordering, x[0] <= 0.5, objective bounds; no boundary fixing. |
| 4, 5 | Contact with each square side, in addition to applicable group-1 restrictions. |
| 6 | Boundary contact and the fixed y-bound table. |
| 7, 8 | Two bottom points, one top point, vertical-edge contact, fixed y-bounds. |
| 9 | Case A: y[0]=0, y[1]<=0.01, y[7]>=0.99, y[8]=1; two selected points near each vertical edge; fixed y-bounds. No x[1]>=x[0] restriction. |
| 10, 11, 12 | Two selected points on each exact edge; bottom x-order and reflection restriction; fixed and relative y-bounds; median-y restriction; edge-sign fixing; left/right half-balance. |

Indices in this table are zero-based, matching the code. The nine-point script
is case A, not the complementary B computation used for Proposition 5.

For n=10,11,12 the additional features are explicit at the start of the execution
block: `use_bottom_order`, `use_relative_y_gaps`, `use_y_median`,
`use_edge_sign_fixing`, and `use_half_balance`. Relative bounds have the original
form `y[i+k-1] >= y[i] + 1/M`. The median restriction uses `y[n//2] >= 0.5`.
Bottom and top triangle-sign variables are fixed consistently with the original
signed-area expression. Half-counts are 5/5, 5/6 or 6/5, and 6/6; a point on
x=0.5 is assigned to one half. These rules do not impose full mirror-paired or
central symmetry of the configuration.

The larger-n restrictions are exploratory and are not asserted to preserve an
unrestricted optimum. This includes the legacy y-bound coefficients and the
inferred n=12 table entry `(9,2)`. The n=9 interpretation likewise depends on
the complementary-class exclusion and the numerical-bound qualifications in
the paper.

## Output and interpretation

Each invocation creates its own directory under `results/`, containing the
executed source, configuration snapshots, solver log, model, parameter file,
machine/software metadata, and JSON/text results. A solution file and a plot
are written when available. Solver output uses relative ASCII filenames inside
that directory, avoiding the previously encountered Unicode log-path failure.

`ObjVal`, `ObjBound`, termination status, solver Runtime, work/node statistics,
and model fingerprint are recorded separately. A time limit without an incumbent
still preserves the available dual bound. The saved decimal coordinates are
checked against the unit square and every triangle area is recomputed with
exact rational arithmetic; no clipping or rounding is used in that check.
Plotting failures do not discard numerical results. Solver time excludes later
plotting and reporting; total elapsed time is also recorded.

The numerical dual bound concerns the implemented model. For the restricted
larger-n models it is not an unrestricted Heilbronn upper bound. A verified
configuration in the unit square supplies a feasible lower bound for the
unrestricted problem independently of those search restrictions. A solver
incumbent objective is not substituted for the recomputed configuration value.

## Historical n=10,11,12 variant

`historical/legacy_n10_n12.py` preserves the author's supplied `final code.txt`
with an explanatory comment header only. It includes the full original code,
including its other approach functions; the archived n=10,11,12 results concern
Approach 1. This is an interactive historical script, not the new entry point.

Its Approach 1 settings are `TimeLimit=86000`, `MIPGap=1e-6`, `MIPGapAbs=0`,
and `Symmetry=2`; Threads is not explicitly set. An input LB of zero is replaced
by `log(n)/n**2`. Its legacy packing functions lack the exactly-one assignment
constraints and therefore do not enforce the intended occupancy restrictions.
These historical behaviors are retained rather than silently presented as the
corrected group-4 formulation.

The short result files provided for these runs did not record their input LB/UB
values. `historical/provenance.json` leaves those inputs unknown. Consequently,
the available materials do not establish an exact reexecution recipe for these
older rows. New results from `final_results.py` must be identified as new runs,
with their own source/configuration snapshots, rather than attached to the old
runtime or bound values.

No solver optimization was performed while preparing this repository section.
