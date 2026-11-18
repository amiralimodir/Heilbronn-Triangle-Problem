# Author formulations and the additional n=7,8 comparison

`author_formulations_historical.py` keeps the author's complete procedural source
in one file: all formulations, helper functions, prompts, and plotting remain
together. Only explanatory header comments were added. The source body is
byte-for-byte identical to the supplied file; `provenance.json` records its hash.

## Current entry point

Run `python formulations.py` to use the same complete procedural source with a
corrected launcher. Approach 3 now accepts `MILP` or `MIQCP` (case-insensitive),
dispatches to the corresponding existing model function, and rejects unknown
choices. Importing this file does not prompt or solve. The plot title is now
"Computed configuration". All mathematical model/helper bodies are unchanged;
only that plot-title constant changes inside a function.

This launcher correction does not standardize the inherited formulations or
settings: the historical horizontal-strip assignments are still incomplete,
Approach 2 retains its active helper behavior, and Approach 2/exact MIQCP still
have no explicit time limit in their original model functions. It therefore
does not establish a controlled baseline comparison or retrofit an MIQCP run
into previously reported MILP results. For the selected final-run protocol use
`../04_final_results/final_results.py`; for corrected four-group comparisons use
`../03_enhancement_comparison/enhancement_comparison.py`.

To retain console output in PowerShell, run from a fresh working directory:

```powershell
python path/to/01_formulations/formulations.py 2>&1 | Tee-Object -FilePath console.log
```

This saves the existing console reporting; it does not add solver attributes
that the original functions did not print. The complete new-run result export
is provided by the dedicated final/group scripts.

## Supported historical comparison

This source was confirmed by the author as the implementation for the additional
Approach 1 versus Approach 3 MILP comparison at n=7 and n=8. The input sequence is:

```console
python author_formulations_historical.py
```

| Prompt | Approach 1 | Approach 3 |
|---|---|---|
| n | 7 or 8 | Same point count |
| Upper bound? | 0.5 | 0.5 |
| Lower bound? | 0 | 0 |
| Approach? | 1 | 3 |
| H | not requested | 10 |
| MIQCP or MILP? | not requested | MILP |

The code uses `H` for the binary expansion depth called `P` in the paper.
Entering LB=0 **does not apply a zero lower-bound input**: it is changed to
log(n)/n^2, namely 0.039712452021537004 for n=7 and 0.032491274088747434 for n=8.
Both approaches use it in the signed-area inequalities; Approach 1 additionally
sets the objective variable's lower bound to this number.

Both have TimeLimit=86400. Approach 1 uses MIPGap=1e-6 and MIPGapAbs=0; Approach 3
keeps the Gurobi default stopping-gap settings (MIPGap=1e-4, MIPGapAbs=1e-10 for
the recorded Gurobi 11 family). Do not make these settings identical and then
describe the results as the original experiment. The selected geometric subset
includes y-ordering, two bottom points, one top point, fixed y-bounds and product
ordering. The original horizontal-strip function lacks mandatory assignments,
so its indicators can all be zero. This historical behavior is preserved here.

The supplied result summaries report:

| n | Approach 1 numerical dual bound | Time (s) | Approach 3 numerical dual bound | Time (s) |
|---|---:|---:|---:|---:|
| 7 | 0.08385946233900 | 9.557921648025513 | 0.08399940673136 | 46.02381491661072 |
| 8 | 0.07237702163144 | 67.14493179321289 | 0.07252825572868 | 245.07035613059998 |

These are archived reported results, not fresh runs. All four summaries report
OPTIMAL within their stopping tolerance and coincident reported incumbent/dual
values. The Approach 3 model is a relaxation; its incumbent is not the minimum
triangle area of the returned coordinates. Floating-point equality of reported
incumbent and bound is not an exact-arithmetic certificate.

## Historical limitations that remain visible

The code contains functions for Approach 2 and exact-MIQCP Approach 3, but the
original launcher sends both MILP and non-MILP selections to the MILP function.
Consequently, this file must not be cited as having reproduced an exact-MIQCP
experiment. The Approach 2 branch also inherits active helper behavior and has
no explicit one-day time limit here. This is not a verified source for the
separate unenhanced baseline table. Those experiments need their own archived
sources and settings; they must not be reconstructed by silently toggling this
file and attaching old timings.

The 2026-10-30 corrected-rerun package changes assignments and aligns stopping
tolerances. That is a new protocol and is not substituted for this source.

The original launcher prompts interactively, writes a fixed log/image name in
the working directory, and shows the plot. Run each case from a fresh working
directory if retaining its outputs. Some original print statements call the
incumbent an "Optimal value"; interpret results with the status, log, numerical
tolerances and independent coordinate check. The original script does not
consistently export a finite dual bound on a time limit. It is preserved as
historical evidence, not presented as the final-run logging implementation.

Dependencies are Python, gurobipy with a valid solver license, and matplotlib.
Python 3.8 syntax was checked without optimizing. The exact patch version and
effective parameters of any new run should be archived; the summary files alone
do not establish a complete parameter/version record for these four old runs.
