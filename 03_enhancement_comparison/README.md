# Enhancement-group comparison

This directory corresponds to the enhancement-comparison experiments. All model
functions remain together in `enhancement_comparison.py`: the procedural model,
constraint insertion order, and numerical output processing are preserved from
the archived comparison source. The 48 formerly separate scripts differed in
their instance/group flags; the command-line selection reproduces those flags.

Importing the script does not prompt, create files, or start a solve.
`python enhancement_comparison.py --help` works without Gurobi installed.
The solver is imported only after CLI and numerical-input validation.

## One run

Use an ASCII-only checkout path, such as `C:/Heilbronn`. This source passes
absolute output paths to Gurobi; the recorded Windows Gurobi version can fail
to write logs/model files when a directory contains non-ASCII characters.
The script checks the path before importing Gurobi or creating run outputs.

From this directory:

```console
python enhancement_comparison.py --n 9 --groups 123
```

For a model export without optimization:

```console
python enhancement_comparison.py --n 7 --groups 1234 --build-only
```

`--n` accepts 7, 8, or 9. `--groups` accepts `none`, `1`, `2`, `3`, `4`,
`12`, `13`, `14`, `23`, `24`, `34`, `123`, `124`, `134`, `234`, or `1234`.
The full design is 16 combinations for each of three point counts (48 runs).
Run them sequentially, with no unrelated optimization workload if comparing time.
For example, in PowerShell, from this directory:

```powershell
$groups = 'none','1','2','3','4','12','13','14','23','24','34','123','124','134','234','1234'
foreach ($points in 7,8,9) {
    foreach ($group in $groups) {
        python enhancement_comparison.py --n $points --groups $group
        if ($LASTEXITCODE -ne 0) { throw 'Comparison run failed; inspect its output.' }
    }
}
```

## Recorded comparison profile

The source and both JSON inputs come from the saved n=9, G123 run started on
2026-11-14, using Python 3.8.12 and Gurobi 11.0.1. Its parameters are preserved:
TimeLimit=86400, MIPGap=1e-4, MIPGapAbs=1e-10, Threads=0, Seed=0,
FeasibilityTol=1e-6, IntFeasTol=1e-5, OptimalityTol=1e-6, NonConvex=-1.
Other Gurobi parameters retain their version-specific defaults. Threads=0 is the
automatic setting, not a fixed thread count. This is a comparison profile; the
tuned final-run settings must not be copied into this JSON.

The saved common input file contains:

| n | Configuration lower-bound input | Upper-bound input |
|---|---:|---:|
| 7 | 0.06661976335812875 | 0.125 |
| 8 | 0.0499921020091318 | 0.08386 |
| 9 | 0.02572465439694456 | 0.072377 |

These are recorded numerical inputs, not newly derived certificates. This whole
input file was archived with the n=9 run; the repository does not infer that every
outstanding comparison has finished or that all run archives are available.

## What each group changes

* G1: objective bounds, ordering/product ordering, and reflection restriction.
* G2: the instance-specific boundary restrictions and associated product equalities.
* G3: the fixed, instance-specific y-bounds from the implementation; they are not
  recomputed when the lower-bound input changes, and no dynamic callback is used.
* G4: horizontal strips, vertical strips, and squares, with the corrected
  exactly-one assignment constraints. The final cell width is retained.

Ordering is also included whenever G2 or G3 needs it. With G1 disabled, the
objective bounds and signed-area floors use [0,0.5], even when G3 or G4 uses a
positive lower-bound input geometrically. With no groups, or only G2, no positive
configuration lower-bound input is applied. G4 requires a positive input above
1e-6 for its grid construction.

For n=9, G2 activates the implemented near-boundary case-A restrictions. The
numerical bound concerns that implemented model; unrestricted interpretation
depends on the coverage argument and on the validity of all added restrictions.
The fixed G3 bounds retain the qualifications described in the manuscript.

## Outputs and scope

Each execution writes its own timestamped `results/` directory with the source,
input snapshots, model MPS, parameter file, log, run configuration, final status,
incumbent, dual bound, gap, saved coordinates, and independently recomputed exact
decimal-coordinate minimum area. Plotting is optional in the source; failures
there do not discard numerical results. `--build-only` outputs are not solved runs.

The archived tuned n=9 run using Gurobi 11.0.3, MIPGap=1e-6, MIPGapAbs=0 and
Symmetry=2 is a distinct experiment. It is not the default comparison profile.
The old three-group table cannot be relabelled as this four-group experiment.

`provenance.json` records source hashes and the packaging checks. All model
function ASTs match the archived source, and all 48 flag combinations match the
original manifest. No long optimization was run while packaging these files.
