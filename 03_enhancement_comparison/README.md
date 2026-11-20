# Enhancement-group comparison

The default command runs all 16 group combinations for n=7,8,9 sequentially
(48 runs). Filters select a subset:

```console
python 03_enhancement_comparison/enhancement_comparison.py --list
python 03_enhancement_comparison/enhancement_comparison.py --check
python 03_enhancement_comparison/enhancement_comparison.py
python 03_enhancement_comparison/enhancement_comparison.py --n 9 --groups 123
python 03_enhancement_comparison/enhancement_comparison.py --n 7 8 --groups none 1 2 12
python 03_enhancement_comparison/enhancement_comparison.py --n 7 --groups 1234 --build-only
```

`--groups` accepts `none`, `1`, `2`, `3`, `4`, `12`, `13`, `14`, `23`, `24`,
`34`, `123`, `124`, `134`, `234`, or `1234`. The complete design has 48 runs.
Use an ASCII-only checkout path, such as `C:/Heilbronn`, for reliable Gurobi
log/model export on Windows.

| Group | Constraints |
|---|---|
| G1 | Objective bounds, ordering/product ordering, reflection restriction |
| G2 | Instance-specific boundary restrictions |
| G3 | Fixed bounds on ordered y-coordinates |
| G4 | Horizontal strips, vertical strips, and square occupancy with mandatory assignments |

Ordering is also activated when G2 or G3 requires it. With G1 off, objective
bounds and signed-area floors use [0,0.5]. G4 still uses the common positive
lower bound to construct its partitions, and requires an input above 1e-6.
G3 uses the fixed paper table independently of the initialization LB.

`bounds_by_n.json` contains the recorded common comparison inputs:

| n | LB | UB |
|---|---:|---:|
| 7 | 0.066619 | 0.125 |
| 8 | 0.049992 | 0.083859 |
| 9 | 0.025724 | 0.072377 |

These are the actual run inputs reported in the manuscript, not rounded
display versions of more precise inputs.

`solver_settings.json` uses Gurobi 11.0.1, TimeLimit=86400, MIPGap=1e-4,
MIPGapAbs=1e-10, Threads=0, Seed=0, NonConvex=-1, and the recorded default
feasibility tolerances. Other parameters use that solver version's defaults.
Keep these comparison settings separate from the tuned final-run settings.

Each run creates a directory under `results/` containing its inputs, parameters,
model, log, status, primal/dual bounds, gap, solution, and exact coordinate-area
recomputation. `--build-only` exports the model without optimization.
For timing comparisons, run sequentially with the same hardware and software.
Batch runs also save `summary.csv` and `summary.json`, with the dual bound,
solver runtime, configuration area, and result path for each completed case.
They stop if a case fails; completed results remain saved.
`--list` and `--check` do not require Gurobi and do not optimize.

