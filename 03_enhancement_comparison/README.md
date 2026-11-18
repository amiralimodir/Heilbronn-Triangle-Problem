# Enhancement-group comparison

One complete model file supports n=7,8,9 and all 16 group combinations:

```console
python 03_enhancement_comparison/enhancement_comparison.py --n 9 --groups 123
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
bounds and signed-area floors use [0,0.5]. G3/G4 can still use the common lower
bound geometrically. G4 requires a lower-bound input above 1e-6. Changing the
input LB does not recompute the fixed G3 table.

`bounds_by_n.json` contains the recorded common comparison inputs:

| n | LB | UB |
|---|---:|---:|
| 7 | 0.06661976335812875 | 0.125 |
| 8 | 0.0499921020091318 | 0.08386 |
| 9 | 0.02572465439694456 | 0.072377 |

`solver_settings.json` uses Gurobi 11.0.1, TimeLimit=86400, MIPGap=1e-4,
MIPGapAbs=1e-10, Threads=0, Seed=0, NonConvex=-1, and the recorded default
feasibility tolerances. Other parameters use that solver version's defaults.
Keep these comparison settings separate from the tuned final-run settings.

Each run creates a directory under `results/` containing its inputs, parameters,
model, log, status, primal/dual bounds, gap, solution, and exact coordinate-area
recomputation. `--build-only` exports the model without optimization.
For timing comparisons, run sequentially with the same hardware and software.
The complete result collection is still being finalized.

