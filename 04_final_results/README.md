# Final runs

`final_results.py` contains the full Approach-1 model for n=3…12.
G1–G3 are enabled; G4 is disabled. With no filters, all ten instances run
sequentially, each with its own one-day solver limit.

1. Review the supplied `lower_bound` for each intended n in
   `bounds_by_n.json`; these are the author's chosen run inputs.
2. Check the upper bound and its source in the same file.
3. Validate and launch the selected instance:

```console
python 04_final_results/final_results.py --list
python 04_final_results/final_results.py --check
python 04_final_results/final_results.py
python 04_final_results/final_results.py --n 9 --check
python 04_final_results/final_results.py --n 9
python 04_final_results/final_results.py --n 10 11 12
```

`--check` validates inputs without loading Gurobi. `--build-only` constructs and
exports the model without optimizing. A zero LB stays zero in this script.

The reference solver is Gurobi 11.0.3. `solver_settings.json` specifies
TimeLimit=86400, MIPGap=1e-6, MIPGapAbs=0, Threads=32, Seed=0, Symmetry=2,
ConcurrentMIP=1, NonConvex=-1, FeasibilityTol=1e-6, IntFeasTol=1e-5, and
OptimalityTol=1e-6.

| n | Instance-specific restrictions |
|---|---|
| 3 | Ordering, reflection, and objective bounds |
| 4,5 | Contact with each square side |
| 6 | Boundary contact and fixed y-bounds |
| 7,8 | Two bottom points, one top point, vertical-edge contact, fixed y-bounds |
| 9 | Case A near-boundary restrictions with epsilon=0.01; fixed y-bounds |
| 10,11,12 | Two selected points on each exact edge, bottom ordering, reflection, fixed/relative y-bounds, median-y restriction, edge-sign fixing, and left/right half-balance |

The n=10…12 additional restrictions are unproved. Their dual bounds apply to the
implemented restricted models. Verified coordinates in the square still provide
feasible lower bounds for the unrestricted problem.

Each run writes a separate `results/` directory with configuration, source,
model, parameter, log, solution, environment, status, objective, dual bound,
gap, and exact minimum-area outputs. Runtime and total elapsed time are recorded
separately. Batches additionally save `summary.csv` and `summary.json` with a
row for each completed instance, and stop on a failed case without discarding
earlier outputs. The selected profile is for these runs; existing numerical results
must retain the input and parameter settings under which they were obtained.

`--check` checks the format and range of the supplied bounds; it does not prove
their mathematical validity. Upper-bound inputs are used exactly as supplied.
The n=9,...,12 entries identify the author's numerical caps without attributing
them to an inapplicable rounding or literature argument. A bound obtained for a
restricted predecessor model is not, by itself, an upper bound for the
unrestricted next instance.

