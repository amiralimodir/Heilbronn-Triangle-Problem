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

The reference solver is Gurobi 11.0.1. `solver_settings.json` specifies
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

For n=10…12, the combined additional restrictions have not been proved to
preserve an unrestricted optimum. The dual bounds apply to the implemented
restricted models. Verified coordinates in the square still provide feasible
lower bounds for the unrestricted problem.

For n=10,11,12, `y_bounds` applies the fixed table entries and the relative
inequalities `y[i+r-1] - y[i] >= 1/kappa` for every `i=0,...,n-r`, using
the following `(r,kappa)` pairs:

| n | (r,kappa) |
|---|---|
| 10 | (3,11), (4,11), (5,5), (6,3), (7,2), (8,2) |
| 11 | (3,14), (4,14), (5,6), (6,4), (7,3), (8,2) |
| 12 | (3,16), (4,16), (5,6), (6,4), (7,3), (8,3), (9,2) |

All code indices are zero-based. The extra median constraint is
`y[n//2] >= 0.5`. The first two points have y=0 and the last two have y=1;
two points are selected on each vertical edge. A corner can serve both a
horizontal and a vertical edge, so these constraints do not require eight
distinct boundary points. Bottom ordering is `x[0] <= x[1]` and reflection
uses `x[0] <= 0.5`. Bottom and top triangle signs are fixed to one, with the
top sign convention corresponding to `x[n-2] >= x[n-1]`. Half assignments
are 5/5 for n=10, either 5/6 or 6/5 for n=11, and 6/6 for n=12; a point on
x=0.5 is assigned to exactly one half. No coordinate-pair mirror symmetry
is imposed.

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

