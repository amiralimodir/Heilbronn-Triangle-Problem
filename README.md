# Heilbronn Triangle Problem

Code and available computational evidence for **Solving the Heilbronn Triangle
Problem using Global Optimization Methods**, by Amirali Modir, Amir Monji,
and Burak Kocuk.

This repository follows the paper's formulations, solver comparison,
enhancement comparison, final experiments, auxiliary computations, and
coordinate appendix. Each optimization entry point contains its complete
procedural model and helper functions. Model construction has not been split
into a shared library or rewritten as a new framework.

## Start here

| Paper component | Directory | Entry point |
|---|---|---|
| Approaches 1, 2, and 3; matched formulation comparison | [01_formulations](01_formulations/) | `formulations.py`; archived comparison source alongside it |
| BARON, SCIP, and Couenne comparison | [02_solver_comparison](02_solver_comparison/) | `solver_comparison_original.py`, `rerun_missing_duals.py` |
| Four enhancement groups, 16 combinations per n | [03_enhancement_comparison](03_enhancement_comparison/) | `enhancement_comparison.py` |
| Selected G1–G3 final-run recipe, n=3…12 | [04_final_results](04_final_results/) | `final_results.py` |
| Proposition 5 and auxiliary-source status | [05_auxiliary](05_auxiliary/) | `proposition5_case_b.py` |
| Independently check every triangle area | [06_validation](06_validation/) | `verify_areas.py` |
| Available logs, models, inputs, solutions, coordinates | [data](data/) | Recorded evidence with provenance |

**First, verify the appendix without installing any solver:**

```console
python 06_validation/verify_areas.py --output validation_results/appendix
```

Python 3.8 or later is sufficient for this check. It interprets finite decimal
coordinates as exact rational numbers, checks the unit square, evaluates all
triangles, and rounds valid configuration lower bounds downward. It does not
verify solver dual bounds or prove unrestricted optimality.

## Optimization environments

The recorded comparison environment uses Python 3.8.12, Gurobi 11.0.1 and,
for the Pyomo comparison, Pyomo 6.8.2. The selected final recipe references
Gurobi 11.0.3. Create separate environments so changing the final-run solver
does not silently change the comparison environment:

```console
python -m pip install -r requirements-comparison.txt
```

or, in the final-run environment:

```console
python -m pip install -r requirements-final.txt
```

SCIP 10.1.0, Couenne 0.5.8, and BARON 22.4.20 are external executables in the
reported Pyomo experiments. Configure their paths as explained in
[02_solver_comparison/README.md](02_solver_comparison/README.md).
Solver licenses and binaries are not distributed here.
See [docs/ENVIRONMENTS.md](docs/ENVIRONMENTS.md) for what is known and what was
not recorded. These requirements are installation profiles, not recovered
complete historical environment locks.

## Keep the experimental protocols distinct

**Enhancement comparison:**

```console
python 03_enhancement_comparison/enhancement_comparison.py --n 9 --groups 123
```

This uses the saved comparison parameters, including a 86400-second limit,
MIPGap=1e-4 and MIPGapAbs=1e-10. The recorded input file is included. Groups
without objective-bound tightening still use objective bounds [0,0.5], as in
the original design. See its README before running the complete 48-case design.

**Selected final recipe:**

Set the desired lower bound for each n in
`04_final_results/bounds_by_n.json`; `null` deliberately requires an explicit
input. Then use:

```console
python 04_final_results/final_results.py --n 9 --check
python 04_final_results/final_results.py --n 9
```

This uses G1–G3, disables G4, and applies the specified larger-n restrictions
for n=10,11,12. Its reference profile uses a 86400-second limit, MIPGap=1e-6,
MIPGapAbs=0, Symmetry=2, Threads=32 and Seed=0. It is a selected recipe for
new runs, **not** a claim that all existing paper rows were produced with
these settings. In particular, the saved n=10…12 excerpts use an older
86000-second implementation, retained in `04_final_results/historical/`.

Run comparisons sequentially on the same intended hardware and environment.
Record the actual configuration and solver output for each run. No claim is
made that a different machine, solver version, thread count, or bound input
reproduces a particular wall-clock time.

## Evidence and remaining material

The repository includes the manuscript-v15 printed coordinates for n=6…12,
their exact independent evaluations, saved n=10…12 console excerpts, two
distinct n=9 model/log/configuration archives, available formulation excerpts,
and the partial multi-solver results. The two n=9 archives are a comparison
run and a tuning run; neither is relabeled as the paper's 633.331-second row.

The full enhancement-result collection and corrected missing Couenne/SCIP
dual-bound results are still being completed. The author will provide the
final max-random and y-bound auxiliary source files separately. Missing
historical inputs/logs are identified explicitly in
[docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

Numerical dual bounds concern the implemented models and solver tolerances.
The n=10…12 additional restrictions are unproved; their restricted-model
dual bounds are not unrestricted upper-bound certificates. Explicit
configurations in the square can still supply feasible lower bounds.

The pre-revision repository is preserved unchanged under
[legacy/pre_revision](legacy/pre_revision/), with its original commit recorded
there. It is historical material, not the default reproduction workflow.

