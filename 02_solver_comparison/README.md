# Solver comparison

This directory corresponds to **Computational Experiments / Solver Comparison**.
Each script contains the complete procedural Pyomo model. Model construction is
not delegated to a separate helper package.

## Files and evidence

| File | Role |
| --- | --- |
| `solver_comparison_original.py` | Original submitted Pyomo experiment with a descriptive header and help-only argument handling. SCIP and Couenne: n=3,...,7; BARON: n=7. Preserved for provenance, including its original dual-bound extraction limitation. |
| `rerun_missing_duals.py` | Same mathematical model with corrected bound/status recovery and per-run output files. By default, only Couenne n=3,...,7 and SCIP n=7 are repeated. |
| `source_as_received/solver_comparison_pyomo.txt` | Byte-for-byte supplied source. |
| `source_provenance.json` | SHA-256 of the received source. |

The previously supplied JSON did not contain usable finite dual bounds for those
six cases. `Infinity` in that export must not be read as evidence that the solver
could not produce a finite bound: the reader did not reliably recover it.
Corrected rerun results are **pending**; this repository does not manufacture
their values or infer optimality from a good incumbent.

## Recorded software environment

The author reported Python **3.8.12**, Pyomo **6.8.2**, SCIP **10.1.0**,
Couenne **0.5.8**, and BARON **22.4.20**. The Gurobi installation in the
comparison environment was reported as **11.0.1**. New reruns record their
actual environment separately. An installed version is not automatically proof
of the version used for every historical experiment.

## Mathematical model and settings

- Original Approach-1 Pyomo model: x,y,w in [0,1], signed areas in [-0.5,0.5],
  triangle sign binaries, bilinear product equalities, two absolute-value
  disjunction inequalities per triangle.
- Input lower bound **0**, upper bound **0.5**; no automatic logarithmic
  substitution, no geometric enhancement groups.
- No added sign-strengthening inequalities from the enhanced Gurobi scripts.
- Requested limit **86400 seconds**: SCIP `limits/time`, Couenne `time_limit`
  in its `couenne.opt`, BARON `MaxTime`.
- Feasibility, optimality-gap, thread, and other solver settings are left at the
  corresponding solver's defaults. The final-result Gurobi tuning is not applied.

Elapsed solve-call time and the requested solver limit are separate quantities;
startup, shutdown, solver time accounting, and interfaces can produce a recorded
elapsed time beyond the requested limit.

## Run only the outstanding cases

Install Pyomo and the solver executables. BARON requires its own installation and
license; no solver executable or license is distributed here. The original
Windows executable locations are editable constants at the start of
`rerun_missing_duals.py`:

```python
SCIP_EXECUTABLE = r"C:\Program Files\SCIPOptSuite 10.1.0\bin\scip.exe"
COUENNE_EXECUTABLE = r"C:\Program Files\Couenne\ampl.mswin64\couenne.exe"
```

From the repository root:

```console
python 02_solver_comparison/rerun_missing_duals.py
```

Both entry scripts support `--help` before importing Pyomo; this cannot start a
solve. Unknown arguments are rejected. Their existing main guards prevent solves
when the files are imported.

The `RUNS_TO_REPEAT` list selects the cases. Remove entries to run a subset.
`OUTPUT_ROOT` defaults to `Path.home() / "heilbronn_solver_reruns"`; use a writable
ASCII-only directory to avoid the file-path issue encountered on the servers.
The script checkpoints after every case and stops before the next case if bound
recovery or the solve fails. Keep that protection enabled when collecting data.

Each run retains source, model/interface files, raw solver log, options,
environment information, status, available objective, available dual bound, and
elapsed time. Failed extraction is represented as a missing value with a reason,
not a fabricated numerical bound. Coordinates are saved when available but may
violate the square by numerical tolerances; validate them before claiming a
feasible configuration value.

## Bound recovery details

Couenne's internal minimization problem has objective `-z`. Its final global
**lower** bound therefore becomes the maximization dual **upper** bound after
negating it. The parser deliberately does not treat a node's local bound or an
incumbent as the final global bound. The exported value includes an allowance
for the final printed digit and outward floating-point conversion; the directly
printed value and provenance are also retained. This allowance addresses printed
precision, not all numerical error in the optimization algorithm.

For SCIP the final log's dual bound is used when available, with documented
interface fallbacks. Retained logs make the extraction auditable. These are
floating-point solver bounds, not independently validated exact certificates.

The original script is included for provenance rather than as the recommended
exporter. Do not merge its old missing-bound entries over successful corrected
rerun records. BARON n=7 already has a reported finite dual bound and is excluded
from the default missing-bound rerun list.
