# Environments and numerical settings

The authors reported a Dell workstation with two Intel Xeon Gold 6248R CPUs,
48 physical cores / 96 hardware threads in total, and 256 GB RAM. Machine
names are omitted from the public run metadata; CPU and solver information
remain in the logs. Archived operating-system builds differ, so the two
saved n=9 runs must not be treated as a controlled one-parameter experiment.

| Experiment | Python / interface | Solver | Nondefault or explicit settings |
|---|---|---|---|
| Saved enhancement comparison profile | Python 3.8.12 | Gurobi 11.0.1 | TimeLimit=86400, MIPGap=1e-4, MIPGapAbs=1e-10, Threads=0, Seed=0, default feasibility tolerances explicitly recorded, NonConvex=-1 |
| Saved n=9 tuning trial | Python 3.8.12 | Gurobi 11.0.3 | TimeLimit=86400, MIPGap=1e-6, MIPGapAbs=0, Symmetry=2; Threads not explicitly set |
| Selected new final recipe | Python 3.8-compatible source | Gurobi 11.0.3 reference | TimeLimit=86400, MIPGap=1e-6, MIPGapAbs=0, Threads=32, Seed=0, Symmetry=2, ConcurrentMIP=1, feasibility tolerances pinned in JSON |
| Saved n=10…12 excerpt/source pair | Historical runtime version not recovered from excerpts | Gurobi, version not recorded in excerpt | Supplied source sets TimeLimit=86000, MIPGap=1e-6, MIPGapAbs=0, Symmetry=2; actual LB/UB inputs absent |
| Multi-solver comparison | Python 3.8.12, Pyomo 6.8.2 | SCIP 10.1.0, Couenne 0.5.8, BARON 22.4.20 | See the preserved procedural source and section README for solver-specific options |

Gurobi feasibility, integrality, and dual-feasibility defaults in the recorded
profile are 1e-6, 1e-5, and 1e-6, respectively. The complete explicit parameter
JSONs and available parameter files are authoritative for each run. Unlisted
parameters use the defaults of that run's solver version. Threads=0 means
automatic selection; it is not equivalent to a fixed 32-thread setting.

The old matched A1/A3 comparison deliberately retains its historical settings:
input LB=0 is converted to log(n)/n²; A1 requests a tighter stopping gap than
A3; A3 is the MILP relaxation at P=10. Do not overwrite these settings with
the new final profile when reporting the old comparison.

For a new execution, save package versions using:

```console
python tools/report_environment.py --output environment.json
python -m pip freeze > installed_packages.txt
```

The portable scripts have been checked for Python 3.8 syntax. Historical
NumPy/Matplotlib versions and full operating-system/software images were not
saved; no complete bit-for-bit historical environment is claimed.

