# Heilbronn Triangle Problem

Code for **Solving the Heilbronn Triangle Problem using Global Optimization
Methods**, by Amirali Modir, Amir Monji, and Burak Kocuk.

Each section has one main Python file containing its complete procedural model.

| Experiment | Main file |
|---|---|
| Approaches 1, 2, and 3 | [01_formulations/formulations.py](01_formulations/formulations.py) |
| Solver comparison | [02_solver_comparison/solver_comparison.py](02_solver_comparison/solver_comparison.py) |
| Enhancement-group comparison | [03_enhancement_comparison/enhancement_comparison.py](03_enhancement_comparison/enhancement_comparison.py) |
| Final runs, n=3…12 | [04_final_results/final_results.py](04_final_results/final_results.py) |
| Proposition 5, Case B | [05_auxiliary/proposition5_case_b.py](05_auxiliary/proposition5_case_b.py) |
| Independent area verification | [06_validation/verify_areas.py](06_validation/verify_areas.py) |

## Setup

Use Python 3.8 or later. Install the appropriate dependencies in separate environments:

```console
python -m pip install -r requirements-comparison.txt
```

```console
python -m pip install -r requirements-final.txt
```

The comparison uses Gurobi 11.0.1; the selected final-run profile uses Gurobi
11.0.3. The Pyomo comparison uses Pyomo 6.8.2, SCIP 10.1.0, Couenne 0.5.8,
and BARON 22.4.20. External solver executables and licenses must be installed
separately. Each section's README explains its inputs and outputs.

## Main commands

```console
python 03_enhancement_comparison/enhancement_comparison.py --n 9 --groups 123
python 04_final_results/final_results.py --n 9 --check
python 04_final_results/final_results.py --n 9
python 05_auxiliary/proposition5_case_b.py
python 06_validation/verify_areas.py --output validation_results/appendix
```

Before a final run, enter its feasible lower bound in
[04_final_results/bounds_by_n.json](04_final_results/bounds_by_n.json).
Null entries require an explicit input. The comparison's recorded bounds and
parameters are kept in its own directory; do not replace them with final-run
settings.

Final runs enable G1–G3 and disable G4. For n=10,11,12, the additional unproved
restrictions are enabled; numerical dual bounds then concern those restricted
models. Coordinate verification supplies feasible lower bounds independently.

Available coordinates and run outputs are under [data](data/).
The remaining enhancement results, corrected Couenne/SCIP results, and final
max-random/y-bound source files are listed in
[docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

