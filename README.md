# Heilbronn Triangle Problem

Code for **Solving the Heilbronn Triangle Problem using Global Optimization
Methods**, by Amirali Modir, Amir Monji, and Burak Kocuk.

Each section has one main Python file containing its complete procedural model.

| Experiment | Main file |
|---|---|
| Formulation tables: 15 unenhanced + 4 shared-enhancement cases | [01_formulations/formulations.py](01_formulations/formulations.py) |
| Solver table: 4 solvers × n=3…7, 20 cases | [02_solver_comparison/solver_comparison.py](02_solver_comparison/solver_comparison.py) |
| Enhancement table: 16 combinations × n=7,8,9, 48 cases | [03_enhancement_comparison/enhancement_comparison.py](03_enhancement_comparison/enhancement_comparison.py) |
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
python 01_formulations/formulations.py --list
python 02_solver_comparison/solver_comparison.py --list
python 03_enhancement_comparison/enhancement_comparison.py --list
python 04_final_results/final_results.py --list
```

Each of these four scripts runs its complete experiment set when called
without arguments. Use `--help` for filters, or `--list` to inspect the selected
cases without loading a solver. All batches run sequentially and save separate
results for each case. Examples:

```console
python 02_solver_comparison/solver_comparison.py --solvers gurobi baron --n 3 4 5 6 7
python 03_enhancement_comparison/enhancement_comparison.py --n 9 --groups 123
python 04_final_results/final_results.py --n 9 --check
python 04_final_results/final_results.py --n 9
python 05_auxiliary/proposition5_case_b.py --check
python 05_auxiliary/proposition5_case_b.py
```

The selected final-run LB/UB inputs are in
[04_final_results/bounds_by_n.json](04_final_results/bounds_by_n.json).
The comparison's recorded bounds and
parameters are kept in its own directory; do not replace them with final-run
settings.

Final runs enable G1–G3 and disable G4. For n=10,11,12, the additional unproved
restrictions are enabled; numerical dual bounds then concern those restricted
models. Coordinate verification supplies feasible lower bounds independently.

Results are generated when experiments run; no bundled `data/` directory is
required. The area verifier accepts saved coordinates or run outputs; see its
[instructions](06_validation/README.md).

The auxiliary max-random and y-bound source files are still awaiting the
authors' final versions. See [reproducibility notes](docs/REPRODUCIBILITY.md)
for the distinction between executable experiments and completed numerical results.

