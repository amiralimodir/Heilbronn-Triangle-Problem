# Formulations

Run the complete interactive model file:

```console
python 01_formulations/formulations.py
```

Enter the point count, upper bound, lower bound, and approach number. Approach 3
also asks for the expansion depth `H` (called `P` in the paper) and `MILP` or
`MIQCP`. The selection calls the corresponding model function.

For the n=7,8 Approach-1/Approach-3 comparison:

| Input | Value |
|---|---|
| n | 7 or 8 |
| UB | 0.5 |
| LB entered | 0 |
| Approach | 1 or 3 |
| Approach-3 depth | 10 |
| Approach-3 model | MILP |

In this file, entering LB=0 changes it to `log(n)/n**2`. Approach 1 uses
TimeLimit=86400, MIPGap=1e-6, and MIPGapAbs=0. Approach 3 MILP uses the same
time limit with Gurobi's default gap settings. These settings are part of the
comparison and should not be replaced by the final-run profile.

The Approach-2 and exact-MIQCP functions have no explicit time limit. The
horizontal-strip function in this comparison does not force point assignments;
use the dedicated enhancement-comparison script for the complete G4 model.

Console output reports the computed solution, and the figure is written in the
working directory. Run each case from a separate directory to retain its files.
An Approach-3 MILP incumbent is a relaxation value; recompute the actual minimum
triangle area from its coordinates before interpreting it as a feasible bound.

