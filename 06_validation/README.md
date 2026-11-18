# Independent area verification

Python 3.8+ is sufficient; no solver or third-party package is required.

```console
python 06_validation/verify_areas.py --output validation_results/appendix
python 06_validation/verify_areas.py data/coordinates/saved_n10_n12.csv --output validation_results/saved_runs
```

Input CSV columns are `n,i,x,y`. Every finite decimal is treated as an exact
rational number. The script checks the unit square, evaluates every triangle,
reports the exact minimum and a minimizing triple, and rounds feasible lower
bounds down to seven decimal places. Coordinates are never clipped.

The default input is the n=6…12 coordinate table from manuscript version 15.
The outputs include a JSON summary and all triangle areas. These checks verify
configurations, not global optimality or numerical solver dual bounds.

