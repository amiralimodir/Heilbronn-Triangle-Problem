# Independent area verification

Python 3.8+ is sufficient; no solver or third-party package is required.
Supply your saved coordinates or the results produced by an experiment:

```console
python 06_validation/verify_areas.py coordinates.csv
python 06_validation/verify_areas.py 04_final_results/results --output validation_results/final
python 06_validation/verify_areas.py path/to/result.json --output validation_results/one_run
```

CSV columns are `n,i,x,y`, with point indices `i=1,...,n`. JSON inputs contain
`x` and `y` coordinate arrays; individual and aggregate experiment results are
accepted. A directory is searched recursively for `result.json` files. If none
exist, its JSON and CSV files are examined. Runs with null coordinates are
skipped. Multiple configurations with the same n are checked separately.

Every finite decimal is treated as an exact rational number. The script checks
the unit square, evaluates every triangle, reports the exact minimum and a
minimizing triple, and rounds feasible lower bounds down to seven decimal
places. Coordinates are never clipped. Points outside the square produce no
feasible lower bound and a nonzero exit code.

Outputs are `validation.json` and `all_triangle_areas.csv`, with a source field
identifying each configuration. These checks verify configurations, not global
optimality or numerical solver dual bounds. No bundled data file is required.
