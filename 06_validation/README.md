# Exact configuration-area verification

No optimizer or third-party package is required; Python 3.8+ suffices.

From the repository root:

```console
python 06_validation/verify_areas.py --output validation_results/appendix
python 06_validation/verify_areas.py data/coordinates/saved_n10_n12.csv --output validation_results/saved_runs
python tools/check_repository.py
```

The CSV columns are `n,i,x,y`. Each finite decimal is interpreted as an exact
rational number. For every configuration the script checks every coordinate
against [0,1], evaluates every triangle, identifies a minimizing triple, and
produces the exact minimum and a downward-rounded seven-decimal lower bound.
Out-of-square coordinates produce no feasible lower bound; they are not clipped.

The checks validate the explicit configurations. They do not establish global
optimality, validate the solver's dual bound, or prove the added restrictions.

