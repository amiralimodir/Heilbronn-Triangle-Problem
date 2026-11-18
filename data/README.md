# Results and coordinates

- `coordinates/appendix_v15.csv`: the printed coordinates for n=6…12.
- `coordinates/appendix_expected.json`: the exact minimum areas stated in the appendix.
- `coordinates/saved_n10_n12.csv`: full decimal coordinates from the saved runs.
- `coordinates/verified_appendix/`: independently computed exact triangle areas.
- `recorded_runs/group_n9_G123_default/`: the available n=9 enhancement-comparison run, including its model, solution, parameters, and log.
- `recorded_runs/tuning_n9_G123/`: the separate n=9 run with tighter solver settings.
- `recorded_runs/paper_n10/`, `paper_n11/`, `paper_n12/`: saved output and result summaries.
- `recorded_runs/formulation_excerpts/`: available formulation-comparison outputs.
- `recorded_runs/solver_comparison_available.json`: available solver results; missing dual bounds are null.
- `reported_auxiliary/proposition5_case_b.json`: the reported Case-B result and confirmed inputs.

Each run retains its own settings and elapsed time. The saved n=10…12 runs use
a 86000-second limit; the current final-run configuration uses 86400 seconds.
Missing input values or solver versions remain null. Neither saved n=9 directory
is labeled as the paper's 633.331-second run.

Validate the coordinates without a solver:

```console
python 06_validation/verify_areas.py
```
