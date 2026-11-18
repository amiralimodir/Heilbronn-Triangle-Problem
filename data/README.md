# Available evidence

`coordinates/appendix_v15.csv` contains exactly the printed decimal coordinates
from manuscript version 15 for n=6…12. `appendix_expected.json` records the
manuscript's stated exact minima. The verifier checks these independently.
`figure_v15.csv` preserves the different, longer decimals used in its drawings;
`saved_n10_n12.csv` comes from the saved run excerpts. A figure or rounded
appendix configuration is not silently substituted for a solver solution.
In version 15, the n=10 figure still contains the older run's coordinates;
the updated appendix and `saved_n10_n12.csv` contain the newer run. The figure
CSV deliberately documents the manuscript as received rather than concealing
this remaining manuscript inconsistency.

`recorded_runs/group_n9_G123_default/` is the saved Gurobi 11.0.1 default-profile
comparison run. `recorded_runs/tuning_n9_G123/` is the distinct Gurobi 11.0.3
tuning run. Both include source, model, solution, inputs, parameters, log and
status/bound metadata. Neither is asserted to be the paper's main n=9 row.

`recorded_runs/paper_n10/`, `paper_n11/`, and `paper_n12/` preserve console
excerpts and structured summaries. Their elapsed times are 5495.305539369583,
86001.00913786888, and 86000.80664873123 seconds. The supplied historical source
sets a time limit of 86000 seconds; these outputs are not results of the new
86400-second final profile. Unknown LB/UB inputs and unrecorded versions remain
null rather than being filled from a different experiment.

`recorded_runs/solver_comparison_available.json` normalizes the latest supplied
multi-solver data. Unrecovered infinite dual values are null with an explanatory
status. This is incomplete evidence; see the targeted rerun instructions.
`reported_auxiliary/proposition5_case_b.json` transcribes the author's reported
console result and input clarification, not an independently recovered log.

Original and published hashes for the two complete archives are in
`recorded_runs/archive_integrity.json`. Only local machine/path fields and,
where present, personal academic-license identifiers are omitted. Model,
parameter, coordinate, and numerical solver information is not changed.

