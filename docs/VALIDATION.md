# Packaging checks

The repository preparation performed the following checks without launching
optimization experiments:

- The final model/helper function syntax trees and order match all ten
  prepared n=3…12 source files; all instance-specific flags match.
- The enhancement-comparison model/helper function syntax trees match the
  archived source; all 48 n/group flag combinations match the earlier manifest.
- Historical source copies are retained with their hashes; the current general
  formulation selector correction is documented separately.
- The Proposition 5 constraint construction and original helper functions are
  preserved, with the author's effective-LB=0 clarification explicitly recorded.
- Active scripts parse using Python 3.8 grammar. Entry-point/help checks were
  performed without starting a solver job; historical files are archival inputs.
- Every one of the 700 triangles for the seven appendix configurations was
  evaluated exactly and agrees with the manuscript's printed-coordinate minima.
- Available archive hashes and the distinction between comparison and final
  parameter profiles are checked by `tools/check_repository.py`.

These are source-preservation, data-integrity, and license-free execution checks.
They are not fresh optimization runs, solver performance benchmarks, or a
validation of floating-point dual bounds. Numerical runtimes remain the saved
experimental observations identified in the archive.
