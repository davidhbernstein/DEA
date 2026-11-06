## Update submission

This is an update of `DEA` from 1.0.0, on CRAN since September 2026, to 1.1.0.

It adds fifteen exported functions, listed in `NEWS.md`. The main ones are
inference (a test of returns to scale, m out of n subsampling intervals, a
two-stage efficiency regression), panel data with the Malmquist index, weight
restrictions, and models for environmental, categorical and undesirable
variables. It also hardens the linear-programming layer against numerical
failures that previously left estimates silently missing.

**One behaviour change:** `dea_rts_test()`, which is new since 1.0.0, now uses
subsampling by default rather than the bootstrap. It did not exist in the CRAN
version, so no CRAN user's code changes. The reason is recorded in `NEWS.md`.

The authorship is unchanged from 1.0.0. Zuleyka Diaz-Martinez and Jose
Fernandez-Menendez, the authors of the archived `DEA` 0.1-2, are authors of
this package. No code is shared with 0.1-2, and there is still no compiled code
and no third-party copyright to declare.

## Test environments

* local macOS (Apple M4), R 4.5.2: `R CMD check --as-cran`
* local Windows 11, R 4.6.1: `R CMD check --as-cran` (without the vignette:
  no pandoc on that machine)
* GitHub Actions, five-platform matrix (`.github/workflows/R-CMD-check.yaml`)
* the declared floor, R >= 4.0.0, by a token-level scan of every R file
  (`tools/check_r_version_floor.R`, not shipped)

win-builder (R-devel) has not been run for this version yet.

## R CMD check results

0 errors | 0 warnings | 2 notes, both from the local machine rather than the
package:

* "checking for future file timestamps ... unable to verify current time":
  the time server could not be reached.
* "checking HTML version of manual": the local HTML Tidy is too old to run.

## Reverse dependencies

None: no CRAN package depends on, imports or suggests `DEA`.
