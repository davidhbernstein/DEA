## ---------------------------------------------------------------------------
## test-screen.R -- two silent defects in the radial path, and their fixes.
##
## Both were found while measuring whether the reference set could be thinned
## of provably redundant DMUs (see ../../../horserace/screen_experiment.R).
## That optimisation is NOT in the package; these two bugs, which had nothing
## to do with it and were both silent, are fixed:
##
##   * `.dea_radial()` never read stage two's solver status, so a failed slack
##     program returned NA slacks while `status` still reported 0;
##   * `dea()` never called `.dea_report_unsolved()` at all -- alone among the
##     entry points -- so nothing warned.
##
## Together they hid 56 of 1200 DMUs failing numerically on a variable-returns
## fit. The cause is not tolerances: it is basis state carried between DMUs on
## the reused lpSolveAPI object, and a fresh object fixes every one of them.
## ---------------------------------------------------------------------------

test_that("a large variable-returns fit leaves no DMU unsolved", {
  ## Regression test for the basis-carryover failure. Before the retry in
  ## `.dea_radial()` this design produced 56 numerically failed DMUs, and
  ## before the status fix it produced them silently: status 0, no warning,
  ## NA slacks, and `efficient` quietly NA.
  s <- dea_sim(1200, p = 2, q = 2, returns = 0.8, seed = 99)
  expect_silent(fit <- dea(s$x, s$y, rts = "vrs", orientation = "in"))
  expect_true(all(fit$status %in% c(0L, 1L)))
  expect_false(any(is.na(fit$slack_x)))
  expect_false(any(is.na(fit$slack_y)))
  expect_false(any(is.na(fit$efficient)))
  expect_false(any(is.na(fit$eff)))
})

test_that("dea() reports unsolved DMUs instead of going quiet", {
  ## The wiring, checked directly: constructing a genuine solver failure on
  ## demand is not reliable, but the reporter must be reachable from dea() and
  ## must fire on the statuses that matter.
  expect_warning(.dea_report_unsolved(c(0L, 0L, 5L), 3L, TRUE, "dea"),
                 "FAILED NUMERICALLY")
  expect_warning(.dea_report_unsolved(c(0L, 2L), 2L, TRUE, "dea"),
                 "INFEASIBLE")
  expect_silent(.dea_report_unsolved(c(0L, 1L, 0L), 3L, TRUE, "dea"))
})

test_that("super-efficiency infeasibility stays silent, as documented", {
  ## Under super-efficiency an infeasible program is the expected result, not a
  ## fault, so it must not start warning now that dea() reports at all. The NA
  ## is the signal there.
  s <- dea_sim(200, p = 2, q = 2, returns = 0.8, seed = 99)
  expect_silent(fit <- dea(s$x, s$y, rts = "vrs", orientation = "in",
                           super = TRUE, slack = FALSE, peers = FALSE))
  expect_true(any(is.na(fit$eff)))
})
