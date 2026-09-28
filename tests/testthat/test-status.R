## ---------------------------------------------------------------------------
## test-status.R -- unsolved DMUs must say they are unsolved.
##
## Three defects, all the same shape: a DMU that could not be scored came back
## as NA while `status` reported 0 ("optimal") and nothing warned, so a failure
## was indistinguishable from an answer through the documented field meant to
## tell them apart. Two were found in the radial LP path; the third, in free
## disposal, was reported by GitHub Copilot on issue #1 and confirmed here.
##
## The first two were found while measuring whether the reference set could be
## thinned of provably redundant DMUs (see ../../../horserace/screen_experiment.R).
## That optimisation is NOT in the package; these bugs, which had nothing to do
## with it, are:
##
##   * `.dea_radial()` never read stage two's solver status, so a failed slack
##     program returned NA slacks while `status` still reported 0;
##   * `dea()` never called `.dea_report_unsolved()` at all -- alone among the
##     entry points -- so nothing warned.
##
##   * `.dea_fdh()` returned `integer(n)` as its status vector -- 0 for every
##     DMU including the ones with no dominating peer.
##
## The first two together hid 56 of 1200 DMUs failing numerically on a
## variable-returns fit. The cause is not tolerances: it is basis state carried
## between DMUs on the reused lpSolveAPI object, and a fresh object fixes every
## one of them.
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

test_that("free disposal reports an undominated DMU as infeasible", {
  ## Copilot, issue #1. Confirmed against the LP paths on identical data: both
  ## `rts = "vrs"` and dea_sbm() returned status 2 and warned, while
  ## `rts = "fdh"` returned NA with status 0 and said nothing.
  ##
  ## Input-oriented FDH needs a reference DMU whose outputs weakly dominate the
  ## evaluated DMU's. Here the second evaluated DMU's outputs exceed every
  ## reference DMU's, so nothing dominates it and the score does not exist.
  ## One input, one output, and enough reference DMUs that the n < 3(p+q)
  ## dimensionality warning does not fire -- the assertions below are about
  ## infeasibility and should not have to filter another warning out.
  xr <- matrix(seq_len(8), ncol = 1)
  yr <- matrix(seq_len(8), ncol = 1)
  x  <- matrix(c(2, 5),  ncol = 1)
  y  <- matrix(c(2, 99), ncol = 1)

  expect_warning(fit <- dea(x, y, rts = "fdh", orientation = "in",
                            xref = xr, yref = yr),
                 "INFEASIBLE")
  expect_true(is.na(fit$eff[2]))
  expect_identical(fit$status, c(0L, 2L))
  expect_true(is.na(fit$sum_lambda[2]))

  ## The same data through the LP path, which has always been right. FDH must
  ## agree with it about WHICH DMU failed and what to call the failure.
  expect_warning(lp <- dea(x, y, rts = "vrs", orientation = "in",
                           xref = xr, yref = yr), "INFEASIBLE")
  expect_identical(fit$status, lp$status)

  ## And the wording is specialised: free disposal needs one dominating DMU,
  ## not a combination of them.
  expect_warning(dea(x, y, rts = "fdh", orientation = "in", xref = xr, yref = yr),
                 "a single reference DMU")
})

test_that("a self-referenced free-disposal fit never reports infeasible", {
  ## The fix must not start warning on ordinary use. It cannot: every DMU
  ## dominates itself, so a self-referenced FDH fit has a peer for every unit
  ## and the new status can only fire against an external reference set.
  s <- dea_sim(200, p = 2, q = 2, returns = 0.8, seed = 5)
  for (ori in c("in", "out")) {
    expect_silent(fit <- dea(s$x, s$y, rts = "fdh", orientation = ori))
    expect_true(all(fit$status == 0L), info = ori)
    expect_false(any(is.na(fit$eff)), info = ori)
  }
})
