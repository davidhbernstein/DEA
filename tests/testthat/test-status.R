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

## ---------------------------------------------------------------------------
## The fourth of the same shape, found while checking complementary slackness:
## the MULTIPLIER sweep had neither the retry nor the reporting. A dual that
## failed returned NA weights while `status` -- which describes the envelopment
## program -- went on saying 0. Only super-efficiency fits were affected;
## ordinary fits show no dual failures even at n = 800.
## ---------------------------------------------------------------------------

test_that("a reused LP object no longer loses multipliers under super-efficiency", {
  ## Four DMUs here returned lpSolve status 2 on the reused object and solved
  ## to optimality on a fresh one. Infeasibility is a contradiction for this
  ## program: the envelopment program is feasible and bounded for every DMU
  ## under crs, so strong duality makes its dual feasible and bounded too.
  d <- dea_sim(40, p = 2, q = 2, returns = 0.9, seed = 5)
  for (ori in c("in", "out")) {
    f <- suppressWarnings(dea(d$x, d$y, rts = "crs", orientation = ori,
                              super = TRUE, slack = FALSE, multipliers = TRUE))
    expect_true(all(f$status == 0L), info = ori)
    expect_true(all(f$mult_status %in% c(0L, 1L)), info = ori)
    expect_false(any(is.na(f$v)), info = ori)
    expect_false(any(is.na(f$u)), info = ori)
  }
})

test_that("every surviving dual failure pairs with an infeasible envelopment program", {
  ## Under vrs some duals still do not solve, and that is the DOCUMENTED
  ## super-efficiency case rather than a fault: the envelopment program is
  ## infeasible, the score is NA, and the dual has no bounded optimum. What
  ## must not happen is a dual failing where the envelopment program solved.
  d <- dea_sim(40, p = 2, q = 2, returns = 0.9, seed = 5)
  for (ori in c("in", "out")) {
    f <- suppressWarnings(dea(d$x, d$y, rts = "vrs", orientation = ori,
                              super = TRUE, slack = FALSE, multipliers = TRUE))
    bad <- !f$mult_status %in% c(0L, 1L)
    expect_true(all(f$status[bad] == 2L), info = ori)
  }
  ## And that documented case must stay silent -- warning about it would be
  ## noise about an answer the user has already been given as NA.
  expect_silent(dea(d$x, d$y, rts = "crs", super = TRUE, slack = FALSE,
                    multipliers = TRUE))
})

test_that("the multiplier reporter warns on a dual failure the duality theorem forbids", {
  ## Driven directly, because a solver inconsistency cannot be produced on
  ## demand from data. status 2 on the dual where the envelopment program
  ## solved (0) contradicts strong duality and must be reported, under
  ## super-efficiency as much as without it.
  rep <- DEA:::.dea_report_mult_unsolved
  expect_warning(rep(c(0L, 2L), c(0L, 0L), 2L, super = FALSE), "NA MULTIPLIERS")
  expect_warning(rep(c(0L, 2L), c(0L, 0L), 2L, super = TRUE),  "NA MULTIPLIERS")
  ## The documented pairing is masked only under super-efficiency, where it is
  ## the expected result; without it, the same pairing is still a fault.
  expect_silent(rep(c(0L, 3L), c(0L, 2L), 2L, super = TRUE))
  expect_warning(rep(c(0L, 3L), c(0L, 2L), 2L, super = FALSE), "NA MULTIPLIERS")
  ## Whatever the solver calls the unbounded dual, the pairing is what counts.
  expect_silent(rep(c(0L, 5L), c(0L, 2L), 2L, super = TRUE))
  expect_silent(rep(c(0L, 1L), c(0L, 0L), 2L, super = FALSE))
})

## ---------------------------------------------------------------------------
## The fifth of the same shape: dea_ddf() recorded `status` and reported
## nothing, so a failed program returned beta = NA in silence. Found by
## auditing every solve site for the three things the earlier four were missing
## -- status recorded, retried on a fresh program, and reported.
## ---------------------------------------------------------------------------

test_that("dea_ddf() reports an unsolved DMU instead of returning NA quietly", {
  ## direction = "in" holds g_y at zero, so outputs cannot expand: a DMU whose
  ## output exceeds everything the reference set can produce has no feasible
  ## point. Reachable from ordinary code, not a contrived program.
  r <- dea_sim(30, p = 1, q = 1, seed = 2)
  xe <- matrix(c(1.5, 1.5), ncol = 1)
  ye <- matrix(c(0.5, max(r$y) * 3), ncol = 1)

  expect_warning(f <- dea_ddf(xe, ye, direction = "in", rts = "vrs",
                              xref = r$x, yref = r$y),
                 "INFEASIBLE")
  expect_identical(f$status, c(0L, 2L))
  expect_true(is.na(f$beta[2]))

  ## The contrast that made it a bug rather than a choice: dea() warns on the
  ## same data, and beta = 1 - theta ties the two answers together.
  expect_warning(g <- dea(xe, ye, rts = "vrs", orientation = "in",
                          xref = r$x, yref = r$y), "INFEASIBLE")
  expect_identical(g$status, f$status)
  expect_equal(unname(f$beta[1]), 1 - unname(g$eff[1]), tolerance = 1e-9)
})

## ---------------------------------------------------------------------------
## The retry guard itself. Every sweep reuses one linear program and rewrites
## only the evaluated DMU's part of it, and lpSolveAPI carries basis state on
## that object; .lp_solve_retry() rebuilds once when a solve gives up for a
## reason that is not an answer about the data.
##
## The logic is tested directly rather than only through a sweep, because the
## branch is rare by construction: a guard whose first real execution is the
## day it matters has never been run.
## ---------------------------------------------------------------------------

## Three one-variable programs with known, deterministic statuses: 0, 2 and 3.
tiny_lp <- function(kind) {
  lp <- lpSolveAPI::make.lp(1L, 1L)
  lpSolveAPI::lp.control(lp, sense = if (kind == "unbounded") "max" else "min",
                         verbose = "neutral")
  lpSolveAPI::set.column(lp, 1L, 1, 1L)
  lpSolveAPI::set.objfn(lp, 1)
  lpSolveAPI::set.constr.type(lp, ">=", 1L)
  if (kind == "unbounded") {
    lpSolveAPI::set.rhs(lp, 1, 1L)
    lpSolveAPI::set.bounds(lp, lower = 0, upper = DEA:::.LP_INF, columns = 1L)
  } else if (kind == "infeasible") {
    lpSolveAPI::set.rhs(lp, 5, 1L)
    lpSolveAPI::set.bounds(lp, lower = 0, upper = 1, columns = 1L)
  } else {
    lpSolveAPI::set.rhs(lp, 2, 1L)
    lpSolveAPI::set.bounds(lp, lower = 0, upper = 10, columns = 1L)
  }
  lp
}

test_that("the retry rebuilds only when the status is not an answer", {
  retry <- DEA:::.lp_solve_retry
  ## Confirm the fixtures really do produce the statuses the test relies on.
  expect_identical(solve(tiny_lp("ok")),         0L)
  expect_identical(solve(tiny_lp("infeasible")), 2L)
  expect_identical(solve(tiny_lp("unbounded")),  3L)

  ## A solved program is not rebuilt.
  calls <- 0L
  fresh <- function() { calls <<- calls + 1L; tiny_lp("ok") }
  r <- retry(tiny_lp("ok"), fresh)
  expect_identical(r$status, 0L)
  expect_identical(calls, 0L)

  ## Nor is an INFEASIBLE one, by default: that is a statement about the data,
  ## and re-solving cannot change it. This is the case the guard must not
  ## touch -- under super-efficiency it is the expected result.
  calls <- 0L
  r <- retry(tiny_lp("infeasible"), fresh)
  expect_identical(r$status, 2L)
  expect_identical(calls, 0L)

  ## Anything else is rebuilt once, and the CALLER GETS THE FRESH PROGRAM --
  ## reading results off the original would report the failed solve's state.
  calls <- 0L
  r <- retry(tiny_lp("unbounded"), fresh)
  expect_identical(r$status, 0L)          ## the rebuilt one solves
  expect_identical(calls, 1L)
  expect_equal(lpSolveAPI::get.objective(r$lp), 2)

  ## `final` is what lets the multiplier program invert the rule: there an
  ## infeasible dual contradicts strong duality and IS retried, while an
  ## unbounded dual is the answer and is not.
  calls <- 0L
  r <- retry(tiny_lp("infeasible"), fresh, final = c(0L, 1L, 3L))
  expect_identical(calls, 1L)
  expect_identical(r$status, 0L)
  calls <- 0L
  r <- retry(tiny_lp("unbounded"), fresh, final = c(0L, 1L, 3L))
  expect_identical(calls, 0L)
  expect_identical(r$status, 3L)
})

test_that("the oriented slacks-based measure recovers a DMU the reuse loses", {
  ## n = 150, seed 12: DMU 108 gives lpSolve status 5 on the reused program.
  ## Verified as recoverable AND correct -- solved independently under four
  ## lpSolve scaling modes (none/geometric/curtisreid/extreme) rho is exactly
  ## 1 in all four, and the package now reports exactly 1.
  s <- dea_sim(150, p = 3, q = 2, returns = 0.85, seed = 12)

  ## First, the failure is real: replicate the pre-fix sweep, which reuses one
  ## program and never rebuilds. If this stops failing the regression test
  ## below has quietly stopped testing anything.
  sc <- DEA:::.dea_scale(s$x, s$y, TRUE)
  X <- sc$X; Y <- sc$Y; n <- 150L; p <- 3L; q <- 2L
  lp <- DEA:::.lp_slack_build(X, Y, "vrs")$lp
  st <- integer(n)
  for (o in seq_len(n)) {
    xo <- pmax(X[o, ], DEA:::.DEA_CONSTANTS$MIN_POS)
    ob <- numeric(n + p + q); ob[n + seq_len(p)] <- 1/(p * xo)
    lpSolveAPI::set.objfn(lp, ob)
    lpSolveAPI::set.rhs(lp, c(X[o, ], Y[o, ]), seq_len(p + q))
    st[o] <- solve(lp)
  }
  expect_identical(which(!st %in% c(0L, 1L)), 108L)

  ## With the guard the DMU is solved, not lost, and silently so.
  expect_silent(f <- dea_sbm(s$x, s$y, rts = "vrs", orientation = "in"))
  expect_true(all(f$status %in% c(0L, 1L)))
  expect_false(any(is.na(f$eff)))
  expect_equal(unname(f$eff[108]), 1, tolerance = 1e-9)
})
