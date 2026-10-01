## ---------------------------------------------------------------------------
## Non-discretionary inputs and outputs.
##
## The LP change is small and in one place, so what is checked is that the
## program really has the constraints it claims -- read off the solved lambda
## rather than off the code -- plus the ordering theorem, which is the thing
## most likely to be stated backwards.
## ---------------------------------------------------------------------------

nd_fixture <- function(n = 60, seed = 31) {
  set.seed(seed)
  teachers <- runif(n, 10, 40)
  pupils <- runif(n, 150, 600)
  x <- cbind(teachers = teachers, pupils = pupils)
  y <- matrix((teachers^0.3 * pupils^0.6) * exp(-abs(rnorm(n, 0, 0.2))), n, 1,
              dimnames = list(NULL, "attainment"))
  list(x = x, y = y)
}

test_that("the fixed row is not scaled and the free row is", {
  ## Read the constraints off the solved program, which is the only way to
  ## check that the theta column and the right-hand side were BOTH changed --
  ## editing one without the other would impose sum(lambda x) <= 0.
  f <- nd_fixture()
  r <- dea_nd(f$x, f$y, nd_x = "pupils")
  peer <- r$lambda %*% f$x
  expect_lt(max(peer[, "pupils"] - f$x[, "pupils"]), 1e-7)        # <= x_io
  expect_lt(max(peer[, "teachers"] - r$eff * f$x[, "teachers"]), 1e-7)  # <= theta x_io
  ## and the fixed row is NOT satisfied at the tighter, scaled level for all --
  ## if it were, the mask would be doing nothing.
  expect_gt(max(peer[, "pupils"] - r$eff * f$x[, "pupils"]), 1e-3)
  ## outputs are still met
  expect_gt(min((r$lambda %*% f$y)[, 1] - f$y[, 1]), -1e-7)
})

test_that("fixing a variable can only lower theta, never raise it", {
  ## The ordering theorem, and the thing most likely to be remembered
  ## backwards: for theta < 1 the fixed row is LOOSER than the scaled one, so
  ## the feasible set grows and the minimum falls.
  for (sd in 1:6) {
    set.seed(sd)
    n <- 50; p <- 3; q <- 2
    X <- matrix(runif(n * p, 2, 20), n, p)
    Y <- matrix(runif(n * q, 1, 10), n, q)
    for (rts in c("vrs", "crs", "nirs", "ndrs")) {
      a <- suppressWarnings(dea(X, Y, rts = rts, slack = FALSE))
      b <- suppressWarnings(dea_nd(X, Y, nd_x = c(FALSE, TRUE, TRUE),
                                   rts = rts, slack = FALSE))
      expect_lt(max(b$eff - a$eff), 1e-8)
      ## output orientation is the mirror
      ao <- suppressWarnings(dea(X, Y, rts = rts, orientation = "out", slack = FALSE))
      bo <- suppressWarnings(dea_nd(X, Y, nd_y = c(TRUE, FALSE), rts = rts,
                                    orientation = "out", slack = FALSE))
      expect_lt(max(ao$eff - bo$eff), 1e-8)
    }
  }
})

test_that("with nothing marked it is exactly dea()", {
  f <- nd_fixture()
  for (ori in c("in", "out")) {
    a <- dea(f$x, f$y, orientation = ori)
    b <- suppressWarnings(dea_nd(f$x, f$y, orientation = ori))
    expect_equal(unname(b$eff), unname(a$eff))
    expect_equal(unname(b$slack_x), unname(a$slack_x))
    expect_equal(unname(b$slack_y), unname(a$slack_y))
    expect_equal(unname(b$lambda), unname(a$lambda))
  }
  expect_warning(dea_nd(f$x, f$y), "ordinary dea\\(\\) fit")
})

test_that("a fixed slack is kept in the constraints and dropped from the objective", {
  ## Banker and Morey's own rule, and it applies in BOTH orientations because
  ## the orientation decides which side stage one scales, not which slacks
  ## count at stage two.
  f <- nd_fixture()
  free <- dea(f$x, f$y)
  fix <- dea_nd(f$x, f$y, nd_x = "pupils")
  ## The program still balances: the slack exists.
  peer <- fix$lambda %*% f$x
  expect_equal(unname(peer[, "pupils"] + fix$slack_x[, "pupils"]),
               unname(f$x[, "pupils"]), tolerance = 1e-6)
  ## Pareto-Koopmans is judged over the discretionary coordinates only, so a
  ## DMU with residual slack on `pupils` alone is still efficient.
  only_fixed <- fix$eff == 1 & fix$slack_x[, "pupils"] > 1e-6 &
    fix$slack_x[, "teachers"] <= 1e-8
  if (any(only_fixed)) expect_true(all(fix$efficient[only_fixed]))
  expect_false(identical(unname(free$slack_x), unname(fix$slack_x)))
})

test_that("the impossible configurations are refused", {
  f <- nd_fixture(40)
  expect_error(dea_nd(f$x, f$y, nd_x = c(TRUE, TRUE)), "nothing to minimise")
  expect_error(dea_nd(f$x, f$y, nd_y = TRUE, orientation = "out"),
               "nothing to maximise")
  expect_error(dea_nd(f$x, f$y, nd_x = "pupil"), "no such column")
  expect_error(dea_nd(f$x, f$y, nd_x = c(TRUE, FALSE, TRUE)), "one entry per input")
  expect_error(dea_nd(f$x, f$y, nd_x = 3), "whole numbers in 1:2")
  expect_error(dea_nd(f$x, f$y, nd_x = list(1)), "logical mask")
})

test_that("the three ways of naming a variable agree, and nd_y is announced", {
  f <- nd_fixture(40)
  a <- dea_nd(f$x, f$y, nd_x = "pupils")
  b <- dea_nd(f$x, f$y, nd_x = 2L)
  c3 <- dea_nd(f$x, f$y, nd_x = c(FALSE, TRUE))
  expect_equal(a$eff, b$eff); expect_equal(a$eff, c3$eff)
  expect_equal(a$nd_x, c(teachers = FALSE, pupils = TRUE))

  ## Marking an output non-discretionary is right much less often than marking
  ## an input, so it says so every time rather than once in the help.
  expect_message(dea_nd(f$x, f$y, nd_y = TRUE), "genuinely exogenous")
  expect_output(print(a), "Non-discretionary")
  expect_output(print(a), "pupils")
  expect_s3_class(a, "dea")
})

test_that("it composes with an external reference set", {
  f <- nd_fixture(60)
  g <- nd_fixture(15, seed = 77)
  ## Scoring against someone else's technology can be infeasible, as for every
  ## other estimator here, and it is reported rather than hidden.
  r <- suppressWarnings(dea_nd(g$x, g$y, nd_x = "pupils", xref = f$x, yref = f$y))
  expect_equal(r$n, 15L)
  expect_equal(r$nref, 60L)
  ok <- is.finite(r$eff)
  expect_true(any(ok))
  peer <- r$lambda[ok, , drop = FALSE] %*% f$x
  expect_lt(max(peer[, "pupils"] - g$x[ok, "pupils"]), 1e-7)
})
