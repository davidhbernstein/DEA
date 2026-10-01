## ---------------------------------------------------------------------------
## The multiplicative / piecewise Cobb-Douglas model.
##
## No CRAN package implements this, so there is no reference to pin against.
## What is checked instead is of the first and third kinds in CLAUDE.md: the
## two properties that define the model (units invariance, and a geometric
## rather than arithmetic peer target) and recovery of a truth -- on a frontier
## that IS Cobb-Douglas, the envelope is a hyperplane in log space and the model
## should be able to sit on it exactly.
## ---------------------------------------------------------------------------

mult_fixture <- function(n = 40, seed = 3) {
  set.seed(seed)
  x <- matrix(runif(n * 2, 2, 10), n, 2,
              dimnames = list(NULL, c("labour", "capital")))
  y <- matrix((x[, 1] * x[, 2])^0.4 * exp(-abs(rnorm(n, 0, 0.25))), n, 1,
              dimnames = list(NULL, "output"))
  list(x = x, y = y)
}

test_that("the model is units invariant under vrs and not otherwise", {
  ## This is what the 1983 paper is for, and it is the one property that makes
  ## the score mean anything. Multiplying a column by c translates it in log
  ## space, so this IS the additive model's translation invariance.
  f <- mult_fixture()
  base <- dea_mult(f$x, f$y)$eff
  for (cc in c(10, 1e4)) {
    x2 <- f$x; x2[, 1] <- x2[, 1] * cc
    expect_equal(dea_mult(x2, f$y)$eff, base, tolerance = 1e-9)
  }
  expect_equal(dea_mult(f$x, f$y * 7)$eff, base, tolerance = 1e-9)

  ## Without the convexity row it moves, and by a lot -- the log-space cone is
  ## anchored at x = y = 1, which is not a point the data knows about.
  for (rts in c("crs", "nirs")) {
    b <- suppressWarnings(dea_mult(f$x, f$y, rts = rts)$eff)
    x2 <- f$x; x2[, 1] <- x2[, 1] * 1000
    a <- suppressWarnings(dea_mult(x2, f$y, rts = rts)$eff)
    expect_gt(max(abs(a - b)), 0.05)
  }
})

test_that("the peer target is a weighted geometric mean, not an arithmetic one", {
  f <- mult_fixture()
  m <- dea_mult(f$x, f$y)
  ## target = prod_j x_ij ^ lambda_j, which is what distinguishes this model.
  expect_equal(unname(f$x * m$ratio_x), unname(exp(m$lambda %*% log(f$x))),
               tolerance = 1e-8)
  expect_equal(unname(f$y * m$ratio_y), unname(exp(m$lambda %*% log(f$y))),
               tolerance = 1e-8)
  ## And it is NOT the arithmetic mean -- otherwise the model adds nothing.
  expect_gt(max(abs(exp(m$lambda %*% log(f$x)) - m$lambda %*% f$x)), 0.1)
  ## lambda is a convex combination under vrs.
  expect_equal(unname(rowSums(m$lambda)), rep(1, nrow(f$x)), tolerance = 1e-8)
})

test_that("the ratios, log slacks and caller-unit slacks agree with each other", {
  f <- mult_fixture()
  m <- dea_mult(f$x, f$y)
  expect_equal(m$ratio_x, exp(-m$log_slack_x))
  expect_equal(m$ratio_y, exp(m$log_slack_y))
  expect_equal(unname(m$slack_x), unname(f$x * (1 - m$ratio_x)))
  expect_equal(unname(m$slack_y), unname(f$y * (m$ratio_y - 1)))
  ## eff is exp of minus the total log slack, which is the product of the
  ## input factors times the reciprocal of the output factors.
  expect_equal(unname(m$eff),
               unname(exp(-(rowSums(m$log_slack_x) + rowSums(m$log_slack_y)))),
               tolerance = 1e-10)
  expect_true(all(m$eff > 0 & m$eff <= 1 + 1e-9))
  expect_true(all(m$ratio_x <= 1 + 1e-9))
  expect_true(all(m$ratio_y >= 1 - 1e-9))
  expect_equal(m$efficient, abs(m$eff - 1) < 1e-9)
})

test_that("it projects exactly onto a frontier that really is Cobb-Douglas", {
  ## dea_sim(frontier = "cobb") and this fixture both have f(x) = prod x_j^b, so
  ## in log space the true frontier is a HYPERPLANE. A piecewise-log-linear
  ## envelope built from points on it IS it, so every projection must land on
  ## it exactly. No piecewise-linear envelope can do that, and the additive
  ## model on the same data is carried along to show the difference is the
  ## model rather than noise.
  ##
  ## The reference set is given explicitly and lies entirely on the surface. On
  ## a self-referenced sample some peers are the inefficient points themselves,
  ## the envelope is then below the true surface, and the property does not
  ## hold -- which is a fact about the estimator, not a failure of the model.
  set.seed(11)
  xr <- matrix(runif(40 * 2, 2, 10), 40, 2)
  yr <- matrix((xr[, 1] * xr[, 2])^0.45, ncol = 1)
  set.seed(12)
  xe <- matrix(runif(10 * 2, 3.5, 8), 10, 2)
  ye <- matrix((xe[, 1] * xe[, 2])^0.45 * runif(10, 0.4, 0.9), ncol = 1)
  f <- function(X) (X[, 1] * X[, 2])^0.45

  m <- dea_mult(xe, ye, xref = xr, yref = yr)
  tx <- xe * m$ratio_x; ty <- as.numeric(ye * m$ratio_y)
  expect_lt(max(abs(ty - f(tx))), 1e-8)
  expect_true(all(m$eff < 1))              ## every evaluated point is inside

  a <- dea_add(xe, ye, measure = "unweighted", xref = xr, yref = yr,
               scaling = FALSE)
  ax <- xe - a$slack_x; ay <- as.numeric(ye + a$slack_y)
  expect_gt(max(abs(ay - f(ax))), 1e-3)

  ## A point ON the surface is efficient, whichever way it is approached.
  mm <- dea_mult(xr, yr)
  expect_true(all(mm$efficient))
  expect_equal(unname(mm$eff), rep(1, 40), tolerance = 1e-9)
})

test_that("zero and negative data are refused, because the model takes logs", {
  f <- mult_fixture(30)
  x0 <- f$x; x0[1, 1] <- 0
  expect_error(dea_mult(x0, f$y), "strictly positive")
  y0 <- f$y; y0[2, 1] <- 0
  expect_error(dea_mult(f$x, y0), "strictly positive")
  xn <- f$x; xn[1, 1] <- -1
  expect_error(dea_mult(xn, f$y), "Negative values")
})

test_that("non-vrs warns, and the warning is the only notice there is", {
  f <- mult_fixture(30)
  expect_warning(dea_mult(f$x, f$y, rts = "crs"), "units invariant")
  expect_warning(dea_mult(f$x, f$y, rts = "ndrs"), "1982 model")
  expect_silent(dea_mult(f$x, f$y))
  expect_silent(dea_mult(f$x, f$y, rts = "vrs"))
})

test_that("an external reference set works and the methods run", {
  f <- mult_fixture(40)
  g <- mult_fixture(12, seed = 99)
  ## Scoring against someone else's technology can be infeasible -- there may
  ## be no combination of reference DMUs dominating the point -- and that is
  ## reported rather than hidden, as it is for dea_add().
  expect_warning(m <- dea_mult(g$x, g$y, xref = f$x, yref = f$y), "INFEASIBLE")
  expect_true(sum(is.na(m$eff)) > 0L)
  expect_equal(m$n, 12L)
  expect_equal(m$nref, 40L)
  expect_equal(dim(m$lambda), c(12L, 40L))
  expect_false(m$self_ref)

  full <- dea_mult(f$x, f$y)
  expect_output(print(full), "piecewise Cobb-Douglas")
  expect_output(print(full), "units invariant")
  expect_output(summary(full), "GEOMETRIC")
  expect_equal(nobs(full), 40L)
  expect_equal(efficiency(full), full$eff)
})
