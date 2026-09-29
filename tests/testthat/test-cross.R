## Cross-efficiency. The substantive claim being tested is the BRACKET: the
## benevolent and aggressive secondary goals must straddle whatever an
## arbitrary choice among alternate optima would have produced, and the gap
## between them must be visible rather than numerical -- if it were negligible
## the secondary goal would not be worth having.

test_that("the diagonal of the cross matrix is the ordinary DEA score", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  f <- dea_cross(d$x, d$y, secondary = "none")
  own <- suppressWarnings(dea(d$x, d$y, rts = "crs", orientation = "in",
                              slack = FALSE))$eff
  expect_equal(unname(diag(f$cross_matrix)), unname(own), tolerance = 1e-8)
  expect_equal(unname(f$own), unname(own), tolerance = 1e-12)
})

test_that("cross-efficiency never exceeds a DMU's own score", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  for (sec in c("none", "benevolent", "aggressive")) {
    f <- dea_cross(d$x, d$y, secondary = sec)
    expect_true(all(f$eff <= f$own + 1e-8), info = sec)
    expect_true(all(f$eff > 0), info = sec)
    expect_true(all(f$maverick >= -1e-8), info = sec)
  }
})

test_that("benevolent brackets aggressive from above, by a visible margin", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  b <- dea_cross(d$x, d$y, secondary = "benevolent")
  a <- dea_cross(d$x, d$y, secondary = "aggressive")
  expect_true(all(b$eff >= a$eff - 1e-8))
  ## Not a rounding difference: if this ever collapses to zero the secondary
  ## programs have stopped doing anything and the test should fail loudly.
  expect_gt(mean(b$eff - a$eff), 1e-3)
  ## The secondary goal must not disturb the DMU's OWN score, which it holds
  ## fixed by construction.
  expect_equal(b$own, a$own, tolerance = 1e-9)
})

test_that("cross-efficiency breaks ties that the DEA score cannot", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  f <- dea_cross(d$x, d$y, secondary = "benevolent")
  n_eff <- sum(f$own >= 1 - 1e-9)
  skip_if(n_eff < 2, "no ties to break on this sample")
  expect_gt(n_eff, 1)
  expect_identical(length(unique(round(f$eff, 8))), length(f$eff))
})

test_that("excluding the self-appraisal lowers every score it can change", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  a <- dea_cross(d$x, d$y, secondary = "none", self = TRUE)
  b <- dea_cross(d$x, d$y, secondary = "none", self = FALSE)
  ## The diagonal is the most generous entry in its column, so dropping it can
  ## only pull the mean down.
  expect_true(all(b$eff <= a$eff + 1e-9))
})

test_that("a secondary goal under variable returns is refused, with a reason", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  expect_error(dea_cross(d$x, d$y, rts = "vrs", secondary = "benevolent"),
               "not bounded by 1")
  ## But descriptive cross-efficiency under vrs is allowed.
  f <- dea_cross(d$x, d$y, rts = "vrs", secondary = "none")
  expect_s3_class(f, "dea_cross")
  expect_false(is.null(f$u0))
})

test_that("an external rater panel scores outsiders, without a maverick index", {
  ref <- toy(50, p = 2, q = 2, returns = 0.9, seed = 5)
  ev  <- toy(12, p = 2, q = 2, returns = 0.9, seed = 6)
  f <- dea_cross(ev$x, ev$y, secondary = "benevolent",
                 xref = ref$x, yref = ref$y)
  expect_identical(dim(f$cross_matrix), c(50L, 12L))
  expect_identical(length(f$eff), 12L)
  expect_null(f$own)
  expect_null(f$maverick)
  expect_error(efficiency(f, "maverick"), "needs each rated DMU")
})

test_that("the methods report what they should", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  f <- dea_cross(d$x, d$y, secondary = "aggressive")
  expect_output(print(f), "Cross-efficiency")
  expect_output(print(f), "aggressive")
  expect_output(print(dea_cross(d$x, d$y, secondary = "none")), "ARBITRARY")
  expect_identical(nobs(f), 30L)
  expect_equal(efficiency(f, "cross"), f$eff)
  m <- multipliers(f)
  expect_identical(dim(m), c(30L, 4L))
  expect_true(all(m >= -1e-9))
})

test_that("a DMU with no usable appraisal is NA, not NaN, and says so", {
  ## Every input zero means v'x is zero whatever v is, so the cross-efficiency
  ## ratio has no denominator and every rater's appraisal of this DMU is
  ## missing. That used to surface as eff = NaN and spread = -Inf, explained
  ## only by base R's "no non-missing arguments to max" -- which tells the
  ## reader nothing they can act on.
  set.seed(3)
  xr <- cbind(runif(12, 1, 3), runif(12, 1, 3)); yr <- cbind(runif(12, 1, 3))
  xe <- rbind(c(2, 2), c(0, 0)); ye <- rbind(2, 1)

  expect_warning(f <- dea_cross(xe, ye, xref = xr, yref = yr, rts = "crs",
                                secondary = "benevolent"),
                 "NO usable appraisal")
  expect_true(is.na(f$eff[2]))
  expect_true(is.na(f$spread[2]))
  expect_false(is.nan(f$eff[2]))          ## NA, and specifically not NaN
  expect_true(is.finite(f$eff[1]))        ## the ordinary DMU is untouched
  expect_true(is.finite(f$spread[1]))
})

test_that("ordinary cross-efficiency is unchanged by the NA handling", {
  ## The guard must not perturb the normal path: every DMU is appraised by the
  ## whole panel, so no warning and no missing entries.
  s <- dea_sim(20, p = 2, q = 1, seed = 4)
  expect_silent(f <- dea_cross(s$x, s$y, rts = "crs", secondary = "benevolent"))
  expect_false(any(is.na(f$eff)))
  expect_true(all(is.finite(f$spread)))
  ## The mean over a full column is exactly colMeans, which is what the old
  ## code computed and what the new code must still compute.
  expect_equal(unname(f$eff), unname(colMeans(f$cross_matrix)), tolerance = 1e-12)
})
