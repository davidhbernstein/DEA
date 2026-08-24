test_that("the one-sided directions reproduce the radial models exactly", {
  d <- toy(n = 50)
  for (rts in c("crs", "vrs", "nirs", "ndrs")) {
    bi <- dea_ddf(d$x, d$y, direction = "in",  rts = rts)
    ti <- dea(d$x, d$y, orientation = "in",  rts = rts, slack = FALSE)
    expect_equal(unname(bi$beta), unname(1 - ti$eff), tolerance = 1e-8, info = rts)

    bo <- dea_ddf(d$x, d$y, direction = "out", rts = rts)
    to <- dea(d$x, d$y, orientation = "out", rts = rts, slack = FALSE)
    expect_equal(unname(bo$beta), unname(to$eff - 1), tolerance = 1e-8, info = rts)
  }
})

test_that("beta is non-negative for observed DMUs and zero on the frontier", {
  d <- toy(n = 50)
  f <- dea_ddf(d$x, d$y, direction = "both", rts = "vrs")
  expect_true(all(f$beta >= -1e-9))
  expect_true(any(f$beta == 0))
  ## Every DMU is in its own technology, so beta = 0 is always feasible and the
  ## solution can never be negative here.
  expect_true(all(f$efficient == (f$beta == 0)))
})

test_that("the directional model accepts data the others cannot", {
  d <- toy(n = 40, p = 1, q = 2)
  y <- d$y
  y[1:3, 1] <- 0            ## a zero output
  expect_error(dea_sbm(d$x, y), "positive")
  f <- dea_ddf(d$x, y, direction = "unit", rts = "vrs")
  expect_true(all(is.finite(f$beta)))

  ## And negative values, under a fixed direction.
  yneg <- d$y; yneg[1:3, 1] <- -0.2
  f2 <- dea_ddf(d$x, yneg, direction = "unit", rts = "vrs")
  expect_true(all(is.finite(f2$beta)))
  expect_error(dea(d$x, yneg), "Negative values")
})

test_that("a fixed direction is invariant to the units of the data", {
  d <- toy(n = 40, p = 2, q = 1)
  a <- dea_ddf(d$x, d$y, direction = c(1, 1, 1), rts = "vrs")
  ## Rescaling an input and rescaling the direction with it must leave beta
  ## alone; rescaling only the data must not.
  b <- dea_ddf(cbind(d$x[, 1] * 10, d$x[, 2]), d$y,
               direction = c(10, 1, 1), rts = "vrs")
  expect_equal(unname(a$beta), unname(b$beta), tolerance = 1e-7)
})

test_that("direction arguments are validated", {
  d <- toy(n = 40, p = 2, q = 1)
  expect_error(dea_ddf(d$x, d$y, direction = c(1, 1)), "length p \\+ q")
  expect_error(dea_ddf(d$x, d$y, direction = matrix(1, 2, 3)), "must be")
  expect_error(dea_ddf(d$x, d$y, direction = "sideways"), "must be one of")
  expect_warning(dea_ddf(d$x, d$y, direction = c(-1, 1, 1)), "negative entries")
})

test_that("efficiency(type = \"score\") refuses an unconvertible direction", {
  d <- toy(n = 40)
  f <- dea_ddf(d$x, d$y, direction = "both", rts = "vrs")
  expect_error(efficiency(f, "score"), "ADDITIVE")
  expect_equal(efficiency(f, "natural"), f$beta)

  fi <- dea_ddf(d$x, d$y, direction = "in", rts = "vrs")
  ti <- dea(d$x, d$y, orientation = "in", rts = "vrs", slack = FALSE)
  expect_equal(unname(efficiency(fi, "score")), unname(ti$eff), tolerance = 1e-8)
})
