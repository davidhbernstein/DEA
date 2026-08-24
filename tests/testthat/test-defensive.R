test_that("bad data is refused with an actionable message", {
  d <- toy(n = 40, p = 2, q = 1)

  xna <- d$x; xna[1, 1] <- NA
  expect_error(dea(xna, d$y), "missing or non-finite")
  xinf <- d$x; xinf[1, 1] <- Inf
  expect_error(dea(xinf, d$y), "missing or non-finite")

  expect_error(dea(d$x[1:10, ], d$y), "different numbers of DMUs")

  xneg <- d$x; xneg[1, 1] <- -1
  expect_error(dea(xneg, d$y), "Negative values")
  expect_error(dea_sbm(xneg, d$y), "Negative values")

  y0 <- d$y; y0[1, 1] <- 0
  expect_error(dea(d$x, y0), "no output at all")

  expect_error(dea("x1", d$y, data = NULL), "needs `data`")
  s <- dea_sim(40, p = 2, q = 1, seed = 1)
  expect_error(dea(c("x1", "nope"), "y1", data = s$data), "not found in `data`")
})

test_that("a sample too small for its dimensions warns rather than lying", {
  s <- dea_sim(8, p = 2, q = 2, seed = 1)
  expect_warning(dea(s$x, s$y, rts = "vrs"), "efficient by dimension alone")
  ## And the warning is right: nearly everything comes back efficient.
  fit <- suppressWarnings(dea(s$x, s$y, rts = "vrs"))
  expect_gt(mean(fit$eff == 1), 0.5)
})

test_that("zeros are handled where they are defined and refused where they are not", {
  s <- dea_sim(40, p = 2, q = 1, seed = 2)
  x0 <- s$x; x0[1:3, 1] <- 0
  ## A radial input-oriented program tolerates a zero input.
  fit <- dea(x0, s$y, rts = "vrs", orientation = "in")
  expect_true(all(is.finite(fit$eff)))
  ## The slacks-based measure divides by it and says so.
  expect_error(dea_sbm(x0, s$y), "strictly positive")
})

test_that("degenerate samples do not crash", {
  ## Every DMU identical: all efficient, no slack, one peer set.
  x <- matrix(rep(2, 20), ncol = 1); y <- matrix(rep(3, 20), ncol = 1)
  fit <- suppressWarnings(dea(x, y, rts = "vrs"))
  expect_true(all(fit$eff == 1))
  expect_true(all(fit$efficient))
  ## A single reference DMU.
  one <- suppressWarnings(dea(x, y, rts = "crs",
                              xref = matrix(1, 1, 1), yref = matrix(1, 1, 1)))
  expect_true(all(is.finite(one$eff)))
})

test_that("argument validation covers the switches", {
  d <- toy(n = 40)
  expect_error(dea(d$x, d$y, orientation = "sideways"), "must be one of")
  expect_error(dea(d$x, d$y, rts = "fdh", super = TRUE), "not implemented for")
  expect_error(dea_sbm(d$x, d$y, rts = "fdh"), "must be one of")
  expect_error(dea_sbm(d$x, d$y, orientation = "both"), "must be one of")
})
