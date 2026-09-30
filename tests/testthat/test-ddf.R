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

## ---------------------------------------------------------------------------
## The range directional model (Portela, Thanassoulis and Simpson 2004).
## ---------------------------------------------------------------------------

test_that("the directional model matches Benchmarking on every direction and technology", {
  ## The DDF had no cross-package reference until 1.0.3, which mattered here
  ## more than elsewhere: the direction is an ARGUMENT, so a misread convention
  ## shows up in no other model. The design is one where nirs and ndrs each
  ## differ from both vrs and crs; see tools/make_reference_values.R.
  d <- toy(40, 2, 2, 12, returns = 0.8)
  ref <- ref_values("ddf_Benchmarking.csv")
  for (dir in unique(ref$direction)) for (rr in unique(ref$rts)) {
    want <- ref$beta[ref$direction == dir & ref$rts == rr]
    got <- suppressWarnings(
      dea_ddf(d$x, d$y, direction = dir, rts = rr, peers = FALSE)$beta)
    expect_equal(unname(got), want, tolerance = 1e-7,
                 info = paste("ddf", dir, rr))
  }
})

test_that("the range direction is the distance to the reference set's ideal point", {
  d <- toy(30, 2, 1, 5)
  f <- dea_ddf(d$x, d$y, direction = "range", rts = "vrs")

  ## g is exactly (x_o - min x, max y - y_o), and is reported.
  expect_equal(unname(f$g[, 1:2]),
               unname(sweep(d$x, 2L, apply(d$x, 2L, min), "-")))
  expect_equal(unname(f$g[, 3, drop = FALSE]),
               unname(-sweep(d$y, 2L, apply(d$y, 2L, max), "-")))

  ## and passing that matrix explicitly must give the identical answer -- the
  ## shorthand is allowed to be convenient, not to be a different program.
  G <- f$g
  expect_identical(
    unname(suppressWarnings(dea_ddf(d$x, d$y, direction = G, rts = "vrs")$beta)),
    unname(f$beta))
})

test_that("beta lies in [0, 1] under vrs, which is what the range direction buys", {
  ## Under vrs, X'lambda is a convex combination of the reference inputs and so
  ## cannot fall below their minimum; feasibility then forces beta <= 1. So the
  ## score reads as "the fraction of the distance to the ideal point that this
  ## DMU could travel", which a plain directional beta does not.
  for (sd in 1:6) {
    d <- toy(30, 2, 2, sd)
    b <- dea_ddf(d$x, d$y, direction = "range", rts = "vrs")$beta
    expect_true(all(b >= -1e-9 & b <= 1 + 1e-9), info = paste("seed", sd))
  }
})

test_that("beta is translation invariant under vrs, and is NOT under crs", {
  ## This is the property the model exists for: it is what lets negative data
  ## be handled without translating it first. The direction is a DIFFERENCE of
  ## two values in a column, so a shift leaves it alone, and the shift moves
  ## both sides of the constraint by t * sum(lambda) -- equal only when
  ## sum(lambda) = 1.
  ##
  ## The crs half is not decoration. Where a claim is "this holds because of
  ## the vrs row", a test that only shows it holding under vrs would pass just
  ## as well if the invariance came from somewhere else entirely.
  d <- toy(35, 2, 1, 4)
  shift_x <- c(-5, 3); shift_y <- -2
  x2 <- sweep(d$x, 2L, shift_x, "+"); y2 <- d$y + shift_y
  expect_true(any(x2 < 0) && any(y2 < 0))

  b_vrs <- dea_ddf(d$x, d$y, direction = "range", rts = "vrs")$beta
  expect_equal(unname(dea_ddf(x2, y2, direction = "range", rts = "vrs")$beta),
               unname(b_vrs), tolerance = 1e-8)

  b_crs <- dea_ddf(d$x, d$y, direction = "range", rts = "crs")$beta
  expect_gt(max(abs(dea_ddf(x2, y2, direction = "range", rts = "crs")$beta -
                    b_crs)), 0.01)
})

test_that("the range direction is units invariant", {
  d <- toy(30, 2, 1, 6)
  b <- dea_ddf(d$x, d$y, direction = "range", rts = "vrs")$beta
  x2 <- sweep(d$x, 2L, c(1000, 0.01), "*"); y2 <- d$y * 7
  expect_equal(unname(dea_ddf(x2, y2, direction = "range", rts = "vrs")$beta),
               unname(b), tolerance = 1e-8)
})

test_that("a DMU that is the ideal point is answered without a program", {
  ## Best in every input AND every output: its direction is zero, so the
  ## program would be UNBOUNDED rather than infeasible -- beta appears in no
  ## constraint. beta = 0 is the answer, and peers are NA because nothing was
  ## solved.
  x <- matrix(c(1, 3, 4, 5, 6, 7, 8), ncol = 1)
  y <- matrix(c(9, 4, 6, 5, 7, 3, 6), ncol = 1)
  expect_warning(f <- dea_ddf(x, y, direction = "range", rts = "vrs"),
                 "ideal point")
  expect_equal(unname(f$beta[1]), 0)
  expect_true(f$efficient[1])
  expect_true(all(is.na(f$lambda[1, ])))
  expect_false(anyNA(f$lambda[-1, ]))
  expect_identical(f$status[1], 0L)

  ## Every other direction still refuses a zero direction, because for those it
  ## really is a mistake rather than a corner of the data.
  expect_error(dea_ddf(x, matrix(0, nrow(x), 1), direction = "out"), "all zero")
})

test_that("a DMU outside the reference set's range has that coordinate clamped", {
  ## Only reachable with a separate xref/yref: self-referenced, a DMU cannot be
  ## better than the minimum of a set it belongs to. A negative component would
  ## ask an input that already beats the whole technology to GROW.
  d <- toy(30, 1, 1, 9)
  keep <- d$x[, 1] > min(d$x[, 1])
  w <- character()
  f <- withCallingHandlers(
    dea_ddf(d$x, d$y, direction = "range", rts = "vrs",
            xref = d$x[keep, , drop = FALSE], yref = d$y[keep, , drop = FALSE]),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("outside the reference set", w)))
  expect_true(all(f$g >= 0))

  ## And the two facts travel together on this design: the DMU that sits below
  ## the reference set's minimum input is also outside the technology the
  ## reference set spans, so its program is infeasible. The clamp keeps the
  ## direction meaningful; it does not manufacture a score.
  expect_true(any(grepl("INFEASIBLE", w)))
  expect_true(anyNA(f$beta))
})

test_that("dea_profit accepts the range direction and reports why it cannot at the ideal point", {
  d <- toy(25, 2, 1, 3)
  w <- matrix(1, nrow(d$x), 2); r <- matrix(1, nrow(d$y), 1)
  f <- dea_profit(d$x, d$y, w, r, direction = "range", rts = "vrs")
  expect_equal(unname(f$eff), unname(f$tech + f$alloc), tolerance = 1e-8)

  x <- matrix(c(1, 3, 4, 5, 6, 7, 8), ncol = 1)
  y <- matrix(c(9, 4, 6, 5, 7, 3, 6), ncol = 1)
  expect_error(
    suppressWarnings(dea_profit(x, y, matrix(1, nrow(x), 1), matrix(1, nrow(x), 1),
                                direction = "range", rts = "vrs")),
    "zero or negative")
})
