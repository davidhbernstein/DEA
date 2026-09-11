## Cost, revenue and profit efficiency.
##
## The identity that matters is the decomposition. If it holds EXACTLY -- not
## to 1e-6 but to the solver's own precision -- then the technical factor and
## the price factor were computed against the same frontier, which is the one
## thing that is easy to get wrong and impossible to see in the numbers.

pw <- function(n, k, lo = 1, hi = 3, seed = 1) {
  set.seed(seed); matrix(stats::runif(n * k, lo, hi), n, k)
}

test_that("cost efficiency factors exactly into technical x allocative", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  for (rts in c("crs", "vrs", "nirs", "ndrs")) {
    f <- dea_cost(d$x, d$y, pw(30, 2), rts = rts)
    expect_equal(unname(f$eff), unname(f$technical * f$allocative),
                 tolerance = 1e-12, info = rts)
    expect_true(all(f$eff > 0 & f$eff <= 1 + 1e-9), info = rts)
    expect_true(all(f$allocative > 0 & f$allocative <= 1 + 1e-9), info = rts)
    ## Cost efficiency can never exceed technical efficiency: the cost program
    ## is free to substitute as well as to contract.
    expect_true(all(f$eff <= f$technical + 1e-9), info = rts)
    expect_true(all(f$optimal <= f$observed + 1e-9), info = rts)
  }
})

test_that("revenue efficiency factors the same way", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  for (rts in c("crs", "vrs")) {
    f <- dea_revenue(d$x, d$y, pw(30, 2, 2, 5), rts = rts)
    expect_equal(unname(f$eff), unname(f$technical * f$allocative),
                 tolerance = 1e-12, info = rts)
    expect_true(all(f$eff > 0 & f$eff <= 1 + 1e-9), info = rts)
    expect_true(all(f$optimal >= f$observed - 1e-9), info = rts)
  }
})

test_that("the Nerlovian profit gap splits into technical plus allocative", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  f <- dea_profit(d$x, d$y, pw(30, 2), pw(30, 2, 2, 5, seed = 2))
  expect_equal(unname(f$eff), unname(f$technical + f$allocative),
               tolerance = 1e-12)
  ## Both parts are losses, so both are non-negative and BIGGER IS WORSE.
  expect_true(all(f$eff >= -1e-9))
  expect_true(all(f$technical >= -1e-9))
  expect_true(all(f$allocative >= -1e-9))
  expect_true(all(f$profit_max >= f$profit_obs - 1e-9))
})

test_that("scaling changes nothing in the price models", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  W <- pw(30, 2)
  a <- dea_cost(d$x, d$y, W, rts = "vrs", scaling = TRUE)
  b <- dea_cost(d$x, d$y, W, rts = "vrs", scaling = FALSE)
  expect_equal(a$eff, b$eff, tolerance = 1e-9)
  expect_equal(a$optimal, b$optimal, tolerance = 1e-8)
  expect_equal(a$optimal_q, b$optimal_q, tolerance = 1e-8)
})

test_that("the price models match Benchmarking's recorded optima", {
  d <- toy(25, p = 2, q = 2, returns = 0.9, seed = 12)
  n <- nrow(d$x)
  W <- matrix(rep(c(1.2, 0.8), each = n), n, 2)
  P <- matrix(rep(c(2.0, 3.1), each = n), n, 2)
  ref <- ref_values("price_Benchmarking.csv")
  for (rts in c("vrs", "crs", "nirs", "ndrs")) {
    r <- ref[ref$rts == rts, ]; r <- r[order(r$dmu), ]
    expect_equal(unname(dea_cost(d$x, d$y, W, rts = rts)$optimal),
                 r$cost_opt, tolerance = 1e-8, info = rts)
    expect_equal(unname(dea_revenue(d$x, d$y, P, rts = rts)$optimal),
                 r$revenue_opt, tolerance = 1e-8, info = rts)
  }
  r <- ref[ref$rts == "vrs", ]; r <- r[order(r$dmu), ]
  expect_equal(unname(dea_profit(d$x, d$y, W, P, rts = "vrs")$profit_max),
               r$profit_max_vrs, tolerance = 1e-8)
})

test_that("profit is refused under constant and non-decreasing returns", {
  d <- toy(25, p = 2, q = 1)
  expect_error(dea_profit(d$x, d$y, pw(25, 2), pw(25, 1, 2, 5), rts = "crs"),
               "unbounded")
  expect_error(dea_profit(d$x, d$y, pw(25, 2), pw(25, 1, 2, 5), rts = "ndrs"),
               "unbounded")
  ## The alias too -- "irs" is ndrs by another name.
  expect_error(dea_profit(d$x, d$y, pw(25, 2), pw(25, 1, 2, 5), rts = "irs"),
               "unbounded")
})

test_that("a length-p price vector is one price list shared by every DMU", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  a <- dea_cost(d$x, d$y, c(1.2, 0.8))
  b <- dea_cost(d$x, d$y, matrix(rep(c(1.2, 0.8), each = 30), 30, 2))
  expect_equal(a$eff, b$eff, tolerance = 1e-12)
})

test_that("ambiguous and impossible price arguments are refused", {
  ## n == p makes "one list for everyone" and "one price each" the same shape.
  d <- toy(2, p = 2, q = 1)
  expect_error(suppressWarnings(dea_cost(d$x, d$y, c(1, 2))), "ambiguous")
  e <- toy(30, p = 2, q = 2, returns = 0.9)
  expect_error(dea_cost(e$x, e$y, c(1, 2, 3)), "neither")
  expect_error(dea_cost(e$x, e$y, c(1, -2)), "non-positive")
  expect_error(dea_cost(e$x, e$y, c(1, 0)), "non-positive")
})

test_that("prices can be named columns of `data`", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  df <- data.frame(d$x, d$y, w1 = 1.2, w2 = 0.8)
  a <- dea_cost(c("x1", "x2"), c("y1", "y2"), c("w1", "w2"), data = df)
  b <- dea_cost(d$x, d$y, c(1.2, 0.8))
  expect_equal(unname(a$eff), unname(b$eff), tolerance = 1e-12)
})

test_that("the price models take an external reference set", {
  ref <- toy(50, p = 2, q = 2, returns = 0.9, seed = 5)
  ev  <- toy(12, p = 2, q = 2, returns = 0.9, seed = 6)
  ## Infeasibility here is CORRECT, not a failure. The cost program needs a
  ## convex combination of the reference DMUs to match the evaluated DMU's
  ## output, and under variable returns an outsider producing more than any
  ## reference DMU has none -- the same situation dea_sbm() warns about. The
  ## honest answer is NA, so what is tested is that it warns and returns NA
  ## rather than silently reporting a number.
  expect_warning(f <- dea_cost(ev$x, ev$y, pw(12, 2), rts = "vrs",
                               xref = ref$x, yref = ref$y),
                 "INFEASIBLE")
  expect_identical(length(f$eff), 12L)
  expect_identical(dim(f$lambda), c(12L, 50L))
  expect_true(any(is.na(f$eff)))
  expect_identical(sum(is.na(f$eff)), sum(!f$status %in% c(0L, 1L)))
  ok <- is.finite(f$eff)
  expect_equal(unname(f$eff[ok]), unname((f$technical * f$allocative)[ok]),
               tolerance = 1e-12)

  ## Under constant returns the cone always reaches, so nothing is infeasible.
  g <- dea_cost(ev$x, ev$y, pw(12, 2), rts = "crs",
                xref = ref$x, yref = ref$y)
  expect_true(all(is.finite(g$eff)))
})

test_that("the methods work and fitted() refuses the profit model", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  f <- dea_cost(d$x, d$y, pw(30, 2))
  expect_output(print(f), "Cost efficiency")
  expect_output(print(f), "technical  x  allocative")
  expect_identical(dim(fitted(f)), dim(d$x))
  expect_identical(nobs(f), 30L)
  expect_equal(efficiency(f, "technical"), f$technical)
  expect_true(all(c("dmu", "peer", "lambda") %in% names(peers(f))))
  g <- dea_profit(d$x, d$y, pw(30, 2), pw(30, 2, 2, 5, seed = 2))
  expect_output(print(g), "Nerlovian")
  expect_error(fitted(g), "moves inputs AND outputs")
})
