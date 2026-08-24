test_that("dea_sim() produces a technology matching its own stated truth", {
  s <- dea_sim(200, p = 2, q = 3, returns = 0.8, seed = 1)
  expect_s3_class(s, "dea_sim")
  expect_equal(dim(s$x), c(200L, 2L))
  expect_equal(dim(s$y), c(200L, 3L))
  ## ||y|| = f(x) exp(-u), so phi = f(x)/||y|| = exp(u) exactly.
  expect_equal(sqrt(rowSums(s$y^2)), s$frontier * exp(-s$u), tolerance = 1e-12)
  expect_equal(s$phi, exp(s$u), tolerance = 1e-12)
  expect_equal(s$theta, exp(-s$u / 0.8), tolerance = 1e-12)
  expect_true(all(s$theta > 0 & s$theta <= 1))
  expect_true(all(s$phi >= 1))
  ## E[u] is what was asked for.
  expect_equal(mean(s$u), 0.3, tolerance = 0.1)
})

test_that("a seed reproduces, and restores the caller's stream", {
  a <- dea_sim(50, seed = 4)
  b <- dea_sim(50, seed = 4)
  expect_equal(a$x, b$x)
  expect_equal(a$u, b$u)
  set.seed(11); before <- runif(1)
  set.seed(11); invisible(dea_sim(50, seed = 4)); after <- runif(1)
  expect_equal(before, after)
})

test_that("a constant-returns fit recovers the truth on a constant-returns design", {
  s <- dea_sim(600, p = 1, q = 1, returns = 1, seed = 6)
  fit <- dea(s$x, s$y, rts = "crs", orientation = "in", slack = FALSE)
  ## The frontier here is a ray, so the estimator is essentially exact.
  expect_lt(mean(fit$eff - s$theta), 0.01)
  expect_gt(cor(fit$eff, s$theta), 0.999)
  ## And the bias is one-signed: the estimated frontier is inside the true one.
  expect_true(all(fit$eff >= s$theta - 1e-8))
})

test_that("the estimator is biased toward 1, always", {
  s <- dea_sim(100, p = 2, q = 1, returns = 0.9, seed = 12)
  fo <- dea(s$x, s$y, rts = "vrs", orientation = "out", slack = FALSE)
  ## phi_hat <= phi: the estimated frontier cannot lie above the true one.
  expect_true(all(fo$eff <= s$phi + 1e-8))
  expect_lt(mean(fo$eff - s$phi), 0)
})

test_that("dea_rate() states the published rates", {
  ## Kneip, Park and Simar (1998); Park, Simar and Weiner (2000).
  expect_equal(dea_rate(1, 1, "vrs"), -4 / 3)
  expect_equal(dea_rate(2, 1, "vrs"), -1)
  expect_equal(dea_rate(1, 1, "crs"), -2)
  expect_equal(dea_rate(2, 2, "crs"), -1)
  expect_equal(dea_rate(1, 1, "fdh"), -1)
  expect_equal(dea_rate(4, 4, "vrs"), -4 / 9)
  ## The estimator's own rate is half the mean-squared-error rate.
  expect_equal(dea_rate(2, 1, "vrs", "estimator"), -1 / 2)
  ## The curse of dimensionality, as a monotone statement.
  d <- vapply(1:8, function(k) dea_rate(k, k, "vrs"), numeric(1))
  expect_true(all(diff(d) > 0))
})

test_that("x_range moves the input support and nothing else", {
  a <- dea_sim(200, p = 2, q = 1, x_range = c(1, 2), seed = 1)
  b <- dea_sim(200, p = 2, q = 1, x_range = c(0.3, 2), seed = 1)
  expect_true(all(a$x >= 1 & a$x <= 2))
  expect_true(all(b$x >= 0.3 & b$x <= 2))
  expect_lt(min(b$x), 1)
  ## The inefficiency draws are independent of x, so widening the support does
  ## not change them -- which is what makes the two designs comparable.
  expect_equal(a$u, b$u)
  ## And the stated truth still matches the realized data.
  expect_equal(sqrt(rowSums(b$y^2)), b$frontier * exp(-b$u), tolerance = 1e-12)
})

test_that("widening the support downward makes theta estimable", {
  ## On the default support a variable-returns fit cannot reach theta, because
  ## the projection leaves the sample. Given room underneath, it can.
  ev <- function(m, lo, hi, seed = 99) {
    set.seed(seed)
    x <- matrix(runif(m, lo, hi), m, 1)
    u <- rexp(m, 1 / 0.3)
    list(x = x, y = matrix(x^0.8 * exp(-u), m, 1), theta = exp(-u / 0.8))
  }
  E <- ev(60, 1.3, 1.7)
  mse <- function(xlo, n) {
    v <- numeric(12)
    for (b in seq_len(12)) {
      s <- dea_sim(n, p = 1, q = 1, returns = 0.8, x_range = c(xlo, 2), seed = 500 + b)
      f <- suppressWarnings(dea(E$x, E$y, rts = "vrs", orientation = "in",
                                slack = FALSE, peers = FALSE,
                                xref = s$x, yref = s$y))
      v[b] <- mean((f$eff - E$theta)^2)
    }
    mean(v)
  }
  narrow <- c(mse(1.0, 200), mse(1.0, 800))
  wide   <- c(mse(0.2, 200), mse(0.2, 800))
  ## On the narrow support the error barely moves; on the wide one it falls.
  expect_gt(narrow[2] / narrow[1], 0.8)
  expect_lt(wide[2] / wide[1], 0.6)
  expect_lt(wide[2], narrow[2] / 5)
})

test_that("dea_sim() validates its design arguments", {
  expect_error(dea_sim(50, returns = 1.5), "must lie in")
  expect_error(dea_sim(50, returns = 0), "must lie in")
  expect_error(dea_sim(0), "positive number of DMUs")
  expect_error(dea_sim(50, p = 0), "at least 1")
  expect_error(dea_sim(50, mean_ineff = -1), "positive")
  expect_error(dea_sim(50, x_range = c(2, 1)), "increasing positive")
  expect_error(dea_sim(50, x_range = c(0, 2)), "increasing positive")
  expect_error(dea_sim(50, x_range = 1), "increasing positive")
})
