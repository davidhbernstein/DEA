## ---------------------------------------------------------------------------
## The two-stage procedure of Simar and Wilson (2007), and the truncated
## regression it is built on.
## ---------------------------------------------------------------------------

## A technology whose TRUE Shephard distances are a known function of z. This
## is the only check here that can catch a misunderstanding rather than a
## coding error, which is why it exists in the same spirit as dea_sim().
sw_dgp <- function(n, beta, sigma = 0.15, p = 2, seed = 1) {
  set.seed(seed)
  X <- matrix(stats::runif(n * p, 1, 10), n, p)
  z <- cbind(1, stats::runif(n), stats::rnorm(n))
  delta <- DEA:::.rtnorm_left(as.numeric(z %*% beta), sigma, 1)
  f <- apply(X, 1L, function(r) prod(r)^(0.9 / p))
  list(x = X, y = matrix(f / delta, ncol = 1L), z = z[, -1L, drop = FALSE],
       delta = delta)
}

test_that("the truncated regression matches the reference implementations", {
  ## .truncreg() is written here because neither truncreg nor npsf is a
  ## dependency of this package and neither will become one. The pinned values
  ## are what they computed; see tools/make_reference_values.R.
  ref <- ref_values("truncreg_reference.csv")
  set.seed(7)
  n <- 400
  z1 <- stats::runif(n); z2 <- stats::rnorm(n)
  y <- 1.2 + 0.8 * z1 - 0.5 * z2 + stats::rnorm(n, 0, 0.6)
  keep <- y > 1
  Z <- cbind(`(Intercept)` = 1, z1 = z1[keep], z2 = z2[keep])
  got <- DEA:::.truncreg(y[keep], Z, ll = 1)
  expect_true(got$converged)
  expect_equal(unname(c(got$beta, got$sigma)), ref$estimate, tolerance = 1e-5)
  expect_equal(got$logLik, ref$logLik[1], tolerance = 1e-6)
})

test_that("the truncated likelihood survives a sample far from the boundary", {
  ## 1 - Phi(a) underflows to zero past about a = 8, which on this scale is an
  ## ordinary DMU a long way from the frontier, and a likelihood written in
  ## levels returns Inf there rather than a number. The tail is computed in
  ## logs for that reason.
  set.seed(3)
  y <- 50 + stats::rnorm(200, 0, 0.5)
  Z <- cbind(1, stats::rnorm(200))
  f <- DEA:::.truncreg(y, Z, ll = 1)
  expect_true(f$converged)
  expect_true(is.finite(f$logLik))
  ## far from the truncation point it must agree with ordinary least squares
  expect_equal(unname(f$beta), unname(stats::lm.fit(Z, y)$coefficients),
               tolerance = 1e-4)
})

test_that("truncated draws stay above the truncation point and have the right mean", {
  set.seed(5)
  mu <- rep(0.5, 20000)
  d <- DEA:::.rtnorm_left(mu, 0.4, 1)
  expect_true(all(d > 1))
  ## E[X | X > lo] = mu + sigma * lambda, the inverse Mills ratio
  a <- (1 - 0.5) / 0.4
  want <- 0.5 + 0.4 * stats::dnorm(a) / stats::pnorm(a, lower.tail = FALSE)
  expect_equal(mean(d), want, tolerance = 0.01)

  ## and a draw whose mean is far above the point must not collapse: computing
  ## p + (1 - p) * u in levels loses every digit of u here.
  far <- DEA:::.rtnorm_left(rep(40, 5000), 1, 1)
  expect_gt(stats::sd(far), 0.9)
  expect_lt(abs(mean(far) - 40), 0.1)
})

test_that("both algorithms recover a beta the design was given", {
  ## delta_hat is biased toward 1, so a regression on it inherits that shift --
  ## and the shift lands on the INTERCEPT, because it is roughly common across
  ## DMUs. Algorithm #2 exists to remove it. Both halves are asserted: #1 comes
  ## out below the truth, and #2 comes out nearer than #1.
  d <- sw_dgp(150, c(1.05, 0.5, -0.25), seed = 111)
  fit <- dea(d$x, d$y, rts = "vrs", orientation = "out", slack = FALSE)
  expect_lt(mean(fit$eff), mean(d$delta))      # the bias is real and signed

  r1 <- dea_reg(fit, d$z, algorithm = "sw1", L2 = 200, seed = 1,
                progress = FALSE)
  r2 <- dea_reg(fit, d$z, algorithm = "sw2", L1 = 40, L2 = 200, seed = 1,
                progress = FALSE)
  truth <- c(1.05, 0.5, -0.25)

  expect_lt(r1$coefficients[[1]], truth[1])
  expect_lt(abs(r2$coefficients[[1]] - truth[1]),
            abs(r1$coefficients[[1]] - truth[1]))

  ## The slopes are not where the bias goes, and both algorithms should find
  ## them. Judged against the fit's OWN bootstrap standard error rather than
  ## a tolerance picked to make one seed pass: at n = 150 a slope estimate has
  ## a real sampling spread, and a fixed tolerance either hides that or turns
  ## an ordinary draw into a failure.
  for (r in list(r1, r2)) for (j in 2:3) {
    expect_lt(abs(r$coefficients[[j]] - truth[j]), 3 * r$table$se[j])
  }
})

test_that("the bias correction moves the scores the way the estimator is wrong", {
  ## delta_hat under-states the distance to the true frontier, so a correction
  ## that worked must move it UP -- for every DMU, since the DEA bias is
  ## one-signed.
  d <- sw_dgp(120, c(1.05, 0.5, -0.25), seed = 12)
  fit <- dea(d$x, d$y, rts = "vrs", orientation = "out", slack = FALSE)
  r <- dea_reg(fit, d$z, algorithm = "sw2", L1 = 40, L2 = 100, seed = 2,
               progress = FALSE)
  expect_true(all(r$delta_bc >= r$delta - 1e-8))
  expect_gt(mean(r$delta_bc - r$delta), 0)
  ## and toward the truth on average
  expect_lt(abs(mean(r$delta_bc) - mean(d$delta)),
            abs(mean(r$delta) - mean(d$delta)))
  expect_null(dea_reg(fit, d$z, algorithm = "sw1", L2 = 50, seed = 2,
                      progress = FALSE)$delta_bc)
})

test_that("the regression is in Shephard's scale, so both orientations agree", {
  ## delta = 1/theta input-oriented and phi output-oriented, so the same
  ## technology fitted either way puts the same quantity on the left-hand side.
  ## Under crs the two are reciprocal exactly, so the coefficients must match.
  d <- sw_dgp(100, c(1.05, 0.5, -0.25), seed = 21)
  a <- dea(d$x, d$y, rts = "crs", orientation = "in", slack = FALSE)
  b <- dea(d$x, d$y, rts = "crs", orientation = "out", slack = FALSE)
  expect_equal(unname(DEA:::.reg_delta(a)), unname(DEA:::.reg_delta(b)),
               tolerance = 1e-7)
  ra <- dea_reg(a, d$z, algorithm = "sw1", L2 = 100, seed = 3, progress = FALSE)
  rb <- dea_reg(b, d$z, algorithm = "sw1", L2 = 100, seed = 3, progress = FALSE)
  expect_equal(ra$coefficients, rb$coefficients, tolerance = 1e-6)
})

test_that("an intercept is added unless one was supplied", {
  d <- sw_dgp(80, c(1.05, 0.5, -0.25), seed = 31)
  fit <- dea(d$x, d$y, rts = "vrs", orientation = "out", slack = FALSE)
  a <- dea_reg(fit, d$z, algorithm = "sw1", L2 = 50, seed = 4, progress = FALSE)
  expect_identical(a$z_names[1], "(Intercept)")
  expect_length(a$coefficients, 4L)      # intercept, two slopes, sigma
  b <- dea_reg(fit, cbind(1, d$z), algorithm = "sw1", L2 = 50, seed = 4,
               progress = FALSE)
  expect_length(b$coefficients, 4L)
  expect_equal(unname(a$coefficients), unname(b$coefficients), tolerance = 1e-6)
})

test_that("column names and data frames select the same covariates", {
  d <- sw_dgp(80, c(1.05, 0.5, -0.25), seed = 33)
  fit <- dea(d$x, d$y, rts = "vrs", orientation = "out", slack = FALSE)
  df <- data.frame(a = d$z[, 1], b = d$z[, 2])
  m <- dea_reg(fit, d$z, algorithm = "sw1", L2 = 40, seed = 5, progress = FALSE)
  n <- dea_reg(fit, c("a", "b"), data = df, algorithm = "sw1", L2 = 40,
               seed = 5, progress = FALSE)
  o <- dea_reg(fit, ~ a + b, data = df, algorithm = "sw1", L2 = 40, seed = 5,
               progress = FALSE)
  expect_equal(unname(m$coefficients), unname(n$coefficients))
  expect_equal(unname(m$coefficients), unname(o$coefficients))
})

test_that("dea_reg refuses fits it has no procedure for", {
  d <- sw_dgp(60, c(1.05, 0.5, -0.25), seed = 41)
  expect_error(dea_reg(dea_sbm(d$x, d$y), d$z), "radial")
  expect_error(dea_reg(dea(d$x, d$y, rts = "fdh", slack = FALSE), d$z),
               "free disposal")
  expect_error(dea_reg(dea(d$x, d$y, super = TRUE, slack = FALSE), d$z),
               "super-efficiency")
  fit <- dea(d$x, d$y, rts = "vrs", orientation = "out", slack = FALSE)
  expect_error(dea_reg(fit, d$z[-1, ]), "rows but the fit has")
  expect_error(dea_reg(fit, d$z, alpha = 0), "in \\(0, 1\\)")
  ## a sample with too few inefficient DMUs cannot support the regression
  eff_only <- dea(d$x, d$x, rts = "crs", slack = FALSE)
  expect_error(dea_reg(eff_only, d$z), "inefficient")
})

test_that("print, summary and the accessors work", {
  d <- sw_dgp(80, c(1.05, 0.5, -0.25), seed = 51)
  fit <- dea(d$x, d$y, rts = "vrs", orientation = "out", slack = FALSE)
  r <- dea_reg(fit, d$z, algorithm = "sw2", L1 = 20, L2 = 60, seed = 6,
               progress = FALSE)
  expect_output(print(r), "Algorithm #2")
  expect_output(print(r), "LESS efficiency")
  expect_output(summary(r), "bias correction")
  expect_equal(nobs(r), r$n_used)
  expect_equal(coef(r), r$coefficients)
  expect_true(all(r$table$lower <= r$table$estimate + 1e-8 |
                  is.na(r$table$lower)))
})
