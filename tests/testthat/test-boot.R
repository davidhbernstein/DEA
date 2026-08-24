test_that("the bootstrap corrects a bias in the right direction", {
  d <- toy(n = 30, p = 1, q = 1, seed = 7)
  fit <- dea(d$x, d$y, rts = "vrs", orientation = "in")
  b <- dea_boot(fit, B = 60, seed = 1, progress = FALSE)

  expect_s3_class(b, "dea_boot")
  expect_equal(dim(b$boot), c(30L, 60L))
  ## The estimated frontier is inside the true one, so theta_hat is biased
  ## TOWARD 1 and the bootstrap replications sit above the original scores.
  expect_true(mean(b$table$bias) > 0)
  expect_true(all(b$table$bias_corrected <= b$table$eff + 1e-9))
  ## The interval brackets the corrected estimate by construction.
  expect_true(all(b$table$ci_lower <= b$table$bias_corrected + 1e-8))
  expect_true(all(b$table$ci_upper >= b$table$bias_corrected - 1e-8))
  expect_true(all(b$table$se >= 0))
})

test_that("a seed makes the bootstrap reproducible and leaves the RNG alone", {
  d <- toy(n = 25, p = 1, q = 1, seed = 8)
  fit <- dea(d$x, d$y, rts = "vrs")
  set.seed(99); before <- runif(1)
  set.seed(99)
  b1 <- dea_boot(fit, B = 20, seed = 42, progress = FALSE)
  after <- runif(1)
  b2 <- dea_boot(fit, B = 20, seed = 42, progress = FALSE)
  expect_equal(b1$table$bias, b2$table$bias)
  ## The caller's stream is where they left it, not where the bootstrap left it.
  expect_equal(before, after)
})

test_that("dea_boot() refuses fits it is not defined for", {
  d <- toy(n = 30, p = 1, q = 1, seed = 7)
  expect_error(dea_boot(dea(d$x, d$y, rts = "fdh"), B = 5), "Jeong and Simar")
  expect_error(dea_boot(dea(d$x, d$y, super = TRUE), B = 5), "super-efficiency")
  expect_error(dea_boot(dea_sbm(d$x, d$y), B = 5), "sbm fit")
  expect_error(dea_boot(dea(d$x, d$y, xref = d$x, yref = d$y), B = 5),
               "reference set IS the sample")
  expect_error(dea_boot(fit <- dea(d$x, d$y), B = 1), "at least 2")
})

test_that("output orientation reflects on the correct side of 1", {
  d <- toy(n = 30, p = 1, q = 1, seed = 7)
  fit <- dea(d$x, d$y, rts = "vrs", orientation = "out")
  b <- dea_boot(fit, B = 40, seed = 2, progress = FALSE)
  ## phi >= 1, so the bias runs the other way: phi_hat understates phi.
  expect_true(mean(b$table$bias) < 0)
  expect_true(all(b$table$bias_corrected >= b$table$eff - 1e-9))
  ## Bootstrap scores may cross 1: the pseudo-technology is a shrunken copy
  ## of the estimated one, so an original observation need not lie inside it.
  ## That tail is the interval, not an error -- clamping it would truncate the
  ## very thing being measured.
  expect_true(any(b$boot < 1))
  expect_true(all(is.finite(b$boot)))
})
