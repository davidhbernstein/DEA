test_that("charnes1981 has the shape and provenance its help page claims", {
  data(charnes1981, envir = environment())
  expect_s3_class(charnes1981, "data.frame")
  expect_equal(dim(charnes1981), c(70L, 10L))
  expect_named(charnes1981, c("site", "pft", paste0("x", 1:5), paste0("y", 1:3)))
  expect_type(charnes1981$pft, "logical")
  expect_equal(sum(charnes1981$pft), 49L)
  ## The two groups are contiguous, rows 1-49 and 50-70, as distributed.
  expect_equal(range(which(charnes1981$pft)), c(1L, 49L))
  expect_equal(range(which(!charnes1981$pft)), c(50L, 70L))
  expect_length(unique(charnes1981$site), 38L)
  ## Every measurement is strictly positive, so the whole data set is usable by
  ## every model in the package including the slacks-based measure.
  num <- as.matrix(charnes1981[, paste0(rep(c("x", "y"), c(5, 3)), c(1:5, 1:3))])
  expect_true(all(is.finite(num)))
  expect_true(all(num > 0))
})

test_that("charnes1981 matches the copies already on CRAN", {
  ## The whole claim of the help page's Source section, asserted. If either
  ## upstream copy is ever corrected, this fails and the discrepancy gets
  ## looked at rather than quietly persisting.
  skip_if_not_installed("Benchmarking")
  data(charnes1981, envir = environment())
  e <- new.env(); utils::data("charnes1981", package = "Benchmarking", envir = e)
  up <- get("charnes1981", envir = e)
  for (v in c(paste0("x", 1:5), paste0("y", 1:3))) {
    expect_equal(charnes1981[[v]], up[[v]], tolerance = 0, info = v)
  }
  expect_equal(charnes1981$pft, up$pft == 1L)
})

test_that("the original CCR analysis runs and gives a sensible answer", {
  data(charnes1981, envir = environment())
  x <- charnes1981[, paste0("x", 1:5)]
  y <- charnes1981[, paste0("y", 1:3)]

  ccr <- dea(x, y, rts = "crs", orientation = "in")
  expect_length(ccr$eff, 70L)
  expect_true(all(ccr$eff > 0 & ccr$eff <= 1))
  ## Eight dimensions on seventy units: a large efficient set is expected, and
  ## a small one would mean something is wrong.
  expect_gt(sum(ccr$eff == 1), 10L)
  expect_lt(sum(ccr$eff == 1), 45L)

  ## The three technologies nest on real data as they do on simulated data.
  bcc <- dea(x, y, rts = "vrs", orientation = "in", slack = FALSE)
  expect_true(all(ccr$eff <= bcc$eff + 1e-9))

  r <- dea_rts(x, y, orientation = "in")
  expect_true(all(r$table$scale_eff > 0 & r$table$scale_eff <= 1 + 1e-9))

  ## And it is a real data set, so the slacks-based measure works on it too.
  rho <- dea_sbm(x, y, rts = "vrs")
  expect_true(all(rho$eff > 0 & rho$eff <= 1 + 1e-9))
})
