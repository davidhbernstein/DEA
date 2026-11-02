## ---------------------------------------------------------------------------
## The Simar-Wilson (2002) test of returns to scale.
##
## The statistic is deterministic and is pinned against rDEA and npsf; the
## p-value is not, so what is checked about it is its CONTRACT -- that it is a
## valid Monte Carlo p-value, that it moves the right way, and that the
## left-tailed convention holds in both orientations.
## ---------------------------------------------------------------------------

test_that("the observed statistic matches rDEA and npsf", {
  ref <- ref_values("rts_test_reference.csv")
  for (i in seq_len(nrow(ref))) {
    d <- toy(ref$n[i], 2, 1, ref$seed[i], returns = ref$returns[i])
    r <- suppressWarnings(dea_rts_test(d$x, d$y, orientation = ref$orientation[i],
                                       h0 = ref$h0[i], B = 2, seed = 1,
                                       progress = FALSE))
    expect_equal(r$statistic[["ratio_of_means"]], ref$ratio_of_means[i],
                 tolerance = 1e-9)
    expect_equal(r$statistic[["mean_of_ratios"]], ref$mean_of_ratios[i],
                 tolerance = 1e-9)
  }
})

test_that("the statistic is the aggregation of dea_rts()'s own scale efficiency", {
  ## mean_of_ratios is literally mean(scale_eff), in both orientations, and the
  ## per-DMU column is the same vector. If these ever drift apart, one of the
  ## two functions has changed its mind about what scale efficiency is.
  d <- toy(40, 2, 1, 5, returns = 0.7)
  for (ori in c("in", "out")) {
    g <- dea_rts(d$x, d$y, orientation = ori)$table
    r <- suppressWarnings(dea_rts_test(d$x, d$y, orientation = ori, B = 2,
                                       seed = 1, progress = FALSE))
    expect_equal(r$table$scale_eff, g$scale_eff, tolerance = 1e-10)
    expect_equal(r$statistic[["mean_of_ratios"]], mean(g$scale_eff),
                 tolerance = 1e-10)
    ## And the statistic is on the (0, 1] scale whichever orientation it is.
    expect_lte(r$statistic[["ratio_of_means"]], 1 + 1e-10)
    expect_lte(r$statistic[["mean_of_ratios"]], 1 + 1e-10)
  }
})

test_that("the nirs null sits between the crs null and vrs", {
  ## crs is a restriction of nirs is a restriction of vrs, so the scores are
  ## ordered and so are the statistics. This is an identity of the programs and
  ## does not depend on the data.
  d <- toy(40, 2, 1, 17, returns = 0.6)
  a <- suppressWarnings(dea_rts_test(d$x, d$y, h0 = "crs",  B = 2, seed = 1, progress = FALSE))
  b <- suppressWarnings(dea_rts_test(d$x, d$y, h0 = "nirs", B = 2, seed = 1, progress = FALSE))
  expect_true(all(a$eff_h0 <= b$eff_h0 + 1e-9))
  expect_true(all(b$eff_h0 <= b$eff_vrs + 1e-9))
  expect_lte(a$statistic[["ratio_of_means"]], b$statistic[["ratio_of_means"]] + 1e-9)
  expect_lte(b$statistic[["ratio_of_means"]], 1 + 1e-9)
})

test_that("the p-value is a valid Monte Carlo p-value", {
  d <- toy(30, 2, 1, 3, returns = 0.8)
  r <- suppressWarnings(dea_rts_test(d$x, d$y, orientation = "out", B = 39,
                                     seed = 4, progress = FALSE))
  ## (count + 1)/(B + 1): never 0, never above 1, and at least 1/(B+1).
  expect_true(all(r$p_value >= 1 / (r$B + 1) - 1e-12))
  expect_true(all(r$p_value <= 1))
  expect_true(all(r$table$p_value >= 1 / (r$B + 1) - 1e-12))
  expect_equal(dim(r$boot), c(39L, 2L))
  expect_equal(dim(r$boot_ratio), c(30L, 39L))
  expect_equal(unname(r$reject), unname(r$p_value < r$alpha))
  expect_equal(nobs(r), 30L)
  ## The p-value really is the lower-tail fraction of the recorded draws.
  expect_equal(r$p_value[["ratio_of_means"]],
               (sum(r$boot[, 1] <= r$statistic[["ratio_of_means"]]) + 1) / (r$B + 1))
})

test_that("the seed makes it reproducible and does not disturb the stream", {
  d <- toy(30, 2, 1, 6, returns = 0.8)
  f <- function(seed) suppressWarnings(
    dea_rts_test(d$x, d$y, orientation = "out", B = 19, seed = seed, progress = FALSE))
  expect_equal(f(11)$p_value, f(11)$p_value)
  expect_false(isTRUE(all.equal(f(11)$boot[, 1], f(12)$boot[, 1])))
  set.seed(99); before <- runif(1)
  set.seed(99); invisible(f(5)); after <- runif(1)
  expect_equal(before, after)
})

test_that("a true cone is not rejected output-oriented, and a DRS design is", {
  ## The only check here that can catch a MISUNDERSTANDING rather than a coding
  ## error: dea_sim(returns = 1) is a cone, so H0 is exactly true, and
  ## returns = 0.5 is as far from one as this generator goes.
  cone <- dea_sim(80, p = 2, q = 1, returns = 1.0, seed = 11)
  drs  <- dea_sim(80, p = 2, q = 1, returns = 0.5, seed = 11)
  a <- suppressWarnings(dea_rts_test(cone$x, cone$y, orientation = "out",
                                     B = 199, seed = 2, progress = FALSE))
  b <- suppressWarnings(dea_rts_test(drs$x, drs$y, orientation = "out",
                                     B = 199, seed = 2, progress = FALSE))
  expect_gt(a$p_value[["mean_of_ratios"]], 0.05)
  expect_lt(b$p_value[["mean_of_ratios"]], 0.05)
  ## And the statistic itself is nearer 1 under the cone than under the DRS
  ## technology, which is the direction the whole test rests on.
  expect_gt(a$statistic[["mean_of_ratios"]], b$statistic[["mean_of_ratios"]])
})

test_that("the support gate fires where the vrs program provably cannot follow", {
  ## Input-oriented on a cone: a third to a half of the sample has a crs
  ## projection below the smallest observed input, so the statistic is partly a
  ## constant and the test is invalid. The warning is the only notice a caller
  ## gets, so it is asserted rather than assumed.
  cone <- dea_sim(60, p = 2, q = 1, returns = 1.0, seed = 11)
  w <- character()
  r <- withCallingHandlers(
    dea_rts_test(cone$x, cone$y, orientation = "in", B = 19, seed = 1,
                 progress = FALSE),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("provably cannot reach", w)))
  expect_true(any(grepl("resolution", w)))         ## the small-B notice
  expect_true(mean(!r$reachable) > 0.25)

  ## The flag is exact, not a rule of thumb: reproduce it from the definition.
  X <- as.matrix(cone$x)
  expect_equal(r$reachable,
               !apply(sweep(X * r$eff_h0, 2L, apply(X, 2L, min), "<"), 1L, any))

  ## Output-oriented on the same data, hardly anything is flagged.
  o <- suppressWarnings(dea_rts_test(cone$x, cone$y, orientation = "out",
                                     B = 19, seed = 1, progress = FALSE))
  expect_lt(mean(!o$reachable), 0.25)
})

test_that("arguments are validated and technology aliases are honoured", {
  d <- toy(30, 2, 1, 8, returns = 0.8)
  expect_error(dea_rts_test(d$x, d$y, B = 1), "at least 2")
  expect_error(dea_rts_test(d$x, d$y, alpha = 0), "\\(0, 1\\)")
  expect_error(dea_rts_test(d$x, d$y, alpha = 1), "\\(0, 1\\)")
  expect_error(dea_rts_test(d$x, d$y, h0 = "vrs"), "h0")
  expect_error(dea_rts_test(d$x, d$y, orientation = "both"), "orientation")
  expect_error(dea_rts_test(d$x, d$y, "crs"), "`data` must be")
  expect_warning(dea_rts_test(d$x, d$y, orientation = "out", B = 5,
                              seed = 1, progress = FALSE), "resolution")

  ## "drs" is Benchmarking's spelling of sum(lambda) <= 1 and must select the
  ## same technology that "nirs" does, not raise an error.
  a <- suppressWarnings(dea_rts_test(d$x, d$y, h0 = "nirs", B = 2, seed = 1, progress = FALSE))
  b <- suppressWarnings(dea_rts_test(d$x, d$y, h0 = "DRS",  B = 2, seed = 1, progress = FALSE))
  expect_equal(a$statistic, b$statistic)
  expect_equal(b$h0, "nirs")
})

test_that("the simultaneous level is Sidak and the print methods run", {
  d <- toy(30, 2, 1, 9, returns = 0.7)
  r <- suppressWarnings(dea_rts_test(d$x, d$y, orientation = "out", B = 39,
                                     alpha = 0.05, seed = 1, progress = FALSE))
  expect_equal(r$alpha_individual, 1 - (1 - 0.05)^(1 / 30))
  expect_equal(r$table$scale_efficient, r$table$p_value > r$alpha_individual)
  expect_output(print(r), "Simar-Wilson test of returns to scale")
  expect_output(print(r), "rejects for SMALL values")
  expect_output(summary(r), "per-DMU scale efficiency")
})

test_that("scaling the columns leaves the test alone", {
  ## Radial efficiency is units invariant, so the statistic must be too --
  ## including when the caller turns the internal rescaling off.
  d <- toy(30, 2, 1, 12, returns = 0.75)
  a <- suppressWarnings(dea_rts_test(d$x, d$y, orientation = "out", B = 9,
                                     seed = 3, progress = FALSE, scaling = TRUE))
  b <- suppressWarnings(dea_rts_test(d$x, d$y, orientation = "out", B = 9,
                                     seed = 3, progress = FALSE, scaling = FALSE))
  expect_equal(a$statistic, b$statistic, tolerance = 1e-9)
  x2 <- d$x; x2[, 1] <- x2[, 1] * 1000
  c3 <- suppressWarnings(dea_rts_test(x2, d$y, orientation = "out", B = 9,
                                      seed = 3, progress = FALSE))
  expect_equal(a$statistic, c3$statistic, tolerance = 1e-9)
})

## ---------------------------------------------------------------------------
## method = "subsample": the m out of n test of Simar and Wilson (2011).
## ---------------------------------------------------------------------------

test_that("subsampling leaves the observed statistic alone and is a valid test", {
  d <- toy(40, 2, 1, 21, returns = 0.8)
  a <- suppressWarnings(dea_rts_test(d$x, d$y, orientation = "out", B = 39,
                                     seed = 2, progress = FALSE))
  b <- suppressWarnings(dea_rts_test(d$x, d$y, orientation = "out", B = 39,
                                     method = "subsample", m = 15, seed = 2,
                                     progress = FALSE))
  ## Only the null distribution changes; the statistic and the table's scores
  ## are the same functional of the same data.
  expect_equal(a$statistic, b$statistic)
  expect_equal(a$table$scale_eff, b$table$scale_eff)
  expect_identical(b$method, "subsample")
  expect_identical(b$m, 15L)

  ## kappa is the VRS rate whatever h0 is: the statistic is a contrast, and a
  ## contrast converges no faster than its slower term.
  expect_equal(b$kappa, 2 / (2 + 1 + 1))
  ## tau is n^kappa (1 - S), with no recentring.
  expect_equal(b$tau, 40^b$kappa * (1 - b$statistic))
  expect_equal(b$boot_tau, 15^b$kappa * (1 - b$boot))

  ## UPPER tail of tau -- large tau is small S -- with the observed value
  ## counted as a draw.
  expect_equal(b$p_value[["mean_of_ratios"]],
               (sum(b$boot_tau[, 2] >= b$tau[["mean_of_ratios"]]) + 1) / (b$B + 1))
  expect_true(all(b$p_value >= 1 / (b$B + 1) - 1e-12 & b$p_value <= 1))
  expect_equal(unname(b$reject), unname(b$p_value < b$alpha))

  ## No per-DMU tests, and they say so rather than inventing a number.
  expect_true(all(is.na(b$table$p_value)))
  expect_true(all(is.na(b$table$scale_efficient)))
  expect_output(print(b), "m out of n subsampling")
  expect_output(print(b), "not available under subsampling")
  expect_output(print(a), "smoothed bootstrap")
})

test_that("subsample statistics are the statistic of the subsample itself", {
  ## A draw must be the same functional of a smaller sample -- the subsample is
  ## both the evaluated and the reference set. Rebuild one draw by hand.
  d <- toy(30, 2, 1, 22, returns = 0.7)
  b <- suppressWarnings(dea_rts_test(d$x, d$y, orientation = "out", B = 3,
                                     method = "subsample", m = 12, seed = 5,
                                     progress = FALSE))
  set.seed(5)
  idx <- sample.int(30, 12, replace = FALSE)
  g <- dea_rts(d$x[idx, , drop = FALSE], d$y[idx, , drop = FALSE],
               orientation = "out")$table
  expect_equal(unname(b$boot[1, "mean_of_ratios"]), mean(g$scale_eff),
               tolerance = 1e-8)
})

test_that("the subsample size is chosen by minimum volatility of the critical value", {
  d <- toy(60, 2, 1, 23, returns = 0.8)
  b <- suppressWarnings(dea_rts_test(d$x, d$y, orientation = "out", B = 39,
                                     method = "subsample", B_select = 20,
                                     seed = 3, progress = FALSE))
  expect_true(b$m %in% b$m_grid)
  expect_equal(rownames(b$critical_path), as.character(b$m_grid))
  ## Only full-window candidates are eligible, so the chosen m is never one of
  ## the `window` grid points at either end.
  k <- length(b$m_grid); w <- 3L
  expect_true(match(b$m, b$m_grid) > w && match(b$m, b$m_grid) <= k - w)
  i <- (w + 1L):(k - w)
  expect_equal(b$m, b$m_grid[i[which.min(b$volatility[i])]])
  expect_output(print(b), "minimum volatility")

  ## Reproducible, including the selection step.
  b2 <- suppressWarnings(dea_rts_test(d$x, d$y, orientation = "out", B = 39,
                                      method = "subsample", B_select = 20,
                                      seed = 3, progress = FALSE))
  expect_equal(b$p_value, b2$p_value)
  expect_identical(b$m, b2$m)
})

test_that("subsampling arguments are validated", {
  d <- toy(30, 2, 1, 24, returns = 0.8)
  f <- function(...) dea_rts_test(d$x, d$y, orientation = "out", B = 200,
                                  progress = FALSE, ...)
  expect_error(f(method = "jackknife"), "method")
  expect_error(f(m = 10), "only to method")
  expect_error(f(method = "subsample", m = 2), "3 <= m < n")
  expect_error(f(method = "subsample", m = 30), "3 <= m < n")
  expect_error(f(method = "subsample", m = NA_real_), "3 <= m < n")
  expect_error(f(method = "subsample", m_grid = c(5, 8)), "at least 3")
})

test_that("subsampling rejects a far-from-cone technology and not a cone", {
  ## The external check: dea_sim() makes H0 exactly true at returns = 1. One
  ## dataset each way, so this is a smoke test of direction, not of size --
  ## horserace/rts_subsample_experiment.R measures size and power.
  cone <- dea_sim(100, p = 2, q = 1, returns = 1.0, seed = 3)
  drs  <- dea_sim(100, p = 2, q = 1, returns = 0.5, seed = 3)
  f <- function(d) suppressWarnings(dea_rts_test(d$x, d$y, orientation = "out",
                                                 method = "subsample", m = 25,
                                                 B = 199, seed = 1,
                                                 progress = FALSE))
  expect_gt(f(cone)$p_value[["mean_of_ratios"]], 0.05)
  expect_lt(f(drs)$p_value[["mean_of_ratios"]], 0.05)
})
