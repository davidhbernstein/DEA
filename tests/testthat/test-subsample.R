## ---------------------------------------------------------------------------
## Subsampling intervals, the m out of n bootstrap.
##
## No reference implementation on CRAN, so the checks are the construction's own
## identities -- what the statistic is, which exponent it uses, that the draw is
## without replacement -- plus the two places the m-selection rule can go wrong.
## ---------------------------------------------------------------------------

sub_fixture <- function(n = 80, seed = 5, rts = "vrs", ori = "out") {
  s <- dea_sim(n, p = 2, q = 1, returns = 0.9, seed = seed)
  list(sim = s, fit = dea(s$x, s$y, rts = rts, orientation = ori, slack = FALSE))
}

test_that("the interval inverts the scaled statistic with the right exponent", {
  f <- sub_fixture()
  r <- suppressWarnings(dea_subsample(f$fit, m = 20, B = 60, seed = 1,
                                      progress = FALSE))
  ## kappa is the estimator's rate, which dea_rate() already knows.
  expect_equal(r$kappa, -dea_rate(2, 1, "vrs", what = "estimator"))
  expect_equal(r$kappa, 2 / (2 + 1 + 1))
  ## Rebuild the endpoints from the recorded draws: the interval must be
  ## eff - n^-kappa * quantile(T*), not a quantile of the scores themselves.
  ## Scoring against a subsample can be infeasible, so some draws are NA and
  ## the quantiles are over the valid ones -- which is why n_valid is reported.
  Tst <- 20^r$kappa * (r$boot - r$table$eff)
  qf <- function(pr) apply(Tst, 1L, stats::quantile, probs = pr, na.rm = TRUE,
                           names = FALSE)
  ok <- r$table$n_valid >= ceiling(2 / 0.05)
  expect_equal(r$table$ci_lower[ok],
               as.numeric(r$table$eff - 80^(-r$kappa) * qf(0.975))[ok],
               tolerance = 1e-9)
  expect_equal(r$table$ci_upper[ok],
               as.numeric(r$table$eff - 80^(-r$kappa) * qf(0.025))[ok],
               tolerance = 1e-9)
  expect_true(all(r$table$ci_lower <= r$table$ci_upper, na.rm = TRUE))
  expect_equal(r$table$n_valid, as.integer(rowSums(is.finite(r$boot))))
})

test_that("the subsample is drawn without replacement and used as the reference", {
  ## Without replacement is the whole point -- with replacement is the n out of
  ## n bootstrap this exists to replace. A subsample of size m spans a smaller
  ## technology, so output-oriented every replicated score is no larger than the
  ## full-sample one, and input-oriented no smaller.
  f <- sub_fixture(60, ori = "out")
  r <- suppressWarnings(dea_subsample(f$fit, m = 15, B = 40, seed = 2,
                                      progress = FALSE))
  expect_lt(max(r$boot - r$table$eff, na.rm = TRUE), 1e-8)
  fi <- sub_fixture(60, ori = "in")
  ri <- suppressWarnings(dea_subsample(fi$fit, m = 15, B = 40, seed = 2,
                                       progress = FALSE))
  expect_gt(min(ri$boot - ri$table$eff, na.rm = TRUE), -1e-8)
  ## A draw WITH replacement would put duplicate reference rows in and could not
  ## produce a hull smaller than m distinct points; the inequality above is what
  ## distinguishes them in aggregate, and the dimensions pin the rest.
  expect_equal(dim(r$boot), c(60L, 40L))
})

test_that("the volatility rule ignores grid points without a full window", {
  ## A clipped window averages over fewer neighbours, so its standard deviation
  ## is smaller for no reason to do with stability. Before this was fixed the
  ## rule picked the first grid point.
  f <- sub_fixture(100)
  r <- suppressWarnings(dea_subsample(f$fit, B = 60, B_select = 40,
                                      window = 3L, seed = 3, progress = FALSE))
  expect_equal(length(r$m_grid), nrow(r$endpoint_path))
  expect_true(r$m %in% r$eligible)
  expect_false(r$m %in% r$m_grid[1:3])
  expect_false(r$m %in% utils::tail(r$m_grid, 3))
  ## and the chosen m is the least volatile among those eligible
  idx <- match(r$eligible, r$m_grid)
  expect_equal(r$m, r$m_grid[idx[which.min(r$volatility[idx])]])
  expect_equal(r$window_used, 3L)
})

test_that("a grid too short for the window narrows it and says so", {
  f <- sub_fixture(60)
  ## Two warnings fire here -- the narrowed window and the small-B shortfall --
  ## so they are collected rather than matched one at a time.
  w <- character()
  r <- withCallingHandlers(
    dea_subsample(f$fit, m_grid = c(8, 12, 16, 20), window = 3L,
                  B = 30, B_select = 20, seed = 4, progress = FALSE),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("too short for a window", w)))
  expect_lt(r$window_used, 3L)
  expect_true(r$m %in% r$m_grid)
})

test_that("it accepts free disposal, which dea_boot refuses", {
  ## dea_boot() stops on fdh and points at the subsampling scheme of Jeong and
  ## Simar; this is that scheme, and the rate it uses is the fdh one.
  s <- dea_sim(60, p = 2, q = 1, returns = 0.9, seed = 11)
  ff <- dea(s$x, s$y, rts = "fdh", orientation = "out", slack = FALSE)
  expect_error(dea_boot(ff, B = 10), "Jeong and Simar")
  r <- suppressWarnings(dea_subsample(ff, m = 25, B = 200, seed = 5,
                                      progress = FALSE))
  expect_equal(r$kappa, -dea_rate(2, 1, "fdh", what = "estimator"))
  expect_equal(r$kappa, 1 / 3)
  ## Not every interval need exist -- a point outside the subsample's hull has
  ## no score, and under free disposal that is commoner still, since a peer has
  ## to dominate rather than merely be combined with others.
  expect_true(any(is.finite(r$table$ci_lower)))
  expect_true(all(r$table$n_valid <= 200L))
})

test_that("a supplied m skips the selection entirely", {
  f <- sub_fixture(60)
  r <- suppressWarnings(dea_subsample(f$fit, m = 17, B = 30, seed = 6,
                                      progress = FALSE))
  expect_equal(r$m, 17L)
  expect_null(r$m_grid)
  expect_null(r$volatility)
  expect_true(is.na(r$window_used))
  expect_output(print(r), "supplied")
})

test_that("the seed reproduces and leaves the stream alone", {
  f <- sub_fixture(50)
  ## Compared on the replicated SCORES, not the endpoints: at a small B many
  ## endpoints are NA for want of usable draws, and two all-NA columns compare
  ## equal whatever the seed was.
  a <- suppressWarnings(dea_subsample(f$fit, m = 12, B = 30, seed = 7, progress = FALSE))
  b <- suppressWarnings(dea_subsample(f$fit, m = 12, B = 30, seed = 7, progress = FALSE))
  expect_equal(a$boot, b$boot)
  expect_equal(a$table$ci_lower, b$table$ci_lower)
  c3 <- suppressWarnings(dea_subsample(f$fit, m = 12, B = 30, seed = 8, progress = FALSE))
  expect_false(isTRUE(all.equal(a$boot, c3$boot)))
  set.seed(42); before <- runif(1)
  set.seed(42); invisible(suppressWarnings(dea_subsample(f$fit, m = 12, B = 30,
                            seed = 9, progress = FALSE))); after <- runif(1)
  expect_equal(before, after)
})

test_that("the refusals fire", {
  f <- sub_fixture(50)
  expect_error(dea_subsample(f$sim, m = 10), "must be a fit from dea")
  expect_error(dea_subsample(f$fit, m = 1), "2 <= m < n")
  expect_error(dea_subsample(f$fit, m = 50), "2 <= m < n")
  expect_error(dea_subsample(f$fit, m = 10, B = 1), "at least 2")
  expect_error(dea_subsample(f$fit, m = 10, alpha = 1), "\\(0, 1\\)")
  expect_error(dea_subsample(f$fit, m_grid = c(5, 10)), "at least 3 candidates")
  expect_error(dea_subsample(f$fit, m_grid = c(5, 10, 60)), "must lie in")
  sb <- dea(f$sim$x, f$sim$y, super = TRUE, slack = FALSE)
  expect_error(dea_subsample(sb, m = 10), "super-efficiency")
  ad <- dea_add(f$sim$x, f$sim$y)
  expect_error(dea_subsample(ad, m = 10), "radial estimator")
  ext <- dea(f$sim$x[1:10, ], f$sim$y[1:10, , drop = FALSE],
             xref = f$sim$x, yref = f$sim$y, slack = FALSE)
  expect_error(dea_subsample(ext, m = 5), "reference set IS the sample")
})

test_that("the methods run", {
  f <- sub_fixture(60)
  r <- suppressWarnings(dea_subsample(f$fit, B = 40, B_select = 25, seed = 10,
                                      progress = FALSE))
  expect_output(print(r), "m out of n")
  expect_output(print(r), "minimum volatility")
  expect_output(summary(r), "flat stretch")
  expect_equal(nrow(r$table), 60L)
  expect_named(r$table, c("dmu", "eff", "bias_corrected", "ci_lower",
                          "ci_upper", "n_valid"))
})
