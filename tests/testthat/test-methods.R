test_that("the extractors return what they say they do", {
  d <- toy(n = 40, p = 2, q = 1)
  fit <- dea(d$x, d$y, rts = "vrs", orientation = "in")

  expect_equal(nobs(fit), 40L)
  expect_equal(efficiency(fit), fit$eff)
  expect_equal(efficiency(fit, "score"), fit$eff)

  fo <- dea(d$x, d$y, rts = "vrs", orientation = "out")
  expect_equal(unname(efficiency(fo, "score")), unname(1 / fo$eff))
  expect_true(all(efficiency(fo, "score") <= 1 + 1e-9))

  p <- peers(fit)
  expect_named(p, c("dmu", "peer", "lambda"))
  expect_true(all(p$lambda > 0))
  ## Every DMU appears, since every one has at least one peer.
  expect_setequal(unique(p$dmu), fit$dmu)

  s <- slacks(fit)
  expect_equal(dim(s), c(40L, 3L))
  expect_true(all(s >= 0))
})

test_that("fitted() lands on the frontier", {
  d <- toy(n = 40, p = 2, q = 1)
  fit <- dea(d$x, d$y, rts = "vrs", orientation = "in", slack = TRUE)
  f <- fitted(fit)
  expect_equal(dim(f), c(40L, 3L))
  ## Scoring the projections against the same technology must return 1.
  again <- dea(f[, 1:2], f[, 3, drop = FALSE], rts = "vrs", orientation = "in",
               slack = FALSE, xref = d$x, yref = d$y)
  expect_equal(unname(again$eff), rep(1, 40), tolerance = 1e-6)
})

test_that("print, summary and plot run for every model", {
  d <- toy(n = 40, p = 1, q = 1)
  fits <- list(dea(d$x, d$y, rts = "vrs"),
               dea(d$x, d$y, rts = "fdh"),
               dea(d$x, d$y, rts = "crs", super = TRUE),
               dea_sbm(d$x, d$y, rts = "vrs"),
               dea_ddf(d$x, d$y, rts = "vrs"))
  pf <- tempfile(fileext = ".pdf")
  grDevices::pdf(pf)
  on.exit({ grDevices::dev.off(); unlink(pf) }, add = TRUE)
  for (f in fits) {
    expect_output(print(f), "Data envelopment analysis")
    expect_output(summary(f), "efficiency by DMU")
    expect_silent(plot(f))
  }
  ## Two dimensions or more: falls back to the score distribution.
  d2 <- toy(n = 40, p = 2, q = 2)
  expect_silent(plot(dea(d2$x, d2$y, rts = "vrs")))
})

test_that("the frontier plot traces the technology, not the efficient points", {
  d <- toy(n = 40, p = 1, q = 1)
  pf <- tempfile(fileext = ".pdf"); grDevices::pdf(pf)
  on.exit({ grDevices::dev.off(); unlink(pf) }, add = TRUE)
  fit <- dea(d$x, d$y, rts = "vrs")
  fr <- plot(fit)
  ## The traced frontier is non-decreasing (free disposability) and lies weakly
  ## above every observed point at the same input level.
  expect_true(all(diff(fr$y) >= -1e-8))
  for (i in seq_len(nrow(d$x))) {
    j <- which.min(abs(fr$x - d$x[i, 1]))
    expect_gte(fr$y[j] + 1e-6, d$y[i, 1])
  }
})

test_that("extractors refuse what was never computed", {
  d <- toy(n = 40)
  expect_error(peers(dea(d$x, d$y, peers = FALSE)), "peers = FALSE")
  expect_error(slacks(dea(d$x, d$y, slack = FALSE)), "slack = TRUE")
  expect_error(slacks(dea_ddf(d$x, d$y)), "along g")
})

test_that("print() says when the reference set is external", {
  ref <- dea_sim(50, p = 2, q = 1, seed = 41)
  ev  <- dea_sim(12, p = 2, q = 1, seed = 42)
  expect_output(print(dea(ev$x, ev$y, rts = "vrs", xref = ref$x, yref = ref$y)),
                "EXTERNAL reference set of 50")
  ## And says nothing extra for an ordinary fit.
  out <- utils::capture.output(print(dea(ref$x, ref$y, rts = "vrs")))
  expect_false(any(grepl("EXTERNAL", out)))
})

## ---------------------------------------------------------------------------
## plot() on the classes that are not "dea".
##
## Until these methods existed, plot() on such an object did one of two wrong
## things. `dea_price`, `dea_cross` and `dea_sim` all carry components named
## `x` and `y` -- the input and output matrices -- so plot.default() FOUND them
## and silently drew the raw inputs against the raw outputs: a plausible
## picture of nothing to do with the model. `dea_rts` and `dea_boot` have no
## such components, so the same fall-through failed with "'x' is a list, but
## does not have components 'x' and 'y'".
## ---------------------------------------------------------------------------

## Draw to a null device: the assertions are about dispatch and the returned
## data, not about pixels.
draw <- function(expr) {
  grDevices::pdf(NULL)
  on.exit({ grDevices::dev.off() }, add = TRUE)
  force(expr)
}

test_that("plot() dispatches to a real method for every class", {
  skip_if_not_installed("grDevices")
  d  <- toy(n = 30, p = 2, q = 2, returns = 0.9, seed = 5)
  pw <- matrix(seq(1, 2, length.out = 60), 30, 2)

  ## Each must dispatch to OUR method, not to plot.default. The returned value
  ## is the plotted data, which plot.default does not give.
  r <- draw(plot(dea_rts(d$x, d$y)))
  expect_true(is.data.frame(r))
  expect_true(all(c("sum_lambda_crs", "scale_eff", "rts") %in% names(r)))

  cx <- draw(plot(dea_cross(d$x, d$y, rts = "crs")))
  expect_true(is.data.frame(cx))
  expect_true(all(c("cross", "lo", "hi") %in% names(cx)))
  ## The range really does bracket the mean appraisal, which is what the plot
  ## claims by drawing the point on the segment.
  expect_true(all(cx$lo <= cx$cross + 1e-9))
  expect_true(all(cx$hi >= cx$cross - 1e-9))

  pc <- draw(plot(dea_cost(d$x, d$y, pw)))
  expect_true(is.data.frame(pc))
  ## The decomposition the plot draws must be the one the object reports.
  expect_equal(pc$technical * pc$allocative, pc$overall, tolerance = 1e-8)

  sm <- dea_sim(25, p = 1, q = 1, returns = 0.9, seed = 5)
  s1 <- draw(plot(sm))
  expect_true(is.data.frame(s1))                 ## p = q = 1: the frontier plot
  s2 <- draw(plot(dea_sim(25, p = 2, q = 2, returns = 0.9, seed = 5)))
  expect_type(s2, "double")                      ## otherwise: the theta histogram

  fit <- dea(d$x, d$y, rts = "vrs")
  b <- suppressWarnings(dea_boot(fit, B = 20, seed = 1, progress = FALSE))
  bb <- draw(plot(b))
  expect_true(is.data.frame(bb))
  ## Sorted by the bias-corrected estimate, and the interval brackets it --
  ## both are what the caterpillar is asserting visually.
  expect_false(is.unsorted(bb$bias_corrected))
  expect_true(all(bb$ci_lower <= bb$bias_corrected + 1e-8))
  expect_true(all(bb$ci_upper >= bb$bias_corrected - 1e-8))
  expect_identical(nrow(draw(plot(b, sort = FALSE))), nrow(bb))
})

test_that("plot() survives the NA cases these objects are documented to produce", {
  skip_if_not_installed("grDevices")
  ## A rated DMU with all inputs zero receives no usable appraisal, so its
  ## column of the cross matrix is entirely NA. The plot must drop it rather
  ## than fail on an all-NA range.
  set.seed(3)
  xr <- cbind(runif(12, 1, 3), runif(12, 1, 3)); yr <- cbind(runif(12, 1, 3))
  cx <- suppressWarnings(dea_cross(rbind(c(2, 2), c(0, 0)), rbind(2, 1),
                                   xref = xr, yref = yr, rts = "crs",
                                   secondary = "benevolent"))
  out <- draw(plot(cx))
  expect_identical(nrow(out), 1L)          ## the appraisable DMU only
  expect_true(all(is.finite(out$cross)))

  ## Super-efficiency under vrs returns NA for DMUs with no dominating peer.
  d <- toy(n = 30, p = 2, q = 2, returns = 0.9, seed = 5)
  f <- suppressWarnings(dea(d$x, d$y, rts = "vrs", super = TRUE))
  skip_if(!any(is.na(f$eff)), "no infeasible super-efficiency DMU on this sample")
  expect_silent(draw(plot(f)))
})
