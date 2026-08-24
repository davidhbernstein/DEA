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
