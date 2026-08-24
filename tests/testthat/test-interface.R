test_that("x and y can be given four different ways, all equivalent", {
  s <- dea_sim(40, p = 2, q = 1, seed = 1)
  d <- s$data
  ref <- dea(s$x, s$y, rts = "vrs", slack = FALSE)

  by_df   <- dea(d[, 1:2], d[, 3, drop = FALSE], rts = "vrs", slack = FALSE)
  by_name <- dea(c("x1", "x2"), "y1", data = d, rts = "vrs", slack = FALSE)
  by_form <- dea(~ x1 + x2, ~ y1, data = d, rts = "vrs", slack = FALSE)
  for (f in list(by_df, by_name, by_form)) {
    expect_equal(unname(f$eff), unname(ref$eff), tolerance = 1e-10)
  }
  ## A single input as a bare vector.
  s1 <- dea_sim(40, p = 1, q = 1, seed = 1)
  expect_equal(unname(dea(s1$x[, 1], s1$y[, 1], rts = "vrs", slack = FALSE)$eff),
               unname(dea(s1$x, s1$y, rts = "vrs", slack = FALSE)$eff))
})

test_that("row names become DMU names throughout", {
  s <- dea_sim(20, p = 1, q = 1, seed = 1)
  x <- s$x; y <- s$y
  rownames(x) <- rownames(y) <- paste0("firm", 1:20)
  fit <- suppressWarnings(dea(x, y, rts = "vrs"))
  expect_equal(names(fit$eff), paste0("firm", 1:20))
  expect_equal(rownames(fit$lambda), paste0("firm", 1:20))
  expect_equal(colnames(fit$lambda), paste0("firm", 1:20))
  expect_true(all(peers(fit)$peer %in% paste0("firm", 1:20)))
})

test_that("returns-to-scale names are matched loosely and aliased", {
  d <- toy(n = 40)
  e <- function(r) dea(d$x, d$y, rts = r, slack = FALSE)$eff
  expect_equal(e("VRS"), e("vrs"))
  expect_equal(e("Crs"), e("crs"))
  ## Benchmarking's spellings for the two one-sided technologies.
  expect_equal(e("drs"), e("nirs"))
  expect_equal(e("irs"), e("ndrs"))
  expect_error(dea(d$x, d$y, rts = "quadratic"), "must be one of")
})

test_that("an explicit reference set equal to the sample is a no-op", {
  d <- toy(n = 40, p = 2, q = 1)
  for (f in list(dea, dea_sbm, dea_ddf)) {
    a <- suppressWarnings(f(d$x, d$y, rts = "vrs"))
    b <- suppressWarnings(f(d$x, d$y, rts = "vrs", xref = d$x, yref = d$y))
    expect_equal(unname(a$eff), unname(b$eff), tolerance = 1e-10)
  }
  a <- dea(d$x, d$y, rts = "fdh"); b <- dea(d$x, d$y, rts = "fdh", xref = d$x, yref = d$y)
  expect_equal(unname(a$eff), unname(b$eff))
})

test_that("a smaller reference set can only shrink the technology", {
  d <- toy(n = 60, p = 1, q = 1)
  full <- dea(d$x, d$y, rts = "vrs", orientation = "in", slack = FALSE)
  half <- dea(d$x, d$y, rts = "vrs", orientation = "in", slack = FALSE,
              xref = d$x[1:30, , drop = FALSE], yref = d$y[1:30, , drop = FALSE])
  ## Fewer reference DMUs means a frontier no further out, so input-oriented
  ## scores can only rise -- and may exceed 1, since the evaluated DMU need not
  ## be inside the technology it is scored against.
  ok <- is.finite(half$eff)
  expect_true(all(half$eff[ok] >= full$eff[ok] - 1e-9))
  expect_true(any(half$eff[ok] > 1))
  expect_equal(half$nref, 30L)
})

test_that("mismatched reference sets are refused", {
  d <- toy(n = 40, p = 2, q = 1)
  expect_error(dea(d$x, d$y, xref = d$x), "must be given together")
  expect_error(dea(d$x, d$y, xref = d$x[, 1, drop = FALSE], yref = d$y),
               "same variables in the same order")
  expect_error(dea(d$x, d$y, super = TRUE, xref = d$x, yref = d$y),
               "already outside the technology")
})

test_that("every model works when the reference set is a DIFFERENT SIZE", {
  ## n == nref for ordinary DEA, which hides any place the code indexes the
  ## reference columns by the count of evaluated DMUs. Scoring 12 points
  ## against 50 is what exposes it -- and did.
  ref  <- dea_sim(50, p = 2, q = 1, seed = 31)
  ev   <- dea_sim(12, p = 2, q = 1, seed = 32)

  for (mk in list(
      function() dea(ev$x, ev$y, rts = "vrs", xref = ref$x, yref = ref$y),
      function() dea(ev$x, ev$y, rts = "crs", orientation = "out",
                     xref = ref$x, yref = ref$y),
      function() dea(ev$x, ev$y, rts = "fdh", xref = ref$x, yref = ref$y),
      function() dea_sbm(ev$x, ev$y, rts = "vrs", orientation = "none",
                         xref = ref$x, yref = ref$y),
      function() dea_sbm(ev$x, ev$y, rts = "vrs", orientation = "in",
                         xref = ref$x, yref = ref$y),
      function() dea_sbm(ev$x, ev$y, rts = "vrs", orientation = "out",
                         xref = ref$x, yref = ref$y),
      function() dea_ddf(ev$x, ev$y, rts = "vrs", xref = ref$x, yref = ref$y))) {
    f <- suppressWarnings(mk())
    expect_equal(f$n, 12L)
    expect_equal(f$nref, 50L)
    expect_length(f$eff, 12L)
    expect_equal(dim(f$lambda), c(12L, 50L))
    ## Peers are drawn from the reference set, never from the evaluated one.
    expect_setequal(colnames(f$lambda), ref$data |> nrow() |> seq_len() |> as.character())
  }
})
