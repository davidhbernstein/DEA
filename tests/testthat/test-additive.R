test_that("all three measures return a well-formed fit", {
  d <- toy(n = 50, p = 2, q = 1)
  for (m in c("ram", "mip", "unweighted")) {
    f <- dea_add(d$x, d$y, measure = m, rts = "vrs")
    expect_s3_class(f, "dea")
    expect_identical(f$model, "additive")
    expect_identical(f$measure, m)
    expect_true(all(is.finite(f$eff)), info = m)
    expect_length(f$eff, nrow(d$x))
    expect_true(all(f$slack_x >= 0), info = m)
    expect_true(all(f$slack_y >= 0), info = m)
  }
})

test_that("the normalized measures land in [0, 1] and the raw one does not", {
  d <- toy(n = 50, p = 2, q = 1)
  ## RAM is bounded on BOTH sides -- no slack can exceed its variable's range,
  ## so the objective cannot exceed 1. That is what RAM was designed for.
  e <- dea_add(d$x, d$y, measure = "ram", rts = "vrs")$eff
  expect_true(all(e >= -1e-9 & e <= 1 + 1e-9))
  expect_true(any(e == 1))
  ## MIP is capped at 1 but NOT bounded below: an input slack cannot exceed the
  ## DMU's own input, but an output slack has no such ceiling.
  e <- dea_add(d$x, d$y, measure = "mip", rts = "vrs")$eff
  expect_true(all(e <= 1 + 1e-9))
  expect_true(any(e == 1))
  u <- dea_add(d$x, d$y, measure = "unweighted", rts = "vrs")$eff
  expect_true(all(u >= -1e-9))
  expect_true(any(u == 0))
  expect_gt(max(u), 1)     ## a raw total, not a share
})

test_that("the additive efficient set IS the Pareto-Koopmans set", {
  ## This is the model's defining property and it must hold on every scale:
  ## the objective is zero exactly when no slack remains, whatever the weights.
  d <- toy(n = 50, p = 2, q = 1)
  pk  <- unname(dea(d$x, d$y, rts = "vrs", orientation = "in", slack = TRUE)$efficient)
  rho <- unname(dea_sbm(d$x, d$y, rts = "vrs")$eff == 1)
  expect_equal(pk, rho)
  for (m in c("ram", "mip", "unweighted")) {
    expect_equal(unname(dea_add(d$x, d$y, measure = m, rts = "vrs")$efficient), pk,
                 info = m)
  }
})

test_that("RAM and MIP are units invariant; the unweighted model is not", {
  d <- toy(n = 40, p = 2, q = 1)
  x2 <- d$x; x2[, 1] <- x2[, 1] * 1000
  for (m in c("ram", "mip")) {
    expect_equal(unname(dea_add(x2, d$y, measure = m)$eff),
                 unname(dea_add(d$x, d$y, measure = m)$eff),
                 tolerance = 1e-7, info = m)
  }
  ## The documented defect, asserted rather than merely described: the raw
  ## totals move when the units move, so the ORDERING of inefficient DMUs is
  ## not stable. The efficient SET is invariant in exact arithmetic, and stays
  ## so under a moderate rescale once the tolerances are relative.
  xm <- d$x; xm[, 1] <- xm[, 1] * 8
  a <- dea_add(d$x, d$y, measure = "unweighted")
  b <- dea_add(xm,   d$y, measure = "unweighted")
  expect_false(isTRUE(all.equal(unname(a$eff), unname(b$eff))))
  expect_equal(unname(a$efficient), unname(b$efficient))
})

test_that("the unweighted total matches dea.add's recorded answers", {
  d   <- toy(n = 40, p = 2, q = 2, seed = 8)
  ref <- ref_values("additive_Benchmarking.csv")
  for (rts in c("crs", "vrs", "nirs", "ndrs")) {
    ours <- dea_add(d$x, d$y, measure = "unweighted", rts = rts, scaling = FALSE)
    r <- ref[ref$rts == rts, ]; r <- r[order(r$dmu), ]
    expect_equal(as.numeric(ours$eff), r$sum, tolerance = 1e-7, info = rts)
  }
})

test_that("scaling is ignored for the unweighted measure, honoured otherwise", {
  d <- toy(n = 40, p = 2, q = 1)
  ## The unweighted model is DEFINED on the caller's units, so the flag must
  ## not silently change the answer.
  a <- dea_add(d$x, d$y, measure = "unweighted", scaling = TRUE)
  b <- dea_add(d$x, d$y, measure = "unweighted", scaling = FALSE)
  expect_equal(unname(a$eff), unname(b$eff))
  expect_false(a$scaling)
  ## RAM is invariant either way.
  expect_equal(unname(dea_add(d$x, d$y, measure = "ram", scaling = TRUE)$eff),
               unname(dea_add(d$x, d$y, measure = "ram", scaling = FALSE)$eff),
               tolerance = 1e-7)
})

test_that("dea_add() validates and refuses what it cannot compute", {
  d <- toy(n = 40, p = 2, q = 1)
  expect_error(dea_add(d$x, d$y, measure = "additive"), "must be one of")
  expect_error(dea_add(d$x, d$y, rts = "fdh"), "must be one of")
  ## MIP divides by the DMU's own levels; RAM does not.
  x0 <- d$x; x0[1, 1] <- 0
  expect_error(dea_add(x0, d$y, measure = "mip"), "strictly positive")
  expect_true(all(is.finite(dea_add(x0, d$y, measure = "ram")$eff)))
  ## No common ratio scale for the unweighted total.
  expect_error(efficiency(dea_add(d$x, d$y, measure = "unweighted"), "score"),
               "not defined")
  expect_equal(efficiency(dea_add(d$x, d$y, measure = "ram"), "score"),
               dea_add(d$x, d$y, measure = "ram")$eff)
})

test_that("an external reference set of a different size works", {
  ref <- dea_sim(50, p = 2, q = 1, seed = 31)
  ev  <- dea_sim(12, p = 2, q = 1, seed = 32)
  f <- suppressWarnings(dea_add(ev$x, ev$y, measure = "ram",
                                xref = ref$x, yref = ref$y))
  expect_equal(f$n, 12L)
  expect_equal(f$nref, 50L)
  expect_equal(dim(f$lambda), c(12L, 50L))
})

test_that("fitted() lands on the frontier for the additive model", {
  d <- toy(n = 40, p = 2, q = 1)
  f <- dea_add(d$x, d$y, measure = "ram", rts = "vrs")
  pr <- fitted(f)
  expect_equal(dim(pr), c(40L, 3L))
  ## The projection removes exactly the slacks.
  expect_equal(unname(pr[, 1:2]), unname(d$x - f$slack_x), tolerance = 1e-9)
  expect_equal(unname(pr[, 3, drop = FALSE]), unname(d$y + f$slack_y), tolerance = 1e-9)
})


test_that("an extreme column spread warns rather than failing quietly", {
  ## The unweighted measure is solved in the caller's units by design, so a
  ## wide spread conditions the program badly. On this design one DMU fails
  ## numerically outright -- which must be reported as a conditioning problem,
  ## NOT as the point lying outside the technology.
  d <- toy(n = 40, p = 2, q = 1)
  xw <- d$x; xw[, 1] <- xw[, 1] * 1000
  ## Two distinct warnings, and both matter: one says the measure is a poor
  ## choice at this spread, the other says a specific DMU could not be solved.
  ww <- testthat::capture_warnings(f <- dea_add(xw, d$y, measure = "unweighted"))
  expect_true(any(grepl("span a factor", ww)))
  expect_true(any(grepl("FAILED NUMERICALLY", ww)))
  ## And the numerical failure is NOT reported as infeasibility -- that would
  ## send the reader looking at their data instead of at the conditioning.
  expect_false(any(grepl("INFEASIBLE", ww)))
  expect_true(any(is.na(f$eff)))
  ## RAM normalizes the columns, so it is untroubled by the same data.
  expect_silent(g <- dea_add(xw, d$y, measure = "ram"))
  expect_true(all(is.finite(g$eff)))
  expect_equal(unname(g$eff), unname(dea_add(d$x, d$y, measure = "ram")$eff),
               tolerance = 1e-7)
})

test_that("a stray positional argument is refused, not silently absorbed", {
  ## dea_add(x, y, "ram") puts "ram" in `data` and takes the DEFAULT measure,
  ## returning a well-formed answer to a question nobody asked.
  d <- toy(n = 30, p = 2, q = 1)
  expect_error(dea_add(d$x, d$y, "ram"), "must be a data frame")
  expect_error(dea(d$x, d$y, "vrs"), "must be a data frame")
  expect_error(dea_sbm(d$x, d$y, "crs"), "must be a data frame")
  ## A real data frame still works positionally.
  s <- dea_sim(30, p = 2, q = 1, seed = 2)
  expect_s3_class(dea(c("x1", "x2"), "y1", s$data), "dea")
})
