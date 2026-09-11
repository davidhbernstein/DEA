test_that("the slacks-based measure obeys its defining inequalities", {
  d <- toy(n = 50)
  rho  <- dea_sbm(d$x, d$y, rts = "vrs", orientation = "none")
  rhoI <- dea_sbm(d$x, d$y, rts = "vrs", orientation = "in")
  rhoO <- dea_sbm(d$x, d$y, rts = "vrs", orientation = "out")
  th   <- dea(d$x, d$y, rts = "vrs", orientation = "in", slack = TRUE)

  expect_true(all(rho$eff > 0 & rho$eff <= 1 + 1e-9))
  ## The non-oriented measure scores both sides, so it cannot beat either
  ## one-sided measure.
  expect_true(all(rho$eff <= rhoI$eff + 1e-8))
  expect_true(all(rho$eff <= rhoO$eff + 1e-8))
  ## The radial score ignores slack, so it cannot be smaller.
  expect_true(all(rhoI$eff <= th$eff + 1e-8))
  ## rho = 1 exactly on the Pareto-Koopmans efficient set.
  expect_equal(unname(rho$eff == 1), unname(th$efficient))
})

test_that("the Charnes-Cooper linearization returns usable primal quantities", {
  d <- toy(n = 40)
  rho <- dea_sbm(d$x, d$y, rts = "vrs", orientation = "none")
  ## lambda and the slacks come back divided by t, so they must satisfy the
  ## ORIGINAL balance equations, not the transformed ones.
  lhs_x <- rho$lambda %*% d$x + rho$slack_x
  lhs_y <- rho$lambda %*% d$y - rho$slack_y
  expect_equal(unname(lhs_x), unname(d$x), tolerance = 1e-6)
  expect_equal(unname(lhs_y), unname(d$y), tolerance = 1e-6)
  expect_equal(unname(rowSums(rho$lambda)), rep(1, nrow(d$x)), tolerance = 1e-6)
})

test_that("the oriented measures are not degenerate", {
  ## Minimizing the output slack rather than maximizing it makes rho_O equal 1
  ## for every DMU and nothing errors, so this is worth asserting explicitly.
  d <- toy(n = 40)
  rhoO <- dea_sbm(d$x, d$y, rts = "crs", orientation = "out")
  expect_true(sum(rhoO$eff < 1 - 1e-8) > nrow(d$x) / 4)
  rhoI <- dea_sbm(d$x, d$y, rts = "crs", orientation = "in")
  expect_true(sum(rhoI$eff < 1 - 1e-8) > nrow(d$x) / 4)
})

test_that("dea_sbm() refuses data it cannot divide by", {
  d <- toy(n = 40)
  x0 <- d$x; x0[1, 1] <- 0
  expect_error(dea_sbm(x0, d$y), "strictly positive")
  y0 <- d$y; y0[2, 1] <- 0
  expect_error(dea_sbm(d$x, y0), "no output at all|strictly positive")
})

test_that("results agree with DJL's recorded answers", {
  d   <- toy(n = 40, p = 2, q = 2, seed = 4)
  ref <- ref_values("sbm_DJL.csv")
  for (rts in c("crs", "vrs")) {
    for (ori in c("none", "in", "out")) {
      ours <- dea_sbm(d$x, d$y, rts = rts, orientation = ori)
      r <- ref[ref$rts == rts & ref$orientation == ori, ]
      r <- r[order(r$dmu), ]
      expect_equal(as.numeric(ours$eff), r$eff, tolerance = 1e-6,
                   info = paste(rts, ori))
    }
  }
})

test_that("an out-of-sample point outside the technology is NA, not a number", {
  ## The additive program needs a convex combination of reference DMUs to
  ## dominate the evaluated point in every input and every output. A point
  ## outside the estimated technology has none, and the honest answer is NA.
  ## The radial model has no such requirement and scores it below 1.
  set.seed(4)
  ref  <- dea_sim(60, p = 2, q = 1, seed = 21)
  ## The point must sit INSIDE the input range -- otherwise no convex
  ## combination satisfies X'lambda <= x_o and the radial program is infeasible
  ## too, which would test something else. So: take an existing DMU's inputs
  ## and put its output well above what the frontier there can deliver.
  ex <- ref$x[1, , drop = FALSE]
  ey <- ref$y[1, , drop = FALSE] * 3

  expect_warning(
    sb <- dea_sbm(ex, ey, rts = "vrs", orientation = "out",
                  xref = ref$x, yref = ref$y),
    "INFEASIBLE")
  expect_true(is.na(sb$eff))
  expect_equal(unname(sb$status), 2L)

  ## The radial model returns a finite score below 1 for the same point.
  rd <- dea(ex, ey, rts = "vrs", orientation = "out", slack = FALSE,
            xref = ref$x, yref = ref$y)
  expect_true(is.finite(rd$eff))
  expect_lt(rd$eff, 1)

  ## And a self-referenced fit is never infeasible, because every DMU
  ## dominates itself.
  self <- dea_sbm(ref$x, ref$y, rts = "vrs", orientation = "out")
  expect_true(all(is.finite(self$eff)))
})
