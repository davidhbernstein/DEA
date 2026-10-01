## ---------------------------------------------------------------------------
## Undesirable outputs.
##
## There is no cross-package reference for this one -- deaR is the only CRAN
## package with it and it cannot be loaded on the development machine, since it
## imports rgl, which needs libGLU from X11. So what is checked here is of the
## first and third kinds in CLAUDE.md: a property that must hold identically
## (translation invariance, which is the WHOLE justification for the Seiford-Zhu
## route) and a monotonicity that follows from the definition.
## ---------------------------------------------------------------------------

bad_fixture <- function(n = 40, seed = 4) {
  set.seed(seed)
  x <- matrix(runif(n * 2, 2, 10), n, 2, dimnames = list(NULL, c("x1", "x2")))
  g <- matrix((x[, 1] * x[, 2])^0.4 * exp(-abs(rnorm(n, 0, 0.2))), n, 1,
              dimnames = list(NULL, "good"))
  ## Null-joint with the good output: more good, more bad.
  b <- matrix(g[, 1]^1.2 * exp(rnorm(n, 0, 0.15)), n, 1, dimnames = list(NULL, "CO2"))
  list(x = x, y = g, bad = b)
}

test_that("the translate route does not depend on the constant", {
  ## This is the property the route exists for. If it ever fails, the scores
  ## have become a function of an arbitrary number.
  f <- bad_fixture()
  base <- dea_undesirable(f$x, f$y, f$bad)$eff
  for (w in c(max(f$bad) + 1, max(f$bad) + 50, 100, 1000)) {
    expect_equal(dea_undesirable(f$x, f$y, f$bad, trans = w)$eff, base,
                 tolerance = 1e-8)
  }
  ## And it is the same answer as transforming by hand and calling dea().
  ## dea_undesirable() adds refusals and a diagnostic, not a new program.
  hand <- dea(f$x, cbind(f$y, 1000 - f$bad), rts = "vrs", orientation = "in")
  expect_equal(unname(hand$eff), unname(base), tolerance = 1e-8)
})

test_that("a constant large enough to lose the data warns", {
  ## Invariance is exact in algebra and not in double precision, and the failure
  ## at the far end is SILENT -- every program reports an optimum and the scores
  ## are wrong. So the gate is asserted, not assumed.
  f <- bad_fixture()
  base <- dea_undesirable(f$x, f$y, f$bad)$eff
  expect_warning(dea_undesirable(f$x, f$y, f$bad, trans = 1e5),
                 "relative spread")
  expect_silent(dea_undesirable(f$x, f$y, f$bad, trans = 1e3))
  ## The thing the gate is protecting against: a clean solve and a wrong answer.
  r <- suppressWarnings(dea_undesirable(f$x, f$y, f$bad, trans = 1e9))
  expect_true(all(r$status %in% c(0L, 1L)))
  expect_gt(max(abs(r$eff - base)), 0.1)
})

test_that("the translate route refuses every model it is not valid for", {
  f <- bad_fixture()
  for (rts in c("crs", "nirs", "ndrs")) {
    expect_error(dea_undesirable(f$x, f$y, f$bad, rts = rts), "Seiford-Zhu")
  }
  expect_error(dea_undesirable(f$x, f$y, f$bad, orientation = "out"), "Seiford-Zhu")
  ## The input route needs no invariance, so it is available everywhere.
  for (rts in c("vrs", "crs", "nirs", "ndrs")) for (ori in c("in", "out")) {
    r <- dea_undesirable(f$x, f$y, f$bad, method = "input", rts = rts,
                         orientation = ori)
    expect_s3_class(r, "dea_undesirable")
    expect_equal(r$p, 3L)          ## the bad joined the inputs
    expect_equal(r$q, 1L)
  }
})

test_that("the refused combinations really are the ones that move", {
  ## The refusal is not a matter of taste: translating the bad output changes
  ## the answer in exactly the cells that are refused, and in no others.
  f <- bad_fixture()
  ## Both constants must leave the transformed output positive -- dea() refuses
  ## negative data, which is itself part of why this route needs a constant.
  moved <- function(rts, ori) {
    e <- lapply(c(max(f$bad) + 1, 1000), function(w) {
      yy <- cbind(f$y, w - f$bad)
      suppressWarnings(dea(f$x, yy, rts = rts, orientation = ori,
                           slack = FALSE)$eff)
    })
    max(abs(e[[1]] - e[[2]]))
  }
  ## Calibrated against the allowed cell rather than against a guessed
  ## constant: whatever the solver tolerance is, the refused cells must move by
  ## orders of magnitude more than the cell that is invariant.
  inv <- moved("vrs", "in")
  expect_lt(inv, 1e-8)
  for (cell in list(c("vrs", "out"), c("crs", "in"), c("crs", "out"),
                    c("nirs", "in"), c("ndrs", "in"))) {
    expect_gt(moved(cell[1], cell[2]), max(1e-3, 1e4 * inv))
  }
})

test_that("slack_bad is the reducible bad output under both routes", {
  f <- bad_fixture()
  for (m in c("translate", "input")) {
    r <- dea_undesirable(f$x, f$y, f$bad, method = m)
    expect_equal(dim(r$slack_bad), c(40L, 1L))
    expect_equal(colnames(r$slack_bad), "CO2")
    expect_true(all(r$slack_bad >= -1e-9, na.rm = TRUE))
    ## It cannot exceed the DMU's own bad output: the target is b_o - s.
    expect_true(all(r$slack_bad <= f$bad[, 1] + 1e-9, na.rm = TRUE))
  }
  ## Under the translate route the bad's slack is literally the output slack on
  ## the transformed column, which is what makes it the reducible bad.
  r <- dea_undesirable(f$x, f$y, f$bad)
  expect_equal(unname(r$slack_bad[, 1]), unname(r$slack_y[, 2]))
})

test_that("a cleaner DMU is never scored worse, all else equal", {
  ## Direct from the definition under either route, and the one check here that
  ## would catch the bad being pointed the wrong way round.
  f <- bad_fixture(30)
  for (m in c("translate", "input")) {
    base <- dea_undesirable(f$x, f$y, f$bad, method = m)$eff
    b2 <- f$bad; b2[1] <- b2[1] * 0.5      ## DMU 1 pollutes half as much
    alt <- dea_undesirable(f$x, f$y, b2, method = m)$eff
    expect_gte(alt[1], base[1] - 1e-9)
  }
})

test_that("arguments are validated", {
  f <- bad_fixture(30)
  expect_error(dea_undesirable(f$x, f$y, f$bad[-1, , drop = FALSE]),
               "different numbers of DMUs")
  expect_error(dea_undesirable(f$x, f$y, f$bad, trans = c(1, 2)), "length 1")
  expect_error(dea_undesirable(f$x, f$y, f$bad, trans = NA_real_), "finite")
  ## A constant below the largest bad leaves the transformed output negative.
  expect_error(dea_undesirable(f$x, f$y, f$bad, trans = min(f$bad)),
               "non-positive")
  expect_error(dea_undesirable(f$x, f$y, f$bad, method = "weak"), "method")
  expect_output(print(dea_undesirable(f$x, f$y, f$bad)), "Seiford-Zhu")
  expect_output(print(dea_undesirable(f$x, f$y, f$bad, method = "input")),
                "Hailu-Veeman")
  ## The print method must still show the underlying fit.
  expect_output(print(dea_undesirable(f$x, f$y, f$bad)), "theta")
})
