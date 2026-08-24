test_that("scale efficiency lies in (0, 1] and matches its definition", {
  d <- toy(n = 60, p = 1, q = 1, returns = 0.7, seed = 5)
  r <- dea_rts(d$x, d$y, orientation = "in")
  tb <- r$table
  expect_true(all(tb$scale_eff > 0 & tb$scale_eff <= 1 + 1e-9))
  expect_equal(tb$scale_eff, tb$crs / tb$vrs, tolerance = 1e-9)
  ## The three technologies nest, so the three scores do too.
  expect_true(all(tb$crs <= tb$nirs + 1e-9))
  expect_true(all(tb$nirs <= tb$vrs + 1e-9))
})

test_that("the classification is exhaustive and agrees with the sweep", {
  d <- toy(n = 60, p = 1, q = 1, returns = 0.7, seed = 5)
  r <- dea_rts(d$x, d$y)
  expect_true(all(r$table$rts %in% c("irs", "crs", "drs")))
  ## Scale-efficient DMUs are exactly the ones classified "crs".
  expect_equal(r$table$rts == "crs", abs(r$table$scale_eff - 1) < r$tol)
  ## The component scores are the ones dea() itself returns.
  for (rr in c("crs", "vrs", "nirs")) {
    f <- dea(d$x, d$y, rts = rr, orientation = "in", slack = FALSE)
    expect_equal(unname(r$table[[rr]]), unname(f$eff), tolerance = 1e-9)
  }
})

test_that("the classification agrees with the independent sum(lambda) rule", {
  ## Two rules, one derived from comparing three technologies and one read off
  ## a single program's intensity weights. They are logically distinct, so
  ## agreement is real evidence -- and an inverted IRS/DRS test, which is easy
  ## to write and produces perfectly plausible output, drives it to near zero.
  for (returns in c(0.7, 0.9)) {
    d <- toy(n = 60, p = 1, q = 1, returns = returns, seed = 5)
    for (ori in c("in", "out")) {
      t <- dea_rts(d$x, d$y, orientation = ori)$table
      sl <- ifelse(abs(t$crs - t$vrs) < 1e-6, "crs",
             ifelse(t$sum_lambda_crs < 1, "irs", "drs"))
      expect_equal(t$rts, sl, info = paste(returns, ori))
    }
  }
})

test_that("the input-oriented verdict describes the PROJECTED scale", {
  ## Returns to scale belong to the frontier point, not to the raw DMU. Under
  ## an output orientation inputs are held fixed, so the verdict tracks the
  ## DMU's own input level against the most productive scale size; under an
  ## input orientation it tracks the projected level, and the two disagree for
  ## any DMU whose projection crosses the MPSS.
  d <- toy(n = 60, p = 1, q = 1, returns = 0.7, seed = 5)
  x <- d$x[, 1]; y <- d$y[, 1]
  mpss <- x[which.max(y / x)]
  side <- function(v) ifelse(abs(v - mpss) < 1e-6, "crs",
                      ifelse(v < mpss, "irs", "drs"))

  ro <- dea_rts(d$x, d$y, orientation = "out")$table
  expect_equal(ro$rts, side(x))

  ri <- dea_rts(d$x, d$y, orientation = "in")$table
  expect_equal(ri$rts, side(ri$crs * x))

  expect_true(all(ro$scale_eff > 0 & ro$scale_eff <= 1 + 1e-9))
  ## Which means the two tables genuinely differ, and the test says so.
  expect_true(mean(ri$rts != ro$rts) > 0.1)
})

test_that("a constant-returns technology is classified as such", {
  ## Points exactly on a ray: every DMU is at the most productive scale size.
  x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), ncol = 1)
  y <- x * 2
  r <- suppressWarnings(dea_rts(x, y))
  expect_true(all(r$table$rts == "crs"))
  expect_equal(r$table$scale_eff, rep(1, 9), tolerance = 1e-8)
})
