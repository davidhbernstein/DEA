test_that("the constant-returns score matches a hand-worked answer", {
  fit <- suppressWarnings(dea(tiny$x, tiny$y, rts = "crs", orientation = "in"))
  expect_equal(unname(fit$eff), unname(tiny_crs_in), tolerance = 1e-9)
  ## Output orientation is the reciprocal under constant returns, exactly.
  out <- suppressWarnings(dea(tiny$x, tiny$y, rts = "crs", orientation = "out"))
  expect_equal(unname(out$eff), 1 / unname(tiny_crs_in), tolerance = 1e-9)
})

test_that("every technology and orientation returns a well-formed fit", {
  d <- toy()
  for (rts in c("crs", "vrs", "nirs", "ndrs", "fdh")) {
    for (ori in c("in", "out")) {
      fit <- dea(d$x, d$y, rts = rts, orientation = ori)
      expect_s3_class(fit, "dea")
      expect_identical(fit$rts, rts)
      expect_true(all(is.finite(fit$eff)), info = paste(rts, ori))
      if (ori == "in") {
        expect_true(all(fit$eff <= 1 + 1e-9), info = paste(rts, ori))
        expect_true(all(fit$eff > 0), info = paste(rts, ori))
      } else {
        expect_true(all(fit$eff >= 1 - 1e-9), info = paste(rts, ori))
      }
      expect_true(any(fit$eff == 1), info = paste(rts, ori))
      expect_length(fit$eff, nrow(d$x))
    }
  }
})

test_that("the technologies nest in the order the theory says", {
  d <- toy()
  e <- function(r) dea(d$x, d$y, rts = r, orientation = "in", slack = FALSE)$eff
  crs <- e("crs"); vrs <- e("vrs"); nirs <- e("nirs"); ndrs <- e("ndrs"); fdh <- e("fdh")
  ## A larger technology cannot give a larger input-oriented score.
  expect_true(all(crs <= nirs + 1e-9))
  expect_true(all(crs <= ndrs + 1e-9))
  expect_true(all(nirs <= vrs + 1e-9))
  expect_true(all(ndrs <= vrs + 1e-9))
  ## Dropping convexity shrinks the technology again, so FDH is largest.
  expect_true(all(vrs <= fdh + 1e-9))
})

test_that("slacks make the efficiency judgement Pareto-Koopmans", {
  d <- toy(n = 50)
  with_slack <- dea(d$x, d$y, rts = "vrs", slack = TRUE)
  no_slack   <- dea(d$x, d$y, rts = "vrs", slack = FALSE)
  expect_equal(with_slack$eff, no_slack$eff, tolerance = 1e-9)
  ## The radial score is unaffected; the verdict may be stricter.
  expect_true(all(with_slack$efficient <= (with_slack$eff == 1)))
  expect_null(no_slack$slack_x)
  expect_true(all(with_slack$slack_x >= 0))
  expect_true(all(with_slack$slack_y >= 0))
  ## A DMU with any slack cannot be Pareto-Koopmans efficient.
  has <- rowSums(with_slack$slack_x) + rowSums(with_slack$slack_y) > 1e-8
  expect_false(any(with_slack$efficient & has))
})

test_that("scaling changes no answer", {
  d <- toy()
  for (f in list(dea, dea_sbm, dea_ddf)) {
    a <- suppressWarnings(f(d$x, d$y, rts = "vrs", scaling = TRUE))
    b <- suppressWarnings(f(d$x, d$y, rts = "vrs", scaling = FALSE))
    expect_equal(a$eff, b$eff, tolerance = 1e-7)
  }
  ## And neither does the units the caller happens to use.
  s <- dea(d$x * 1000, d$y / 7, rts = "vrs", orientation = "in", slack = FALSE)
  u <- dea(d$x,        d$y,     rts = "vrs", orientation = "in", slack = FALSE)
  expect_equal(s$eff, u$eff, tolerance = 1e-7)
})

test_that("super-efficiency exceeds ordinary efficiency, and is NA when infeasible", {
  d <- toy(n = 40)
  ord <- dea(d$x, d$y, rts = "crs", orientation = "in", slack = FALSE)
  sup <- dea(d$x, d$y, rts = "crs", orientation = "in", super = TRUE)
  ## Removing a DMU from its own reference set can only enlarge its score.
  expect_true(all(sup$eff >= ord$eff - 1e-9))
  ## Inefficient DMUs are unaffected: they were not spanning anything.
  ineff <- ord$eff < 1 - 1e-8
  expect_equal(sup$eff[ineff], ord$eff[ineff], tolerance = 1e-8)
  ## Input-oriented constant returns is always feasible; variable returns is not.
  expect_false(any(is.na(sup$eff)))
  supv <- dea(d$x, d$y, rts = "vrs", orientation = "in", super = TRUE)
  ordv <- dea(d$x, d$y, rts = "vrs", orientation = "in", slack = FALSE)
  ## Only a DMU that was ON the frontier can become infeasible or score above 1;
  ## an inefficient one keeps its score exactly, since it never spanned anything.
  expect_true(all(is.na(supv$eff[ordv$eff == 1]) | supv$eff[ordv$eff == 1] >= 1 - 1e-9))
  expect_false(any(is.na(supv$eff[ordv$eff < 1 - 1e-8])))
})

test_that("results agree with Benchmarking's recorded answers", {
  d   <- toy(n = 60, p = 2, q = 2, seed = 9)
  ref <- ref_values("radial_Benchmarking.csv")
  for (rts in c("crs", "vrs", "nirs", "ndrs", "fdh")) {
    for (ori in c("in", "out")) {
      ours <- dea(d$x, d$y, rts = rts, orientation = ori, slack = FALSE)
      r <- ref[ref$rts == rts & ref$orientation == ori, ]
      r <- r[order(r$dmu), ]
      expect_equal(unname(ours$eff), r$eff, tolerance = 1e-7,
                   info = paste(rts, ori))
    }
  }
  ## Super-efficiency too, including which DMUs come back infeasible. They
  ## report Inf where this package reports NA; what must match is WHICH units.
  ours <- dea(d$x, d$y, rts = "vrs", orientation = "in", super = TRUE)
  s <- ref_values("super_Benchmarking.csv"); s <- s[order(s$dmu), ]
  expect_equal(unname(is.na(ours$eff)), !as.logical(s$feasible))
  ok <- as.logical(s$feasible)
  expect_equal(unname(ours$eff[ok]), s$eff[ok], tolerance = 1e-7)
})
