## The multiplier form. Two families of check: strong duality (the two programs
## must attain the same value) and the sign table on u0, which is the only part
## a reader cannot verify by inspection and the only part that fails silently.

test_that("the multiplier program attains the envelopment value", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  for (rts in c("crs", "vrs", "nirs", "ndrs")) {
    for (ori in c("in", "out")) {
      env <- suppressWarnings(dea(d$x, d$y, rts = rts, orientation = ori,
                                  slack = FALSE))
      mul <- suppressWarnings(dea(d$x, d$y, rts = rts, orientation = ori,
                                  slack = FALSE, multipliers = TRUE))
      u0 <- if (is.null(mul$u0)) 0 else mul$u0
      obj <- if (ori == "in") rowSums(mul$u * d$y) - u0
             else             rowSums(mul$v * d$x) - u0
      expect_equal(unname(obj), unname(env$eff), tolerance = 1e-8,
                   info = paste(rts, ori))
      ## The score itself must not move: multipliers = TRUE adds a program, it
      ## does not change the one that was already there.
      expect_equal(mul$eff, env$eff, tolerance = 1e-12)
    }
  }
})

test_that("the normalization holds in the CALLER's units, not the solver's", {
  ## Scaling divides each column by its mean before solving, so a weight that
  ## satisfies v'x = 1 inside the LP satisfies it outside only after being
  ## divided by the same factor. Getting that wrong is invisible in the score.
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  for (sc in c(TRUE, FALSE)) {
    fin <- suppressWarnings(dea(d$x, d$y, rts = "vrs", orientation = "in",
                                slack = FALSE, scaling = sc, multipliers = TRUE))
    expect_equal(unname(rowSums(fin$v * d$x)), rep(1, nrow(d$x)),
                 tolerance = 1e-8, info = paste("scaling", sc))
    fout <- suppressWarnings(dea(d$x, d$y, rts = "vrs", orientation = "out",
                                 slack = FALSE, scaling = sc, multipliers = TRUE))
    expect_equal(unname(rowSums(fout$u * d$y)), rep(1, nrow(d$y)),
                 tolerance = 1e-8, info = paste("scaling", sc))
  }
})

test_that("scaling leaves the multipliers themselves unchanged", {
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  a <- suppressWarnings(dea(d$x, d$y, rts = "vrs", slack = FALSE,
                            scaling = TRUE,  multipliers = TRUE))
  b <- suppressWarnings(dea(d$x, d$y, rts = "vrs", slack = FALSE,
                            scaling = FALSE, multipliers = TRUE))
  ## Only where the optimum is unique -- an efficient DMU has a whole face of
  ## optimal weights and the two runs may stop at different vertices.
  ineff <- a$eff < 1 - 1e-6
  expect_equal(a$v[ineff, ], b$v[ineff, ], tolerance = 1e-6)
  expect_equal(a$u[ineff, ], b$u[ineff, ], tolerance = 1e-6)
})

## The residual of the technology constraint at DMU o against reference DMU j.
## Input orientation:   u'Y_j - v'X_j - u0 <= 0
## Output orientation:  v'X_j - u'Y_j - u0 >= 0
## Written once because the two tests below both need it, and because getting
## the orientation's sign wrong here would make both of them vacuous.
mult_residual <- function(f, x, y) {
  u0 <- if (is.null(f$u0)) rep(0, nrow(x)) else f$u0
  if (f$orientation == "in") f$u %*% t(y) - f$v %*% t(x) - u0
  else                       f$v %*% t(x) - f$u %*% t(y) - u0
}

test_that("the weights are FEASIBLE in the multiplier program, both orientations", {
  ## This is what makes the weights a supporting price vector rather than an
  ## arbitrary positive one. Note what it does NOT catch: with the u0 sign
  ## table inverted the weights stay feasible to 1e-12, because they are then
  ## a valid price vector for a DIFFERENT technology. The complementary
  ## slackness test below is what has teeth against that.
  d <- toy(30, p = 2, q = 2, returns = 0.9)
  for (rts in c("crs", "vrs", "nirs", "ndrs")) {
    for (ori in c("in", "out")) {
      f <- suppressWarnings(dea(d$x, d$y, rts = rts, orientation = ori,
                                slack = FALSE, multipliers = TRUE))
      R <- mult_residual(f, d$x, d$y)
      if (ori == "in") expect_true(all(R <=  1e-7), info = paste(rts, ori))
      else             expect_true(all(R >= -1e-7), info = paste(rts, ori))
      expect_true(all(f$v >= -1e-9) && all(f$u >= -1e-9), info = paste(rts, ori))
    }
  }
})

test_that("complementary slackness ties the peers to the weights", {
  ## For a linear program, complementary slackness holds between ANY optimal
  ## primal solution and ANY optimal dual solution -- not merely between a
  ## corresponding pair. So although the envelopment and multiplier programs
  ## are solved independently here, a reference DMU carrying positive lambda
  ## MUST have its multiplier constraint exactly tight. Non-uniqueness of the
  ## weights does not weaken this; it is a theorem about the optimal faces.
  ##
  ## WHY IT EARNS ITS PLACE. It is the only check here that catches an
  ## inverted u0 sign table in every case. Feasibility does not catch it at
  ## all, and the value-agreement test above catches it only where the
  ## restriction binds -- which is why that one carries a skip_if. Measured
  ## against a deliberately inverted table, the residual at a positive lambda
  ## goes from 1e-12 to between 0.4 and 1.6 for all four of nirs/ndrs by
  ## in/out.
  d <- dea_sim(40, p = 2, q = 2, returns = 0.9, seed = 5)
  for (rts in c("crs", "vrs", "nirs", "ndrs")) {
    for (ori in c("in", "out")) {
      f <- suppressWarnings(dea(d$x, d$y, rts = rts, orientation = ori,
                                slack = FALSE, multipliers = TRUE))
      R <- mult_residual(f, d$x, d$y)
      pos <- f$lambda > 1e-7
      ## A degenerate sample with no positive lambda anywhere would make this
      ## vacuous; assert there is something to check before checking it.
      expect_true(sum(pos) > 0, info = paste(rts, ori))
      expect_lt(max(abs(R[pos])), 1e-6)
    }
  }
})

test_that("u0 obeys the sign table for nirs and ndrs, in both orientations", {
  ## Checked on data where the restriction BINDS -- nirs/ndrs scores that
  ## differ from the vrs ones. On data where it does not bind, every sign
  ## convention agrees and the test proves nothing. This is how the table in
  ## R/lp.R was fixed in the first place.
  d <- dea_sim(40, p = 2, q = 1, returns = 0.6, ineff = "hnorm",
               mean_ineff = 0.3, x_range = c(1, 12), seed = 103)
  want <- list(nirs = c(`in` = 1, out = -1), ndrs = c(`in` = -1, out = 1))
  for (rts in c("nirs", "ndrs")) {
    for (ori in c("in", "out")) {
      vrs <- suppressWarnings(dea(d$x, d$y, rts = "vrs", orientation = ori,
                                  slack = FALSE))$eff
      f <- suppressWarnings(dea(d$x, d$y, rts = rts, orientation = ori,
                                slack = FALSE, multipliers = TRUE))
      skip_if(sum(abs(f$eff - vrs) > 1e-6) == 0,
              paste(rts, ori, "does not bind on this sample"))
      s <- want[[rts]][[ori]]
      expect_true(all(s * f$u0 >= -1e-7), info = paste(rts, ori))
    }
  }
})

test_that("multipliers agree with Benchmarking wherever the optimum is unique", {
  ## The reference file already carries Benchmarking's weights under THIS
  ## package's names. Theirs are transposed relative to ours: dea(DUAL = TRUE)
  ## returns `ux` for the INPUT weights and `vy` for the OUTPUT weights, the
  ## reverse of the CCR notation. The swap is undone once, in
  ## tools/make_reference_values.R, so it does not have to be re-derived here.
  d <- toy(25, p = 2, q = 2, returns = 0.9, seed = 12)
  m <- suppressWarnings(dea(d$x, d$y, rts = "crs", orientation = "in",
                            slack = FALSE, multipliers = TRUE))
  ref <- ref_values("multipliers_Benchmarking.csv")
  wide <- function(w, k) {
    r <- ref[ref$which == w, ]
    matrix(r$value[order(r$j, r$dmu)], nrow(d$x), k)
  }
  bv <- wide("v", ncol(d$x)); bu <- wide("u", ncol(d$y))

  ineff <- m$eff < 1 - 1e-6
  expect_true(any(ineff))
  expect_equal(unname(m$v[ineff, ]), bv[ineff, ], tolerance = 1e-6)
  expect_equal(unname(m$u[ineff, ]), bu[ineff, ], tolerance = 1e-6)
  ## And where they disagree, the DMU is efficient -- alternate optima, not a
  ## discrepancy. Both vectors give the same score; only the vertex differs.
  differs <- rowSums(abs(unname(m$v) - bv)) > 1e-6
  expect_true(all(m$eff[differs] > 1 - 1e-6))
})

test_that("multipliers are refused for fdh, and explained", {
  d <- toy(20)
  expect_error(dea(d$x, d$y, rts = "fdh", multipliers = TRUE), "not convex")
})

test_that("multipliers() needs a fit that has them", {
  d <- toy(20)
  f <- suppressWarnings(dea(d$x, d$y))
  expect_error(multipliers(f), "multipliers = TRUE")
  g <- suppressWarnings(dea(d$x, d$y, multipliers = TRUE))
  m <- multipliers(g)
  expect_true(is.matrix(m))
  expect_identical(nrow(m), nrow(d$x))
  expect_true("u0" %in% colnames(m))          ## vrs by default
  expect_identical(rownames(m), g$dmu)
})

test_that("super-efficiency has a multiplier form too", {
  d <- toy(30, p = 2, q = 1)
  env <- suppressWarnings(dea(d$x, d$y, rts = "crs", slack = FALSE,
                              super = TRUE))
  mul <- suppressWarnings(dea(d$x, d$y, rts = "crs", slack = FALSE,
                              super = TRUE, multipliers = TRUE))
  ok <- is.finite(env$eff)
  expect_equal(unname(rowSums(mul$u * d$y)[ok]), unname(env$eff[ok]),
               tolerance = 1e-7)
})

test_that("multipliers work against an external reference set", {
  d <- toy(50, p = 2, q = 2, returns = 0.9, seed = 5)
  e <- toy(12, p = 2, q = 2, returns = 0.9, seed = 6)
  env <- suppressWarnings(dea(e$x, e$y, rts = "vrs", slack = FALSE,
                              xref = d$x, yref = d$y))
  mul <- suppressWarnings(dea(e$x, e$y, rts = "vrs", slack = FALSE,
                              xref = d$x, yref = d$y, multipliers = TRUE))
  expect_identical(dim(mul$v), c(12L, 2L))
  expect_equal(unname(rowSums(mul$u * e$y) - mul$u0), unname(env$eff),
               tolerance = 1e-7)
})
