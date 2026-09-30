## ---------------------------------------------------------------------------
## Weight restrictions, assurance regions and the cone-ratio model.
## ---------------------------------------------------------------------------

test_that("with no restrictions the restricted score IS the ordinary score", {
  ## The anchor for everything else here. dea_wr() reads its score off the
  ## MULTIPLIER program where dea() reads it off the envelopment one, so with
  ## the restriction set empty the two must agree by strong duality -- and this
  ## is the last point at which they can be compared at all.
  d <- toy(30, 2, 2, 3)
  for (rr in c("crs", "vrs", "nirs", "ndrs")) for (ori in c("in", "out")) {
    expect_equal(
      unname(dea_wr(d$x, d$y, NULL, rts = rr, orientation = ori)$eff),
      unname(dea(d$x, d$y, rts = rr, orientation = ori, slack = FALSE)$eff),
      tolerance = 1e-7, info = paste(rr, ori))
  }
})

test_that("assurance regions match Benchmarking, on a design where they bind", {
  ## Benchmarking's DUAL takes ratios against the first input and the first
  ## output. Which row is which was pinned by measurement with asymmetric
  ## bounds; see tools/make_reference_values.R.
  d <- toy(30, 2, 2, 8)
  DU <- matrix(c(0.3, 0.8, 1.2, 4.0), nrow = 2L, byrow = TRUE)
  rs <- c(wr_ratio("v", 2, 1, lower = DU[1, 1], upper = DU[1, 2]),
          wr_ratio("u", 2, 1, lower = DU[2, 1], upper = DU[2, 2]))
  ref <- ref_values("weights_Benchmarking.csv")
  for (ori in c("in", "out")) {
    want <- ref$eff[ref$orientation == ori]
    got <- suppressWarnings(
      dea_wr(d$x, d$y, rs, rts = "crs", orientation = ori)$eff)
    expect_equal(unname(got), want, tolerance = 1e-7, info = ori)
    ## and it has to bind, or the agreement is with an unrestricted fit
    free <- dea(d$x, d$y, rts = "crs", orientation = ori, slack = FALSE)$eff
    expect_gt(max(abs(want - unname(free))), 0.02)
  }
})

test_that("a restriction can only make a DMU look worse", {
  ## Restricting the weight set shrinks the feasible region of a max, so the
  ## input-oriented optimum falls; the output-oriented program is a min, so its
  ## optimum rises. Either way the DMU loses ground, and a restriction that
  ## improved a score would mean a row with the wrong sign.
  d <- toy(40, 2, 2, 3)
  rs <- c(wr_ratio("v", 2, 1, lower = 0.5, upper = 2),
          wr_ratio("u", 2, 1, lower = 0.5, upper = 2))
  free_in <- dea(d$x, d$y, rts = "crs", slack = FALSE)$eff
  got_in <- dea_wr(d$x, d$y, rs, rts = "crs", orientation = "in")$eff
  expect_true(all(got_in <= free_in + 1e-9))
  expect_gt(sum(abs(got_in - free_in) > 1e-8), 10)

  free_out <- dea(d$x, d$y, rts = "crs", orientation = "out", slack = FALSE)$eff
  got_out <- dea_wr(d$x, d$y, rs, rts = "crs", orientation = "out")$eff
  expect_true(all(got_out >= free_out - 1e-9))
})

test_that("restrictions are read in the caller's units, whatever the scaling", {
  ## The solver divides every column by its mean, so a weight it reports
  ## applies to x_i / sx_i. A restriction written in the caller's units has to
  ## be divided by the same factors on the way in; getting that wrong would not
  ## error, it would impose a different restriction.
  d <- toy(40, 2, 2, 3)
  rs <- c(wr_ratio("v", 2, 1, lower = 0.5, upper = 2),
          wr_ratio("u", 2, 1, lower = 0.5, upper = 2))
  on  <- dea_wr(d$x, d$y, rs, rts = "crs", scaling = TRUE)$eff
  off <- dea_wr(d$x, d$y, rs, rts = "crs", scaling = FALSE)$eff
  expect_equal(unname(on), unname(off), tolerance = 1e-8)

  ## and a change of units carries the restriction with it: scaling input 1 by
  ## 100 and input 2 by 0.01 multiplies the admissible v[2]/v[1] by 100/0.01.
  x2 <- sweep(d$x, 2L, c(100, 0.01), "*")
  rs2 <- c(wr_ratio("v", 2, 1, lower = 0.5 * 1e4, upper = 2 * 1e4),
           wr_ratio("u", 2, 1, lower = 0.5, upper = 2))
  expect_equal(unname(dea_wr(x2, d$y, rs2, rts = "crs")$eff), unname(on),
               tolerance = 1e-8)
})

test_that("a bound and a ratio say what they claim about the weights", {
  d <- toy(30, 2, 1, 5)
  f <- dea_wr(d$x, d$y, wr_bound("v", 1, lower = 0.05), rts = "crs")
  ok <- is.finite(f$eff)
  expect_true(all(f$v[ok, 1] >= 0.05 - 1e-8))

  g <- dea_wr(d$x, d$y, wr_ratio("v", 2, 1, lower = 0.8, upper = 1.25),
              rts = "crs")
  ok <- is.finite(g$eff) & g$v[, 1] > 1e-10
  r <- g$v[ok, 2] / g$v[ok, 1]
  expect_true(all(r >= 0.8 - 1e-6 & r <= 1.25 + 1e-6))
})

test_that("the cone-ratio model equals its polar written as inequalities", {
  ## Two structurally different routes to one estimand: dea_cone() transforms
  ## the DATA and calls dea(), dea_wr() adds ROWS to the multiplier program.
  ## They are not wrappers for each other, and agreeing is worth something.
  d <- toy(30, 2, 2, 12)
  expect_equal(unname(dea_cone(d$x, d$y, rts = "crs")$eff),
               unname(dea(d$x, d$y, rts = "crs")$eff))

  ## A cone inside the non-negative orthant, so that the cone model's implicit
  ## alpha >= 0 and dea_wr()'s implicit v >= 0 describe the same set.
  A <- rbind(c(1, 0.4), c(0.3, 1))
  B <- rbind(c(1, 0.25), c(0.5, 1))
  co <- dea_cone(d$x, d$y, A, B, rts = "crs", slack = FALSE)$eff
  iA <- solve(t(A)); iB <- solve(t(B))
  rs <- c(wr_linear(v = iA[1, ], type = ">=", rhs = 0),
          wr_linear(v = iA[2, ], type = ">=", rhs = 0),
          wr_linear(u = iB[1, ], type = ">=", rhs = 0),
          wr_linear(u = iB[2, ], type = ">=", rhs = 0))
  expect_equal(unname(co),
               unname(dea_wr(d$x, d$y, rs, rts = "crs")$eff),
               tolerance = 1e-7)
  ## and the cone binds, so the agreement is not two unrestricted fits
  expect_gt(max(abs(co - dea(d$x, d$y, rts = "crs", slack = FALSE)$eff)), 0.05)
})

test_that("contradictory restrictions are infeasible and say which to look at", {
  ## The one place a DEA program here can be genuinely infeasible on ordinary
  ## data. v[2] >= 5 v[1] and v[2] <= v[1] force both to zero, and then the
  ## normalization v'x = 1 cannot hold.
  d <- toy(12, 2, 1, 4)
  rs <- c(wr_ratio("v", 2, 1, lower = 5), wr_ratio("v", 2, 1, upper = 1))
  expect_warning(f <- dea_wr(d$x, d$y, rs, rts = "crs"), "INFEASIBLE")
  expect_true(all(is.na(f$eff)))
  expect_true(all(f$status == 2L))
  ## the same data unrestricted is fine, which is what makes the message's
  ## advice -- look at the restrictions, not the data -- true
  expect_false(anyNA(dea(d$x, d$y, rts = "crs", slack = FALSE)$eff))
})

test_that("a restricted fit refuses to report peers or slacks", {
  ## The trap this whole entry point exists to avoid: the dual of a restricted
  ## multiplier program carries one non-negative residue term per restriction,
  ## sitting exactly where a slack sits. Reporting it as a slack, or the lambda
  ## beside it as a peer weight, describes a projection onto a point that is
  ## not in the technology.
  d <- toy(20, 2, 1, 5)
  f <- dea_wr(d$x, d$y, wr_ratio("v", 2, 1, lower = 0.5, upper = 2),
              rts = "crs")
  expect_error(peers(f), "residue")
  expect_error(slacks(f), "residue")
  expect_null(f$lambda)
  expect_null(f$slack_x)
  expect_output(print(f), "No peers or slacks")
})

test_that("restriction constructors validate what they are given", {
  expect_error(wr_bound("v", 1), "at least one")
  expect_error(wr_bound("v", 1, lower = 2, upper = 1), "above upper")
  expect_error(wr_bound("v", 0, lower = 1), "positive whole number")
  expect_error(wr_ratio("u", 2, 2, lower = 1), "the same index")
  expect_error(wr_ratio("u", 1, 2), "at least one")
  expect_error(wr_linear(), "needs coefficients")
  expect_error(wr_linear(v = c(1, NA)), "finite")
  expect_error(c(wr_bound("v", 1, lower = 1), 3), "combine only with")

  d <- toy(15, 2, 1, 5)
  expect_error(dea_wr(d$x, d$y, list(a = 1), rts = "crs"), "must be built by")
  expect_error(dea_wr(d$x, d$y, wr_linear(v = c(1, 1, 1)), rts = "crs"),
               "3 coefficient")
  expect_error(dea_cone(d$x, d$y, A = matrix(1, 2, 3), rts = "crs"),
               "one column per input")
  expect_warning(dea_cone(d$x, d$y, A = matrix(c(1, 1), 1, 2), rts = "crs"),
                 "1-dimensional subspace")
  expect_equal(length(c(wr_bound("v", 1, lower = 1, upper = 2))$rows), 2L)
  expect_output(print(wr_ratio("v", 2, 1, lower = 1)), "v\\[2\\]/v\\[1\\] >= 1")
})
