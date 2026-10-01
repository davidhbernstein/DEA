## ---------------------------------------------------------------------------
## Weak disposability of undesirable outputs.
##
## No CRAN package to pin against, so the checks are structural: the equality
## row that IS the model, the containment ordering between the two VRS
## formulations, their coincidence under crs, and the data property the whole
## technology rests on.
## ---------------------------------------------------------------------------

weak_fixture <- function(n = 50, seed = 7, coupling = 1.1) {
  set.seed(seed)
  x <- matrix(runif(n * 2, 5, 20), n, 2,
              dimnames = list(NULL, c("capital", "energy")))
  g <- matrix((x[, 1] * x[, 2])^0.4 * exp(-abs(rnorm(n, 0, 0.2))), n, 1,
              dimnames = list(NULL, "product"))
  ## Null-joint by construction: bads rise with goods, so nobody is clean.
  b <- matrix(g[, 1]^coupling * exp(rnorm(n, 0, 0.15)), n, 1,
              dimnames = list(NULL, "CO2"))
  list(x = x, y = g, bad = b)
}

test_that("the bad output row is an equality, which is the whole model", {
  ## Strong disposability would leave a gap here; weak disposability cannot.
  f <- weak_fixture()
  w <- dea_weak(f$x, f$y, f$bad)
  ok <- is.finite(w$eff)
  lhs <- as.numeric(w$lambda[ok, , drop = FALSE] %*% f$bad)
  rhs <- f$bad[ok, 1] - w$eff[ok] * f$bad[ok, 1]
  expect_lt(max(abs(lhs - rhs)), 1e-7)
  ## The good-output row is an inequality and the input row too.
  expect_gt(min(as.numeric(w$lambda[ok, , drop = FALSE] %*% f$y) -
                  (f$y[ok, 1] + w$eff[ok] * f$y[ok, 1])), -1e-7)
  expect_lt(max(as.numeric(w$lambda[ok, , drop = FALSE] %*% f$x[, 1]) -
                  f$x[ok, 1]), 1e-7)
})

test_that("Kuosmanen contains Fare-Grosskopf, and coincides with it under crs", {
  ## mu = 0 is always feasible, so the two-vector technology contains the
  ## one-vector one and its beta can only be the larger. That is an identity,
  ## so a violation is a bug rather than a surprise.
  for (dirn in c("goods_bads", "bads", "all")) {
    for (coup in c(1.1, 0.3)) {
      f <- weak_fixture(40, seed = 12, coupling = coup)
      k <- suppressWarnings(dea_weak(f$x, f$y, f$bad, direction = dirn,
                                     technology = "kuosmanen"))
      fg <- suppressWarnings(dea_weak(f$x, f$y, f$bad, direction = dirn,
                                      technology = "fare"))
      expect_gt(min(k$eff - fg$eff, na.rm = TRUE), -1e-8)
      ## Under crs mu only adds slack free disposal already allows.
      kc <- suppressWarnings(dea_weak(f$x, f$y, f$bad, direction = dirn,
                                      rts = "crs", technology = "kuosmanen"))
      fc <- suppressWarnings(suppressMessages(
        dea_weak(f$x, f$y, f$bad, direction = dirn, rts = "crs",
                 technology = "fare")))
      expect_equal(unname(kc$eff), unname(fc$eff), tolerance = 1e-8)
    }
  }
  expect_message(dea_weak(weak_fixture(40)$x, weak_fixture(40)$y,
                          weak_fixture(40)$bad, rts = "crs", technology = "fare"),
                 "ignored here")
})

test_that("the extra intensity vector is used, and only where it can be", {
  ## If mu were never positive the two formulations would be the same program
  ## and `technology` would be decoration. The `bads` direction is where
  ## scaling the activity back is the attractive move.
  f <- weak_fixture(60, seed = 3, coupling = 0.3)
  k <- suppressWarnings(dea_weak(f$x, f$y, f$bad, direction = "bads"))
  expect_true(any(rowSums(k$mu) > 1e-7, na.rm = TRUE))
  ## Where mu is positive, sum(lambda) is below one: that IS the abatement.
  used <- which(rowSums(k$mu) > 1e-7)
  expect_true(all(rowSums(k$lambda)[used] < 1 - 1e-8))
  ## And the two vectors together are a convex combination.
  expect_equal(unname(rowSums(k$lambda) + rowSums(k$mu)), rep(1, 60),
               tolerance = 1e-7)
  ## The one-vector form has no mu at all.
  expect_null(suppressWarnings(
    dea_weak(f$x, f$y, f$bad, technology = "fare"))$mu)
})

test_that("beta is zero on the frontier and never negative on an observed DMU", {
  f <- weak_fixture()
  w <- dea_weak(f$x, f$y, f$bad)
  expect_true(all(w$eff >= -1e-9, na.rm = TRUE))
  expect_gt(sum(w$efficient), 0)
  expect_equal(w$efficient, is.finite(w$eff) & w$eff == 0)
  ## beta = 0 means its own lambda is a solution.
  eff1 <- which(w$efficient)[1]
  expect_equal(sum(w$lambda[eff1, ]) + sum(w$mu[eff1, ]), 1, tolerance = 1e-7)
})

test_that("null-jointness is checked, not imposed", {
  ## It cannot be imposed by a constraint -- it is a property of the data --
  ## and a clean DMU lets the frontier claim clean production is feasible.
  f <- weak_fixture(40)
  expect_true(dea_weak(f$x, f$y, f$bad)$null_joint)
  b2 <- f$bad; b2[1:3, 1] <- 0
  expect_warning(w <- dea_weak(f$x, f$y, b2), "null-jointness")
  expect_false(w$null_joint)
  expect_equal(w$n_clean, 3L)
  expect_output(print(w), "VIOLATED")
  ## A bad nobody produces is refused outright: its equality row pins beta.
  expect_error(dea_weak(f$x, f$y, matrix(0, 40, 1)), "zero for every DMU")
})

test_that("the direction is validated and its sign convention enforced", {
  f <- weak_fixture(40)
  ## Every component is a MAGNITUDE; the signs are in the program.
  expect_error(dea_weak(f$x, f$y, f$bad, direction = c(0, 0, 1, -1)),
               "negative entries")
  expect_error(dea_weak(f$x, f$y, f$bad, direction = c(0, 0, 0, 0)),
               "all zero")
  expect_error(dea_weak(f$x, f$y, f$bad, direction = c(1, 2)), "length p \\+ q \\+ b")
  expect_error(dea_weak(f$x, f$y, f$bad, direction = "sideways"), "direction")
  expect_error(dea_weak(f$x, f$y, f$bad[-1, , drop = FALSE]),
               "different numbers of DMUs")
  expect_error(dea_weak(f$x, -f$y, f$bad), "Negative values")
  ## A user-supplied direction is honoured.
  w <- dea_weak(f$x, f$y, f$bad, direction = c(0, 0, 1, 1))
  expect_equal(unname(w$direction[1, ]), c(0, 0, 1, 1))
})

test_that("the methods run and rts is restricted to the two that are defined", {
  f <- weak_fixture(40)
  w <- dea_weak(f$x, f$y, f$bad)
  expect_s3_class(w, "dea")
  expect_output(print(w), "Weak disposability")
  expect_output(print(w), "null-jointness")
  expect_equal(nobs(w), 40L)
  expect_error(dea_weak(f$x, f$y, f$bad, rts = "nirs"), "rts")
})
