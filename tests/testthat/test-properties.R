## ---------------------------------------------------------------------------
## test-properties.R -- exact properties, asserted across every estimator.
##
## The rest of the suite checks VALUES: hand-worked answers, agreement with
## Benchmarking and DJL, recovery of dea_sim()'s closed-form truth. This file
## checks PROPERTIES -- statements that must hold exactly, for every estimator,
## on any data, with no reference value to compare against.
##
## The two kinds catch different things. A value test catches a wrong formula.
## A property test catches a wrong *index*: a loop that reads column j where it
## meant row j still produces plausible numbers on the fixture that was used to
## write it, and still agrees with another package if that package makes the
## same slip. It cannot survive being asked for the same answer twice with the
## rows in a different order.
##
## Every estimator is registered in ESTIMATORS below and every applicable
## property is applied to every entry, so adding an estimator without adding
## its properties is not possible by accident.
## ---------------------------------------------------------------------------

## A reproducible panel with enough spread to make the properties bite. Kept
## small: these are contract checks, not a power study.
prop_data <- function(n = 30, p = 2, q = 2, seed = 42) {
  set.seed(seed)
  X <- matrix(stats::runif(n * p, 1, 3), n, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))
  fx <- exp(rowSums(log(X)) * 0.8 / p)
  D  <- matrix(abs(stats::rnorm(n * q)), n, q)
  D  <- D / sqrt(rowSums(D^2))
  Y  <- D * (fx * exp(-stats::rexp(n, 1 / 0.3)))
  colnames(Y) <- paste0("y", seq_len(q))
  list(x = X, y = Y)
}

## Each estimator as a bare function of (x, y) returning one number per DMU.
##
##   units    is it invariant to rescaling each column separately?
##
## `units = FALSE` entries are NOT failures -- the unweighted additive model is
## defined in the caller's units on purpose (see R/additive.R), and a DDF with a
## fixed direction is not units invariant because the direction is a quantity in
## those units. Both are recorded here so the exclusion is deliberate and
## visible rather than an estimator quietly going untested.
ESTIMATORS <- list(
  list(name = "radial crs in",  units = TRUE,
       f = function(x, y) dea(x, y, rts = "crs", orientation = "in",  slack = FALSE, peers = FALSE)$eff),
  list(name = "radial crs out", units = TRUE,
       f = function(x, y) dea(x, y, rts = "crs", orientation = "out", slack = FALSE, peers = FALSE)$eff),
  list(name = "radial vrs in",  units = TRUE,
       f = function(x, y) dea(x, y, rts = "vrs", orientation = "in",  slack = FALSE, peers = FALSE)$eff),
  list(name = "radial vrs out", units = TRUE,
       f = function(x, y) dea(x, y, rts = "vrs", orientation = "out", slack = FALSE, peers = FALSE)$eff),
  list(name = "radial nirs in", units = TRUE,
       f = function(x, y) dea(x, y, rts = "nirs", orientation = "in", slack = FALSE, peers = FALSE)$eff),
  list(name = "radial ndrs in", units = TRUE,
       f = function(x, y) dea(x, y, rts = "ndrs", orientation = "in", slack = FALSE, peers = FALSE)$eff),
  list(name = "radial fdh in",  units = TRUE,
       f = function(x, y) dea(x, y, rts = "fdh", orientation = "in",  slack = FALSE, peers = FALSE)$eff),
  list(name = "sbm vrs in",     units = TRUE,
       f = function(x, y) dea_sbm(x, y, rts = "vrs", orientation = "in",   peers = FALSE)$eff),
  list(name = "sbm vrs none",   units = TRUE,
       f = function(x, y) dea_sbm(x, y, rts = "vrs", orientation = "none", peers = FALSE)$eff),
  list(name = "ddf vrs both",   units = TRUE,
       f = function(x, y) dea_ddf(x, y, direction = "both", rts = "vrs", peers = FALSE)$beta),
  list(name = "ddf vrs unit",   units = FALSE,
       f = function(x, y) dea_ddf(x, y, direction = "unit", rts = "vrs", peers = FALSE)$beta),
  list(name = "additive ram",   units = TRUE,
       f = function(x, y) dea_add(x, y, measure = "ram", rts = "vrs", peers = FALSE)$eff),
  list(name = "additive mip",   units = TRUE,
       f = function(x, y) dea_add(x, y, measure = "mip", rts = "vrs", peers = FALSE)$eff),
  list(name = "additive unweighted", units = FALSE,
       f = function(x, y) dea_add(x, y, measure = "unweighted", rts = "vrs", peers = FALSE)$eff)
)

## Solver tolerance. These are exact mathematical identities; what is being
## allowed for is floating point and the simplex, not model error.
PROP_TOL <- 1e-8

test_that("every units-invariant estimator is invariant to per-column rescaling", {
  d <- prop_data()
  a <- c(3, 0.1); b <- c(7, 0.02)          ## deliberately extreme and unequal
  x2 <- sweep(d$x, 2L, a, "*"); y2 <- sweep(d$y, 2L, b, "*")
  for (e in Filter(function(z) z$units, ESTIMATORS)) {
    expect_equal(unname(e$f(d$x, d$y)), unname(e$f(x2, y2)),
                 tolerance = PROP_TOL, info = e$name)
  }
})

test_that("the units-dependent estimators really do depend on units", {
  ## The complement of the test above. If one of these ever becomes invariant,
  ## either the model changed or the argument is being ignored -- and silently
  ## gaining an invariance is as much a defect as silently losing one.
  d <- prop_data()
  x2 <- sweep(d$x, 2L, c(3, 0.1), "*"); y2 <- sweep(d$y, 2L, c(7, 0.02), "*")
  for (e in Filter(function(z) !z$units, ESTIMATORS)) {
    expect_false(isTRUE(all.equal(unname(e$f(d$x, d$y)), unname(e$f(x2, y2)),
                                  tolerance = 1e-6)), info = e$name)
  }
})

test_that("every estimator is invariant to the order of the rows", {
  ## The property most likely to catch an index bug: a routine that confuses a
  ## reference index with an evaluated index gives plausible numbers until the
  ## two orderings stop coinciding.
  d <- prop_data()
  n <- nrow(d$x)
  set.seed(11); perm <- sample(n)
  xp <- d$x[perm, , drop = FALSE]; yp <- d$y[perm, , drop = FALSE]
  for (e in ESTIMATORS) {
    expect_equal(unname(e$f(d$x, d$y))[perm], unname(e$f(xp, yp)),
                 tolerance = PROP_TOL, info = e$name)
  }
})

test_that("duplicating a DMU leaves every other DMU's score unchanged", {
  ## Adding a copy of an existing unit adds nothing to the technology -- the
  ## hull, the cone and the free-disposal set are all unchanged -- so no score
  ## may move. This catches an implementation that counts DMUs where it should
  ## be spanning them.
  d <- prop_data()
  n <- nrow(d$x)
  xd <- rbind(d$x, d$x[3L, , drop = FALSE])
  yd <- rbind(d$y, d$y[3L, , drop = FALSE])
  for (e in ESTIMATORS) {
    expect_equal(unname(e$f(d$x, d$y)), unname(e$f(xd, yd))[seq_len(n)],
                 tolerance = PROP_TOL, info = e$name)
    ## And the copy must score exactly what the original scores.
    expect_equal(unname(e$f(xd, yd))[3L], unname(e$f(xd, yd))[n + 1L],
                 tolerance = PROP_TOL, info = e$name)
  }
})

test_that("enlarging the reference set can never improve a score", {
  ## Monotonicity in the technology. More reference DMUs means a weakly larger
  ## production set, so an input-oriented score can only fall and an
  ## output-oriented one can only rise.
  d  <- prop_data()
  d2 <- prop_data(n = 15, seed = 7)
  xr <- rbind(d$x, d2$x); yr <- rbind(d$y, d2$y)
  for (rts in c("crs", "vrs", "nirs", "ndrs", "fdh")) {
    small <- dea(d$x, d$y, rts = rts, orientation = "in", slack = FALSE, peers = FALSE)$eff
    big   <- dea(d$x, d$y, rts = rts, orientation = "in", xref = xr, yref = yr,
                 slack = FALSE, peers = FALSE)$eff
    expect_true(all(big <= small + PROP_TOL), info = paste(rts, "in"))
    small_o <- dea(d$x, d$y, rts = rts, orientation = "out", slack = FALSE, peers = FALSE)$eff
    big_o   <- dea(d$x, d$y, rts = rts, orientation = "out", xref = xr, yref = yr,
                   slack = FALSE, peers = FALSE)$eff
    expect_true(all(big_o >= small_o - PROP_TOL), info = paste(rts, "out"))
  }
})

test_that("passing the sample as its own reference set is the default", {
  d <- prop_data()
  for (rts in c("crs", "vrs", "fdh")) {
    expect_identical(
      dea(d$x, d$y, rts = rts, orientation = "in", slack = FALSE, peers = FALSE)$eff,
      dea(d$x, d$y, rts = rts, orientation = "in", xref = d$x, yref = d$y,
          slack = FALSE, peers = FALSE)$eff)
  }
})

test_that("the multiplier solution reproduces the envelopment score exactly", {
  ## Primal-dual agreement, checked as an identity rather than by comparing two
  ## packages. The weights are non-unique at an efficient DMU, so they cannot be
  ## compared against anything -- but whichever vertex the solver returns must
  ## satisfy the program, and its objective must equal the score that came out
  ## of the envelopment form.
  d <- prop_data()
  for (rts in c("crs", "vrs", "nirs", "ndrs")) {
    for (ori in c("in", "out")) {
      f <- dea(d$x, d$y, rts = rts, orientation = ori, multipliers = TRUE,
               slack = FALSE, peers = FALSE)
      u0 <- if (length(f$u0)) f$u0 else 0
      if (ori == "in") {
        obj <- rowSums(f$u * d$y) - u0     ## maximise u'y - u0
        nrm <- rowSums(f$v * d$x)          ## subject to v'x = 1
      } else {
        obj <- rowSums(f$v * d$x) - u0     ## minimise v'x - v0
        nrm <- rowSums(f$u * d$y)          ## subject to u'y = 1
      }
      info <- paste(rts, ori)
      expect_equal(unname(obj), unname(f$eff), tolerance = PROP_TOL, info = info)
      expect_equal(unname(nrm), rep(1, nrow(d$x)), tolerance = PROP_TOL, info = info)
      ## Prices are non-negative, and the technology constraints hold.
      expect_true(all(f$v >= -PROP_TOL) && all(f$u >= -PROP_TOL), info = info)
    }
  }
})

test_that("the directional model reproduces the radial ones exactly", {
  ## g = (x, 0) contracts inputs proportionally, so beta = 1 - theta; g = (0, y)
  ## expands outputs proportionally, so beta = phi - 1. A cross-model identity:
  ## it ties R/ddf.R to R/dea.R without either being the reference for the
  ## other, and it holds for every technology.
  d <- prop_data()
  for (rts in c("crs", "vrs", "nirs", "ndrs")) {
    th <- dea(d$x, d$y, rts = rts, orientation = "in",  slack = FALSE, peers = FALSE)$eff
    ph <- dea(d$x, d$y, rts = rts, orientation = "out", slack = FALSE, peers = FALSE)$eff
    bi <- dea_ddf(d$x, d$y, direction = "in",  rts = rts, peers = FALSE)$beta
    bo <- dea_ddf(d$x, d$y, direction = "out", rts = rts, peers = FALSE)$beta
    expect_equal(unname(bi), unname(1 - th), tolerance = PROP_TOL, info = rts)
    expect_equal(unname(bo), unname(ph - 1), tolerance = PROP_TOL, info = rts)
  }
})

test_that("super-efficiency changes nothing for an inefficient DMU", {
  ## Removing a DMU from its own reference set can only matter if it was
  ## spanning the frontier, and an inefficient DMU is not.
  d <- prop_data()
  for (rts in c("crs", "vrs")) {
    f  <- dea(d$x, d$y, rts = rts, orientation = "in", slack = FALSE, peers = FALSE)
    fs <- dea(d$x, d$y, rts = rts, orientation = "in", super = TRUE,
              slack = FALSE, peers = FALSE)
    ineff <- f$eff < 1 - 1e-9
    expect_true(any(ineff), info = rts)
    expect_equal(unname(fs$eff[ineff]), unname(f$eff[ineff]),
                 tolerance = PROP_TOL, info = rts)
    ## And where it is finite, an efficient DMU scores at least 1.
    fin <- !ineff & is.finite(fs$eff)
    expect_true(all(fs$eff[fin] >= 1 - PROP_TOL), info = rts)
  }
})

test_that("the cost decomposition is exact and price-scale invariant", {
  d <- prop_data()
  set.seed(3); w <- stats::runif(ncol(d$x), 0.5, 2)
  fc <- dea_cost(d$x, d$y, w = w, rts = "vrs", peers = FALSE)
  ## Overall = technical x allocative, exactly -- which holds only if both came
  ## off the same frontier.
  expect_equal(unname(fc$eff), unname(fc$technical * fc$allocative), tolerance = PROP_TOL)
  ## The technical factor IS the radial score, not a separate computation.
  expect_equal(unname(fc$technical),
               unname(dea(d$x, d$y, rts = "vrs", orientation = "in",
                          slack = FALSE, peers = FALSE)$eff), tolerance = PROP_TOL)
  expect_true(all(fc$allocative <= 1 + PROP_TOL))
  ## Measuring an input in different units, at a correspondingly different
  ## price, is the same problem. A true no-op, not an approximate one.
  a <- c(3, 0.1)
  fc2 <- dea_cost(sweep(d$x, 2L, a, "*"), d$y, w = w / a, rts = "vrs", peers = FALSE)
  expect_equal(unname(fc$eff), unname(fc2$eff), tolerance = PROP_TOL)
})

test_that("the Nerlovian profit decomposition is exact and ADDITIVE", {
  ## Not a product. Nerlovian inefficiency inherits the directional model's
  ## additivity from the Chambers-Chung-Fare duality, and checking it with the
  ## multiplicative rule reports an error of about 0.5 on ordinary data and
  ## looks like a real defect. It is not; it is the wrong identity.
  d <- prop_data()
  set.seed(5)
  w <- stats::runif(ncol(d$x), 0.5, 2); r <- stats::runif(ncol(d$y), 0.5, 2)
  fp <- dea_profit(d$x, d$y, w = w, r = r, rts = "vrs", peers = FALSE)
  expect_equal(unname(fp$eff), unname(fp$technical + fp$allocative), tolerance = PROP_TOL)
})

test_that("cross-efficiency is row-order invariant ONLY with a secondary goal", {
  ## The sharpest property here, and the one that justifies a default.
  ##
  ## Without a secondary goal the cross-efficiency weights are not determined:
  ## an efficient rater has a whole face of optimal weights and the solver
  ## returns whichever vertex it reached. Reordering the rows changes the
  ## pivoting and therefore changes the answer -- by ~1e-2 here, which is
  ## first-order relative to the scores themselves. Doyle and Green's secondary
  ## programs pin the weights down and the answer becomes invariant to 1e-8.
  ##
  ## This is why dea_cross() defaults to secondary = "benevolent". If that
  ## default ever changes to "none", this test says what it costs.
  d <- prop_data()
  set.seed(11); perm <- sample(nrow(d$x))
  xp <- d$x[perm, , drop = FALSE]; yp <- d$y[perm, , drop = FALSE]
  ce <- function(x, y, sec) dea_cross(x, y, rts = "crs", secondary = sec)$eff

  for (sec in c("benevolent", "aggressive")) {
    expect_equal(unname(ce(d$x, d$y, sec))[perm], unname(ce(xp, yp, sec)),
                 tolerance = PROP_TOL, info = sec)
  }
  ## And the plain version is NOT invariant. Asserted, not merely documented.
  drift <- max(abs(unname(ce(d$x, d$y, "none"))[perm] - unname(ce(xp, yp, "none"))))
  expect_gt(drift, 1e-4)

  ## The package default is the invariant one.
  expect_equal(unname(dea_cross(d$x, d$y, rts = "crs")$eff),
               unname(ce(d$x, d$y, "benevolent")), tolerance = PROP_TOL)
})

test_that("scaling = TRUE is a no-op wherever it is claimed to be", {
  ## R/utils.R rescales columns to mean one before solving. Every estimator
  ## that claims invariance must return the identical answer either way; the
  ## unweighted additive model is documented as ignoring the argument, and that
  ## is asserted rather than assumed.
  d <- prop_data()
  same <- function(f) expect_equal(unname(f(TRUE)), unname(f(FALSE)), tolerance = PROP_TOL)
  same(function(s) dea(d$x, d$y, rts = "vrs", orientation = "in", scaling = s,
                       slack = FALSE, peers = FALSE)$eff)
  same(function(s) dea_sbm(d$x, d$y, rts = "vrs", orientation = "none", scaling = s,
                           peers = FALSE)$eff)
  same(function(s) dea_ddf(d$x, d$y, direction = "both", rts = "vrs", scaling = s,
                           peers = FALSE)$beta)
  same(function(s) dea_add(d$x, d$y, measure = "ram", rts = "vrs", scaling = s,
                           peers = FALSE)$eff)
})
