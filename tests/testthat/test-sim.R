test_that("dea_sim() produces a technology matching its own stated truth", {
  s <- dea_sim(200, p = 2, q = 3, returns = 0.8, seed = 1)
  expect_s3_class(s, "dea_sim")
  expect_equal(dim(s$x), c(200L, 2L))
  expect_equal(dim(s$y), c(200L, 3L))
  ## ||y|| = f(x) exp(-u), so phi = f(x)/||y|| = exp(u) exactly.
  expect_equal(sqrt(rowSums(s$y^2)), s$frontier * exp(-s$u), tolerance = 1e-12)
  expect_equal(s$phi, exp(s$u), tolerance = 1e-12)
  expect_equal(s$theta, exp(-s$u / 0.8), tolerance = 1e-12)
  expect_true(all(s$theta > 0 & s$theta <= 1))
  expect_true(all(s$phi >= 1))
  ## E[u] is what was asked for.
  expect_equal(mean(s$u), 0.3, tolerance = 0.1)
})

test_that("a seed reproduces, and restores the caller's stream", {
  a <- dea_sim(50, seed = 4)
  b <- dea_sim(50, seed = 4)
  expect_equal(a$x, b$x)
  expect_equal(a$u, b$u)
  set.seed(11); before <- runif(1)
  set.seed(11); invisible(dea_sim(50, seed = 4)); after <- runif(1)
  expect_equal(before, after)
})

test_that("a constant-returns fit recovers the truth on a constant-returns design", {
  s <- dea_sim(600, p = 1, q = 1, returns = 1, seed = 6)
  fit <- dea(s$x, s$y, rts = "crs", orientation = "in", slack = FALSE)
  ## The frontier here is a ray, so the estimator is essentially exact.
  expect_lt(mean(fit$eff - s$theta), 0.01)
  expect_gt(cor(fit$eff, s$theta), 0.999)
  ## And the bias is one-signed: the estimated frontier is inside the true one.
  expect_true(all(fit$eff >= s$theta - 1e-8))
})

test_that("the estimator is biased toward 1, always", {
  s <- dea_sim(100, p = 2, q = 1, returns = 0.9, seed = 12)
  fo <- dea(s$x, s$y, rts = "vrs", orientation = "out", slack = FALSE)
  ## phi_hat <= phi: the estimated frontier cannot lie above the true one.
  expect_true(all(fo$eff <= s$phi + 1e-8))
  expect_lt(mean(fo$eff - s$phi), 0)
})

test_that("dea_rate() states the published rates", {
  ## Kneip, Park and Simar (1998); Park, Simar and Weiner (2000).
  expect_equal(dea_rate(1, 1, "vrs"), -4 / 3)
  expect_equal(dea_rate(2, 1, "vrs"), -1)
  expect_equal(dea_rate(1, 1, "crs"), -2)
  expect_equal(dea_rate(2, 2, "crs"), -1)
  expect_equal(dea_rate(1, 1, "fdh"), -1)
  expect_equal(dea_rate(4, 4, "vrs"), -4 / 9)
  ## The estimator's own rate is half the mean-squared-error rate.
  expect_equal(dea_rate(2, 1, "vrs", "estimator"), -1 / 2)
  ## The curse of dimensionality, as a monotone statement.
  d <- vapply(1:8, function(k) dea_rate(k, k, "vrs"), numeric(1))
  expect_true(all(diff(d) > 0))
})

test_that("x_range moves the input support and nothing else", {
  a <- dea_sim(200, p = 2, q = 1, x_range = c(1, 2), seed = 1)
  b <- dea_sim(200, p = 2, q = 1, x_range = c(0.3, 2), seed = 1)
  expect_true(all(a$x >= 1 & a$x <= 2))
  expect_true(all(b$x >= 0.3 & b$x <= 2))
  expect_lt(min(b$x), 1)
  ## The inefficiency draws are independent of x, so widening the support does
  ## not change them -- which is what makes the two designs comparable.
  expect_equal(a$u, b$u)
  ## And the stated truth still matches the realized data.
  expect_equal(sqrt(rowSums(b$y^2)), b$frontier * exp(-b$u), tolerance = 1e-12)
})

test_that("widening the support downward makes theta estimable", {
  ## On the default support a variable-returns fit cannot reach theta, because
  ## the projection leaves the sample. Given room underneath, it can.
  ev <- function(m, lo, hi, seed = 99) {
    set.seed(seed)
    x <- matrix(runif(m, lo, hi), m, 1)
    u <- rexp(m, 1 / 0.3)
    list(x = x, y = matrix(x^0.8 * exp(-u), m, 1), theta = exp(-u / 0.8))
  }
  E <- ev(60, 1.3, 1.7)
  mse <- function(xlo, n) {
    v <- numeric(12)
    for (b in seq_len(12)) {
      s <- dea_sim(n, p = 1, q = 1, returns = 0.8, x_range = c(xlo, 2), seed = 500 + b)
      f <- suppressWarnings(dea(E$x, E$y, rts = "vrs", orientation = "in",
                                slack = FALSE, peers = FALSE,
                                xref = s$x, yref = s$y))
      v[b] <- mean((f$eff - E$theta)^2)
    }
    mean(v)
  }
  narrow <- c(mse(1.0, 200), mse(1.0, 800))
  wide   <- c(mse(0.2, 200), mse(0.2, 800))
  ## On the narrow support the error barely moves; on the wide one it falls.
  expect_gt(narrow[2] / narrow[1], 0.8)
  expect_lt(wide[2] / wide[1], 0.6)
  expect_lt(wide[2], narrow[2] / 5)
})

test_that("dea_sim() validates its design arguments", {
  expect_error(dea_sim(50, returns = 1.5), "must lie in")
  expect_error(dea_sim(50, returns = 0), "must lie in")
  expect_error(dea_sim(0), "positive number of DMUs")
  expect_error(dea_sim(50, p = 0), "at least 1")
  expect_error(dea_sim(50, mean_ineff = -1), "positive")
  expect_error(dea_sim(50, x_range = c(2, 1)), "increasing positive")
  expect_error(dea_sim(50, x_range = c(0, 2)), "increasing positive")
  expect_error(dea_sim(50, x_range = 1), "increasing positive")
})

## ---------------------------------------------------------------------------
## Frontier geometry.  The whole reason a second shape is cheap is that the
## closed forms need homogeneity of degree r and NOTHING else -- not
## smoothness, not strict concavity, nothing Cobb-Douglas.  That is an argument,
## so it gets a test that does not assume it: the distance is brute-forced from
## its DEFINITION with uniroot and compared against the formula the object
## reports.  If the reasoning were wrong for the non-smooth shape, this is where
## it would show.
## ---------------------------------------------------------------------------

## f as each shape defines it, written out independently of R/sim.R's vectorised
## code so that a shared slip cannot cancel out.
sim_f <- function(shape, par, r, p) {
  switch(shape,
    cobb  = function(x) prod(x)^(r / p),
    ces   = function(x) (mean(x^par))^(r / par),
    facet = { A <- DEA:::.sim_facets(as.integer(par), p)
              function(x) (min(A %*% x) / p)^r })
}

geom_cases <- list(
  list(shape = "cobb",  par = NULL, p = 2, r = 1.0),
  list(shape = "cobb",  par = NULL, p = 3, r = 0.8),
  list(shape = "ces",   par = 0.5,  p = 2, r = 0.9),
  list(shape = "ces",   par = 1,    p = 2, r = 1.0),   # flat frontier
  list(shape = "ces",   par = -2,   p = 3, r = 0.7),   # past Leontief-ward
  list(shape = "facet", par = 2,    p = 2, r = 1.0),
  list(shape = "facet", par = 4,    p = 3, r = 0.85))

test_that("every frontier shape is homogeneous of degree `returns`", {
  for (cs in geom_cases) {
    f <- sim_f(cs$shape, cs$par, cs$r, cs$p)
    s <- dea_sim(40, p = cs$p, q = 1, returns = cs$r, frontier = cs$shape,
                 frontier_par = cs$par, seed = 11)
    set.seed(5)
    for (i in 1:20) {
      x <- s$x[i, ]; t <- 0.3 + 2 * runif(1)
      expect_equal(f(t * x), t^cs$r * f(x), tolerance = 1e-10,
                   info = paste(cs$shape, cs$par))
    }
  }
})

test_that("the closed-form truth survives every shape, brute-forced", {
  for (cs in geom_cases) {
    f <- sim_f(cs$shape, cs$par, cs$r, cs$p)
    s <- dea_sim(60, p = cs$p, q = 2, returns = cs$r, frontier = cs$shape,
                 frontier_par = cs$par, seed = 11)
    ny <- sqrt(rowSums(s$y^2))
    lab <- paste(cs$shape, cs$par)
    ## The generative identity, independent of any distance.
    expect_equal(ny, vapply(seq_len(60), function(i) f(s$x[i, ]), 0) * exp(-s$u),
                 tolerance = 1e-10, info = lab)
    ## theta = the smallest t with f(t x) >= ||y||, solved numerically.
    brute_t <- vapply(seq_len(60), function(i)
      stats::uniroot(function(t) f(t * s$x[i, ]) - ny[i], c(1e-8, 1),
                     tol = 1e-14)$root, 0)
    expect_equal(brute_t, unname(s$theta), tolerance = 1e-8, info = lab)
    ## phi = the largest scaling of y still inside the output set.
    brute_p <- vapply(seq_len(60), function(i) f(s$x[i, ]) / ny[i], 0)
    expect_equal(brute_p, unname(s$phi), tolerance = 1e-10, info = lab)
  }
})

test_that("the shapes agree on the diagonal, so `returns` means one thing", {
  ## f(c, ..., c) = c^r for all three, which is what lets a geometry table be
  ## read down a column: only the shape changes, not the scale elasticity.
  for (p in 2:3) for (r in c(0.7, 1)) {
    fs <- list(sim_f("cobb", NULL, r, p), sim_f("ces", 0.4, r, p),
               sim_f("ces", -3, r, p), sim_f("facet", 2, r, p))
    for (cc in c(1, 1.5, 2)) for (f in fs)
      expect_equal(f(rep(cc, p)), cc^r, tolerance = 1e-10)
  }
})

test_that("the default shape is unchanged by the shape argument existing", {
  ## `frontier` was added to a released function; the Cobb-Douglas path must be
  ## bit-identical or every number recorded before it moves.
  for (cs in list(list(p = 1, q = 1, r = 1), list(p = 2, q = 1, r = 0.9),
                  list(p = 3, q = 2, r = 0.7))) {
    a <- dea_sim(90, p = cs$p, q = cs$q, returns = cs$r, seed = 3)
    b <- dea_sim(90, p = cs$p, q = cs$q, returns = cs$r, seed = 3,
                 frontier = "cobb")
    for (f in c("x", "y", "theta", "phi", "u", "frontier"))
      expect_identical(a[[f]], b[[f]])
    expect_null(a$facet)
  }
})

test_that("a facet design reports which facets actually bind", {
  ## A facet that is never the minimum contributes nothing, so `frontier_par`
  ## would overstate the geometry. p >= 3 is where interior facets can bind.
  s <- dea_sim(400, p = 3, q = 1, frontier = "facet", frontier_par = 4, seed = 2)
  expect_length(s$facet, 400L)
  expect_setequal(unique(s$facet), 1:4)
  expect_equal(dim(s$facet_normals), c(4L, 3L))
  expect_true(all(s$facet_normals > 0))
  expect_equal(unname(rowSums(s$facet_normals)), rep(3, 4), tolerance = 1e-12)
})

test_that("dea_sim() refuses geometries that would misreport themselves", {
  ## Collinear normals at p = 2: only the two endpoints can ever be the
  ## minimum, so m > 2 there is always an overstatement.
  expect_error(dea_sim(50, p = 2, frontier = "facet", frontier_par = 3),
               "only 2 of them can ever be active")
  expect_error(dea_sim(50, p = 1, frontier = "facet"), "at least 2 inputs")
  expect_error(dea_sim(50, p = 2, frontier = "facet", frontier_par = 1),
               "single whole number")
  ## rho > 1 makes the technology non-convex, so the reported theta stops being
  ## the estimand rather than merely being hard to reach.
  expect_error(dea_sim(50, p = 2, frontier = "ces", frontier_par = 1.5),
               "convex hull")
  expect_error(dea_sim(50, p = 2, frontier = "ces", frontier_par = 1e-6),
               "Cobb-Douglas")
  expect_error(dea_sim(50, frontier = "cobb", frontier_par = 2),
               "no meaning")
  expect_error(dea_sim(50, p = 2, frontier = "nope"), "frontier")
})
