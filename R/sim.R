## ---------------------------------------------------------------------------
## dea_sim() -- simulate a technology whose efficiency is known in closed form.
##
## This is the package's testing substrate. Everything else here can be checked
## against another implementation of the same linear program, which catches
## coding errors but not misunderstandings: if two packages make the same
## mistake they agree perfectly.  What cannot be faked is recovering a truth
## that was written down before the estimator ran.
##
## THE DESIGN.  Inputs x are uniform on `x_range`^p, [1, 2] by default.  The
## frontier is a Cobb-Douglas aggregate
##
##     f(x) = prod_j x_j ^ (r/p),        r = `returns`
##
## and the output set at x is the non-negative part of the ball of radius f(x):
##
##     T = { (x, y) : ||y|| <= f(x) }.
##
## This is a genuine production technology in the sense DEA assumes -- f is
## concave for r <= 1 so T is convex, both f and the norm are monotone so T is
## freely disposable, and at r = 1 it is a cone, so a constant-returns fit is
## consistent there and only there.
##
## Efficiency then has a closed form.  With inefficiency u >= 0 drawn
## independently of x and the observed point set at ||y|| = f(x) e^{-u},
##
##     output-oriented Farrell   phi = f(x)/||y|| = e^{u}
##     input-oriented Farrell    theta = e^{-u/r}
##
## the second because scaling every input by theta scales f by theta^r.
##
## WHY THE OUTPUT DIRECTION IS RANDOM.  With q > 1 the direction of y within
## the ball is drawn uniformly and independently of u, so the output MIX
## carries no information about efficiency.  Fixing the mix instead would put
## every DMU on one ray, and a q-output problem on a single ray is a 1-output
## problem wearing q columns -- the estimator would look far better than it is.
##
## WHY THE INPUT SUPPORT IS AN ARGUMENT.  theta and beta are distances to the
## technology as defined for ALL x > 0, and reaching them means moving inputs
## DOWN.  A DEA hull cannot extrapolate below min(x), so an estimand whose
## projection leaves the sample's input support is not identified from that
## sample however large it gets -- the mean squared error hits a floor and the
## convergence rate reads as zero.  Widening `x_range` below the range of
## interest gives the reference sample room underneath, which is what makes
## input-oriented and two-sided estimands testable at all.  Output-oriented
## estimands never need it: they move y up at fixed x.
## ---------------------------------------------------------------------------

dea_sim <- function(n, p = 1L, q = 1L,
                    returns = 1,
                    ineff = c("exp", "hnorm"),
                    mean_ineff = 0.3,
                    x_range = c(1, 2),
                    seed = NULL) {

  ineff <- .match_arg_ci(ineff, c("exp", "hnorm"), "ineff")
  if (!is.numeric(n) || length(n) != 1L || n < 1) {
    stop("`n` must be a single positive number of DMUs.", call. = FALSE)
  }
  n <- as.integer(n); p <- as.integer(p); q <- as.integer(q)
  if (p < 1L || q < 1L) stop("`p` and `q` must be at least 1.", call. = FALSE)
  if (!is.numeric(returns) || length(returns) != 1L || returns <= 0 || returns > 1) {
    stop("`returns` is the elasticity of scale of the true frontier and must ",
         "lie in (0, 1]: 1 gives a constant-returns cone, below 1 a strictly ",
         "concave technology on which a crs fit is INCONSISTENT by ",
         "construction.", call. = FALSE)
  }
  if (!is.numeric(mean_ineff) || length(mean_ineff) != 1L || mean_ineff <= 0) {
    stop("`mean_ineff` must be a single positive number: E[u].", call. = FALSE)
  }
  if (!is.numeric(x_range) || length(x_range) != 2L || any(!is.finite(x_range)) ||
      x_range[1] <= 0 || x_range[2] <= x_range[1]) {
    stop("`x_range` must be two increasing positive numbers, the support of ",
         "each input.", call. = FALSE)
  }

  if (!is.null(seed)) {
    old <- .rng_snapshot(); on.exit(.rng_restore(old), add = TRUE)
    set.seed(seed)
  }

  X <- matrix(stats::runif(n * p, x_range[1], x_range[2]), n, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))
  fx <- exp(rowSums(log(X)) * (returns / p))

  ## u >= 0, with E[u] = mean_ineff under either family, so the two designs are
  ## comparable at the same nominal inefficiency.
  u <- if (ineff == "exp") stats::rexp(n, rate = 1 / mean_ineff)
       else abs(stats::rnorm(n, 0, mean_ineff * sqrt(pi / 2)))

  ## Uniform direction on the non-negative unit sphere.
  D <- matrix(abs(stats::rnorm(n * q)), n, q)
  D <- D / sqrt(rowSums(D^2))
  Y <- D * (fx * exp(-u))
  colnames(Y) <- paste0("y", seq_len(q))

  structure(list(
    x = X, y = Y,
    data = as.data.frame(cbind(X, Y)),
    theta = exp(-u / returns),   ## true input-oriented Farrell efficiency
    phi   = exp(u),              ## true output-oriented Farrell efficiency
    u = u, frontier = fx,
    design = list(n = n, p = p, q = q, returns = returns, ineff = ineff,
                  mean_ineff = mean_ineff, x_range = x_range, seed = seed)
  ), class = "dea_sim")
}

print.dea_sim <- function(x, ...) {
  d <- x$design
  cat("--- Simulated DEA design ---\n")
  cat("n = ", d$n, "   inputs = ", d$p, "   outputs = ", d$q,
      "   elasticity of scale = ", d$returns, "\n", sep = "")
  cat("input support: [", d$x_range[1], ", ", d$x_range[2], "]\n", sep = "")
  cat("inefficiency: ", d$ineff, ", E[u] = ", d$mean_ineff, "\n", sep = "")
  cat("true input efficiency theta:  ",
      paste(round(stats::quantile(x$theta, c(0, .5, 1)), 4), collapse = " / "),
      "   (min / median / max)\n", sep = "")
  ## The rate the DEA estimator can attain here, which is what the design is
  ## usually built to test.  See ?dea_rate.
  cat("attainable MSE slope, vrs: ", round(dea_rate(d$p, d$q, "vrs"), 3),
      "   crs: ", round(dea_rate(d$p, d$q, "crs"), 3), "\n", sep = "")
  invisible(x)
}

## ---------------------------------------------------------------------------
## dea_rate() -- the rate a DEA estimator can attain, as an MSE slope.
##
## DEA is consistent but NOT root-n consistent, and the rate depends on the
## dimension of the problem.  For the input-oriented estimator with p inputs
## and q outputs,
##
##   |theta_hat - theta| = O_p( n^(-2/(p+q+1)) )   variable returns
##                         O_p( n^(-2/(p+q))   )   constant returns
##                         O_p( n^(-1/(p+q))   )   free disposal hull
##
## (Kneip, Park and Simar 1998; Park, Simar and Weiner 2000; Park, Simar and
## Weiner 2000 for FDH.)  Squaring doubles the exponent, so regressing
## log MSE on log n over a grid of sample sizes has slope
##
##   -4/(p+q+1),  -4/(p+q),  -2/(p+q)
##
## respectively.  That number is the honest benchmark to test a DEA
## implementation against, and it is NOT -1: the familiar root-n slope is
## reached only at p + q = 3 under variable returns.  At p + q = 8 the slope is
## -4/9, meaning a hundredfold increase in sample size buys about a factor of
## 8 in MSE where a parametric estimator would buy 100.  That is the curse of
## dimensionality stated as something you can measure.
## ---------------------------------------------------------------------------
dea_rate <- function(p, q, rts = c("vrs", "crs", "nirs", "ndrs", "fdh"),
                     what = c("mse", "estimator")) {
  rts <- .match_arg_ci(rts, c("vrs", "crs", "nirs", "ndrs", "fdh"), "rts")
  what <- .match_arg_ci(what, c("mse", "estimator"), "what")
  d <- p + q
  ## nirs/ndrs impose one inequality where vrs imposes an equality; the rate is
  ## the vrs rate, since the binding constraint set has the same dimension.
  r <- switch(rts,
    crs = 2 / d,
    fdh = 1 / d,
    2 / (d + 1))
  if (what == "estimator") -r else -2 * r
}

## Save and restore the caller's random number stream, so that a seeded call
## does not silently reset the session's RNG the way set.seed() alone would.
.rng_snapshot <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else NULL
}
.rng_restore <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}
