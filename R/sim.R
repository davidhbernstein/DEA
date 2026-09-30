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
                    seed = NULL,
                    ## Appended after `seed` rather than placed next to
                    ## `returns`, where they belong logically: the package is
                    ## released, so inserting an argument mid-signature would
                    ## silently change the meaning of any positional call.
                    frontier = c("cobb", "ces", "facet"),
                    frontier_par = NULL) {

  ineff <- .match_arg_ci(ineff, c("exp", "hnorm"), "ineff")
  frontier <- .match_arg_ci(frontier, c("cobb", "ces", "facet"), "frontier")
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
  fr <- .sim_frontier(X, returns, frontier, frontier_par)
  fx <- fr$fx

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
    ## Which facet each DMU sits under, for frontier = "facet" and NULL
    ## otherwise. Reported rather than discarded because a facet design with
    ## only one active facet is a LINEAR design wearing m normals, and nothing
    ## else in the returned object would reveal that.
    facet = fr$facet, facet_normals = fr$A,
    design = list(n = n, p = p, q = q, returns = returns, ineff = ineff,
                  mean_ineff = mean_ineff, x_range = x_range, seed = seed,
                  shape = frontier, shape_par = fr$par)
  ), class = "dea_sim")
}

## ---------------------------------------------------------------------------
## The frontier shape.
##
## WHY ANOTHER GEOMETRY IS CHEAP, WHICH IS NOT OBVIOUS -- the ROADMAP costed
## this at two days, "most of it deriving and checking the distance functions".
## Almost none of that is needed, because both closed forms above survive ANY
## change of shape provided f stays HOMOGENEOUS OF DEGREE r:
##
##   phi   = f(x) / ||y||                            -- true for any f at all
##   theta = min{ t : f(tx) >= ||y|| }
##         = min{ t : t^r f(x) >= f(x) e^{-u} }  =  e^{-u/r}
##
## The second line uses homogeneity and nothing else: not smoothness, not
## strict concavity, and nothing Cobb-Douglas. So every shape offered here is
## built homogeneous of degree `returns`, the truth needs no re-derivation and
## no root-finding, and `test-sim.R` checks that claim against the estimator
## rather than taking it on trust.
##
## All three agree on the diagonal, f(c, ..., c) = c^r, so `returns` means the
## same thing across shapes and a convergence table can be read down a column.
##
##   cobb   f(x) = (prod_j x_j)^(r/p)
##          Smooth, strictly concave. The original, and still the default.
##
##   ces    f(x) = ((1/p) sum_j x_j^rho)^(r/rho)
##          Nests the other two: rho -> 0 IS cobb (the geometric mean is the
##          limit of the power mean), rho = 1 is an arithmetic mean and so a
##          flat frontier, rho -> -Inf is Leontief min(). So `rho` moves
##          curvature continuously, which is the point of it -- the estimator's
##          difficulty ought to vary with curvature, and a single shape cannot
##          show whether it does.
##
##   facet  f(x) = (min_k a_k'x / p)^r  over m facets
##          Piecewise linear: the shape a DEA hull can represent EXACTLY, so
##          the estimator should do BEST here, and the shape with flat faces,
##          so the optimal face is not a point and the projection is not
##          unique. That last is the geometry the convergence study was
##          missing and the one most likely to find something.
##
## rho > 1 is REFUSED rather than clamped. It makes f convex, so T is not
## convex; a vrs estimator then consistently estimates the convex hull of T
## instead of T. `theta` above would still be the true distance to T, and would
## no longer be what the estimator converges to -- the design would quietly
## score the estimator against the wrong number, which is worse than not
## offering the option.
##
## THE FACET NORMALS ARE DETERMINISTIC, and that is not housekeeping. Drawn
## from the RNG they would move with `seed`, so "the same geometry at n = 100
## and at n = 800" would be false and the convergence slope would be fitted
## across a moving target -- the one thing the study must not do. They come
## from a fixed golden-ratio phase pattern instead. A plain k/m by j/p lattice
## was tried first and is WRONG: with p = 2 and m = 3 it aliases, facets 2 and
## 3 coming out identical, so `m` would overstate the geometry. The irrational
## phase is what prevents that, and `.sim_facets()` checks for coincident rows
## anyway rather than trusting the argument.
##
## Every normal is scaled to sum to p, which puts every facet through the
## diagonal ray. The frontier is then a fan of m flat faces meeting along that
## ridge: each face wins in its own region, so all m are genuinely active, and
## near the diagonal they TIE -- which is the non-unique projection, present by
## construction rather than by luck.
## ---------------------------------------------------------------------------
.sim_frontier <- function(X, returns, frontier, par) {
  p <- ncol(X)

  if (frontier == "cobb") {
    if (!is.null(par)) {
      stop("`frontier_par` has no meaning for frontier = \"cobb\": the ",
           "Cobb-Douglas shape is fixed once `returns` is given. Use ",
           "frontier = \"ces\" for a tunable curvature.", call. = FALSE)
    }
    ## Left exactly as it was written, so the default design is bit-identical
    ## to every result recorded before the shape argument existed.
    return(list(fx = exp(rowSums(log(X)) * (returns / p)),
                par = NULL, facet = NULL, A = NULL))
  }

  if (frontier == "ces") {
    rho <- if (is.null(par)) 0.5 else par
    if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho)) {
      stop("`frontier_par` for frontier = \"ces\" is rho, a single finite ",
           "number.", call. = FALSE)
    }
    if (rho > 1) {
      stop("`frontier_par` = ", rho, " > 1 makes the CES frontier CONVEX, so ",
           "the technology is not convex and a vrs estimator converges to the ",
           "convex hull of it rather than to it. The closed-form `theta` this ",
           "function returns would still be the true distance and would no ",
           "longer be the estimand, so the design would score the estimator ",
           "against the wrong number. rho must be <= 1.", call. = FALSE)
    }
    if (abs(rho) < 1e-3) {
      stop("`frontier_par` = ", rho, " is numerically at the rho -> 0 limit, ",
           "where CES IS Cobb-Douglas; the power-mean formula loses all its ",
           "precision there. Use frontier = \"cobb\", which is that limit in ",
           "closed form.", call. = FALSE)
    }
    fx <- (rowMeans(X^rho))^(returns / rho)
    return(list(fx = fx, par = rho, facet = NULL, A = NULL))
  }

  ## --- facet ----------------------------------------------------------------
  if (p < 2L) {
    stop("frontier = \"facet\" needs at least 2 inputs. With p = 1 every ",
         "facet normal reduces to the same number, so the shape collapses to ",
         "f(x) = x^returns -- which is exactly frontier = \"cobb\" at p = 1. ",
         "Allowing it would put a duplicate row in a geometry comparison and ",
         "make the two shapes look like independent evidence.", call. = FALSE)
  }
  m <- if (is.null(par)) 3L else par
  if (!is.numeric(m) || length(m) != 1L || !is.finite(m) || m < 2 || m != round(m)) {
    stop("`frontier_par` for frontier = \"facet\" is the number of facets: a ",
         "single whole number >= 2. One facet is a flat frontier, for which ",
         "frontier = \"ces\" with rho = 1 is the honest name.", call. = FALSE)
  }
  m <- as.integer(m)
  ## With every normal scaled to sum to p, the normals are COLLINEAR when
  ## p = 2 -- they lie on a line in the sum = p plane -- and a'x is linear in a,
  ## so a minimum over collinear points is always attained at an endpoint. At
  ## p = 2, therefore, exactly two facets can ever be active however many are
  ## constructed, and the rest are dominated everywhere. Measured: m = 3, 4 and
  ## 6 at p = 2 all give 2 active facets, splitting the sample 51/49. Refused
  ## up front rather than warned about, because it is knowable from p and m
  ## alone and is always an overstatement of the geometry.
  if (p == 2L && m > 2L) {
    stop("frontier = \"facet\" with ", m, " facets and p = 2 inputs: only 2 ",
         "of them can ever be active. Every normal is scaled to sum to p, so ",
         "at p = 2 they are collinear, and min_k a_k'x over collinear a_k is ",
         "always attained at an endpoint -- the other ", m - 2L, " facet(s) ",
         "are dominated at every x. Use frontier_par = 2, or p >= 3 where ",
         "interior facets do bind.", call. = FALSE)
  }
  A <- .sim_facets(m, p)
  lin <- X %*% t(A) / p
  k <- apply(lin, 1L, which.min)
  ## The real guard on the geometry is which facets BIND, not whether the
  ## normals are distinct: a facet that is never the minimum contributes
  ## nothing, and `m` would then overstate the shape being studied.
  if (length(unique(k)) < m) {
    warning("frontier = \"facet\" asked for ", m, " facets but only ",
            length(unique(k)), " bind anywhere in this sample, so the ",
            "geometry is simpler than `frontier_par` says. Check $facet.",
            call. = FALSE)
  }
  list(fx = apply(lin, 1L, min)^returns, par = m, facet = k, A = A)
}

## Deterministic facet normals: m rows, p columns, every entry positive and
## every row summing to p.  The phase advances by the golden-ratio conjugate
## across inputs, which is irrational and so cannot come back into step with
## the k/m advance across facets -- see the aliasing note above.
.sim_facets <- function(m, p) {
  g <- (sqrt(5) - 1) / 2
  A <- outer(seq_len(m), seq_len(p),
             function(k, j) 1 + 0.5 * cos(2 * pi * ((k - 1) / m + (j - 1) * g)))
  A <- A * p / rowSums(A)
  ## Coincident rows would mean the design has fewer facets than it says. The
  ## construction is meant to prevent it; this is here so that a silent
  ## reduction is impossible rather than merely unlikely.
  if (m > 1L && min(stats::dist(A)) < 1e-8) {
    stop("the facet construction produced coincident normals at m = ", m,
         ", p = ", p, ", so the design would have fewer facets than ",
         "requested. Pick a different `frontier_par`.", call. = FALSE)
  }
  dimnames(A) <- list(paste0("facet", seq_len(m)), paste0("x", seq_len(p)))
  A
}

print.dea_sim <- function(x, ...) {
  d <- x$design
  cat("--- Simulated DEA design ---\n")
  cat("n = ", d$n, "   inputs = ", d$p, "   outputs = ", d$q,
      "   elasticity of scale = ", d$returns, "\n", sep = "")
  cat("input support: [", d$x_range[1], ", ", d$x_range[2], "]\n", sep = "")
  ## The shape is printed with the number of facets that actually BIND, not the
  ## number requested: those can differ, and the difference is the whole point
  ## of recording $facet.
  shape_txt <- switch(d$shape,
    cobb  = "Cobb-Douglas (smooth, strictly concave)",
    ces   = paste0("CES, rho = ", d$shape_par,
                   if (d$shape_par == 1) " (flat)" else ""),
    facet = paste0("piecewise linear, ", d$shape_par, " facets, ",
                   length(unique(x$facet)), " binding"))
  cat("frontier shape: ", shape_txt, "\n", sep = "")
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
