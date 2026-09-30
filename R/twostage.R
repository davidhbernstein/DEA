## ---------------------------------------------------------------------------
## dea_reg() -- the two-stage procedure of Simar and Wilson (2007).
##
## WHY THIS EXISTS. Regressing DEA scores on covariates is one of the most
## common things done with DEA and one of the most common things done wrong.
## The usual route -- OLS, or Tobit because the scores are bounded -- is
## invalid, and not marginally: the scores share a frontier, so they are
## correlated by construction in a way that depends on the sample, and
## conventional standard errors do not describe their sampling distribution.
## Tobit does not fix it. Censoring is a statement about how a variable is
## OBSERVED; the mass at 1 here is an estimation artefact of a boundary
## estimator, not a censored observation of something else.
##
## Simar and Wilson's answer is a TRUNCATED regression -- the efficient DMUs
## are dropped rather than censored, because delta = 1 carries no information
## about z once it is known the frontier was hit -- plus a bootstrap for the
## standard errors.
##
## THE SCALE IS SHEPHARD'S, NOT FARRELL'S. The regression is in
## delta >= 1, with delta = 1/theta input-oriented and delta = phi output-
## oriented, truncated from the left at 1. So a POSITIVE coefficient means the
## covariate makes a DMU LESS efficient. That sign is the opposite of what a
## reader expects from a regression on theta, it is what the paper reports,
## and print() says so on every fit rather than leaving it to be rediscovered.
##
## ALGORITHM #1 is the single bootstrap: fit the truncated regression on the
## inefficient DMUs, then draw parametrically from the fitted model and refit,
## L2 times, to get the sampling distribution of beta. It treats delta_hat as
## if it were delta.
##
## ALGORITHM #2 adds the step that makes the first one honest: delta_hat is a
## BIASED estimate of delta, biased toward 1, so a regression on delta_hat
## regresses on something systematically shifted. The inner loop builds
## pseudo-samples whose true distances are known, re-estimates against them,
## and bias-corrects delta_hat before the regression is run at all. It costs
## L1 DEA sweeps on top of the L2 regressions and is the procedure the paper
## recommends.
##
## WHAT IS NOT BOOTSTRAPPED. The covariates. z is taken as fixed, exactly as
## in the paper; the resampling is of the error term and, in #2, of the
## frontier. A package that also resampled rows would be answering a different
## question.
## ---------------------------------------------------------------------------

dea_reg <- function(object, z, data = NULL,
                    algorithm = c("sw2", "sw1"),
                    L1 = 100L, L2 = 2000L, alpha = 0.05,
                    seed = NULL, progress = interactive()) {
  call <- match.call()
  algorithm <- .match_arg_ci(algorithm, c("sw2", "sw1"), "algorithm")
  if (!inherits(object, "dea"))
    stop("`object` must be a fit from dea().", call. = FALSE)
  if (!identical(object$model, "radial"))
    stop("`object` is a ", object$model, " fit; the Simar-Wilson two-stage ",
         "procedure is defined for the radial efficiency estimator.",
         call. = FALSE)
  if (isTRUE(object$super))
    stop("A super-efficiency score is not an estimator of a distance to the ",
         "true frontier, so there is nothing here for this procedure to ",
         "regress or to bias-correct.", call. = FALSE)
  if (identical(object$rts, "fdh"))
    stop("The pseudo-sample step assumes the convex estimator; under free ",
         "disposal the relevant resampling scheme is a different one and is ",
         "not implemented here.", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be one number in (0, 1).", call. = FALSE)
  L1 <- as.integer(L1); L2 <- as.integer(L2)
  if (!is.null(seed)) set.seed(seed)

  n <- object$n
  Z <- .reg_z(z, data, n)
  delta <- .reg_delta(object)
  if (anyNA(delta))
    stop("The fit has ", sum(is.na(delta)), " unsolved DMU(s). Regressing on ",
         "a score vector with holes in it would silently drop them; resolve ",
         "them first.", call. = FALSE)

  ## The truncated regression uses the INEFFICIENT DMUs only. An efficient one
  ## contributes delta_hat = 1 exactly, which says the frontier was hit and
  ## nothing about z -- keeping it as an observation at the truncation point is
  ## what makes the Tobit version wrong.
  ineff <- delta > 1 + .DEA_CONSTANTS$TOL_EFF
  if (sum(ineff) < ncol(Z) + 2L)
    stop("Only ", sum(ineff), " DMU(s) are inefficient, against ", ncol(Z),
         " coefficient(s) to estimate. The truncated regression is fitted on ",
         "the inefficient DMUs alone; this sample cannot support it.",
         call. = FALSE)

  fit1 <- .truncreg(delta[ineff], Z[ineff, , drop = FALSE], ll = 1)

  bc <- NULL
  if (identical(algorithm, "sw2")) {
    bc <- .reg_bias_correct(object, Z, delta, fit1, L1, progress)
    ## The corrected scores can land below 1, which the truncated regression
    ## has no room for: delta is a distance to a frontier and 1 is its floor.
    ## Simar and Wilson's step 4 refits on the corrected scores that are still
    ## above it.
    delta_use <- bc$delta_bc
    ineff <- delta_use > 1 + .DEA_CONSTANTS$TOL_EFF
    if (sum(ineff) < ncol(Z) + 2L)
      stop("After bias correction only ", sum(ineff), " DMU(s) remain above ",
           "the truncation point, which cannot support ", ncol(Z),
           " coefficient(s). Algorithm #1 does not need this step.",
           call. = FALSE)
    fit2 <- .truncreg(delta_use[ineff], Z[ineff, , drop = FALSE], ll = 1)
  } else {
    delta_use <- delta
    fit2 <- fit1
  }

  boot <- .reg_boot(delta_use[ineff], Z[ineff, , drop = FALSE], fit2, L2,
                    progress)
  est <- c(fit2$beta, sigma = fit2$sigma)
  ci <- t(apply(boot, 2L, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2),
                na.rm = TRUE))
  colnames(ci) <- c("lower", "upper")
  se <- apply(boot, 2L, stats::sd, na.rm = TRUE)
  bias <- colMeans(boot, na.rm = TRUE) - est

  tab <- data.frame(term = names(est), estimate = as.numeric(est),
                    se = as.numeric(se), lower = ci[, 1], upper = ci[, 2],
                    boot_bias = as.numeric(bias), stringsAsFactors = FALSE)
  rownames(tab) <- NULL

  structure(list(
    table = tab, coefficients = est, ci = ci, boot = boot,
    delta = delta, delta_bc = if (is.null(bc)) NULL else bc$delta_bc,
    bias = if (is.null(bc)) NULL else bc$bias,
    used = ineff, n = n, n_used = sum(ineff),
    algorithm = algorithm, L1 = if (identical(algorithm, "sw2")) L1 else NA_integer_,
    L2 = L2, alpha = alpha, orientation = object$orientation, rts = object$rts,
    stage1 = fit1, stage2 = fit2, z_names = colnames(Z),
    converged = c(stage1 = fit1$converged, stage2 = fit2$converged),
    call = call
  ), class = "dea_reg")
}

## Shephard's distance, which is what the regression is in. Both orientations
## map onto delta >= 1; see the header on why the sign of a coefficient then
## reads backwards from the usual expectation.
.reg_delta <- function(object) {
  e <- as.numeric(object$eff)
  if (identical(object$orientation, "in")) 1 / e else e
}

.reg_z <- function(z, data, n) {
  Z <- .dea_matrix(z, data, "z")
  if (nrow(Z) != n)
    stop("`z` has ", nrow(Z), " rows but the fit has ", n, " DMUs.",
         call. = FALSE)
  ## An intercept is added unless the caller already supplied a constant
  ## column, because a truncated regression without one is a strong and almost
  ## always unintended restriction.
  const <- apply(Z, 2L, function(v) stats::sd(v) < .Machine$double.eps^0.5)
  if (!any(const)) Z <- cbind(`(Intercept)` = 1, Z)
  Z
}

## ---------------------------------------------------------------------------
## The truncated-normal regression, by maximum likelihood.
##
## y_i = z_i'beta + eps_i, eps_i ~ N(0, sigma^2), observed only where
## y_i > ll. Written here rather than taken from a package because the two
## that provide it are DEA packages' dependencies or modelling packages, and
## this one declares neither -- the same reasoning as the reference values in
## tests/. It is checked against truncreg::truncreg and npsf::truncreg in
## tools/make_reference_values.R.
##
## Parameterised in log(sigma) so the optimiser is unconstrained, with the
## analytic gradient supplied: the numerical one is workable but the tail
## terms lose precision exactly where the truncation bites, which is the whole
## sample here.
##
## The tail is computed in logs. 1 - Phi(a) underflows to zero for a beyond
## about 8, which on this scale is an ordinary DMU far from the boundary, and
## the likelihood then returns Inf rather than a number.
## ---------------------------------------------------------------------------
.truncreg <- function(y, Z, ll = 1, start = NULL) {
  k <- ncol(Z)
  nll <- function(par) {
    b <- par[seq_len(k)]; s <- exp(par[k + 1L])
    mu <- as.numeric(Z %*% b)
    r <- (y - mu) / s
    a <- (ll - mu) / s
    ltail <- stats::pnorm(a, lower.tail = FALSE, log.p = TRUE)
    ## The 0.5*log(2*pi) term does not move the optimum and is carried anyway,
    ## so that $logLik is the log-likelihood other packages report rather than
    ## that number plus n*0.9189.
    v <- -par[k + 1L] - 0.5 * r^2 - 0.5 * log(2 * pi) - ltail
    if (any(!is.finite(v))) return(1e10)
    -sum(v)
  }
  gr <- function(par) {
    b <- par[seq_len(k)]; s <- exp(par[k + 1L])
    mu <- as.numeric(Z %*% b)
    r <- (y - mu) / s
    a <- (ll - mu) / s
    lam <- exp(stats::dnorm(a, log = TRUE) -
               stats::pnorm(a, lower.tail = FALSE, log.p = TRUE))
    if (any(!is.finite(lam))) return(rep(0, k + 1L))
    g_b <- as.numeric(crossprod(Z, r - lam)) / s
    g_t <- sum(r^2 - 1 - lam * a)
    -c(g_b, g_t)
  }
  if (is.null(start)) {
    ols <- stats::lm.fit(Z, y)
    start <- c(ols$coefficients, log(max(stats::sd(ols$residuals), 1e-4)))
    start[!is.finite(start)] <- 0
  }
  o <- stats::optim(start, nll, gr, method = "BFGS",
                    control = list(maxit = 500L, reltol = 1e-12))
  beta <- o$par[seq_len(k)]
  names(beta) <- colnames(Z)
  list(beta = beta, sigma = exp(o$par[k + 1L]), logLik = -o$value,
       converged = identical(o$convergence, 0L), ll = ll, par = o$par,
       n = length(y))
}

## Draws from N(mu, sigma^2) conditioned on being above `lo`, by the inverse
## CDF. Done in log space for the same underflow reason as the likelihood: for
## a DMU whose fitted value sits far above the truncation point the naive
## p + (1-p) * u loses every digit of u.
.rtnorm_left <- function(mu, sigma, lo) {
  a <- (lo - mu) / sigma
  lp <- stats::pnorm(a, lower.tail = FALSE, log.p = TRUE)
  u <- stats::runif(length(mu))
  mu + sigma * stats::qnorm(log(u) + lp, lower.tail = FALSE, log.p = TRUE)
}

## The parametric bootstrap of beta: redraw the errors from the fitted model,
## keeping z fixed, and refit. This is steps 5-7 of Algorithm #1 and steps
## 5-7 of #2; they differ only in which fit they start from.
.reg_boot <- function(y, Z, fit, L2, progress) {
  k <- ncol(Z)
  out <- matrix(NA_real_, L2, k + 1L,
                dimnames = list(NULL, c(colnames(Z), "sigma")))
  mu <- as.numeric(Z %*% fit$beta)
  pb <- .reg_progress(progress, L2, "bootstrap")
  for (b in seq_len(L2)) {
    ys <- .rtnorm_left(mu, fit$sigma, fit$ll)
    f <- try(.truncreg(ys, Z, ll = fit$ll, start = fit$par), silent = TRUE)
    if (!inherits(f, "try-error") && f$converged)
      out[b, ] <- c(f$beta, f$sigma)
    pb(b)
  }
  keep <- rowSums(is.na(out)) == 0
  if (!any(keep))
    stop("Every bootstrap replication failed to converge; the fitted model is ",
         "too close to degenerate to resample from.", call. = FALSE)
  if (sum(keep) < L2)
    warning(L2 - sum(keep), " of ", L2, " bootstrap replication(s) did not ",
            "converge and are dropped. The interval is from the ", sum(keep),
            " that did.", call. = FALSE)
  out[keep, , drop = FALSE]
}

## Algorithm #2, step 3: the inner loop that bias-corrects delta_hat.
##
## Each pass draws a delta* from the fitted truncated model, builds a
## pseudo-sample whose true distance is that delta*, and re-estimates the
## ORIGINAL points against the pseudo-technology. Scaling a DMU's inputs by a
## factor scales its input distance by the same factor, which is what makes the
## construction exact rather than approximate.
##
## As in dea_boot(): the pseudo-sample is the TECHNOLOGY and the points
## EVALUATED against it are the originals. Getting that backwards produces a
## number for every DMU and measures nothing.
.reg_bias_correct <- function(object, Z, delta, fit, L1, progress) {
  n <- object$n
  X <- object$x; Y <- object$y
  inp <- identical(object$orientation, "in")
  mu <- as.numeric(Z %*% fit$beta)
  acc <- matrix(NA_real_, n, L1)
  pb <- .reg_progress(progress, L1, "bias correction")
  for (b in seq_len(L1)) {
    ds <- .rtnorm_left(mu, fit$sigma, fit$ll)
    k <- ds / delta
    Xs <- if (inp) X * k else X
    Ys <- if (inp) Y else Y / k
    f <- suppressWarnings(dea(X, Y, rts = object$rts,
                              orientation = object$orientation,
                              slack = FALSE, peers = FALSE,
                              scaling = object$scaling, xref = Xs, yref = Ys))
    acc[, b] <- .reg_delta(f)
    pb(b)
  }
  mstar <- rowMeans(acc, na.rm = TRUE)
  bias <- mstar - delta
  list(delta_bc = delta - bias, bias = bias, draws = acc)
}

.reg_progress <- function(on, total, what) {
  if (!isTRUE(on)) return(function(i) invisible(NULL))
  step <- max(1L, total %/% 20L)
  function(i) {
    if (i %% step == 0L || i == total) {
      cat("\r  ", what, ": ", i, "/", total, sep = "")
      utils::flush.console()
      if (i == total) cat("\n")
    }
    invisible(NULL)
  }
}

print.dea_reg <- function(x, ...) {
  cat("--- Two-stage efficiency regression (Simar and Wilson 2007, ",
      if (x$algorithm == "sw2") "Algorithm #2" else "Algorithm #1", ") ---\n",
      sep = "")
  cat("DMUs: ", x$n, "   used (inefficient): ", x$n_used,
      "   rts: ", x$rts, "   orientation: ", x$orientation, "\n", sep = "")
  cat("bootstrap: L2 = ", x$L2,
      if (x$algorithm == "sw2") paste0(", L1 = ", x$L1) else "",
      "   interval: ", round(100 * (1 - x$alpha)), "%\n", sep = "")
  cat("\n")
  tb <- x$table
  tb[-1L] <- lapply(tb[-1L], function(v) round(v, 5))
  print(tb, row.names = FALSE)
  cat("\nThe response is Shephard's distance, delta >= 1, so a POSITIVE\n",
      "coefficient means the covariate is associated with LESS efficiency.\n",
      sep = "")
  if (!all(x$converged))
    cat("\nWARNING: a stage did not converge; see $converged.\n")
  invisible(x)
}

summary.dea_reg <- function(object, ...) {
  print(object)
  if (!is.null(object$bias)) {
    cat("\nbias correction of the scores (Algorithm #2)\n")
    print(round(summary(object$bias), 5))
    cat("\ndelta_hat vs bias-corrected\n")
    print(round(rbind(delta_hat = summary(object$delta),
                      corrected = summary(object$delta_bc)), 4))
  }
  invisible(object)
}

nobs.dea_reg <- function(object, ...) object$n_used
coef.dea_reg <- function(object, ...) object$coefficients
