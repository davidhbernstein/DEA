## ---------------------------------------------------------------------------
## dea() -- radial (Debreu-Farrell) technical efficiency.
##
## One function covers the whole radial family, because CCR, BCC, NIRS, NDRS
## and FDH are the SAME linear program under five different restrictions on
## sum(lambda).  Splitting them into five exported functions, as several
## packages do, hides that and makes the returns-to-scale comparison -- which
## needs three of them fitted on the same data -- look like three unrelated
## calls instead of one sweep.
##
## THE REFERENCE SET.  `xref`/`yref` separate the DMUs being SCORED from the
## DMUs that SPAN the technology.  They default to each other, which is
## ordinary DEA, and pulling them apart is what makes possible:
##
##   * scoring new units against a frontier estimated from a fixed sample;
##   * comparing a group against another group's frontier (metafrontiers);
##   * the Monte Carlo designs in which the estimator's convergence RATE is
##     measured, which require a fixed interior evaluation point and a growing
##     reference sample -- scoring every sample point against itself mixes
##     interior and boundary DMUs, and the boundary ones converge more slowly
##     and dominate the average.
## ---------------------------------------------------------------------------

dea <- function(x, y, data = NULL,
                rts = c("vrs", "crs", "nirs", "ndrs", "fdh"),
                orientation = c("in", "out"),
                slack = TRUE,
                super = FALSE,
                multipliers = FALSE,
                peers = TRUE,
                scaling = TRUE,
                xref = NULL, yref = NULL, dataref = NULL) {

  call <- match.call()
  rts <- .match_arg_ci(rts, c("vrs", "crs", "nirs", "ndrs", "fdh"), "rts")
  orientation <- .match_arg_ci(orientation, c("in", "out"), "orientation")
  t0 <- proc.time()[["elapsed"]]

  d <- .dea_data(x, y, data, xref, yref, dataref, super, "dea")
  X <- d$X; Y <- d$Y; XR <- d$XR; YR <- d$YR
  n <- nrow(X); nr <- nrow(XR); p <- ncol(X); q <- ncol(Y)

  ## Scale by the REFERENCE column means, and apply the same factors to the
  ## evaluated points: the two sets have to sit in one coordinate system or the
  ## comparison between them is not the one that was asked for.
  sc  <- .dea_scale(XR, YR, scaling)
  XRs <- sc$X; YRs <- sc$Y
  Xs  <- sweep(X, 2L, sc$sx, "/"); Ys <- sweep(Y, 2L, sc$sy, "/")

  if (multipliers && identical(rts, "fdh")) {
    stop("multipliers = TRUE is not available for rts = \"fdh\". The free ",
         "disposal hull is not convex, so it is not the feasible set of a ",
         "linear program and has no supporting price vector -- there is no ",
         "multiplier form to return, rather than one that is merely ",
         "unimplemented.", call. = FALSE)
  }

  out <- if (identical(rts, "fdh")) {
    .dea_fdh(Xs, Ys, XRs, YRs, orientation, super, n, nr, p, q)
  } else {
    .dea_radial(Xs, Ys, XRs, YRs, rts, orientation, super, slack, peers, n, nr, p, q)
  }

  ## The multiplier sweep is a second family of programs over the same
  ## technology, run only when asked for: it doubles the solve time and most
  ## callers want the score, not the prices that support it.
  mult <- if (multipliers) {
    .dea_multipliers(Xs, Ys, XRs, YRs, rts, orientation, super, n, nr, p, q,
                     sc$sx, sc$sy, d, colnames(X), colnames(Y))
  } else NULL

  ## Snap to the boundary. A simplex solver returns 0.9999999998 where the
  ## answer is exactly 1, and leaving that in makes every downstream "is this
  ## DMU efficient" test a lottery.
  eff <- out$eff
  eff[is.finite(eff) & abs(eff - 1) < .DEA_CONSTANTS$TOL_EFF] <- 1
  names(eff) <- d$dmu

  sx <- out$sx; sy <- out$sy
  if (!is.null(sx)) {
    sx <- sweep(sx, 2L, sc$sx, "*")   ## back into the caller's units
    sy <- sweep(sy, 2L, sc$sy, "*")
    sx[abs(sx) < .DEA_CONSTANTS$TOL_SLACK] <- 0
    sy[abs(sy) < .DEA_CONSTANTS$TOL_SLACK] <- 0
    dimnames(sx) <- list(d$dmu, colnames(X)); dimnames(sy) <- list(d$dmu, colnames(Y))
  }
  L <- out$lambda
  if (!is.null(L)) {
    L[abs(L) < .DEA_CONSTANTS$TOL_LAMBDA] <- 0
    dimnames(L) <- list(d$dmu, d$ref)
  }

  structure(list(
    eff = eff, model = "radial", rts = rts, orientation = orientation,
    super = super,
    lambda = L, sum_lambda = out$sum_lambda,
    slack_x = sx, slack_y = sy,
    ## Pareto-Koopmans efficiency: radially on the frontier AND no slack left.
    efficient = if (is.null(sx)) eff == 1 else
                  eff == 1 & rowSums(sx) <= .DEA_CONSTANTS$TOL_SLACK &
                  rowSums(sy) <= .DEA_CONSTANTS$TOL_SLACK,
    status = out$status,
    ## The multiplier form: v on the inputs, u on the outputs, u0 the
    ## returns-to-scale intercept. NULL unless multipliers = TRUE.
    v = mult$v, u = mult$u, u0 = mult$u0, mult_status = mult$status,
    multipliers = multipliers,
    x = X, y = Y, xref = XR, yref = YR, self_ref = d$self,
    dmu = d$dmu, ref = d$ref, n = n, nref = nr, p = p, q = q,
    scaling = scaling, slack = slack,
    total_time = proc.time()[["elapsed"]] - t0,
    call = call
  ), class = "dea")
}

## The LP sweep, in scaled coordinates.  Kept separate from dea() so that
## dea_boot() can call it directly on resampled reference sets without
## re-validating the data B times.
.dea_radial <- function(Xs, Ys, XRs, YRs, rts, orientation, super, slack,
                        peers, n, nr, p, q) {
  B   <- .lp_radial_build(XRs, YRs, rts, orientation)
  eff <- numeric(n); st <- integer(n)
  L   <- if (peers) matrix(0, n, nr) else NULL
  suml <- numeric(n)

  for (o in seq_len(n)) {
    r <- .lp_radial_at(B, Xs, Ys, o, exclude = if (super) o else NULL)
    eff[o] <- r$eff; st[o] <- r$status
    suml[o] <- if (all(is.na(r$lambda))) NA_real_ else sum(r$lambda)
    if (peers) L[o, ] <- r$lambda
  }

  ## Stage two runs from the RADIAL projection, so it needs stage one's theta.
  ## Skipped under super-efficiency: the projection there is onto a frontier
  ## the DMU is not part of, and the residual slack is not a Pareto-Koopmans
  ## statement about the observed technology.
  sx <- sy <- NULL
  if (slack && !super) {
    S  <- .lp_slack_build(XRs, YRs, rts)
    sx <- matrix(0, n, p); sy <- matrix(0, n, q)
    for (o in seq_len(n)) {
      if (!is.finite(eff[o])) { sx[o, ] <- NA_real_; sy[o, ] <- NA_real_; next }
      rhs_x <- if (orientation == "in") eff[o] * Xs[o, ] else Xs[o, ]
      rhs_y <- if (orientation == "in") Ys[o, ] else eff[o] * Ys[o, ]
      z <- .lp_slack_at(S, rhs_x, rhs_y)
      sx[o, ] <- z$sx; sy[o, ] <- z$sy
      ## Stage two's lambda is the one that survives a Pareto-Koopmans
      ## projection, so it replaces stage one's where both exist.
      if (peers && !any(is.na(z$lambda))) L[o, ] <- z$lambda
    }
  }
  list(eff = eff, lambda = L, sx = sx, sy = sy, status = st, sum_lambda = suml)
}

## The multiplier sweep, and the unscaling that has to go with it.
##
## SCALING IS NOT NEUTRAL FOR THE WEIGHTS, even though it is neutral for the
## score. Columns were divided by their reference means before solving, so a
## weight returned by the solver applies to x_i / sx_i and not to x_i; the
## weight the caller asked about is therefore v_i / sx_i. Returning the solver's
## numbers unchanged would give weights that satisfy v'x_o = 1 in units nobody
## supplied. u0 needs no such correction -- it multiplies the constant 1, which
## has no units.
.dea_multipliers <- function(Xs, Ys, XRs, YRs, rts, orientation, super,
                             n, nr, p, q, sx, sy, d, xnames, ynames) {
  M  <- .lp_mult_build(XRs, YRs, rts, orientation)
  V  <- matrix(NA_real_, n, p)
  U  <- matrix(NA_real_, n, q)
  u0 <- rep(NA_real_, n)
  st <- integer(n)
  for (o in seq_len(n)) {
    r <- .lp_mult_at(M, Xs, Ys, o, exclude = if (super) o else NULL)
    V[o, ] <- r$v; U[o, ] <- r$u; u0[o] <- r$u0; st[o] <- r$status
  }
  V <- sweep(V, 2L, sx, "/")
  U <- sweep(U, 2L, sy, "/")
  dimnames(V) <- list(d$dmu, xnames)
  dimnames(U) <- list(d$dmu, ynames)
  names(u0) <- d$dmu
  list(v = V, u = U, u0 = if (M$has0) u0 else NULL, status = st)
}

.dea_fdh <- function(Xs, Ys, XRs, YRs, orientation, super, n, nr, p, q) {
  if (super) {
    stop("super = TRUE is not implemented for rts = \"fdh\". Andersen-Petersen ",
         "super-efficiency removes the DMU from its own reference set; under ",
         "free disposal that leaves a DMU on the boundary of the input space ",
         "with no dominating peer at all, so the score is undefined rather ",
         "than merely infeasible for a large share of any sample.",
         call. = FALSE)
  }
  f <- .fdh_radial(Xs, Ys, XRs, YRs, orientation)
  L <- matrix(0, n, nr)
  ok <- is.finite(f$eff) & !is.na(f$peer)
  L[cbind(which(ok), f$peer[ok])] <- 1
  ## The FDH projection is onto the dominating DMU itself, so the remaining
  ## slack is available in closed form -- no second program needed.
  pk <- ifelse(is.na(f$peer), 1L, f$peer)
  projx <- if (orientation == "in") f$eff * Xs else Xs
  projy <- if (orientation == "in") Ys else f$eff * Ys
  sx <- projx - XRs[pk, , drop = FALSE]
  sy <- YRs[pk, , drop = FALSE] - projy
  sx[!ok, ] <- NA_real_; sy[!ok, ] <- NA_real_
  list(eff = f$eff, lambda = L, sx = sx, sy = sy,
       status = integer(n), sum_lambda = ifelse(ok, 1, NA_real_))
}

## ---------------------------------------------------------------------------
## Shared front end: coerce, validate, and work out whether the evaluated and
## reference sets are the same one.
## ---------------------------------------------------------------------------
.dea_data <- function(x, y, data, xref, yref, dataref, super, model,
                      require_positive = FALSE, allow_negative = FALSE) {
  .dea_check_data_arg(data)
  .dea_check_data_arg(dataref)
  X <- .dea_matrix(x, data, "x")
  Y <- .dea_matrix(y, data, "y")
  self <- is.null(xref) && is.null(yref)
  if (xor(is.null(xref), is.null(yref))) {
    stop("`xref` and `yref` must be given together: a reference technology ",
         "needs both its inputs and its outputs.", call. = FALSE)
  }
  if (self) { XR <- X; YR <- Y } else {
    if (is.null(dataref)) dataref <- data
    XR <- .dea_matrix(xref, dataref, "xref")
    YR <- .dea_matrix(yref, dataref, "yref")
    if (ncol(XR) != ncol(X) || ncol(YR) != ncol(Y)) {
      stop("The reference set has ", ncol(XR), " inputs and ", ncol(YR),
           " outputs, but the evaluated set has ", ncol(X), " and ", ncol(Y),
           ". They must describe the same variables in the same order.",
           call. = FALSE)
    }
  }
  if (super && !self) {
    stop("super = TRUE removes each DMU from its own reference set, which is ",
         "only meaningful when the evaluated and reference sets are the same. ",
         "With an explicit `xref`/`yref` the DMU is already outside the ",
         "technology it is scored against.", call. = FALSE)
  }
  if (!allow_negative) {
    .dea_check(XR, YR, require_positive, model)
    if (!self) .dea_check_eval(X, Y, require_positive, model)
  }
  list(X = X, Y = Y, XR = XR, YR = YR, self = self,
       dmu = .dea_dmu_names(X, Y), ref = .dea_dmu_names(XR, YR))
}
