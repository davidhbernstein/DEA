## ---------------------------------------------------------------------------
## dea_ddf() -- the directional distance function of Chambers, Chung and Fare
## (1996).
##
##   beta_o = max beta  s.t.  X'lambda <= x_o - beta g_x
##                            Y'lambda >= y_o + beta g_y
##
## The radial measures are the special cases g = (x_o, 0) and g = (0, y_o):
## with the first, beta = 1 - theta; with the second, beta = phi - 1.  What the
## general form buys is threefold and each part matters in practice.
##
## 1. INPUTS AND OUTPUTS MOVE TOGETHER.  A radial measure has to choose a side.
##    g = (x_o, y_o) contracts inputs and expands outputs at once, which is the
##    question most efficiency studies actually ask.
##
## 2. ZEROS AND NEGATIVES ARE FINE.  beta is an ADDITIVE contraction, so
##    nothing is divided by a DMU's own level.  A zero output makes the radial
##    output score undefined and the SBM undefined; here it is just a number.
##    Negative values -- net income, net exports, a pollutant measured as an
##    abatement credit -- work under a fixed direction for the same reason.
##
## 3. THE DIRECTION IS AN ASSUMPTION YOU CAN SEE.  A radial model also picks a
##    direction; it just does not say so.
##
## UNITS.  beta is scale invariant when g is proportional to the DMU's own
## data, and NOT when g is fixed -- doubling an input's units halves the beta
## needed to move one unit of g.  That is a property of the estimand, not a
## defect, but it means a fixed direction must be stated in the same units as
## the data.  When `scaling = TRUE` the direction is rescaled with the data, so
## the reported beta is the one the caller's own units imply.
##
## Variable layout: beta (1), lambda (n).  beta is left FREE: for an observed
## DMU beta >= 0 always, since beta = 0 is feasible, but the same program
## evaluated at a point outside the technology has a genuinely negative
## solution and clamping it at zero would report that point as efficient.
## ---------------------------------------------------------------------------

dea_ddf <- function(x, y, data = NULL,
                    direction = c("both", "in", "out", "unit", "mean"),
                    rts = c("vrs", "crs", "nirs", "ndrs"),
                    peers = TRUE,
                    scaling = TRUE,
                    xref = NULL, yref = NULL, dataref = NULL) {

  call <- match.call()
  rts <- .match_arg_ci(rts, c("vrs", "crs", "nirs", "ndrs"), "rts")
  t0 <- proc.time()[["elapsed"]]

  ## Deliberately NOT .dea_check()'s non-negativity rule -- handling data of
  ## either sign is the reason this estimator exists.
  d <- .dea_data(x, y, data, xref, yref, dataref, FALSE, "dea_ddf",
                 allow_negative = TRUE)
  X <- d$X; Y <- d$Y; XR <- d$XR; YR <- d$YR
  n <- nrow(X); nr <- nrow(XR); p <- ncol(X); q <- ncol(Y)
  if (nr < 3 * (p + q)) {
    warning("Only ", nr, " reference DMUs for ", p, " inputs and ", q,
            " outputs; with n < 3(p+q) most DMUs are efficient by dimension ",
            "alone.", call. = FALSE)
  }
  dmu <- d$dmu

  G <- .ddf_direction(direction, X, Y, n, p, q)
  sc  <- .dea_scale(XR, YR, scaling && all(XR > 0) && all(YR > 0))
  XRs <- sc$X; YRs <- sc$Y
  Xs  <- sweep(X, 2L, sc$sx, "/"); Ys <- sweep(Y, 2L, sc$sy, "/")
  ## Rescaling the direction alongside the data is what keeps beta unchanged.
  Gs <- cbind(sweep(G[, seq_len(p), drop = FALSE], 2L, sc$sx, "/"),
              sweep(G[, p + seq_len(q), drop = FALSE], 2L, sc$sy, "/"))

  rr <- .rts_row(rts)
  nrows <- p + q + (!is.null(rr))
  lp <- lpSolveAPI::make.lp(nrows, nr + 1L)
  lpSolveAPI::lp.control(lp, sense = "max", epsel = .DEA_CONSTANTS$LP_EPSEL,
                         verbose = "neutral")
  rts_one <- if (is.null(rr)) numeric(0) else 1
  for (j in seq_len(nr)) lpSolveAPI::set.column(lp, j + 1L, c(XRs[j, ], YRs[j, ], rts_one))
  lpSolveAPI::set.constr.type(lp, rep("<=", p), seq_len(p))
  lpSolveAPI::set.constr.type(lp, rep(">=", q), p + seq_len(q))
  if (!is.null(rr)) {
    lpSolveAPI::set.constr.type(lp, rr, p + q + 1L)
    lpSolveAPI::set.rhs(lp, 1, p + q + 1L)
  }
  lpSolveAPI::set.bounds(lp, lower = -.LP_INF, upper = .LP_INF, columns = 1L)

  beta <- numeric(n); st <- integer(n)
  L <- if (peers) matrix(0, n, nr) else NULL
  for (o in seq_len(n)) {
    lpSolveAPI::set.column(lp, 1L,
      c(1, Gs[o, seq_len(p)], -Gs[o, p + seq_len(q)]),
      c(0L, seq_len(p), p + seq_len(q)))
    lpSolveAPI::set.rhs(lp, c(Xs[o, ], Ys[o, ]), seq_len(p + q))
    st[o] <- solve(lp)
    if (!st[o] %in% c(0L, 1L)) { beta[o] <- NA_real_; next }
    beta[o] <- lpSolveAPI::get.objective(lp)
    if (peers) L[o, ] <- lpSolveAPI::get.variables(lp)[-1L]
  }
  beta[is.finite(beta) & abs(beta) < .DEA_CONSTANTS$TOL_EFF] <- 0
  names(beta) <- dmu
  if (!is.null(L)) { L[abs(L) < .DEA_CONSTANTS$TOL_LAMBDA] <- 0; dimnames(L) <- list(dmu, d$ref) }
  dimnames(G) <- list(dmu, c(colnames(X), colnames(Y)))

  structure(list(
    ## `eff` on this object is beta, and beta measures INEFFICIENCY: 0 is
    ## efficient and larger is worse, the opposite of every other model here.
    ## Keeping the field name shared lets print/summary/peers work unchanged;
    ## the `model` field is what tells them how to read it.
    eff = beta, beta = beta,
    model = "ddf", rts = rts, orientation = "directional",
    direction = if (is.character(direction)) direction[1] else "user",
    g = G, super = FALSE,
    lambda = L, sum_lambda = if (is.null(L)) NULL else rowSums(L),
    slack_x = NULL, slack_y = NULL,
    efficient = beta == 0, status = st,
    x = X, y = Y, xref = XR, yref = YR, self_ref = d$self,
    dmu = dmu, ref = d$ref, n = n, nref = nr, p = p, q = q,
    scaling = scaling, slack = FALSE,
    total_time = proc.time()[["elapsed"]] - t0,
    call = call
  ), class = "dea")
}

## Build the n x (p+q) direction matrix.  A character shorthand expands to one
## row per DMU; a numeric vector is one fixed direction shared by all of them;
## a matrix is taken as given.
.ddf_direction <- function(direction, X, Y, n, p, q) {
  if (is.character(direction)) {
    direction <- .match_arg_ci(direction, c("both", "in", "out", "unit", "mean"),
                               "direction")
    G <- switch(direction,
      both = cbind(X, Y),
      `in` = cbind(X, matrix(0, n, q)),
      out  = cbind(matrix(0, n, p), Y),
      unit = matrix(1, n, p + q),
      mean = matrix(rep(c(colMeans(X), colMeans(Y)), each = n), n, p + q))
  } else if (is.matrix(direction) || is.data.frame(direction)) {
    G <- as.matrix(direction)
    if (nrow(G) != n || ncol(G) != p + q) {
      stop("`direction` given as a matrix must be ", n, " x ", p + q,
           " (one row per DMU, inputs then outputs); got ", nrow(G), " x ",
           ncol(G), ".", call. = FALSE)
    }
    storage.mode(G) <- "double"
  } else if (is.numeric(direction)) {
    if (length(direction) != p + q) {
      stop("`direction` given as a vector must have length p + q = ", p + q,
           " (inputs then outputs); got ", length(direction), ".", call. = FALSE)
    }
    G <- matrix(rep(direction, each = n), n, p + q)
  } else {
    stop("`direction` must be one of \"both\", \"in\", \"out\", \"unit\", ",
         "\"mean\", or a numeric vector/matrix.", call. = FALSE)
  }
  if (any(!is.finite(G))) stop("`direction` contains non-finite values.", call. = FALSE)
  if (any(rowSums(abs(G)) <= 0)) {
    bad <- which(rowSums(abs(G)) <= 0)
    stop("`direction` is all zero for DMU(s) ", paste(utils::head(bad, 5), collapse = ", "),
         ". A zero direction has no distance to measure along; this usually ",
         "means direction = \"out\" was used on a DMU whose outputs are all 0.",
         call. = FALSE)
  }
  if (any(G < 0)) {
    warning("`direction` has negative entries. A direction vector points the ",
            "way improvement lies, so a negative input entry asks for the ",
            "input to GROW; check this is intended.", call. = FALSE)
  }
  G
}
