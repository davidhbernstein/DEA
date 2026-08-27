## ---------------------------------------------------------------------------
## dea_add() -- the additive model, and the weighted measures built on it.
##
##   max  sum_i w_i^- s_i^-  +  sum_r w_r^+ s_r^+
##   s.t. X'lambda + s^- = x_o,  Y'lambda - s^+ = y_o,  returns to scale
##
## THE MODEL IS THE PROJECTION, THE WEIGHTS ARE THE MEASURE.  All three options
## below solve the SAME program over the SAME technology and differ only in
## what the objective weights are.  That is worth being explicit about, because
## the literature presents RAM and MIP as separate "models" when they are one
## model read on three scales.
##
## WHY THE CLASSIC VERSION NEEDS A CAVEAT.  With w = 1 the objective adds a
## slack measured in staff-hours to one measured in euros, so the answer
## depends on the units the data happens to be recorded in -- double an input's
## units and its slack halves, changing which DMU looks worst.  The projection
## is still correct and the efficient SET is still right; it is the ORDERING of
## the inefficient DMUs that is not units invariant.  This is exactly the
## defect Tone's slacks-based measure was designed to fix, and the weighted
## measures below fix it a different way, by normalizing each slack before it
## is added.
##
##   "ram"         Range Adjusted Measure (Cooper, Park and Pastor 1999).
##                 w = 1/((p+q) R), R the range of that variable across the
##                 REFERENCE set. Units invariant, and rho = 1 - objective
##                 lands in [0, 1] with 1 exactly on the efficient set.
##   "mip"         Measure of Inefficiency Proportions (Cooper, Seiford and
##                 Tone). w = 1/((p+q) x_io), i.e. each slack as a fraction of
##                 the DMU's OWN level. Units invariant; needs strictly
##                 positive data for the same reason dea_sbm() does.
##   "unweighted"  the original of Charnes, Cooper, Golany, Seiford and Stutz
##                 (1985). Reported as a raw slack total: 0 is efficient and
##                 larger is worse, and it carries the units of the data.
##
## WHAT IT IS FOR.  The additive model has no orientation: it does not ask by
## what factor inputs could shrink or outputs grow, it asks how far the DMU is
## from the efficient frontier in every coordinate at once.  Its efficient set
## is the Pareto-Koopmans set exactly -- the same set dea_sbm() scores 1 -- so
## it is the natural companion when the question is "which units are
## dominated" rather than "by how much".
## ---------------------------------------------------------------------------

dea_add <- function(x, y, data = NULL,
                    measure = c("ram", "mip", "unweighted"),
                    rts = c("vrs", "crs", "nirs", "ndrs"),
                    peers = TRUE,
                    scaling = TRUE,
                    xref = NULL, yref = NULL, dataref = NULL) {

  call <- match.call()
  measure <- .match_arg_ci(measure, c("ram", "mip", "unweighted"), "measure")
  rts <- .match_arg_ci(rts, c("vrs", "crs", "nirs", "ndrs"), "rts")
  t0 <- proc.time()[["elapsed"]]

  ## MIP divides by each DMU's own levels; RAM and the unweighted model do not.
  d <- .dea_data(x, y, data, xref, yref, dataref, FALSE, "dea_add",
                 require_positive = identical(measure, "mip"))
  X <- d$X; Y <- d$Y; XR <- d$XR; YR <- d$YR
  dmu <- d$dmu
  n <- nrow(X); nr <- nrow(XR); p <- ncol(X); q <- ncol(Y)

  ## The unweighted model is the one case where scaling is NOT neutral: it is
  ## defined on the caller's units and rescaling the columns would silently
  ## change the answer, which is the whole point of the caveat above.
  use_scaling <- scaling && !identical(measure, "unweighted")
  ## The unweighted model is solved in the caller's units, so a wide spread of
  ## column magnitudes conditions the program badly -- and it is exactly the
  ## case where the measure is least meaningful anyway, since the objective is
  ## then dominated by whichever variable happens to be recorded in the largest
  ## units.
  if (!use_scaling && identical(measure, "unweighted")) {
    mag <- c(colMeans(abs(XR)), colMeans(abs(YR)))
    mag <- mag[mag > 0]
    if (length(mag) > 1L && max(mag) / min(mag) > 1e3) {
      warning("Column magnitudes span a factor of ",
              signif(max(mag) / min(mag), 3), ". The unweighted additive ",
              "objective adds raw slacks, so it is dominated by whichever ",
              "variable is recorded in the largest units, and the linear ",
              "program is poorly conditioned at this spread. Use ",
              "measure = \"ram\", which normalizes each slack by its range.",
              call. = FALSE)
    }
  }
  sc  <- .dea_scale(XR, YR, use_scaling)
  XRs <- sc$X; YRs <- sc$Y
  Xs  <- sweep(X, 2L, sc$sx, "/"); Ys <- sweep(Y, 2L, sc$sy, "/")

  res <- .dea_additive(Xs, Ys, XRs, YRs, rts, measure, peers, n, nr, p, q)

  eff <- res$eff
  eff[is.finite(eff) & abs(eff - res$efficient_at) < .DEA_CONSTANTS$TOL_EFF] <-
    res$efficient_at
  names(eff) <- dmu

  sx <- sweep(res$sx, 2L, sc$sx, "*")
  sy <- sweep(res$sy, 2L, sc$sy, "*")
  ## RELATIVE tolerances. An absolute 1e-8 is meaningless on a column measured
  ## in millions and brutal on one measured in thousandths, and the unweighted
  ## measure is solved in the caller's own units precisely so that it cannot be
  ## rescaled -- so the tolerance has to carry the scale instead. Without this,
  ## multiplying one input by 1000 flipped four DMUs in and out of the
  ## efficient set on slack totals of order 1e-9, against a smallest genuine
  ## slack of 0.105.
  tolx <- .DEA_CONSTANTS$TOL_SLACK * pmax(1, colMeans(abs(XR)))
  toly <- .DEA_CONSTANTS$TOL_SLACK * pmax(1, colMeans(abs(YR)))
  sx[abs(sx) < rep(tolx, each = nrow(sx))] <- 0
  sy[abs(sy) < rep(toly, each = nrow(sy))] <- 0
  dimnames(sx) <- list(dmu, colnames(X)); dimnames(sy) <- list(dmu, colnames(Y))
  L <- res$lambda
  if (!is.null(L)) { L[abs(L) < .DEA_CONSTANTS$TOL_LAMBDA] <- 0; dimnames(L) <- list(dmu, d$ref) }

  .dea_report_unsolved(res$status, n, d$self, "The additive program")

  structure(list(
    eff = eff, model = "additive", measure = measure,
    rts = rts, orientation = "none", super = FALSE,
    lambda = L, sum_lambda = if (is.null(L)) NULL else rowSums(L),
    slack_x = sx, slack_y = sy,
    ## Read off the SLACKS rather than the objective. The two agree
    ## mathematically -- the objective is zero exactly when no slack remains --
    ## but the slacks carry per-column tolerances and the objective is a sum
    ## that can hide a negligible term behind a large one.
    efficient = is.finite(eff) & rowSums(sx) == 0 & rowSums(sy) == 0,
    status = res$status,
    x = X, y = Y, xref = XR, yref = YR, self_ref = d$self,
    dmu = dmu, ref = d$ref, n = n, nref = nr, p = p, q = q,
    scaling = use_scaling, slack = TRUE,
    total_time = proc.time()[["elapsed"]] - t0,
    call = call
  ), class = "dea")
}

## The sweep. Only the objective changes between measures, and only MIP's
## objective changes between DMUs -- RAM's weights come from the reference
## set's ranges and the unweighted model's are all 1, so both set the objective
## once and then vary nothing but the right-hand side.
.dea_additive <- function(X, Y, XR, YR, rts, measure, peers, n, nr, p, q) {
  S  <- .lp_slack_build(XR, YR, rts)
  lp <- S$lp

  wx <- wy <- NULL
  if (identical(measure, "ram")) {
    ## Ranges over the REFERENCE set: the technology's spread is what makes a
    ## slack large or small, not the evaluated point's own level.
    rx <- apply(XR, 2L, function(v) diff(range(v)))
    ry <- apply(YR, 2L, function(v) diff(range(v)))
    ## A variable with no spread contributes no scope for improvement; giving
    ## it weight 0 is the limit of 1/R as R -> 0, and avoids a division by zero
    ## that would otherwise make every score -Inf.
    wx <- ifelse(rx > .DEA_CONSTANTS$MIN_POS, 1 / ((p + q) * rx), 0)
    wy <- ifelse(ry > .DEA_CONSTANTS$MIN_POS, 1 / ((p + q) * ry), 0)
    lpSolveAPI::set.objfn(lp, c(rep(0, nr), wx, wy))
  } else if (identical(measure, "unweighted")) {
    lpSolveAPI::set.objfn(lp, c(rep(0, nr), rep(1, p + q)))
  }

  eff <- numeric(n); st <- integer(n)
  sx <- matrix(NA_real_, n, p); sy <- matrix(NA_real_, n, q)
  L  <- if (peers) matrix(0, n, nr) else NULL

  for (o in seq_len(n)) {
    if (identical(measure, "mip")) {
      xo <- pmax(X[o, ], .DEA_CONSTANTS$MIN_POS)
      yo <- pmax(Y[o, ], .DEA_CONSTANTS$MIN_POS)
      lpSolveAPI::set.objfn(lp, c(rep(0, nr), 1/((p + q) * xo), 1/((p + q) * yo)))
    }
    lpSolveAPI::set.rhs(lp, c(X[o, ], Y[o, ]), seq_len(p + q))
    st[o] <- solve(lp)
    if (!st[o] %in% c(0L, 1L)) { eff[o] <- NA_real_; next }
    z <- lpSolveAPI::get.objective(lp)
    ## RAM and MIP are normalized so that the objective is an inefficiency
    ## share; the unweighted model has no such normalization and is reported
    ## as the raw total it is.
    eff[o] <- if (identical(measure, "unweighted")) z else 1 - z
    v <- lpSolveAPI::get.variables(lp)
    if (peers) L[o, ] <- v[seq_len(nr)]
    sx[o, ] <- v[nr + seq_len(p)]
    sy[o, ] <- v[nr + p + seq_len(q)]
  }
  list(eff = eff, lambda = L, sx = sx, sy = sy, status = st,
       efficient_at = if (identical(measure, "unweighted")) 0 else 1)
}
