## ---------------------------------------------------------------------------
## The linear-programming layer.
##
## Every model in this package is n linear programs -- one per DMU -- that
## differ from each other only in a handful of numbers.  The whole design of
## this file follows from that: the constraint matrix is built ONCE per
## (technology, orientation) and then mutated in place, because the columns
## holding the reference technology (the lambda columns) are identical for
## every DMU.  Only the evaluated DMU's own column and right-hand side change.
##
## That is where the speed comes from.  Rebuilding an n-column LP n times is
## O(n^2) work before the solver has done anything at all, and at n = 2000 that
## dominates the solve time.
##
## Row layout, shared by every builder here:
##   rows 1..p            one per input
##   rows p+1..p+q        one per output
##   row  p+q+1           the returns-to-scale row, present unless rts = "crs"
## ---------------------------------------------------------------------------

## lpSolveAPI's stand-in for infinity.  Passing Inf works but round-trips
## through a double comparison inside the C code; 1e30 is the documented value.
.LP_INF <- 1e30

.rts_row <- function(rts) {
  switch(rts,
    crs  = NULL,
    vrs  = "=",
    nirs = "<=",
    ndrs = ">=",
    fdh  = "=",
    stop("Unknown rts '", rts, "'.", call. = FALSE))
}

## ---------------------------------------------------------------------------
## Radial (Debreu-Farrell) programs.
##
## input orientation      min theta  s.t.  X'lambda <= theta x_o
##                                         Y'lambda >= y_o
## output orientation     max phi    s.t.  X'lambda <= x_o
##                                         Y'lambda >= phi y_o
##
## Variable 1 is theta (or phi); variables 2..(n+1) are lambda.
## ---------------------------------------------------------------------------
.lp_radial_build <- function(X, Y, rts, orientation) {
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  rr <- .rts_row(rts)
  nrows <- p + q + (!is.null(rr))
  lp <- lpSolveAPI::make.lp(nrows, n + 1L)
  lpSolveAPI::lp.control(lp, sense = if (orientation == "in") "min" else "max",
                         epsel = .DEA_CONSTANTS$LP_EPSEL, verbose = "neutral")

  ## The lambda columns: input loadings, output loadings, and a 1 in the
  ## returns-to-scale row.  Fixed for the whole sweep.
  rts_one <- if (is.null(rr)) numeric(0) else 1
  for (j in seq_len(n)) {
    lpSolveAPI::set.column(lp, j + 1L, c(X[j, ], Y[j, ], rts_one))
  }

  ## Constraint senses and the parts of the right-hand side that never move.
  lpSolveAPI::set.constr.type(lp, rep("<=", p), seq_len(p))
  lpSolveAPI::set.constr.type(lp, rep(">=", q), p + seq_len(q))
  if (!is.null(rr)) {
    lpSolveAPI::set.constr.type(lp, rr, p + q + 1L)
    lpSolveAPI::set.rhs(lp, 1, p + q + 1L)
  }

  obj <- numeric(n + 1L); obj[1L] <- 1
  lpSolveAPI::set.objfn(lp, obj)

  if (orientation == "in") {
    ## theta enters only the input rows, whose rhs is then 0.
    lpSolveAPI::set.rhs(lp, rep(0, p), seq_len(p))
  } else {
    ## phi enters only the output rows, whose rhs is then 0.
    lpSolveAPI::set.rhs(lp, rep(0, q), p + seq_len(q))
  }

  list(lp = lp, n = n, p = p, q = q, rts = rts, orientation = orientation,
       has_rts_row = !is.null(rr))
}

## Point the built program at DMU `o` and solve it.
##
## `exclude` drops a DMU from its own reference set, which is what makes the
## Andersen-Petersen super-efficiency score.  Under vrs/nirs/ndrs that program
## can be genuinely INFEASIBLE -- a DMU at the boundary of the input space may
## have no other DMU able to dominate it -- and the honest answer there is NA,
## not a large number.  See ?dea, section 'super-efficiency'.
.lp_radial_at <- function(B, X, Y, o, exclude = NULL) {
  lp <- B$lp; p <- B$p; q <- B$q
  ## The 0 index is the objective row. lpSolveAPI's set.column REPLACES the
  ## whole column, objective coefficient included, so omitting it silently
  ## zeroes the objective and every DMU then solves min 0 -- which returns
  ## status 0 (optimal) and a feasible-but-arbitrary theta. Nothing errors;
  ## the scores are simply wrong.
  if (B$orientation == "in") {
    lpSolveAPI::set.column(lp, 1L, c(1, -X[o, ]), c(0L, seq_len(p)))
    lpSolveAPI::set.rhs(lp, Y[o, ], p + seq_len(q))
  } else {
    lpSolveAPI::set.column(lp, 1L, c(1, -Y[o, ]), c(0L, p + seq_len(q)))
    lpSolveAPI::set.rhs(lp, X[o, ], seq_len(p))
  }
  if (!is.null(exclude)) lpSolveAPI::set.bounds(lp, upper = 0, columns = exclude + 1L)
  st <- solve(lp)
  if (!is.null(exclude)) lpSolveAPI::set.bounds(lp, upper = .LP_INF, columns = exclude + 1L)

  ## 0 optimal, 1 suboptimal but usable, 2 infeasible, 3 unbounded.
  if (!st %in% c(0L, 1L)) {
    return(list(status = st, eff = NA_real_, lambda = rep(NA_real_, B$n)))
  }
  v <- lpSolveAPI::get.variables(lp)
  list(status = st, eff = v[1L], lambda = v[-1L])
}

## ---------------------------------------------------------------------------
## Stage two: the slack-maximizing program.
##
## The radial score alone can leave a projected point on a vertical or
## horizontal face of the frontier -- technically on the boundary, but still
## dominated.  Stage two finds the largest remaining slack, which is what makes
## the efficiency judgement Pareto-Koopmans rather than merely radial:
##
##   max  sum(s_minus) + sum(s_plus)
##   s.t. X'lambda + s_minus = theta* x_o
##        Y'lambda - s_plus  = y_o
##
## Only the right-hand side depends on the DMU, so the entire matrix is fixed
## across the sweep -- this stage is close to free relative to stage one.
##
## Variables: lambda (n), s_minus (p), s_plus (q).
## ---------------------------------------------------------------------------
.lp_slack_build <- function(X, Y, rts) {
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  rr <- .rts_row(rts)
  nrows <- p + q + (!is.null(rr))
  lp <- lpSolveAPI::make.lp(nrows, n + p + q)
  lpSolveAPI::lp.control(lp, sense = "max", epsel = .DEA_CONSTANTS$LP_EPSEL,
                         verbose = "neutral")
  rts_one <- if (is.null(rr)) numeric(0) else 1
  for (j in seq_len(n)) {
    lpSolveAPI::set.column(lp, j, c(X[j, ], Y[j, ], rts_one))
  }
  for (i in seq_len(p)) lpSolveAPI::set.column(lp, n + i, 1, i)
  for (r in seq_len(q)) lpSolveAPI::set.column(lp, n + p + r, -1, p + r)

  lpSolveAPI::set.constr.type(lp, rep("=", p + q), seq_len(p + q))
  if (!is.null(rr)) {
    lpSolveAPI::set.constr.type(lp, rr, p + q + 1L)
    lpSolveAPI::set.rhs(lp, 1, p + q + 1L)
  }
  obj <- c(rep(0, n), rep(1, p + q))
  lpSolveAPI::set.objfn(lp, obj)
  list(lp = lp, n = n, p = p, q = q)
}

.lp_slack_at <- function(S, xo, yo) {
  lp <- S$lp; p <- S$p; q <- S$q
  lpSolveAPI::set.rhs(lp, c(xo, yo), seq_len(p + q))
  st <- solve(lp)
  if (!st %in% c(0L, 1L)) {
    return(list(status = st, lambda = rep(NA_real_, S$n),
                sx = rep(NA_real_, p), sy = rep(NA_real_, q)))
  }
  v <- lpSolveAPI::get.variables(lp)
  list(status = st, lambda = v[seq_len(S$n)],
       sx = v[S$n + seq_len(p)], sy = v[S$n + p + seq_len(q)])
}

## ---------------------------------------------------------------------------
## Free disposal hull, by enumeration rather than by mixed-integer programming.
##
## FDH is the same program as VRS with lambda restricted to be binary, and it
## is routinely implemented that way -- but the binary program has a closed
## form, because exactly one lambda is 1 at the optimum.  Input-oriented,
##
##   theta_o = min over j dominating o of  max_i (x_ij / x_io)
##
## where j dominates o if y_j >= y_o componentwise.  That is an O(n^2 (p+q))
## sweep of arithmetic against n branch-and-bound solves, and it is exact
## rather than exact-up-to-a-gap-tolerance.
## ---------------------------------------------------------------------------
.fdh_radial <- function(X, Y, XR, YR, orientation) {
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  eff  <- numeric(n)
  peer <- integer(n)
  for (o in seq_len(n)) {
    if (orientation == "in") {
      ok <- .rowsAllGE(YR, Y[o, ])
      ratio <- do.call(pmax, lapply(seq_len(p), function(i)
        XR[, i] / max(X[o, i], .DEA_CONSTANTS$MIN_POS)))
      ratio[!ok] <- Inf
      if (!any(ok)) { eff[o] <- NA_real_; peer[o] <- NA_integer_; next }
      peer[o] <- which.min(ratio); eff[o] <- ratio[peer[o]]
    } else {
      ok <- .rowsAllLE(XR, X[o, ])
      ratio <- do.call(pmin, lapply(seq_len(q), function(r)
        YR[, r] / max(Y[o, r], .DEA_CONSTANTS$MIN_POS)))
      ratio[!ok] <- -Inf
      if (!any(ok)) { eff[o] <- NA_real_; peer[o] <- NA_integer_; next }
      peer[o] <- which.max(ratio); eff[o] <- ratio[peer[o]]
    }
  }
  list(eff = eff, peer = peer)
}

.rowsAllGE <- function(M, v) {
  ok <- rep(TRUE, nrow(M))
  for (k in seq_len(ncol(M))) ok <- ok & (M[, k] >= v[k] - 1e-12)
  ok
}
.rowsAllLE <- function(M, v) {
  ok <- rep(TRUE, nrow(M))
  for (k in seq_len(ncol(M))) ok <- ok & (M[, k] <= v[k] + 1e-12)
  ok
}
