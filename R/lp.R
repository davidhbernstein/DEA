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

## ---------------------------------------------------------------------------
## The MULTIPLIER (dual) programs.
##
## Every radial score above can be read off either of two programs.  The
## envelopment form asks "what combination of the other DMUs dominates this
## one"; the multiplier form asks "what set of prices would make this DMU look
## as good as possible".  Strong duality makes them attain the same value, so
## the score alone gives no reason to solve both -- but the multiplier form
## returns something the envelopment form does not: the weights themselves,
## v on the inputs and u on the outputs, which are what the DMU would have to
## believe about relative worth to justify its own score.
##
## Input orientation, at DMU o:
##
##   max  u'y_o - u0
##   s.t. v'x_o = 1
##        u'Y_j - v'X_j - u0 <= 0      for every reference DMU j
##        u >= 0, v >= 0
##
## Output orientation:
##
##   min  v'x_o - v0
##   s.t. u'y_o = 1
##        v'X_j - u'Y_j - v0 >= 0
##
## THE SIGN ON u0 IS THE WHOLE RETURNS-TO-SCALE ASSUMPTION, and it is not
## symmetric between the two orientations.  Under the convention above --
## u0 entering the objective and the constraints with the same MINUS sign --
##
##            input oriented      output oriented
##   crs      no u0               no u0
##   vrs      free                free
##   nirs     u0 >= 0             u0 <= 0
##   ndrs     u0 <= 0             u0 >= 0
##
## Getting one of these backwards does not error: it solves a different
## technology's program and returns a plausible number.  The table was fixed by
## checking the optimal value against the envelopment program on samples where
## the restriction actually binds (nirs and ndrs scores that differ from the
## vrs ones), because on data where it does not bind every sign convention
## agrees and the check proves nothing.  test-multipliers.R keeps that.
##
## As above, the matrix is built once: the constraint rows are the reference
## technology and never move.  Only row 1 (which holds the evaluated DMU) and
## the objective change from DMU to DMU.
## ---------------------------------------------------------------------------

## u0's bounds, given the technology and the orientation.  NULL means the model
## has no u0 at all.
.mult_u0_bounds <- function(rts, orientation) {
  if (identical(rts, "crs")) return(NULL)
  if (identical(rts, "vrs")) return(c(-.LP_INF, .LP_INF))
  nonneg <- c(0, .LP_INF)
  nonpos <- c(-.LP_INF, 0)
  switch(rts,
    nirs = if (orientation == "in") nonneg else nonpos,
    ndrs = if (orientation == "in") nonpos else nonneg,
    stop("rts = \"", rts, "\" has no multiplier form here.", call. = FALSE))
}

.lp_mult_build <- function(XR, YR, rts, orientation) {
  nr <- nrow(XR); p <- ncol(XR); q <- ncol(YR)
  b0 <- .mult_u0_bounds(rts, orientation)
  has0 <- !is.null(b0)
  nv <- p + q + has0
  lp <- lpSolveAPI::make.lp(1L + nr, nv)
  lpSolveAPI::lp.control(lp, sense = if (orientation == "in") "max" else "min",
                         epsel = .DEA_CONSTANTS$LP_EPSEL, verbose = "neutral")

  ## Row 1 is the normalization and is rewritten per DMU; rows 2..(nr+1) are
  ## the reference technology and are written once, here.
  if (orientation == "in") {
    for (i in seq_len(p)) lpSolveAPI::set.column(lp, i,      c(0, -XR[, i]))
    for (r in seq_len(q)) lpSolveAPI::set.column(lp, p + r,  c(0,  YR[, r]))
    if (has0) lpSolveAPI::set.column(lp, nv, c(0, rep(-1, nr)))
    lpSolveAPI::set.constr.type(lp, c("=", rep("<=", nr)))
  } else {
    for (i in seq_len(p)) lpSolveAPI::set.column(lp, i,      c(0,  XR[, i]))
    for (r in seq_len(q)) lpSolveAPI::set.column(lp, p + r,  c(0, -YR[, r]))
    if (has0) lpSolveAPI::set.column(lp, nv, c(0, rep(-1, nr)))
    lpSolveAPI::set.constr.type(lp, c("=", rep(">=", nr)))
  }
  lpSolveAPI::set.rhs(lp, c(1, rep(0, nr)))
  if (has0) {
    lpSolveAPI::set.bounds(lp, lower = b0[1L], upper = b0[2L], columns = nv)
  }
  list(lp = lp, nref = nr, p = p, q = q, nv = nv, has0 = has0,
       rts = rts, orientation = orientation)
}

.lp_mult_at <- function(M, X, Y, o, exclude = NULL) {
  lp <- M$lp; p <- M$p; q <- M$q; nv <- M$nv
  ## The full row, not a sparse update: set.row with an index vector zeroes the
  ## entries it is not given, and writing the whole thing is both cheaper to
  ## reason about and no slower at these widths.
  row1 <- numeric(nv)
  if (M$orientation == "in") {
    row1[seq_len(p)] <- X[o, ]
    obj <- c(rep(0, p), Y[o, ], if (M$has0) -1 else NULL)
  } else {
    row1[p + seq_len(q)] <- Y[o, ]
    obj <- c(X[o, ], rep(0, q), if (M$has0) -1 else NULL)
  }
  lpSolveAPI::set.row(lp, 1L, row1)
  lpSolveAPI::set.objfn(lp, obj)

  ## Super-efficiency drops DMU o from its own reference set. In the
  ## envelopment form that is a bound on a column; here it is the CONSTRAINT
  ## contributed by o, relaxed by pushing its right-hand side out of reach.
  if (!is.null(exclude)) {
    lpSolveAPI::set.rhs(lp, if (M$orientation == "in") .LP_INF else -.LP_INF,
                        exclude + 1L)
  }
  st <- solve(lp)
  if (!is.null(exclude)) lpSolveAPI::set.rhs(lp, 0, exclude + 1L)

  if (!st %in% c(0L, 1L)) {
    return(list(status = st, eff = NA_real_, v = rep(NA_real_, p),
                u = rep(NA_real_, q), u0 = NA_real_))
  }
  z <- lpSolveAPI::get.variables(lp)
  list(status = st,
       eff = lpSolveAPI::get.objective(lp),
       v   = z[seq_len(p)],
       u   = z[p + seq_len(q)],
       ## u0 is reported in the sign the TABLE above uses, which is the sign it
       ## carries in the objective -- so that for vrs its sign is the usual
       ## returns-to-scale reading rather than its negative.
       u0  = if (M$has0) z[nv] else NA_real_)
}
