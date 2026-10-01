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
## ---------------------------------------------------------------------------
## Solve an already-aimed program, and rebuild it once if it gives up.
##
## THE FAILURE THIS GUARDS. Every sweep in this package builds ONE linear
## program per (technology, orientation) and rewrites only the evaluated DMU's
## column and right-hand side -- that is where the speed comes from. But
## lpSolveAPI carries basis and factorisation state on that object, and for
## some right-hand sides the inherited state is bad enough that the solve gives
## up with status 5. It is not a tolerance problem: of 56 failures on a
## 1200-DMU variable-returns radial fit, loosening epsel from 1e-12 to 1e-9
## fixed one, every scaling mode fixed at most a quarter and guess.basis()
## fixed one, while REBUILDING THE OBJECT fixed all 56. The state is the cause
## and a fresh object is the remedy.
##
## It is a property of the REUSE rather than of any one program, which is why
## the guard belongs here rather than being re-derived at each sweep. Observed
## so far in the radial stage-two slack program, the multiplier program, the
## additive sweep (one DMU of 1200) and the oriented slacks-based measure
## (seeds 2 and 3 at n = 900 and 1600).
##
## WHAT IS NOT RETRIED. `final` lists the statuses that are answers rather than
## failures. It defaults to infeasible-included, because for an envelopment
## program infeasibility is a statement about the data -- under
## super-efficiency it is the expected one -- and re-solving cannot change it.
## The multiplier program passes a different set: there an infeasible dual
## contradicts strong duality whenever the primal is feasible and bounded, so
## it IS retried, while an unbounded dual is the answer and is not.
##
## `fresh` must return a newly built program already aimed at the same DMU.
## Returning the lp that produced the status lets the caller read its results
## without knowing which of the two it got.
.lp_solve_retry <- function(lp, fresh, final = c(0L, 1L, 2L)) {
  st <- solve(lp)
  if (st %in% final) return(list(status = st, lp = lp))
  lp2 <- fresh()
  st <- solve(lp2)
  if (st %in% final) return(list(status = st, lp = lp2))
  ## STAGE THREE: a fresh object with a different scaling mode. See the second
  ## half of the note above -- rebuilding cannot help a program that fails
  ## identically however often it is rebuilt, and for that class a scaling
  ## change does. Applied to the object `fresh()` returns rather than asking
  ## every builder for a scaling argument: lpSolveAPI applies scaling at solve
  ## time, so setting it after the program is aimed gives the same answer as
  ## setting it before (checked both ways on a failing program).
  for (mode in .DEA_CONSTANTS$LP_SCALING_FALLBACK) {
    lp3 <- fresh()
    invisible(lpSolveAPI::lp.control(lp3, scaling = mode))
    st <- solve(lp3)
    if (st %in% final) return(list(status = st, lp = lp3))
  }
  list(status = st, lp = lp3)
}

## Same three stages for the second-stage slack program, which likewise aims and
## solves in one call. Stage three here is UNEXERCISED: every stage-two failure
## observed so far -- 56 of them on a 1200-DMU variable-returns fit -- was fixed
## by the rebuild alone. It is present because the failure class that needs it
## was found in the stage-one program of the same sweep, and a guard that exists
## only where the bug has already been seen is a guard against the past.
.lp_slack_solve <- function(S, XRs, YRs, rts, rhs_x, rhs_y) {
  ## The rebuilds must carry the same objective mask, or a retry would answer a
  ## different question from the first attempt.
  nd_x <- S$nd_x; nd_y <- S$nd_y
  z <- .lp_slack_at(S, rhs_x, rhs_y)          # the hot path, as above
  if (z$status %in% c(0L, 1L)) return(z)
  ok <- function(z) z$status %in% c(0L, 1L)
  fresh <- function(mode = NULL) {
    S2 <- .lp_slack_build(XRs, YRs, rts, nd_x, nd_y)
    if (!is.null(mode)) invisible(lpSolveAPI::lp.control(S2$lp, scaling = mode))
    S2
  }
  z <- .lp_slack_at(fresh(), rhs_x, rhs_y)
  if (ok(z)) return(z)
  for (mode in .DEA_CONSTANTS$LP_SCALING_FALLBACK) {
    z <- .lp_slack_at(fresh(mode), rhs_x, rhs_y)
    if (ok(z)) return(z)
  }
  z
}

## The radial sweep's own version of the three stages. It needs its own because
## .lp_radial_at() aims and solves in one call, so there is no aimed-but-unsolved
## program for .lp_solve_retry() to take; the RULE is the same and is documented
## once, above.
## `B` is the ONE object built outside the caller's loop and reused across every
## DMU -- that reuse is where this package's speed comes from, so the first
## attempt must use it and only a failure may rebuild.
.lp_radial_solve <- function(B, XRs, YRs, rts, orientation, Xs, Ys, o, exclude,
                             fixed = NULL) {
  ## The first attempt is the hot path -- it runs n times per fit and succeeds
  ## essentially always -- so it is written without the helper closures below
  ## rather than constructing them once per DMU for nothing. That is tidiness,
  ## NOT a measured win: wrapping .lp_radial_at() in this function costs 0.983x
  ## on an 800-DMU sweep and 0.987x end to end, both inside noise. An earlier
  ## note here claimed 5.6%, which was an artifact of always timing the old
  ## package first -- whichever of the two ran second measured slower, in both
  ## orders.
  r <- .lp_radial_at(B, Xs, Ys, o, exclude = exclude, fixed = fixed)
  if (r$status %in% c(0L, 1L, 2L)) return(r)
  ok <- function(r) r$status %in% c(0L, 1L, 2L)
  fresh <- function(mode = NULL) {
    B2 <- .lp_radial_build(XRs, YRs, rts, orientation)
    if (!is.null(mode)) invisible(lpSolveAPI::lp.control(B2$lp, scaling = mode))
    B2
  }
  r <- .lp_radial_at(fresh(), Xs, Ys, o, exclude = exclude, fixed = fixed)
  if (ok(r)) return(r)
  for (mode in .DEA_CONSTANTS$LP_SCALING_FALLBACK) {
    r <- .lp_radial_at(fresh(mode), Xs, Ys, o, exclude = exclude, fixed = fixed)
    if (ok(r)) return(r)
  }
  r
}

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
## `fixed` marks NON-DISCRETIONARY rows on the oriented side: a logical of
## length p input-oriented, of length q output-oriented. Those rows keep
## constraining the comparison set but theta does not multiply them, so the row
## reads sum(lambda_j x_ij) <= x_io rather than <= theta x_io. In the program
## that is two changes together -- a zero in the theta column and the DMU's own
## level moved into the right-hand side, which the build left at 0 -- and doing
## only the first would quietly impose x_ij <= 0. NULL is the ordinary model and
## touches nothing.
.lp_radial_at <- function(B, X, Y, o, exclude = NULL, fixed = NULL) {
  lp <- B$lp; p <- B$p; q <- B$q
  ## The 0 index is the objective row. lpSolveAPI's set.column REPLACES the
  ## whole column, objective coefficient included, so omitting it silently
  ## zeroes the objective and every DMU then solves min 0 -- which returns
  ## status 0 (optimal) and a feasible-but-arbitrary theta. Nothing errors;
  ## the scores are simply wrong.
  if (B$orientation == "in") {
    co <- -X[o, ]
    if (!is.null(fixed)) co[fixed] <- 0
    lpSolveAPI::set.column(lp, 1L, c(1, co), c(0L, seq_len(p)))
    lpSolveAPI::set.rhs(lp, Y[o, ], p + seq_len(q))
    ## Written every time rather than only for the fixed rows: the object is
    ## reused across DMUs, so a level left behind by the previous one would be
    ## applied to this one.
    if (!is.null(fixed)) {
      rx <- numeric(p); rx[fixed] <- X[o, fixed]
      lpSolveAPI::set.rhs(lp, rx, seq_len(p))
    }
  } else {
    co <- -Y[o, ]
    if (!is.null(fixed)) co[fixed] <- 0
    lpSolveAPI::set.column(lp, 1L, c(1, co), c(0L, p + seq_len(q)))
    lpSolveAPI::set.rhs(lp, X[o, ], seq_len(p))
    if (!is.null(fixed)) {
      ry <- numeric(q); ry[fixed] <- Y[o, fixed]
      lpSolveAPI::set.rhs(lp, ry, p + seq_len(q))
    }
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
## `nd_x`/`nd_y` drop a slack from the OBJECTIVE while leaving it in the
## constraints, which is what makes a variable non-discretionary at the second
## stage: the DMU gets no credit for a reduction it cannot make, but the row
## still has to balance. Both sides are maskable in both orientations -- the
## orientation decides which side stage ONE scales, not which slacks count.
.lp_slack_build <- function(X, Y, rts, nd_x = NULL, nd_y = NULL) {
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
  if (!is.null(nd_x)) obj[n + which(nd_x)] <- 0
  if (!is.null(nd_y)) obj[n + p + which(nd_y)] <- 0
  lpSolveAPI::set.objfn(lp, obj)
  list(lp = lp, n = n, p = p, q = q, nd_x = nd_x, nd_y = nd_y)
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

## `restrict`, when given, is the weight-restriction block: a list with `mat`
## (k x (p+q), inputs then outputs, already in the SOLVER's units), `type` and
## `rhs`. Its rows go after the technology rows, so nothing that indexes a
## technology row by `j + 1L` has to change -- including the super-efficiency
## relaxation in .lp_mult_at(). u0 never appears in a restriction: a bound on
## the returns-to-scale intercept is a statement about the technology, not
## about relative worth, and the rts argument is where that belongs.
.lp_mult_build <- function(XR, YR, rts, orientation, restrict = NULL) {
  nr <- nrow(XR); p <- ncol(XR); q <- ncol(YR)
  b0 <- .mult_u0_bounds(rts, orientation)
  has0 <- !is.null(b0)
  nv <- p + q + has0
  k <- if (is.null(restrict)) 0L else nrow(restrict$mat)
  R <- if (k) restrict$mat else matrix(0, 0L, p + q)
  lp <- lpSolveAPI::make.lp(1L + nr + k, nv)
  lpSolveAPI::lp.control(lp, sense = if (orientation == "in") "max" else "min",
                         epsel = .DEA_CONSTANTS$LP_EPSEL, verbose = "neutral")

  ## Row 1 is the normalization and is rewritten per DMU; rows 2..(nr+1) are
  ## the reference technology and are written once, here; the rest, if any, are
  ## the weight restrictions and never move either.
  if (orientation == "in") {
    for (i in seq_len(p)) lpSolveAPI::set.column(lp, i,      c(0, -XR[, i], R[, i]))
    for (r in seq_len(q)) lpSolveAPI::set.column(lp, p + r,  c(0,  YR[, r], R[, p + r]))
    if (has0) lpSolveAPI::set.column(lp, nv, c(0, rep(-1, nr), rep(0, k)))
    lpSolveAPI::set.constr.type(lp, c("=", rep("<=", nr)), seq_len(1L + nr))
  } else {
    for (i in seq_len(p)) lpSolveAPI::set.column(lp, i,      c(0,  XR[, i], R[, i]))
    for (r in seq_len(q)) lpSolveAPI::set.column(lp, p + r,  c(0, -YR[, r], R[, p + r]))
    if (has0) lpSolveAPI::set.column(lp, nv, c(0, rep(-1, nr), rep(0, k)))
    lpSolveAPI::set.constr.type(lp, c("=", rep(">=", nr)), seq_len(1L + nr))
  }
  if (k) {
    lpSolveAPI::set.constr.type(lp, restrict$type, 1L + nr + seq_len(k))
    lpSolveAPI::set.rhs(lp, restrict$rhs, 1L + nr + seq_len(k))
  }
  lpSolveAPI::set.rhs(lp, c(1, rep(0, nr)), seq_len(1L + nr))
  if (has0) {
    lpSolveAPI::set.bounds(lp, lower = b0[1L], upper = b0[2L], columns = nv)
  }
  list(lp = lp, nref = nr, p = p, q = q, nv = nv, has0 = has0,
       nrestrict = k, rts = rts, orientation = orientation)
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
