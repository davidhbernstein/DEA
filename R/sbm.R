## ---------------------------------------------------------------------------
## dea_sbm() -- the slacks-based measure of Tone (2001).
##
## The radial score answers "by what common factor can every input be cut?".
## That factor stops at the first binding input, so a DMU sitting on a vertical
## face of the frontier scores 1 while still wasting an input outright.  The
## SBM answers "what fraction of each input is being wasted?", averaged over
## inputs, and is therefore equal to 1 only for a Pareto-Koopmans efficient
## DMU.  It is also units invariant and monotone in every slack, which the
## radial score plus a slack report is not.
##
##   rho = [1 - (1/p) sum_i s_i^- / x_io] / [1 + (1/q) sum_r s_r^+ / y_ro]
##
## THE FRACTIONAL PROGRAM.  The non-oriented measure is a ratio of two linear
## forms, which is not an LP.  Tone linearizes it by the Charnes-Cooper change
## of variable: introduce t > 0, set Lambda = t*lambda, S = t*s, and force the
## denominator to 1 with a normalization row.  The resulting LP is exact -- not
## an approximation -- and the original quantities come back as lambda =
## Lambda/t, s = S/t.  t > 0 is not imposed: with t = 0 every Lambda and S is
## driven to 0 by the balance rows and the normalization row then reads 0 = 1,
## so the LP itself rules it out.
##
## The two ORIENTED measures need no such trick.  Dropping one side of the
## ratio leaves a linear objective, so input- and output-oriented SBM are
## ordinary LPs over the same constraint set as the radial model's stage two.
##
## Variable layout, non-oriented:  t (1), Lambda (n), S^- (p), S^+ (q)
## Variable layout, oriented:      lambda (n), s^- (p), s^+ (q)
##
## ---------------------------------------------------------------------------
## THE ADDITIVE PROGRAM NEEDS THE POINT TO BE INSIDE THE TECHNOLOGY.
##
## This never bites in ordinary use, and always bites with an external
## reference set, so it is worth stating precisely.  The balance rows are
## EQUALITIES with non-negative slacks,
##     X'lambda + s^- = x_o    and    Y'lambda - s^+ = y_o,
## which together say that some convex combination of reference DMUs weakly
## dominates the evaluated point in every input AND every output.  For a DMU
## scored against its own sample that is free -- it dominates itself.  For a
## point scored against someone else's frontier it can simply be false, and
## then the program is INFEASIBLE and the honest answer is NA.
##
## The radial model has no such problem: its output program maximizes phi
## subject to Y'lambda >= phi y_o, and phi is free to fall below 1, so a point
## outside the technology gets a meaningful score under 1 rather than nothing.
##
## Measured on the convergence harness's design -- 200 interior evaluation
## points against a 200-DMU reference sample, 2 inputs and 2 outputs -- 17 of
## the 200 came back infeasible, every one of them with ZERO dominating
## reference DMUs and a radial phi below 1.  The infeasibility was correct in
## every case.  dea_sbm() warns and reports NA rather than substituting a
## number, and the warning names this cause.
## ---------------------------------------------------------------------------

dea_sbm <- function(x, y, data = NULL,
                    rts = c("vrs", "crs", "nirs", "ndrs"),
                    orientation = c("none", "in", "out"),
                    peers = TRUE,
                    scaling = TRUE,
                    xref = NULL, yref = NULL, dataref = NULL) {

  call <- match.call()
  rts <- .match_arg_ci(rts, c("vrs", "crs", "nirs", "ndrs"), "rts")
  orientation <- .match_arg_ci(orientation, c("none", "in", "out"), "orientation")
  t0 <- proc.time()[["elapsed"]]

  d <- .dea_data(x, y, data, xref, yref, dataref, FALSE, "dea_sbm",
                 require_positive = TRUE)
  X <- d$X; Y <- d$Y; XR <- d$XR; YR <- d$YR
  dmu <- d$dmu
  n <- nrow(X); nr <- nrow(XR); p <- ncol(X); q <- ncol(Y)

  sc  <- .dea_scale(XR, YR, scaling)
  XRs <- sc$X; YRs <- sc$Y
  Xs  <- sweep(X, 2L, sc$sx, "/"); Ys <- sweep(Y, 2L, sc$sy, "/")

  res <- if (orientation == "none")
           .sbm_nonoriented(Xs, Ys, XRs, YRs, rts, peers, n, nr, p, q)
         else
           .sbm_oriented(Xs, Ys, XRs, YRs, rts, orientation, peers, n, nr, p, q)

  eff <- res$eff
  eff[is.finite(eff) & abs(eff - 1) < .DEA_CONSTANTS$TOL_EFF] <- 1
  names(eff) <- dmu

  ## Infeasibility here has one overwhelmingly likely cause, and saying so is
  ## more use than reporting a solver code.
  n_inf <- sum(res$status == 2L)
  if (n_inf > 0L) {
    warning(n_inf, " of ", n, " DMU(s) gave an INFEASIBLE program and are ",
            "reported as NA. The slacks-based measure needs some convex ",
            "combination of the reference DMUs to weakly dominate the ",
            "evaluated point in every input and every output; a point lying ",
            "outside the estimated technology has no such combination. ",
            if (!d$self)
              "This is expected when scoring against an external `xref`/`yref`: use dea() instead, whose radial score falls below 1 for such a point rather than failing."
            else
              "Seeing this on a self-referenced fit is unusual and may indicate near-duplicate DMUs or severe scaling; check `status`.",
            call. = FALSE)
  }

  sx <- sweep(res$sx, 2L, sc$sx, "*")
  sy <- sweep(res$sy, 2L, sc$sy, "*")
  sx[abs(sx) < .DEA_CONSTANTS$TOL_SLACK] <- 0
  sy[abs(sy) < .DEA_CONSTANTS$TOL_SLACK] <- 0
  dimnames(sx) <- list(dmu, colnames(X)); dimnames(sy) <- list(dmu, colnames(Y))
  L <- res$lambda
  if (!is.null(L)) { L[abs(L) < .DEA_CONSTANTS$TOL_LAMBDA] <- 0; dimnames(L) <- list(dmu, d$ref) }

  structure(list(
    eff = eff, model = "sbm", rts = rts, orientation = orientation, super = FALSE,
    lambda = L, sum_lambda = if (is.null(L)) NULL else rowSums(L),
    slack_x = sx, slack_y = sy,
    ## For the SBM the two notions coincide by construction: rho = 1 if and
    ## only if every slack is zero.
    efficient = eff == 1,
    status = res$status,
    x = X, y = Y, xref = XR, yref = YR, self_ref = d$self,
    dmu = dmu, ref = d$ref, n = n, nref = nr, p = p, q = q,
    scaling = scaling, slack = TRUE,
    total_time = proc.time()[["elapsed"]] - t0,
    call = call
  ), class = "dea")
}

## --- non-oriented: the Charnes-Cooper linearization -------------------------
.sbm_nonoriented <- function(X, Y, XR, YR, rts, peers, n, nr, p, q) {
  rr <- .rts_row(rts)
  ## row 1        normalization      t + (1/q) sum S^+_r / y_ro = 1
  ## rows 2..p+1  input balance      sum_j Lambda_j x_ij + S^-_i - t x_io = 0
  ## rows ..+q    output balance     sum_j Lambda_j y_rj - S^+_r - t y_ro = 0
  ## last row     returns to scale   sum_j Lambda_j - t {=,<=,>=} 0
  nrows <- 1L + p + q + (!is.null(rr))
  ncols <- 1L + nr + p + q
  ri <- 1L + seq_len(p); ro <- 1L + p + seq_len(q); rr_i <- nrows

  lp <- lpSolveAPI::make.lp(nrows, ncols)
  lpSolveAPI::lp.control(lp, sense = "min", epsel = .DEA_CONSTANTS$LP_EPSEL,
                         verbose = "neutral")
  rts_one <- if (is.null(rr)) numeric(0) else 1
  for (j in seq_len(nr)) {
    lpSolveAPI::set.column(lp, 1L + j, c(XR[j, ], YR[j, ], rts_one),
                           c(ri, ro, if (is.null(rr)) NULL else rr_i))
  }
  for (i in seq_len(p)) lpSolveAPI::set.column(lp, 1L + nr + i, 1, ri[i])
  for (r in seq_len(q)) lpSolveAPI::set.column(lp, 1L + nr + p + r, -1, ro[r])

  lpSolveAPI::set.constr.type(lp, rep("=", 1L + p + q), seq_len(1L + p + q))
  lpSolveAPI::set.rhs(lp, c(1, rep(0, p + q)), seq_len(1L + p + q))
  if (!is.null(rr)) {
    lpSolveAPI::set.constr.type(lp, rr, rr_i)
    lpSolveAPI::set.rhs(lp, 0, rr_i)
  }

  eff <- numeric(n); st <- integer(n)
  sx <- matrix(NA_real_, n, p); sy <- matrix(NA_real_, n, q)
  L  <- if (peers) matrix(0, n, nr) else NULL

  for (o in seq_len(n)) {
    xo <- pmax(X[o, ], .DEA_CONSTANTS$MIN_POS)
    yo <- pmax(Y[o, ], .DEA_CONSTANTS$MIN_POS)
    ## The t column carries the DMU's own levels, and the objective row.
    lpSolveAPI::set.column(lp, 1L,
      c(1, 1, -xo, -yo, if (is.null(rr)) NULL else -1),
      c(0L, 1L, ri, ro, if (is.null(rr)) NULL else rr_i))
    ## S^- appears only in the objective (-1/(p x_io)) and its own balance row.
    for (i in seq_len(p))
      lpSolveAPI::set.column(lp, 1L + nr + i, c(-1/(p * xo[i]), 1), c(0L, ri[i]))
    ## S^+ appears in the normalization row (1/(q y_ro)) and its balance row.
    for (r in seq_len(q))
      lpSolveAPI::set.column(lp, 1L + nr + p + r, c(1/(q * yo[r]), -1), c(1L, ro[r]))

    st[o] <- solve(lp)
    if (!st[o] %in% c(0L, 1L)) { eff[o] <- NA_real_; next }
    eff[o] <- lpSolveAPI::get.objective(lp)
    v  <- lpSolveAPI::get.variables(lp)
    tt <- max(v[1L], .DEA_CONSTANTS$MIN_POS)
    if (peers) L[o, ] <- v[1L + seq_len(nr)] / tt
    sx[o, ] <- v[1L + nr + seq_len(p)] / tt
    sy[o, ] <- v[1L + nr + p + seq_len(q)] / tt
  }
  list(eff = eff, lambda = L, sx = sx, sy = sy, status = st)
}

## --- oriented: already linear ----------------------------------------------
## Input:  rho_I = 1 - z,        z = max (1/p) sum s^-_i / x_io
## Output: rho_O = 1 / (1 + z),  z = max (1/q) sum s^+_r / y_ro
## Both over  X'lambda + s^- = x_o,  Y'lambda - s^+ = y_o,  returns to scale.
##
## BOTH PROGRAMS MAXIMIZE.  rho is decreasing in the slack it scores, so the
## measure is minimized by making that slack as large as the technology allows.
## Minimizing the output slack instead -- the natural-looking mirror of the
## input program's "min rho_I" -- is degenerate: lambda_o = 1 sets every slack
## to zero, so the program returns rho_O = 1 for every DMU in the sample and
## nothing errors. The un-scored slack on the other side stays in the model as
## a free non-negative variable, which is what makes the oriented measure
## one-sided: rho_I sees only input waste and rho_O only forgone output.
.sbm_oriented <- function(X, Y, XR, YR, rts, orientation, peers, n, nr, p, q) {
  S <- .lp_slack_build(XR, YR, rts)        ## same constraint set as stage two
  lp <- S$lp

  eff <- numeric(n); st <- integer(n)
  sx <- matrix(NA_real_, n, p); sy <- matrix(NA_real_, n, q)
  L  <- if (peers) matrix(0, n, nr) else NULL

  for (o in seq_len(n)) {
    xo <- pmax(X[o, ], .DEA_CONSTANTS$MIN_POS)
    yo <- pmax(Y[o, ], .DEA_CONSTANTS$MIN_POS)
    ## nr, not n: the objective spans the REFERENCE columns plus the slacks.
    ## The two are equal for ordinary DEA, which is why this indexing error
    ## survived every self-referenced test and only surfaced the first time a
    ## point was scored against a differently sized reference set.
    obj <- numeric(nr + p + q)
    if (orientation == "in") obj[nr + seq_len(p)]     <- 1/(p * xo)
    else                     obj[nr + p + seq_len(q)] <- 1/(q * yo)
    lpSolveAPI::set.objfn(lp, obj)
    lpSolveAPI::set.rhs(lp, c(X[o, ], Y[o, ]), seq_len(p + q))

    st[o] <- solve(lp)
    if (!st[o] %in% c(0L, 1L)) { eff[o] <- NA_real_; next }
    z <- lpSolveAPI::get.objective(lp)
    eff[o] <- if (orientation == "in") 1 - z else 1/(1 + z)
    v <- lpSolveAPI::get.variables(lp)
    if (peers) L[o, ] <- v[seq_len(nr)]
    sx[o, ] <- v[nr + seq_len(p)]
    sy[o, ] <- v[nr + p + seq_len(q)]
  }
  list(eff = eff, lambda = L, sx = sx, sy = sy, status = st)
}
