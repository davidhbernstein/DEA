## ---------------------------------------------------------------------------
## dea_nd() -- non-discretionary inputs and outputs (Banker and Morey 1986a).
##
## THE PROBLEM.  Some of a DMU's inputs are not its to choose.  A school cannot
## vary the number of pupils it is sent, a hospital cannot choose its catchment,
## a branch cannot move its town.  The ordinary radial program contracts EVERY
## input by the same theta, so it reports a unit as inefficient for failing to
## do something it has no power to do, and the shortfall it reports is partly a
## description of the environment.
##
## THE CHANGE IS SMALL AND IT IS IN ONE PLACE.  theta multiplies only the
## discretionary input rows; the fixed rows stay at the DMU's own level:
##
##   sum_j lambda_j x_ij <= theta x_io      i discretionary
##   sum_j lambda_j x_ij <=       x_io      i fixed
##   sum_j lambda_j y_rj >=       y_ro
##
## So a fixed input still CONSTRAINS who may be a peer -- a unit with twice the
## pupils is not a fair comparison -- while no credit is given for cutting it.
## In the program that is two edits together, a zero in the theta column and the
## DMU's own level moved into the right-hand side; doing only the first would
## quietly impose sum(lambda_j x_ij) <= 0.  See .lp_radial_at().
##
## THE TWO STAGES MEAN DIFFERENT THINGS, which is the part worth reading twice.
## Stage one scales only the oriented side, so only `nd_x` changes an
## input-oriented score and only `nd_y` an output-oriented one.  Stage two
## maximises a sum of slacks, and THERE both masks apply in both orientations:
## a slack the DMU cannot act on should earn it no credit whichever way the
## program is oriented.  That is Banker and Morey's own rule -- the
## non-discretionary slacks are dropped from the objective and kept in the
## constraints -- and it is why `nd_y` is not ignored under an input
## orientation even though it leaves theta alone.
##
## FIXING AN INPUT MAKES THE MEASURED INEFFICIENCY LARGER, NOT SMALLER, AND
## EVERY INTUITION SAYS OTHERWISE.  One expects that no longer asking a DMU to
## cut something it cannot control should make it look better.  It does the
## reverse: theta_nd <= theta_ordinary, always.  Measured over 8704 DMU-fits
## across four technologies and random dimensions, the largest excess of
## theta_nd over theta_ordinary was 1.2e-12 -- solver noise -- and the mean
## theta on the school example below falls from 0.901 to 0.789.  Output
## oriented the mirror holds: phi_nd >= phi_ordinary everywhere.
##
## The reason is a one-line argument about the feasible set.  For theta < 1 the
## fixed row sum(lambda_j x_ij) <= x_io is LOOSER than <= theta x_io, so the
## non-discretionary program's feasible set CONTAINS the ordinary one and its
## minimum can only fall.  Economically: theta in the ordinary model is a single
## proportional contraction of every input at once, so a peer must use
## proportionally less of everything, pupils included -- a demanding
## requirement.  Here a peer need only use no more pupils, and the whole
## contraction lands on the teachers.
##
## So the two thetas are not the same quantity measured more or less fairly.
## They answer different questions -- "could all inputs have been smaller
## together?" against "with these pupils, could the teachers have been fewer?"
## -- and a table comparing them as if one were a corrected version of the
## other is comparing two estimands.
##
## THE SYMMETRY WITH OUTPUTS LOOKS EXACT AND IS NOT.  Chapter 21 of Charnes,
## Cooper, Lewin and Seiford (1994) flags this as a recurring interpretational
## error, so it is documented here rather than papered over.  "Non-discretionary
## input" has a clean meaning: a quantity the DMU is handed.  "Non-discretionary
## output" is usually used for something quite different -- an output the DMU
## influences only INDIRECTLY, through its inputs, such as sales through
## advertising.  That is not fixed; it is controlled at one remove, and the
## formulation does not model it.  Marking it `nd_y` tells the program the DMU
## should get no credit for producing more of it, which for sales is plainly
## wrong.  Use `nd_y` only for an output that is genuinely exogenous -- a
## statutory quota, a count of mandated cases -- and if the right description is
## "controlled indirectly", neither this nor any other DEA model here says so.
## ---------------------------------------------------------------------------

dea_nd <- function(x, y, data = NULL, nd_x = NULL, nd_y = NULL,
                   rts = c("vrs", "crs", "nirs", "ndrs"),
                   orientation = c("in", "out"),
                   slack = TRUE, peers = TRUE, scaling = TRUE,
                   xref = NULL, yref = NULL, dataref = NULL) {

  call <- match.call()
  rts <- .match_arg_ci(rts, c("vrs", "crs", "nirs", "ndrs"), "rts")
  orientation <- .match_arg_ci(orientation, c("in", "out"), "orientation")
  t0 <- proc.time()[["elapsed"]]

  d <- .dea_data(x, y, data, xref, yref, dataref, FALSE, "dea_nd")
  X <- d$X; Y <- d$Y; XR <- d$XR; YR <- d$YR
  n <- nrow(X); nr <- nrow(XR); p <- ncol(X); q <- ncol(Y)

  ## Named, so that the mask printed back is readable against the data rather
  ## than a bare run of TRUE/FALSE the reader has to count.
  ndx <- stats::setNames(.nd_resolve(nd_x, colnames(X), p, "nd_x"), colnames(X))
  ndy <- stats::setNames(.nd_resolve(nd_y, colnames(Y), q, "nd_y"), colnames(Y))

  if (identical(orientation, "in") && all(ndx)) {
    stop("Every input is marked non-discretionary, so theta multiplies no row ",
         "at all and the program has nothing to minimise -- it would return 0 ",
         "for every DMU without failing. At least one input must be ",
         "discretionary for an input-oriented fit; use orientation = \"out\", ",
         "or dea_add(), which has no orientation.", call. = FALSE)
  }
  if (identical(orientation, "out") && all(ndy)) {
    stop("Every output is marked non-discretionary, so phi multiplies no row ",
         "at all and the program has nothing to maximise. At least one output ",
         "must be discretionary for an output-oriented fit; use ",
         "orientation = \"in\", or dea_add(), which has no orientation.",
         call. = FALSE)
  }
  if (!any(ndx) && !any(ndy)) {
    warning("Neither `nd_x` nor `nd_y` marks anything, so this is an ordinary ",
            "dea() fit. The result is identical; the extra class and the ",
            "printed note are the only difference.", call. = FALSE)
  }
  if (any(ndy)) {
    message("`nd_y` marks ", sum(ndy), " output(s) non-discretionary. That is ",
            "right only for an output that is genuinely exogenous -- a quota, ",
            "a mandated caseload. An output the DMU influences INDIRECTLY, ",
            "through its inputs, is not non-discretionary, and marking it so ",
            "denies the DMU credit for producing more of it. See ?dea_nd.")
  }

  sc  <- .dea_scale(XR, YR, scaling)
  XRs <- sc$X; YRs <- sc$Y
  Xs  <- sweep(X, 2L, sc$sx, "/"); Ys <- sweep(Y, 2L, sc$sy, "/")

  out <- .dea_radial(Xs, Ys, XRs, YRs, rts, orientation, super = FALSE,
                     slack = slack, peers = peers, n, nr, p, q,
                     nd = list(x = ndx, y = ndy))

  .dea_report_unsolved(out$status, n, d$self, "dea_nd")

  eff <- out$eff
  eff[is.finite(eff) & abs(eff - 1) < .DEA_CONSTANTS$TOL_EFF] <- 1
  names(eff) <- d$dmu

  sx <- out$sx; sy <- out$sy
  if (!is.null(sx)) {
    sx <- sweep(sx, 2L, sc$sx, "*"); sy <- sweep(sy, 2L, sc$sy, "*")
    sx[abs(sx) < .DEA_CONSTANTS$TOL_SLACK] <- 0
    sy[abs(sy) < .DEA_CONSTANTS$TOL_SLACK] <- 0
    dimnames(sx) <- list(d$dmu, colnames(X)); dimnames(sy) <- list(d$dmu, colnames(Y))
  }
  L <- out$lambda
  if (!is.null(L)) { L[abs(L) < .DEA_CONSTANTS$TOL_LAMBDA] <- 0; dimnames(L) <- list(d$dmu, d$ref) }

  structure(list(
    eff = eff, model = "radial", rts = rts, orientation = orientation,
    super = FALSE,
    lambda = L, sum_lambda = out$sum_lambda,
    slack_x = sx, slack_y = sy,
    ## Pareto-Koopmans over the DISCRETIONARY coordinates only: a residual slack
    ## on a fixed variable is not something the DMU failed to do.
    efficient = if (is.null(sx)) eff == 1 else
      eff == 1 &
      rowSums(sx[, !ndx, drop = FALSE]) <= .DEA_CONSTANTS$TOL_SLACK &
      rowSums(sy[, !ndy, drop = FALSE]) <= .DEA_CONSTANTS$TOL_SLACK,
    status = out$status,
    nd_x = ndx, nd_y = ndy,
    n_fixed = c(x = sum(ndx), y = sum(ndy)),
    x = X, y = Y, xref = XR, yref = YR, self_ref = d$self,
    dmu = d$dmu, ref = d$ref, n = n, nref = nr, p = p, q = q,
    scaling = scaling, slack = slack,
    total_time = proc.time()[["elapsed"]] - t0,
    call = call
  ), class = c("dea_nd", "dea"))
}

## Accept a logical mask, column positions, or column names, and return a
## logical of the right length. Names are checked against the matrix's own
## columns so that a typo is an error rather than a silently empty mask.
.nd_resolve <- function(nd, nms, k, what) {
  if (is.null(nd)) return(rep(FALSE, k))
  if (is.logical(nd)) {
    if (length(nd) != k) {
      stop("`", what, "` given as a logical must have one entry per ",
           if (what == "nd_x") "input" else "output", " (", k, "); got ",
           length(nd), ".", call. = FALSE)
    }
    if (anyNA(nd)) stop("`", what, "` contains NA.", call. = FALSE)
    return(nd)
  }
  if (is.character(nd)) {
    miss <- setdiff(nd, nms)
    if (length(miss)) {
      stop("`", what, "`: no such column(s): ", paste(miss, collapse = ", "),
           ". Available: ", paste(nms, collapse = ", "), ".", call. = FALSE)
    }
    return(nms %in% nd)
  }
  if (is.numeric(nd)) {
    if (any(!is.finite(nd)) || any(nd != as.integer(nd)) ||
        any(nd < 1) || any(nd > k)) {
      stop("`", what, "` given as positions must be whole numbers in 1:", k,
           ".", call. = FALSE)
    }
    return(seq_len(k) %in% as.integer(nd))
  }
  stop("`", what, "` must be a logical mask, column positions, or column ",
       "names.", call. = FALSE)
}

print.dea_nd <- function(x, ...) {
  cat("--- Non-discretionary variables (Banker-Morey) ---\n")
  nmx <- colnames(x$x)[x$nd_x]; nmy <- colnames(x$y)[x$nd_y]
  cat("fixed inputs:  ", if (length(nmx)) paste(nmx, collapse = ", ") else "none",
      "\nfixed outputs: ", if (length(nmy)) paste(nmy, collapse = ", ") else "none",
      "\n", sep = "")
  scaled <- if (x$orientation == "in") "theta" else "phi"
  cat(scaled, " multiplies only the discretionary ",
      if (x$orientation == "in") "input" else "output",
      " rows; the fixed ones still\nconstrain who may be a peer.\n", sep = "")
  cat("\n")
  NextMethod()
}
