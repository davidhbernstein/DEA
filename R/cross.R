## ---------------------------------------------------------------------------
## dea_cross() -- cross-efficiency (Sexton, Silkman and Hogan 1986).
##
## Ordinary DEA lets every DMU choose the weights that flatter it most, which
## is what makes the scores generous and the efficient set large: at n close to
## p + q almost everything scores 1 and the ranking is empty.  Cross-efficiency
## answers that by making each DMU be judged by everyone else's weights too.
##
##   E[d, o] = (u_d' y_o) / (v_d' x_o)
##
## is DMU o's efficiency priced at DMU d's own optimal weights.  The diagonal
## E[d, d] is d's ordinary DEA score, the column mean is o's cross-efficiency,
## and the ranking that comes out of it is a peer appraisal rather than a
## self-appraisal.  It very rarely ties.
##
## THE WEIGHTS ARE NOT UNIQUE, AND THAT IS THE WHOLE DIFFICULTY.  An efficient
## DMU generally has a whole face of alternate optimal weight vectors, every
## one of them giving it the same score of 1 but giving OTHER DMUs quite
## different cross-efficiencies -- so a cross-efficiency computed from whatever
## vertex the solver happened to stop at is not reproducible across solvers,
## let alone across packages.  This is not a rounding problem; the spread is
## first-order.
##
## Doyle and Green (1994) fix it by choosing among the alternate optima with a
## stated secondary goal, holding the DMU's own score at theta_d and then
##
##   benevolent   maximizing the average cross-efficiency of everyone else
##   aggressive   minimizing it
##
## which brackets the answer instead of pretending to a single one.  Reporting
## both is the honest use: if a ranking survives from aggressive to benevolent
## it is a property of the data, and if it does not, it was a property of the
## solver.  `secondary = "none"` is available and is what most software does,
## but it is not the default here.
## ---------------------------------------------------------------------------

dea_cross <- function(x, y, data = NULL,
                      secondary = c("benevolent", "aggressive", "none"),
                      rts = c("crs", "vrs"),
                      self = TRUE,
                      scaling = TRUE,
                      xref = NULL, yref = NULL, dataref = NULL) {

  call <- match.call()
  secondary <- .match_arg_ci(secondary, c("benevolent", "aggressive", "none"),
                             "secondary")
  rts <- .match_arg_ci(rts, c("crs", "vrs"), "rts")
  t0 <- proc.time()[["elapsed"]]

  if (rts == "vrs" && secondary != "none") {
    stop("secondary = \"", secondary, "\" is defined here only for ",
         "rts = \"crs\". Under variable returns the cross-efficiency ratio ",
         "carries the intercept u0, is not bounded by 1 and can be negative, ",
         "so \"maximize everyone else's average\" is optimizing a quantity ",
         "that is not an efficiency. Use rts = \"crs\", or ",
         "secondary = \"none\" and read the result as descriptive.",
         call. = FALSE)
  }

  d <- .dea_data(x, y, data, xref, yref, dataref, FALSE, "dea_cross")
  X <- d$X; Y <- d$Y; XR <- d$XR; YR <- d$YR
  n <- nrow(X); nr <- nrow(XR); p <- ncol(X); q <- ncol(Y)

  ## The RATERS are the reference DMUs and the RATED are the evaluated ones.
  ## They coincide in ordinary use, and separating them is what lets a fixed
  ## panel of raters score newcomers.
  fit <- dea(XR, YR, rts = rts, orientation = "in", slack = FALSE,
             peers = FALSE, scaling = scaling, multipliers = TRUE)
  V <- fit$v; U <- fit$u; u0 <- fit$u0 %||% rep(0, nr)

  if (secondary != "none") {
    dg <- .cross_secondary(XR, YR, fit$eff, nr, p, q, scaling,
                           benevolent = identical(secondary, "benevolent"))
    V <- dg$v; U <- dg$u
  }

  ## E is nr x n: one ROW per rater, one COLUMN per rated DMU.
  denom <- V %*% t(X)
  numer <- U %*% t(Y) - u0
  bad <- abs(denom) < .DEA_CONSTANTS$MIN_POS
  denom[bad] <- NA_real_
  E <- numer / denom
  dimnames(E) <- list(d$ref, d$dmu)

  ## A DMU's own appraisal is one of the n it receives. Dropping it (Doyle and
  ## Green's "maverick" convention) is available because with small n the
  ## diagonal, which is always the most generous entry in its column, pulls the
  ## average up by 1/n of the whole gap.
  keep <- matrix(TRUE, nr, n)
  if (!self && d$self) diag(keep) <- FALSE
  Em <- E; Em[!keep] <- NA_real_
  cross <- colMeans(Em, na.rm = TRUE)
  spread <- apply(Em, 2L, function(z) diff(range(z, na.rm = TRUE)))
  names(cross) <- names(spread) <- d$dmu

  structure(list(
    eff = cross, model = "cross", rts = rts, secondary = secondary,
    self = self,
    cross_matrix = E, own = if (d$self) fit$eff else NULL, spread = spread,
    v = V, u = U, u0 = fit$u0,
    ## Doyle and Green's maverick index: how much a DMU's self-appraisal
    ## exceeds the appraisal it gets from everyone else. Large means a DMU that
    ## looks efficient only under weights nobody else would choose. It needs a
    ## self-appraisal to compare against, so it is NULL when the rated DMUs are
    ## not among the raters -- an outsider has no diagonal entry.
    maverick = if (d$self) (fit$eff - cross) / cross else NULL,
    x = X, y = Y, xref = XR, yref = YR, self_ref = d$self,
    dmu = d$dmu, ref = d$ref, n = n, nref = nr, p = p, q = q,
    scaling = scaling,
    total_time = proc.time()[["elapsed"]] - t0,
    call = call
  ), class = "dea_cross")
}

## ---------------------------------------------------------------------------
## The Doyle-Green secondary program, for rater d:
##
##   max / min  u' sum_{o != d} y_o
##   s.t.       v' sum_{o != d} x_o = 1
##              u'y_d - theta_d v'x_d = 0        <- d's own score is held fixed
##              u'Y_j - v'X_j <= 0               for every j
##              u, v >= 0
##
## Rows 1 and 2 both move with d, and row 2's coefficients move with theta_d as
## well, so only the technology rows are built once here -- still the majority
## of the matrix.
## ---------------------------------------------------------------------------
.cross_secondary <- function(XR, YR, theta, nr, p, q, scaling, benevolent) {
  sc  <- .dea_scale(XR, YR, scaling)
  Xs  <- sc$X; Ys <- sc$Y
  sumx <- colSums(Xs); sumy <- colSums(Ys)

  nv <- p + q
  lp <- lpSolveAPI::make.lp(2L + nr, nv)
  lpSolveAPI::lp.control(lp, sense = if (benevolent) "max" else "min",
                         epsel = .DEA_CONSTANTS$LP_EPSEL, verbose = "neutral")
  for (i in seq_len(p)) lpSolveAPI::set.column(lp, i,     c(0, 0, -Xs[, i]))
  for (r in seq_len(q)) lpSolveAPI::set.column(lp, p + r, c(0, 0,  Ys[, r]))
  lpSolveAPI::set.constr.type(lp, c("=", "=", rep("<=", nr)))
  lpSolveAPI::set.rhs(lp, c(1, 0, rep(0, nr)))

  V <- matrix(NA_real_, nr, p); U <- matrix(NA_real_, nr, q)
  for (dd in seq_len(nr)) {
    lpSolveAPI::set.row(lp, 1L, c(sumx - Xs[dd, ], rep(0, q)))
    lpSolveAPI::set.row(lp, 2L, c(-theta[dd] * Xs[dd, ], Ys[dd, ]))
    lpSolveAPI::set.objfn(lp, c(rep(0, p), sumy - Ys[dd, ]))
    st <- solve(lp)
    if (!st %in% c(0L, 1L)) next
    z <- lpSolveAPI::get.variables(lp)
    V[dd, ] <- z[seq_len(p)]; U[dd, ] <- z[p + seq_len(q)]
  }
  ## Back into the caller's units, exactly as in .dea_multipliers().
  list(v = sweep(V, 2L, sc$sx, "/"), u = sweep(U, 2L, sc$sy, "/"))
}
