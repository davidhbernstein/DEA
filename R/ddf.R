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
                    direction = c("both", "in", "out", "unit", "mean", "range"),
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

  G <- .ddf_direction(direction, X, Y, XR, YR, n, p, q)
  sc  <- .dea_scale(XR, YR, scaling && all(XR > 0) && all(YR > 0))
  XRs <- sc$X; YRs <- sc$Y
  Xs  <- sweep(X, 2L, sc$sx, "/"); Ys <- sweep(Y, 2L, sc$sy, "/")
  ## Rescaling the direction alongside the data is what keeps beta unchanged.
  Gs <- cbind(sweep(G[, seq_len(p), drop = FALSE], 2L, sc$sx, "/"),
              sweep(G[, p + seq_len(q), drop = FALSE], 2L, sc$sy, "/"))

  ## Built as a function rather than inline so that it can be REBUILT. The
  ## object is reused across all n DMUs and only the evaluated DMU's column and
  ## right-hand side are rewritten -- that is where the speed comes from -- but
  ## lpSolveAPI carries basis and factorisation state with it, and for some
  ## right-hand sides that inherited state is bad enough that the solve gives
  ## up. See the long note in .dea_radial(); a fresh object is the remedy.
  rr <- .rts_row(rts)
  build <- function() {
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
    lp
  }
  lp <- build()

  ## Point the program at DMU o. Kept separate from the loop so the retry can
  ## repeat it exactly rather than approximately.
  aim <- function(lp, o) {
    lpSolveAPI::set.column(lp, 1L,
      c(1, Gs[o, seq_len(p)], -Gs[o, p + seq_len(q)]),
      c(0L, seq_len(p), p + seq_len(q)))
    lpSolveAPI::set.rhs(lp, c(Xs[o, ], Ys[o, ]), seq_len(p + q))
  }

  ## Under direction = "range" a DMU can BE the ideal point, in which case its
  ## direction is zero and the program is unbounded rather than infeasible --
  ## beta appears in no constraint. The answer is beta = 0 and needs no LP: the
  ## DMU dominates the whole reference set, so there is nothing to improve
  ## toward. Peers are left NA rather than faked, because no program was solved
  ## and a lambda of all zeros would violate the vrs row.
  ideal <- attr(G, "ideal")
  if (is.null(ideal)) ideal <- logical(n)
  if (any(ideal)) {
    warning(sum(ideal), " DMU(s) attain the best value in EVERY input and ",
            "output of the reference set (first: ",
            paste(utils::head(which(ideal), 5), collapse = ", "),
            "), so they are the range direction's ideal point. Their beta is 0 ",
            "with no program to solve, and their peers are reported as NA.",
            call. = FALSE)
  }

  beta <- numeric(n); st <- integer(n)
  L <- if (peers) matrix(0, n, nr) else NULL
  for (o in seq_len(n)) {
    if (ideal[o]) {
      beta[o] <- 0; st[o] <- 0L
      if (peers) L[o, ] <- NA_real_
      next
    }
    aim(lp, o)
    st[o] <- solve(lp)
    ## An INFEASIBLE program (2) is not retried: that is an answer about the
    ## data, not a conditioning failure. Anything else is worth a fresh object.
    if (!st[o] %in% c(0L, 1L, 2L)) {
      lp2 <- build(); aim(lp2, o)
      st[o] <- solve(lp2)
      if (st[o] %in% c(0L, 1L)) {
        beta[o] <- lpSolveAPI::get.objective(lp2)
        if (peers) L[o, ] <- lpSolveAPI::get.variables(lp2)[-1L]
        next
      }
    }
    if (!st[o] %in% c(0L, 1L)) { beta[o] <- NA_real_; next }
    beta[o] <- lpSolveAPI::get.objective(lp)
    if (peers) L[o, ] <- lpSolveAPI::get.variables(lp)[-1L]
  }

  ## Report DMUs that did not solve. dea_ddf() recorded `status` and then said
  ## nothing about it, so a failed program returned beta = NA silently -- the
  ## same hole that dea() had before 1.0.2 and free disposal had before 1.0.3.
  ## Reachable from ordinary code: direction = "in" holds g_y at zero, so a DMU
  ## whose outputs exceed everything the reference set can produce has no
  ## feasible point at all.
  .dea_report_unsolved(st, n, d$self, "dea_ddf",
                       requires = "some convex combination of the reference DMUs")
  beta[is.finite(beta) & abs(beta) < .DEA_CONSTANTS$TOL_EFF] <- 0
  names(beta) <- dmu
  if (!is.null(L)) {
    L[!is.na(L) & abs(L) < .DEA_CONSTANTS$TOL_LAMBDA] <- 0
    dimnames(L) <- list(dmu, d$ref)
  }
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
##
## "range" is the odd one out and deliberately so: it is the only shorthand
## taken over the REFERENCE set rather than the evaluated set.  That is not a
## stylistic choice.  The range direction points at the technology's ideal
## point, and it is the reference DMUs that span the technology, so under vrs
## the bound X'lambda >= min(XR) gives beta <= 1 for free -- the property the
## model exists for.  Measured against the evaluated set it would not hold
## whenever the two sets differ.  "mean" keeps using the evaluated set because
## it is only a scale and nothing depends on which set it comes from.
.ddf_direction <- function(direction, X, Y, XR, YR, n, p, q) {
  if (is.character(direction)) {
    direction <- .match_arg_ci(direction,
                               c("both", "in", "out", "unit", "mean", "range"),
                               "direction")
    G <- switch(direction,
      both = cbind(X, Y),
      `in` = cbind(X, matrix(0, n, q)),
      out  = cbind(matrix(0, n, p), Y),
      unit = matrix(1, n, p + q),
      mean = matrix(rep(c(colMeans(X), colMeans(Y)), each = n), n, p + q),
      range = .ddf_range_direction(X, Y, XR, YR, n, p, q))
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
  ## An all-zero row is fatal in general -- there is no distance to measure
  ## along -- but under "range" it is not an error at all, it is the ideal
  ## point, and dea_ddf() answers those DMUs without an LP. See
  ## .ddf_range_direction().
  ideal <- attr(G, "ideal")
  zero <- rowSums(abs(G)) <= 0
  if (!is.null(ideal)) zero <- zero & !ideal
  if (any(zero)) {
    bad <- which(zero)
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

## The range direction of Portela, Thanassoulis and Simpson (2004):
##
##   g_o = (x_o - min XR,  max YR - y_o)
##
## the vector from the DMU to the technology's IDEAL POINT, the corner of the
## reference set's bounding box that is best in every coordinate at once. Two
## properties follow, and both are worth being precise about because both are
## vrs properties and neither survives a cone.
##
## BETA LIES IN [0, 1]. Under vrs, X'lambda is a convex combination of the
## reference inputs, so X'lambda >= min XR componentwise; feasibility then
## forces x_o - beta g_x >= min XR, i.e. beta <= 1, and the same argument on
## the output side. So beta reads directly as "the fraction of the distance to
## the ideal point that this DMU could travel". Under crs, nirs or ndrs the
## convex-combination bound is gone -- lambda can scale a reference DMU past
## the box -- and beta can exceed 1. That is not a bug in the model, it is the
## bound being a statement about the technology rather than about the
## direction, and ?dea_ddf says so.
##
## AND BETA IS TRANSLATION INVARIANT, which is the reason the model exists.
## Adding a constant to an input column leaves g unchanged (it is a difference
## of two values in that column) and shifts both sides of the constraint by
## t * sum(lambda) -- equal only when sum(lambda) = 1. So translation
## invariance is exactly the vrs case too, and it is what lets negative data be
## handled without translating it at all.
.ddf_range_direction <- function(X, Y, XR, YR, n, p, q) {
  lo <- apply(XR, 2L, min)
  hi <- apply(YR, 2L, max)
  gx <- sweep(X, 2L, lo, "-")
  gy <- -sweep(Y, 2L, hi, "-")

  ## An evaluated DMU can sit outside the reference set's box when xref/yref
  ## are supplied separately; self-referenced, this never fires. A negative
  ## component would ask an input that is already better than anything in the
  ## technology to GROW, so it is clamped to zero -- that coordinate simply has
  ## no room to improve toward. Clamping is announced, because it changes the
  ## estimand for those DMUs and silence would hide that.
  out <- rowSums(gx < 0) + rowSums(gy < 0)
  if (any(out > 0)) {
    who <- which(out > 0)
    warning(length(who), " DMU(s) lie outside the reference set's range in at ",
            "least one coordinate (first: ", paste(utils::head(who, 5), collapse = ", "),
            "). Those coordinates are already better than anything the ",
            "technology contains, so their direction components are clamped to ",
            "zero rather than pointing backwards.", call. = FALSE)
    gx[gx < 0] <- 0
    gy[gy < 0] <- 0
  }

  G <- cbind(gx, gy)
  ## A DMU that attains the minimum in EVERY input and the maximum in EVERY
  ## output is the ideal point itself. Its direction is zero, and the program
  ## would be unbounded rather than infeasible -- beta does not appear in any
  ## constraint. It is also, by construction, the answer: there is nowhere to
  ## improve toward, so beta = 0. dea_ddf() fills those in without solving.
  attr(G, "ideal") <- rowSums(abs(G)) <= 0
  G
}
