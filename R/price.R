## ---------------------------------------------------------------------------
## Efficiency when prices are known: cost, revenue and profit.
##
## Every other estimator in this package measures distance to the frontier
## without asking what anything is worth.  That is the right default -- prices
## are usually unavailable and often unreliable -- but when they ARE known they
## answer a question technical efficiency cannot: a DMU can sit exactly on the
## frontier and still be using the wrong INPUT MIX for the prices it faces.
##
## The decomposition is the point:
##
##   cost efficiency  =  technical efficiency  x  allocative efficiency
##
## the first factor asking "could this DMU produce the same output with less of
## everything", the second "could it produce it more cheaply by substituting".
## Both are in (0, 1] and a DMU is cost efficient only if it is both.
##
## THE PROGRAMS REDUCE.  The textbook cost program carries the optimal input
## vector as a free variable,
##
##   min w'x   s.t.  X'lambda <= x,  Y'lambda >= y_o,  (rts),  lambda >= 0
##
## but with w >= 0 the objective pushes x down onto X'lambda at every optimum,
## so x can be substituted out and the program solved in lambda alone:
##
##   min (w'X')lambda   s.t.  Y'lambda >= y_o,  (rts)
##
## which is p rows and p columns smaller, and -- the part that matters for the
## sweep -- leaves a constraint matrix that does not depend on the DMU at all.
## Only the objective (through w_o) and the right-hand side (through y_o) move,
## so the build-once-mutate design of the rest of the package applies unchanged.
## The optimal input vector is recovered afterwards as x* = X'lambda*.
## ---------------------------------------------------------------------------

dea_cost <- function(x, y, w, data = NULL,
                     rts = c("vrs", "crs", "nirs", "ndrs"),
                     peers = TRUE,
                     scaling = TRUE,
                     xref = NULL, yref = NULL, dataref = NULL) {
  .dea_price_fit("cost", x, y, w, NULL, data, rts, peers, scaling,
                 xref, yref, dataref, match.call())
}

dea_revenue <- function(x, y, r, data = NULL,
                        rts = c("vrs", "crs", "nirs", "ndrs"),
                        peers = TRUE,
                        scaling = TRUE,
                        xref = NULL, yref = NULL, dataref = NULL) {
  .dea_price_fit("revenue", x, y, NULL, r, data, rts, peers, scaling,
                 xref, yref, dataref, match.call())
}

## ---------------------------------------------------------------------------
## Profit is not a ratio, and pretending otherwise is the usual mistake.
##
## Cost and revenue efficiency are ratios because both numerator and
## denominator are non-negative.  Observed profit is not: it is routinely zero
## and sometimes negative, so Pi_o / Pi* is undefined, unbounded, or negative
## exactly where the question is most interesting.  The measure used here is
## therefore the NERLOVIAN one -- the profit gap, normalized by the value of a
## direction --
##
##   Nerlovian inefficiency  =  (Pi* - Pi_o) / (w'g_x + r'g_y)
##
## which is 0 for a profit-maximizing DMU, larger for a worse one, and in the
## same units as the directional distance function.  That is not a coincidence:
## the Chambers-Chung-Fare (1996) duality splits it exactly,
##
##   Nerlovian  =  beta (technical)  +  allocative
##
## with beta the directional distance in the same direction g, and the
## remainder the loss from facing these prices with the wrong mix.  Both parts
## are >= 0 and BIGGER IS WORSE, the reverse of the cost and revenue measures
## and for the same reason dea_ddf()'s beta runs that way.
##
## CONSTANT AND NON-DECREASING RETURNS ARE REFUSED, and this is not a gap in
## the implementation.  Maximum profit over a cone is unbounded the moment any
## reference DMU turns a profit at the evaluated DMU's prices -- double it,
## double the profit -- so Pi* is either 0 or infinite and the measure carries
## no information. Profit maximization is a variable-returns question.
## ---------------------------------------------------------------------------
dea_profit <- function(x, y, w, r, data = NULL,
                       direction = c("both", "in", "out", "unit", "mean"),
                       rts = c("vrs", "nirs"),
                       peers = TRUE,
                       scaling = TRUE,
                       xref = NULL, yref = NULL, dataref = NULL) {
  call <- match.call()
  ## Catch the two refused technologies BEFORE the match, so that asking for
  ## one gets the reason rather than a list of the two that remain.
  if (!missing(rts) && length(rts) == 1L) {
    a <- tolower(as.character(rts))
    if (a %in% names(.RTS_ALIAS)) a <- .RTS_ALIAS[[a]]
    if (a %in% c("crs", "ndrs")) {
      stop("dea_profit() has no answer under rts = \"", a, "\". The ",
           "technology is a cone in the profit-increasing direction, so as ",
           "soon as one reference DMU is profitable at this DMU's prices the ",
           "maximum profit is unbounded -- scaling that DMU up scales its ",
           "profit up with it. Maximum profit is a variable-returns concept; ",
           "use rts = \"vrs\" or \"nirs\".", call. = FALSE)
    }
  }
  rts <- .match_arg_ci(rts, c("vrs", "nirs"), "rts")
  t0 <- proc.time()[["elapsed"]]

  d <- .dea_data(x, y, data, xref, yref, dataref, FALSE, "dea_profit",
                 allow_negative = TRUE)
  X <- d$X; Y <- d$Y; XR <- d$XR; YR <- d$YR
  n <- nrow(X); nr <- nrow(XR); p <- ncol(X); q <- ncol(Y)
  W <- .dea_prices(w, data, n, p, "w", colnames(X))
  R <- .dea_prices(r, data, n, q, "r", colnames(Y))
  G <- .ddf_direction(direction, X, Y, n, p, q)

  ## max (r'Y' - w'X')lambda subject to the returns-to-scale row alone: the
  ## input and output constraints are what x and y were substituted out of.
  rr <- .rts_row(rts)
  lp <- lpSolveAPI::make.lp(1L, nr)
  lpSolveAPI::lp.control(lp, sense = "max", epsel = .DEA_CONSTANTS$LP_EPSEL,
                         verbose = "neutral")
  for (j in seq_len(nr)) lpSolveAPI::set.column(lp, j, 1, 1L)
  lpSolveAPI::set.constr.type(lp, rr, 1L)
  lpSolveAPI::set.rhs(lp, 1, 1L)

  profit_max <- numeric(n); st <- integer(n)
  L <- if (peers) matrix(0, n, nr) else NULL
  for (o in seq_len(n)) {
    lpSolveAPI::set.objfn(lp, as.numeric(YR %*% R[o, ] - XR %*% W[o, ]))
    st[o] <- solve(lp)
    if (!st[o] %in% c(0L, 1L)) { profit_max[o] <- NA_real_; next }
    profit_max[o] <- lpSolveAPI::get.objective(lp)
    if (peers) L[o, ] <- lpSolveAPI::get.variables(lp)
  }
  if (any(!st %in% c(0L, 1L))) {
    .dea_report_unsolved(st, n, d$self, "dea_profit")
  }

  profit_obs <- rowSums(R * Y) - rowSums(W * X)
  norm <- rowSums(W * G[, seq_len(p), drop = FALSE]) +
          rowSums(R * G[, p + seq_len(q), drop = FALSE])
  if (any(norm <= 0)) {
    stop("The direction has zero or negative value at ", sum(norm <= 0),
         " DMU(s): w'g_x + r'g_y must be strictly positive, since it is what ",
         "the profit gap is measured in. Check the prices and the direction.",
         call. = FALSE)
  }
  nerl <- (profit_max - profit_obs) / norm

  tech <- suppressWarnings(dea_ddf(X, Y, direction = direction, rts = rts,
                                   peers = FALSE, scaling = scaling,
                                   xref = if (d$self) NULL else XR,
                                   yref = if (d$self) NULL else YR))$eff
  ## Snap: the two programs reach the same frontier point for a technically
  ## efficient DMU, and a residual of 1e-13 should read as zero allocative
  ## loss rather than as a tiny real one.
  alloc <- nerl - tech
  alloc[is.finite(alloc) & abs(alloc) < .DEA_CONSTANTS$TOL_EFF] <- 0
  nerl[is.finite(nerl) & abs(nerl) < .DEA_CONSTANTS$TOL_EFF] <- 0

  if (!is.null(L)) {
    L[abs(L) < .DEA_CONSTANTS$TOL_LAMBDA] <- 0
    dimnames(L) <- list(d$dmu, d$ref)
  }
  names(nerl) <- names(tech) <- names(alloc) <- d$dmu
  names(profit_max) <- names(profit_obs) <- d$dmu

  structure(list(
    eff = nerl, model = "profit", rts = rts, direction = direction,
    technical = tech, allocative = alloc,
    profit_max = profit_max, profit_obs = profit_obs, normalization = norm,
    lambda = L, status = st,
    x = X, y = Y, w = W, r = R, xref = XR, yref = YR, self_ref = d$self,
    dmu = d$dmu, ref = d$ref, n = n, nref = nr, p = p, q = q,
    scaling = scaling,
    total_time = proc.time()[["elapsed"]] - t0,
    call = call
  ), class = "dea_price")
}

## ---------------------------------------------------------------------------
## The shared cost/revenue engine.  The two models are mirror images -- minimize
## the value of the inputs holding output fixed, or maximize the value of the
## outputs holding input fixed -- so they are one function with the sense, the
## constrained side and the ratio's orientation flipped.
## ---------------------------------------------------------------------------
.dea_price_fit <- function(model, x, y, w, r, data, rts, peers, scaling,
                           xref, yref, dataref, call) {
  rts <- .match_arg_ci(rts, c("vrs", "crs", "nirs", "ndrs"), "rts")
  t0 <- proc.time()[["elapsed"]]
  cost <- identical(model, "cost")

  d <- .dea_data(x, y, data, xref, yref, dataref, FALSE, paste0("dea_", model))
  X <- d$X; Y <- d$Y; XR <- d$XR; YR <- d$YR
  n <- nrow(X); nr <- nrow(XR); p <- ncol(X); q <- ncol(Y)
  P <- if (cost) .dea_prices(w, data, n, p, "w", colnames(X))
       else       .dea_prices(r, data, n, q, "r", colnames(Y))

  ## Scaling is a genuine no-op here, not merely a harmless one: dividing a
  ## column by s and multiplying its price by s leaves every cost and every
  ## revenue numerically identical, so only the conditioning changes.
  sc  <- .dea_scale(XR, YR, scaling)
  XRs <- sc$X; YRs <- sc$Y
  Xs  <- sweep(X, 2L, sc$sx, "/"); Ys <- sweep(Y, 2L, sc$sy, "/")
  Ps  <- if (cost) sweep(P, 2L, sc$sx, "*") else sweep(P, 2L, sc$sy, "*")

  ## The constrained side: outputs must be matched for a cost program, inputs
  ## must not be exceeded for a revenue one.
  CM  <- if (cost) YRs else XRs           ## constraint columns
  ce  <- if (cost) Ys  else Xs            ## and the evaluated levels
  k   <- ncol(CM)
  rr  <- .rts_row(rts)
  nrows <- k + (!is.null(rr))
  lp <- lpSolveAPI::make.lp(nrows, nr)
  lpSolveAPI::lp.control(lp, sense = if (cost) "min" else "max",
                         epsel = .DEA_CONSTANTS$LP_EPSEL, verbose = "neutral")
  rts_one <- if (is.null(rr)) numeric(0) else 1
  for (j in seq_len(nr)) lpSolveAPI::set.column(lp, j, c(CM[j, ], rts_one))
  lpSolveAPI::set.constr.type(lp, rep(if (cost) ">=" else "<=", k), seq_len(k))
  if (!is.null(rr)) {
    lpSolveAPI::set.constr.type(lp, rr, k + 1L)
    lpSolveAPI::set.rhs(lp, 1, k + 1L)
  }

  VS  <- if (cost) XRs else YRs           ## the side being valued
  opt <- numeric(n); st <- integer(n)
  L   <- if (peers) matrix(0, n, nr) else NULL
  OPTQ <- matrix(NA_real_, n, if (cost) p else q)
  for (o in seq_len(n)) {
    lpSolveAPI::set.objfn(lp, as.numeric(VS %*% Ps[o, ]))
    lpSolveAPI::set.rhs(lp, ce[o, ], seq_len(k))
    st[o] <- solve(lp)
    if (!st[o] %in% c(0L, 1L)) { opt[o] <- NA_real_; next }
    opt[o] <- lpSolveAPI::get.objective(lp)
    lam <- lpSolveAPI::get.variables(lp)
    if (peers) L[o, ] <- lam
    OPTQ[o, ] <- as.numeric(crossprod(VS, lam))   ## x* (or y*), still scaled
  }
  if (any(!st %in% c(0L, 1L))) {
    .dea_report_unsolved(st, n, d$self, paste0("dea_", model))
  }

  obs <- rowSums(P * (if (cost) X else Y))
  ## Ratios, both oriented so that 1 is best and smaller is worse.
  eff <- if (cost) opt / obs else obs / opt
  eff[is.finite(eff) & abs(eff - 1) < .DEA_CONSTANTS$TOL_EFF] <- 1

  ## The technical factor is the radial score on the SAME side, so that the
  ## product decomposition is exact rather than approximate.
  tech <- suppressWarnings(dea(X, Y, rts = rts,
                               orientation = if (cost) "in" else "out",
                               slack = FALSE, peers = FALSE, scaling = scaling,
                               xref = if (d$self) NULL else XR,
                               yref = if (d$self) NULL else YR))$eff
  if (!cost) tech <- 1 / tech             ## phi >= 1 onto (0, 1]
  alloc <- eff / tech
  alloc[is.finite(alloc) & abs(alloc - 1) < .DEA_CONSTANTS$TOL_EFF] <- 1

  OPTQ <- sweep(OPTQ, 2L, if (cost) sc$sx else sc$sy, "*")
  dimnames(OPTQ) <- list(d$dmu, colnames(if (cost) X else Y))
  if (!is.null(L)) {
    L[abs(L) < .DEA_CONSTANTS$TOL_LAMBDA] <- 0
    dimnames(L) <- list(d$dmu, d$ref)
  }
  names(eff) <- names(tech) <- names(alloc) <- names(opt) <- names(obs) <- d$dmu

  structure(list(
    eff = eff, model = model, rts = rts,
    technical = tech, allocative = alloc,
    optimal = opt, observed = obs, optimal_q = OPTQ,
    lambda = L, status = st,
    x = X, y = Y, prices = P, xref = XR, yref = YR, self_ref = d$self,
    dmu = d$dmu, ref = d$ref, n = n, nref = nr, p = p, q = q,
    scaling = scaling,
    total_time = proc.time()[["elapsed"]] - t0,
    call = call
  ), class = "dea_price")
}

## ---------------------------------------------------------------------------
## Price coercion.  Deliberately NOT .dea_matrix(): a bare numeric vector means
## something different for prices than it does for data.  Given p inputs, a
## length-p vector is ONE price list faced by every DMU -- the common case, and
## a competitive market is the usual assumption -- whereas .dea_matrix() would
## read it as n observations on a single variable.  Guessing between the two
## when n happens to equal p is not possible, so that case is refused.
## ---------------------------------------------------------------------------
.dea_prices <- function(arg, data, n, k, what, varnames) {
  if (is.null(arg)) {
    stop("`", what, "` is required: ", if (what == "w") "input" else "output",
         " prices cannot be defaulted.", call. = FALSE)
  }
  if (inherits(arg, "formula") || is.character(arg)) {
    M <- .dea_matrix(arg, data, what)
  } else if (is.data.frame(arg) || is.matrix(arg)) {
    M <- .dea_num(as.matrix(arg), what)
  } else if (is.numeric(arg) && is.vector(arg)) {
    if (length(arg) == k && (n != k || k == 1L)) {
      M <- matrix(arg, nrow = n, ncol = k, byrow = TRUE)
    } else if (length(arg) == n && k == 1L) {
      M <- matrix(arg, ncol = 1L)
    } else if (length(arg) == k && n == k) {
      stop("`", what, "` is a vector of length ", k, ", and with ", n,
           " DMUs and ", k, " ", if (what == "w") "inputs" else "outputs",
           " that is ambiguous -- one price list shared by every DMU, or one ",
           "DMU-specific price for a single variable? Pass a ", n, " x ", k,
           " matrix to say which.", call. = FALSE)
    } else {
      stop("`", what, "` has length ", length(arg), ", which is neither ", k,
           " (one price list shared by every DMU) nor a ", n, " x ", k,
           " matrix (one list per DMU).", call. = FALSE)
    }
    M <- .dea_num(M, what)
  } else {
    stop("`", what, "` must be a numeric vector, matrix or data frame, or -- ",
         "with `data` -- column names or a one-sided formula.", call. = FALSE)
  }
  if (nrow(M) == 1L && n > 1L) M <- M[rep(1L, n), , drop = FALSE]
  if (nrow(M) != n || ncol(M) != k) {
    stop("`", what, "` is ", nrow(M), " x ", ncol(M), " but must be ", n,
         " x ", k, ": one row per DMU, one column per ",
         if (what == "w") "input" else "output", ".", call. = FALSE)
  }
  if (any(M <= 0)) {
    stop("`", what, "` contains ", sum(M <= 0), " non-positive price(s). The ",
         "reduction that substitutes the optimal quantity out of the program ",
         "relies on the objective pushing it against its bound, which needs ",
         "strictly positive prices.", call. = FALSE)
  }
  colnames(M) <- varnames
  M
}
