## ---------------------------------------------------------------------------
## dea_weak() -- the weak-disposability technology for undesirable outputs, and
## the directional distance function on it (Fare and Grosskopf; Chung, Fare and
## Grosskopf 1997).
##
## WHAT THE OTHER TWO ROUTES CANNOT DO.  dea_undesirable() offers the Seiford-Zhu
## translation and the bad-as-input route, and both are honest rearrangements of
## an ordinary program.  Neither encodes the one thing everybody actually
## believes about pollution: that you cannot simply stop producing it.  Under
## free disposal a DMU may always cut its bad output at no cost, so the model
## will cheerfully prescribe less pollution and the same product.
##
## WEAK DISPOSABILITY IS AN EQUALITY ROW.  The bads are not freely disposable:
## they can be reduced only by scaling the whole activity back.
##
##   sum_j lambda_j x_ij <= x_io - beta g_xi      inputs, freely disposable
##   sum_j lambda_j y_rj >= y_ro + beta g_yr      goods,  freely disposable
##   sum_j lambda_j b_fj  = b_fo - beta g_bf      BADS,   weakly disposable
##
## The direction answers the question that matters: with g = (0, y_o, b_o), beta
## is the proportion by which good output could rise AND bad output fall at the
## same time.  beta = 0 is on the frontier and larger is worse, as in dea_ddf().
##
## NULL-JOINTNESS IS A PROPERTY OF THE DATA, NOT A CONSTRAINT.  The technology
## only means what it is supposed to mean if no good output can be had without
## some bad: a DMU with positive goods and zero bads says the opposite, and one
## such row makes the frontier claim that clean production is already feasible.
## It cannot be imposed by a constraint, so it is CHECKED and reported.
##
## THE VRS FORMULATION IS DISPUTED, AND THE DISPUTE IS SETTLED HERE BY
## MEASUREMENT RATHER THAN BY CITATION.  Fare and Grosskopf's variable-returns
## version applies one abatement factor to the whole sample, which in linear
## form is the single-intensity program above with sum(lambda) = 1.  Kuosmanen
## (2005) argued that this is not the minimal technology satisfying weak
## disposability and variable returns, and gave a two-vector form in which each
## DMU carries its own abatement factor:
##
##   sum_j (lambda_j + mu_j) x_ij <= x_io - beta g_xi
##   sum_j lambda_j y_rj          >= y_ro + beta g_yr
##   sum_j lambda_j b_fj           = b_fo - beta g_bf
##   sum_j (lambda_j + mu_j)       = 1
##
## so that sum(lambda) may fall below 1 while the input side still averages to
## one.  That is strictly more freedom -- mu = 0 is always available -- so the
## Kuosmanen technology CONTAINS the Fare-Grosskopf one and its beta can only be
## the larger.  That much is an identity.  HOW MUCH larger is not, and the
## answer depends almost entirely on the DIRECTION, which neither paper
## discusses.  50 replications at n = 60, share of DMUs where mu is positive and
## the mean and largest gap in beta:
##
##   direction   coupling   mu used   mean gap   max gap
##   goods_bads  prop          6.6%    0.00091     0.113
##   goods_bads  steep         0.7%    0.00017     0.107
##   bads        prop         16.6%    0.00435     0.203
##   bads        flat         30.0%    0.01347     0.393
##   all         flat          2.1%    0.00039     0.097
##
## On the standard goods-and-bads direction the second vector is used for a few
## per cent of DMUs and the mean gap is a thousandth -- but the LARGEST gap is
## 0.1 even there, so the choice bites a few units hard rather than everyone a
## little.  It matters most when the direction asks for bads to be cut and
## nothing else, which is exactly when scaling the activity back is the
## attractive move and mu is what permits it.
##
## What it does NOT much change is who is on the frontier: the efficient set
## differs between the two formulations for 0.2% to 0.7% of DMUs.  So the
## formulation is a question about the SIZE of measured inefficiency far more
## than about its classification.
##
## UNDER CRS THEY COINCIDE, and that is checked rather than assumed: with no
## convexity row mu only adds slack to input rows that free disposal already
## allows, so it is redundant.  `technology` is therefore ignored under crs and
## says so.
##
## The default is Kuosmanen's, because it is the one that satisfies minimal
## extrapolation; Fare and Grosskopf's is kept for reproducing published results,
## the same reason dea_mult() keeps the 1982 model.
## ---------------------------------------------------------------------------

dea_weak <- function(x, y, bad, data = NULL,
                     direction = c("goods_bads", "all", "bads"),
                     rts = c("vrs", "crs"),
                     technology = c("kuosmanen", "fare"),
                     peers = TRUE, scaling = TRUE) {

  call <- match.call()
  .dea_check_data_arg(data)
  rts <- .match_arg_ci(rts, c("vrs", "crs"), "rts")
  technology <- .match_arg_ci(technology, c("kuosmanen", "fare"), "technology")
  t0 <- proc.time()[["elapsed"]]

  X <- .dea_matrix(x, data, "x")
  Y <- .dea_matrix(y, data, "y")
  B <- .dea_matrix(bad, data, "bad")
  if (ncol(B) == 1L && identical(colnames(B), "bad1")) colnames(B) <- "bad"
  if (nrow(B) != nrow(X) || nrow(Y) != nrow(X)) {
    stop("`x`, `y` and `bad` describe different numbers of DMUs (", nrow(X),
         ", ", nrow(Y), ", ", nrow(B), ").", call. = FALSE)
  }
  if (any(X < 0) || any(Y < 0) || any(B < 0)) {
    stop("Negative values in `x`, `y` or `bad`. The weak-disposability ",
         "technology is defined on non-negative quantities.", call. = FALSE)
  }
  n <- nrow(X); p <- ncol(X); q <- ncol(Y); nb <- ncol(B)
  dmu <- .dea_dmu_names(X, Y)
  d <- p + q + nb
  if (n < 3 * d) {
    warning("Only ", n, " DMUs for ", p, " inputs, ", q, " good and ", nb,
            " bad output(s); with n < 3(p+q+b) = ", 3 * d, " most DMUs are ",
            "efficient by dimension alone.", call. = FALSE)
  }

  if (identical(rts, "crs") && identical(technology, "fare")) {
    message("Under rts = \"crs\" the two formulations are the same program -- ",
            "with no convexity row, Kuosmanen's second intensity vector only ",
            "adds slack that free disposal already allows. `technology` is ",
            "ignored here.")
  }

  ## The unusable-column check goes FIRST. An all-zero bad column also makes
  ## every DMU look like a null-jointness violation, so running that check
  ## first would bury a fatal error under a warning about its symptom -- the
  ## same ordering lesson as dea_panel()'s period count before its duplicates.
  if (any(colSums(B) <= 0)) {
    stop("Bad output column(s) ",
         paste(colnames(B)[colSums(B) <= 0], collapse = ", "),
         " are zero for every DMU. A bad nobody produces has an equality row ",
         "that forces beta to 0 for everyone; drop the column.", call. = FALSE)
  }
  nj <- .weak_nulljoint(Y, B)
  if (nj$n_clean > 0L) {
    warning(nj$n_clean, " DMU(s) produce positive good output with ZERO ",
            "bad output (first: ", paste(utils::head(nj$which, 5), collapse = ", "),
            "). That contradicts null-jointness, which is the assumption the ",
            "whole technology rests on: it says clean production is already ",
            "observed, and the frontier will say so too. Null-jointness cannot ",
            "be imposed by a constraint -- it is a property of the data -- so ",
            "this is reported rather than fixed. See ?dea_weak.", call. = FALSE)
  }

  G <- .weak_direction(direction, X, Y, B, n, p, q, nb)

  ## Scale all three blocks, and the direction with them, so beta is unchanged.
  pos <- scaling && all(X > 0) && all(Y > 0) && all(B > 0)
  sx <- if (pos) colMeans(X) else rep(1, p)
  sy <- if (pos) colMeans(Y) else rep(1, q)
  sb <- if (pos) colMeans(B) else rep(1, nb)
  Xs <- sweep(X, 2L, sx, "/"); Ys <- sweep(Y, 2L, sy, "/"); Bs <- sweep(B, 2L, sb, "/")
  Gs <- cbind(sweep(G[, seq_len(p), drop = FALSE], 2L, sx, "/"),
              sweep(G[, p + seq_len(q), drop = FALSE], 2L, sy, "/"),
              sweep(G[, p + q + seq_len(nb), drop = FALSE], 2L, sb, "/"))

  two <- identical(technology, "kuosmanen") && identical(rts, "vrs")
  res <- .weak_solve(Xs, Ys, Bs, Gs, rts, two, peers, n, p, q, nb)

  beta <- res$beta
  beta[is.finite(beta) & abs(beta) < .DEA_CONSTANTS$TOL_EFF] <- 0
  names(beta) <- dmu
  L <- res$lambda
  if (!is.null(L)) { L[abs(L) < .DEA_CONSTANTS$TOL_LAMBDA] <- 0; dimnames(L) <- list(dmu, dmu) }
  M <- res$mu
  if (!is.null(M)) { M[abs(M) < .DEA_CONSTANTS$TOL_LAMBDA] <- 0; dimnames(M) <- list(dmu, dmu) }

  .dea_report_unsolved(res$status, n, TRUE, "dea_weak",
                       requires = "a combination of reference DMUs matching the bad outputs exactly")

  structure(list(
    eff = beta, model = "weak", rts = rts, orientation = "none",
    technology = technology, super = FALSE,
    lambda = L, mu = M, sum_lambda = if (is.null(L)) NULL else rowSums(L),
    direction = G, efficient = is.finite(beta) & beta == 0,
    status = res$status,
    bad = B, nbad = nb, null_joint = nj$n_clean == 0L, n_clean = nj$n_clean,
    x = X, y = Y, xref = X, yref = Y, self_ref = TRUE,
    dmu = dmu, ref = dmu, n = n, nref = n, p = p, q = q,
    scaling = pos, slack = FALSE,
    total_time = proc.time()[["elapsed"]] - t0,
    call = call
  ), class = c("dea_weak", "dea"))
}

## Null-jointness: no good output without some bad. A DMU with positive goods
## and no bads at all is the violation that matters; it is the row that lets the
## frontier claim clean production is already feasible.
.weak_nulljoint <- function(Y, B) {
  clean <- rowSums(B) <= 0 & rowSums(Y) > 0
  list(n_clean = sum(clean), which = which(clean))
}

.weak_direction <- function(direction, X, Y, B, n, p, q, nb) {
  if (is.character(direction) || missing(direction)) {
    direction <- .match_arg_ci(direction, c("goods_bads", "all", "bads"),
                               "direction")
    G <- switch(direction,
      ## The Chung-Fare-Grosskopf default: more good and less bad, together.
      goods_bads = cbind(matrix(0, n, p), Y, B),
      all        = cbind(X, Y, B),
      bads       = cbind(matrix(0, n, p), matrix(0, n, q), B))
  } else if (is.matrix(direction) || is.data.frame(direction)) {
    G <- as.matrix(direction)
    if (nrow(G) != n || ncol(G) != p + q + nb) {
      stop("`direction` as a matrix must be ", n, " x ", p + q + nb,
           " (inputs, then good outputs, then bads); got ", nrow(G), " x ",
           ncol(G), ".", call. = FALSE)
    }
    storage.mode(G) <- "double"
  } else if (is.numeric(direction)) {
    if (length(direction) != p + q + nb) {
      stop("`direction` as a vector must have length p + q + b = ", p + q + nb,
           "; got ", length(direction), ".", call. = FALSE)
    }
    G <- matrix(rep(direction, each = n), n, p + q + nb)
  } else {
    stop("`direction` must be \"goods_bads\", \"all\", \"bads\", or numeric.",
         call. = FALSE)
  }
  if (any(!is.finite(G))) stop("`direction` contains non-finite values.", call. = FALSE)
  if (any(G < 0)) {
    stop("`direction` has negative entries. Every component here is a ",
         "MAGNITUDE: the signs are already in the program, which raises good ",
         "output and lowers inputs and bads. A negative entry would reverse ",
         "one of those silently.", call. = FALSE)
  }
  bad_zero <- rowSums(abs(G)) <= 0
  if (any(bad_zero)) {
    stop("`direction` is all zero for DMU(s) ",
         paste(utils::head(which(bad_zero), 5), collapse = ", "),
         "; there is no distance to measure along.", call. = FALSE)
  }
  G
}

## The program. Columns: beta (1), lambda (n), and mu (n) when the two-vector
## Kuosmanen form is in use. Rows: inputs, goods, bads, then the convexity row
## under vrs.
.weak_solve <- function(X, Y, B, G, rts, two, peers, n, p, q, nb) {
  vrs <- identical(rts, "vrs")
  nrows <- p + q + nb + vrs
  ncols <- 1L + n + (if (two) n else 0L)

  build <- function() {
    lp <- lpSolveAPI::make.lp(nrows, ncols)
    lpSolveAPI::lp.control(lp, sense = "max", epsel = .DEA_CONSTANTS$LP_EPSEL,
                           verbose = "neutral")
    one <- if (vrs) 1 else numeric(0)
    for (j in seq_len(n)) {
      ## lambda enters every row: inputs, goods, bads and the convexity row.
      lpSolveAPI::set.column(lp, j + 1L, c(X[j, ], Y[j, ], B[j, ], one))
      ## mu enters ONLY the input rows and the convexity row. That is the whole
      ## of Kuosmanen's correction: it lets sum(lambda) fall below 1 while the
      ## input side still averages to one.
      if (two) {
        idx <- c(seq_len(p), if (vrs) nrows else integer(0))
        lpSolveAPI::set.column(lp, 1L + n + j, c(X[j, ], one), idx)
      }
    }
    lpSolveAPI::set.constr.type(lp, rep("<=", p), seq_len(p))
    lpSolveAPI::set.constr.type(lp, rep(">=", q), p + seq_len(q))
    ## The equality row IS weak disposability.
    lpSolveAPI::set.constr.type(lp, rep("=", nb), p + q + seq_len(nb))
    if (vrs) {
      lpSolveAPI::set.constr.type(lp, "=", nrows)
      lpSolveAPI::set.rhs(lp, 1, nrows)
    }
    obj <- numeric(ncols); obj[1L] <- 1
    lpSolveAPI::set.objfn(lp, obj)
    ## beta is left free: an evaluated point outside the technology has a
    ## genuinely negative solution, and clamping it would report that point as
    ## efficient. Same reasoning as dea_ddf().
    lpSolveAPI::set.bounds(lp, lower = -.LP_INF, upper = .LP_INF, columns = 1L)
    lp
  }

  aim <- function(lp, o) {
    ## beta's own column: +g_x on the input rows (they tighten as beta grows),
    ## -g_y on the good rows, +g_b on the bad rows.
    lpSolveAPI::set.column(lp, 1L,
      c(1, G[o, seq_len(p)], -G[o, p + seq_len(q)], G[o, p + q + seq_len(nb)]),
      c(0L, seq_len(p + q + nb)))
    lpSolveAPI::set.rhs(lp, c(X[o, ], Y[o, ], B[o, ]), seq_len(p + q + nb))
  }

  lp <- build()
  beta <- numeric(n); st <- integer(n)
  L <- if (peers) matrix(0, n, n) else NULL
  M <- if (peers && two) matrix(0, n, n) else NULL

  for (o in seq_len(n)) {
    aim(lp, o)
    r <- .lp_solve_retry(lp, fresh = function() { l2 <- build(); aim(l2, o); l2 })
    st[o] <- r$status
    if (!r$status %in% c(0L, 1L)) {
      beta[o] <- NA_real_
      if (peers) L[o, ] <- NA_real_
      if (!is.null(M)) M[o, ] <- NA_real_
      next
    }
    v <- lpSolveAPI::get.variables(r$lp)
    beta[o] <- v[1L]
    if (peers) L[o, ] <- v[1L + seq_len(n)]
    if (!is.null(M)) M[o, ] <- v[1L + n + seq_len(n)]
  }
  list(beta = beta, lambda = L, mu = M, status = st)
}

print.dea_weak <- function(x, ...) {
  cat("--- Weak disposability of undesirable outputs ---\n")
  cat("technology:  ", toupper(x$rts),
      if (identical(x$rts, "vrs")) paste0(", ", x$technology, " formulation") else "",
      "\n", sep = "")
  cat("DMUs: ", x$n, "   inputs: ", x$p, "   good outputs: ", x$q,
      "   bad outputs: ", x$nbad, "   (", format(round(x$total_time, 3)),
      " sec)\n", sep = "")
  cat("bad output(s): ", paste(colnames(x$bad), collapse = ", "), "\n", sep = "")
  cat("\nbeta: the proportion by which good output could rise and bad output\n",
      "fall together; 0 is on the frontier and larger is worse\n", sep = "")
  print(round(summary(x$eff[is.finite(x$eff)]), 4))
  cat("\non the frontier (beta = 0): ", sum(x$efficient, na.rm = TRUE), " of ",
      x$n, "\n", sep = "")
  cat("null-jointness: ",
      if (x$null_joint) "holds -- no DMU produces goods without bads"
      else paste0("VIOLATED by ", x$n_clean, " DMU(s); see ?dea_weak"), "\n", sep = "")
  invisible(x)
}
