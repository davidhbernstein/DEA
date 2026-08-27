## ---------------------------------------------------------------------------
## Methods and extractors for class "dea".
##
## dea(), dea_sbm() and dea_ddf() all return the same class, exactly as the
## models they implement are all distances to the same estimated technology.
## The `model` field is what the methods below dispatch on internally, so a
## script that fits one and prints, plots or extracts from it does not need to
## change when it fits another.
## ---------------------------------------------------------------------------

print.dea <- function(x, ...) {
  cat("--- Data envelopment analysis ---\n")
  cat("model:       ", .dea_model_label(x), "\n", sep = "")
  cat("technology:  ", toupper(x$rts), ", ", .dea_orient_label(x), "\n", sep = "")
  cat("DMUs: ", x$n, "   inputs: ", x$p, "   outputs: ", x$q,
      "   (", format(round(x$total_time, 3)), " sec)\n", sep = "")
  ## Only worth a line when the two sets differ -- but then it is essential,
  ## because the scores are no longer bounded by 1 and a reader who has not
  ## noticed the external reference set will think something is broken.
  if (!isTRUE(x$self_ref)) {
    cat("scored against an EXTERNAL reference set of ", x$nref, " DMUs; ",
        "scores are not bounded by 1\n", sep = "")
  }
  if (isTRUE(x$super)) cat("Andersen-Petersen super-efficiency\n")

  e <- x$eff
  bad <- sum(!is.finite(e))
  cat("\n", .dea_measure_label(x), "\n", sep = "")
  print(round(summary(e[is.finite(e)]), 4))
  cat("\nefficient DMUs: ", sum(x$efficient, na.rm = TRUE), " of ", x$n, sep = "")
  if (x$model == "radial" && !is.null(x$slack_x)) {
    onfront <- sum(e == 1, na.rm = TRUE)
    if (onfront > sum(x$efficient, na.rm = TRUE)) {
      cat("  (", onfront, " score 1 radially, but ",
          onfront - sum(x$efficient, na.rm = TRUE),
          " of those still carry slack)", sep = "")
    }
  }
  cat("\n")
  if (bad) {
    cat("NOT SOLVED: ", bad, " DMU(s). ",
        if (isTRUE(x$super))
          "Super-efficiency under a non-constant-returns technology is genuinely infeasible for DMUs at the edge of the input space; NA is the correct answer there, not a large number."
        else "See $status for the lpSolve return codes.", "\n", sep = "")
  }
  invisible(x)
}

summary.dea <- function(object, ...) {
  print(object)
  n <- object$n
  cat("\nefficiency by DMU (first ", min(10L, n), " of ", n, ")\n", sep = "")
  tb <- data.frame(dmu = object$dmu, eff = round(object$eff, 4),
                   row.names = NULL, stringsAsFactors = FALSE)
  if (!is.null(object$slack_x)) {
    tb$slack <- round(rowSums(object$slack_x) + rowSums(object$slack_y), 6)
  }
  if (!is.null(object$sum_lambda)) tb$sum_lambda <- round(object$sum_lambda, 4)
  print(utils::head(tb, 10))
  if (!is.null(object$lambda)) {
    cnt <- colSums(object$lambda > 0, na.rm = TRUE)
    top <- sort(cnt[cnt > 0], decreasing = TRUE)
    cat("\nmost frequently referenced peers\n")
    print(utils::head(top, 5))
    cat("\nPeer sets are not unique when the linear program has multiple\n",
        "optima, which is common. Treat them as one valid benchmark set,\n",
        "not as the benchmark set.\n", sep = "")
  }
  invisible(object)
}

.dea_orient_label <- function(x) {
  if (identical(x$model, "ddf")) return(paste0("direction \"", x$direction, "\""))
  switch(x$orientation,
    `in` = "input orientation", out = "output orientation",
    none = "non-oriented", x$orientation)
}

.dea_model_label <- function(x) switch(x$model,
  radial   = "radial (Debreu-Farrell)",
  sbm      = "slacks-based measure (Tone 2001)",
  ddf      = "directional distance (Chambers-Chung-Fare 1996)",
  additive = paste0("additive (Charnes et al. 1985), ",
                    switch(x$measure, ram = "range adjusted",
                           mip = "inefficiency proportions",
                           unweighted = "unweighted"), " measure"),
  x$model)

.dea_measure_label <- function(x) switch(x$model,
  radial = if (x$orientation == "in") "theta (input-oriented, 1 = on the frontier)"
           else "phi (output-oriented, 1 = on the frontier)",
  sbm    = "rho (1 = Pareto-Koopmans efficient)",
  ddf    = "beta (INEFFICIENCY: 0 = on the frontier, larger is worse)",
  additive = if (identical(x$measure, "unweighted"))
      "total slack (INEFFICIENCY: 0 = efficient; carries the units of the data)"
    else if (identical(x$measure, "ram")) "rho_RAM (1 = Pareto-Koopmans efficient)"
    else "rho_MIP (1 = Pareto-Koopmans efficient)",
  "efficiency")

nobs.dea <- function(object, ...) object$n

## ---------------------------------------------------------------------------
## efficiency() -- the scores, optionally on one common scale.
##
## Four models report on four scales: theta in (0,1], phi in [1,Inf), rho in
## (0,1] and beta in [0,Inf).  type = "score" maps all of them onto (0,1] with
## 1 meaning efficient, which is what a table comparing models needs.
##
## It REFUSES to do so for a directional model whose direction is not purely
## input or purely output.  beta is an additive measure -- units of g -- and
## there is no ratio it corresponds to; inventing one such as (1-beta)/(1+beta)
## would produce a number that looks comparable and is not.
## ---------------------------------------------------------------------------
efficiency <- function(object, ...) UseMethod("efficiency")

efficiency.dea <- function(object, type = c("natural", "score"), ...) {
  type <- .match_arg_ci(type, c("natural", "score"), "type")
  if (type == "natural") return(object$eff)
  switch(object$model,
    radial = if (object$orientation == "in") object$eff else 1 / object$eff,
    sbm    = object$eff,
    additive = {
      if (identical(object$measure, "unweighted")) {
        stop("type = \"score\" is not defined for the unweighted additive ",
             "model. Its objective adds slacks measured in different units, ",
             "so it is a total in the data's own units rather than a ratio, ",
             "and it is not even invariant to those units -- see ?dea_add. ",
             "Use measure = \"ram\" or \"mip\", both of which normalize the ",
             "slacks and return a score in [0, 1].", call. = FALSE)
      }
      object$eff
    },
    ddf    = {
      if (identical(object$direction, "in"))  return(1 - object$eff)
      if (identical(object$direction, "out")) return(1 / (1 + object$eff))
      stop("type = \"score\" is not defined for a directional model with ",
           "direction \"", object$direction, "\". beta is an ADDITIVE distance ",
           "measured in units of g, and only the purely input (g = (x, 0)) and ",
           "purely output (g = (0, y)) directions correspond to a Farrell ",
           "ratio -- 1 - beta and 1/(1 + beta) respectively. Use ",
           "type = \"natural\" and report beta, or refit with ",
           "direction = \"in\" / \"out\".", call. = FALSE)
    },
    object$eff)
}

## Peers and their weights, as a tidy long table: one row per (DMU, peer) pair
## with a non-negligible lambda.
peers <- function(object, ...) UseMethod("peers")

## The default is .DEA_CONSTANTS$TOL_LAMBDA written out: a promise referring
## to an internal object cannot be documented without a codoc mismatch, and a
## user reading args(peers.dea) should see the number.
peers.dea <- function(object, threshold = 1e-8, ...) {
  if (is.null(object$lambda)) {
    stop("This fit was made with peers = FALSE, so the lambda matrix was not ",
         "kept. Refit with peers = TRUE.", call. = FALSE)
  }
  idx <- which(object$lambda > threshold, arr.ind = TRUE)
  if (!nrow(idx)) {
    return(data.frame(dmu = character(0), peer = character(0),
                      lambda = numeric(0), stringsAsFactors = FALSE))
  }
  idx <- idx[order(idx[, 1L], -object$lambda[idx]), , drop = FALSE]
  data.frame(dmu    = object$dmu[idx[, 1L]],
             peer   = object$dmu[idx[, 2L]],
             lambda = object$lambda[idx],
             row.names = NULL, stringsAsFactors = FALSE)
}

slacks <- function(object, ...) UseMethod("slacks")

slacks.dea <- function(object, ...) {
  if (is.null(object$slack_x)) {
    stop("This fit carries no slacks. ",
         if (identical(object$model, "ddf"))
           "The directional model projects onto the frontier along g and does not run a second-stage slack program."
         else "Refit dea() with slack = TRUE.", call. = FALSE)
  }
  out <- cbind(object$slack_x, object$slack_y)
  colnames(out) <- c(paste0("sx_", colnames(object$slack_x)),
                     paste0("sy_", colnames(object$slack_y)))
  out
}

## The frontier projection: where each DMU would sit if it were efficient.
## Radial contraction/expansion first, then the remaining slack, so the result
## is a Pareto-Koopmans point whenever slacks were computed.
fitted.dea <- function(object, ...) {
  X <- object$x; Y <- object$y
  if (identical(object$model, "ddf")) {
    px <- X - object$eff * object$g[, seq_len(object$p), drop = FALSE]
    py <- Y + object$eff * object$g[, object$p + seq_len(object$q), drop = FALSE]
  } else if (identical(object$model, "additive")) {
    px <- X - object$slack_x; py <- Y + object$slack_y
  } else if (identical(object$model, "radial")) {
    px <- if (object$orientation == "in") X * object$eff else X
    py <- if (object$orientation == "in") Y else Y * object$eff
    if (!is.null(object$slack_x)) { px <- px - object$slack_x; py <- py + object$slack_y }
  } else {
    px <- X - object$slack_x; py <- Y + object$slack_y
  }
  out <- cbind(px, py)
  colnames(out) <- c(colnames(X), colnames(Y))
  rownames(out) <- object$dmu
  out
}

## ---------------------------------------------------------------------------
## plot.dea()
##
## With one input and one output, draw the estimated frontier.  It is traced by
## EVALUATING the fitted technology on a grid of input levels rather than by
## joining up the efficient DMUs, which means the same code draws it correctly
## for every returns-to-scale assumption, including FDH's staircase, and puts
## the horizontal and vertical faces where they actually are.
##
## Otherwise, plot the distribution of scores -- with more than two dimensions
## there is nothing honest to draw of the frontier itself.
## ---------------------------------------------------------------------------
plot.dea <- function(x, ngrid = 400, ...) {
  if (x$p == 1L && x$q == 1L) return(.dea_plot_frontier(x, ngrid, ...))
  .dea_plot_scores(x, ...)
}

.dea_plot_frontier <- function(x, ngrid, ...) {
  X <- x$x; Y <- x$y
  sc <- .dea_scale(X, Y, x$scaling)
  xg <- seq(min(X), max(X), length.out = ngrid)
  Xg <- matrix(xg / sc$sx[1], ncol = 1L)
  Yg <- matrix(1 / sc$sy[1],  nrow = ngrid, ncol = 1L)
  rts <- if (identical(x$rts, "fdh")) "fdh" else x$rts
  if (identical(rts, "fdh")) {
    ## Free disposal: the frontier at any x is the largest output among DMUs
    ## with no more input, which is a step function.
    yg <- vapply(xg, function(v) { ok <- X[, 1L] <= v + 1e-12
                                   if (any(ok)) max(Y[ok, 1L]) else NA_real_ }, numeric(1))
  } else {
    B <- .lp_radial_build(sc$X, sc$Y, rts, "out")
    yg <- vapply(seq_len(ngrid), function(i) .lp_radial_at(B, Xg, Yg, i)$eff, numeric(1))
    yg <- yg * sc$sy[1]                     ## phi at y = 1 IS the frontier level
  }
  graphics::plot(X[, 1L], Y[, 1L], xlab = colnames(X)[1L], ylab = colnames(Y)[1L],
                 pch = 16, col = ifelse(x$efficient, "black", "grey55"),
                 main = paste0(toupper(x$rts), " frontier"), ...)
  graphics::lines(xg, yg, lwd = 2, col = "steelblue")
  graphics::legend("bottomright", bty = "n", pch = c(16, 16, NA), lty = c(NA, NA, 1),
                   lwd = c(NA, NA, 2), col = c("black", "grey55", "steelblue"),
                   legend = c("efficient", "inefficient", "estimated frontier"))
  invisible(list(x = xg, y = yg))
}

.dea_plot_scores <- function(x, ...) {
  e <- x$eff[is.finite(x$eff)]
  graphics::hist(e, breaks = "FD", col = "grey85", border = "white",
                 xlab = .dea_measure_label(x),
                 main = paste0(.dea_model_label(x), ", ", toupper(x$rts)), ...)
  graphics::abline(v = if (identical(x$model, "ddf")) 0 else 1,
                   lwd = 2, col = "steelblue")
  invisible(e)
}
