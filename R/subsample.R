## ---------------------------------------------------------------------------
## dea_subsample() -- confidence intervals by the m out of n bootstrap
## (Simar and Wilson 2011; Jeong and Simar 2006; Kneip, Simar and Wilson 2008).
##
## WHY A SECOND INTERVAL, WHEN dea_boot() ALREADY GIVES ONE.  Because the one it
## gives is not known to be valid.  The smoothed homogeneous bootstrap of Simar
## and Wilson (1998) resamples all n observations, and for a boundary estimator
## that is not consistent: the DEA score's limiting distribution is driven by how
## many observations land near the frontier, and an n-out-of-n resample
## reproduces the sample's own boundary rather than the population's.  Simar and
## Wilson (2011) set this out and propose subsampling instead -- draw m < n
## WITHOUT replacement, with m growing more slowly than n.
##
## This package now has a measured reason to care rather than a citation.
## dea_rts_test() implements the 2002 full-sample bootstrap, and on a design
## where the null is exactly true its power against returns = 0.8 is 2.0% --
## BELOW its own size -- with rDEA agreeing on the same datasets.  See
## horserace/README.md.
##
## THE MACHINERY.  Write kappa for the estimator's convergence exponent, which
## dea_rate() already knows: 2/(p+q+1) under vrs, 2/(p+q) under crs, 1/(p+q)
## under free disposal.  Then n^kappa (theta_n - theta) has a limiting
## distribution, and the subsample statistic
##
##     T*_m = m^kappa (theta*_m - theta_n)
##
## estimates it -- where theta*_m scores the SAME observation against a
## subsample of size m as the reference set.  Inverting gives
##
##     [theta_n - n^-kappa q_(1-alpha/2),  theta_n - n^-kappa q_(alpha/2)]
##
## and the formula needs no orientation branch: input-oriented theta is biased
## up and theta*_m is larger still, output-oriented phi is biased down and
## phi*_m smaller still, and the signs carry themselves.
##
## IT ALSO COVERS THE CASE dea_boot() REFUSES.  Free disposal has no smoothed
## homogeneous bootstrap -- the hull is not convex -- and dea_boot() says so and
## stops.  Subsampling is exactly the scheme Jeong and Simar (2006) propose for
## it, so rts = "fdh" is accepted here.
##
## CHOOSING m IS THE WHOLE DIFFICULTY.  The theory asks for m -> infinity with
## m/n -> 0 and says nothing about which m to use at the n in front of you, and
## the interval is genuinely sensitive to it.  The minimum-volatility rule of
## Politis, Romano and Wolf is used: compute the interval across a grid of m,
## and take the m whose neighbourhood of the grid moves the endpoints least, on
## the argument that a stable stretch is where the asymptotics have bitten.
##
## WHAT SUBSAMPLING BUYS IS PAID FOR IN INFEASIBILITY, AND THE BILL FALLS ON THE
## DMUs YOU CARE ABOUT.  Scoring an original observation against a subsample can
## have no solution at all: under variable returns the point may lie outside the
## smaller hull, and then there is no answer rather than a large one.  That is
## not an edge case.  Measured on dea_sim(returns = 0.9), output-oriented vrs,
## the share of replications with no solution:
##
##            m = n^0.4   m = n^0.6   m = n^0.8
##   n = 60     31.3%       15.7%        6.3%     all DMUs
##              58.1%       39.2%       22.1%     DMUs the fit calls efficient
##              26.5%       11.5%        3.5%     the interior
##   n = 150    27.3%       13.3%        5.0%     all DMUs
##              59.6%       42.7%       24.9%     efficient
##
## and the worst single DMU lost 63% to 97% of its draws.  So the interval for a
## frontier DMU -- the one an analyst actually wants -- is built from the fewest
## observations, and a smaller m makes it worse at exactly the rate that makes
## the asymptotics better.  This is the real cost of the method and it is not
## something the papers dwell on.
##
## So `n_valid` is reported per DMU, an interval needs at least ceiling(2/alpha)
## valid draws to be reported at all -- below that the alpha/2 quantile is
## extrapolation into a tail with nothing in it rather than an order statistic --
## and anything short of that comes back NA with a warning naming the count.
##
## ONE CHOICE IS MINE AND THE METHOD DOES NOT MAKE IT.  Politis, Romano and Wolf
## are choosing m for ONE parameter; here there are n intervals at once, and a
## per-DMU choice would be n noisy selections from the same small B.  So the
## volatility is computed on the MEAN endpoints across DMUs and one m is used for
## all of them.  `m_grid`, `endpoint_path` and `volatility` are returned so the
## choice can be inspected and overridden rather than trusted.
## ---------------------------------------------------------------------------

dea_subsample <- function(object, m = NULL, B = 2000, alpha = 0.05,
                          m_grid = NULL, B_select = 200L, window = 3L,
                          seed = NULL, ncores = 1L, progress = interactive()) {

  call <- match.call()
  if (!inherits(object, "dea")) {
    stop("`object` must be a fit from dea().", call. = FALSE)
  }
  if (!identical(object$model, "radial")) {
    stop("`object` is a ", object$model, " fit; dea_subsample() is written for ",
         "the radial estimator, whose convergence rate dea_rate() knows. The ",
         "non-radial measures have rates of their own and are not covered.",
         call. = FALSE)
  }
  if (isTRUE(object$super)) {
    stop("The Andersen-Petersen super-efficiency score is not an estimator of ",
         "a distance to the true frontier, so there is no sampling ",
         "distribution here to approximate.", call. = FALSE)
  }
  if (!isTRUE(object$self_ref)) {
    stop("Subsampling resamples the reference set, so it needs a fit whose ",
         "reference set IS the sample. This fit was scored against an explicit ",
         "`xref`/`yref`.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || B < 2) {
    stop("`B` must be a single number of replications, at least 2.", call. = FALSE)
  }
  B <- as.integer(B)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be in (0, 1).", call. = FALSE)
  }
  eff <- as.numeric(object$eff)
  if (any(!is.finite(eff))) {
    stop(sum(!is.finite(eff)), " of the supplied fit's scores are non-finite; ",
         "there is nothing to build an interval around for those DMUs.",
         call. = FALSE)
  }

  t0 <- proc.time()[["elapsed"]]
  X <- object$x; Y <- object$y
  n <- object$n; p <- object$p; q <- object$q
  rts <- object$rts; ori <- object$orientation
  kappa <- -dea_rate(p, q, rts, what = "estimator")

  if (!is.null(seed)) {
    old <- .rng_snapshot(); on.exit(.rng_restore(old), add = TRUE)
    set.seed(seed)
  }

  ## Scaling once, outside every loop: a subsample is a subset of the SAME
  ## columns, so the column means stay valid and the programs stay conditioned.
  sc <- .dea_scale(X, Y, object$scaling)
  Xs <- sc$X; Ys <- sc$Y

  ## One replication: score every original observation against a subsample of
  ## size m drawn WITHOUT replacement. Without replacement is the whole point --
  ## with replacement is the n-out-of-n bootstrap this exists to replace.
  one_rep <- function(mm) function(b) {
    idx <- sample.int(n, mm, replace = FALSE)
    .dea_radial(Xs, Ys, Xs[idx, , drop = FALSE], Ys[idx, , drop = FALSE],
                rts, ori, super = FALSE, slack = FALSE, peers = FALSE,
                n, mm, p, q)$eff
  }
  ## An interval needs enough valid draws for its quantile to be an order
  ## statistic rather than an extrapolation: at alpha = 0.05 the 0.025 quantile
  ## needs 40 points before any observation sits below it.
  ## `min_valid` differs between the two uses on purpose. The reported interval
  ## needs enough draws for its quantile to be real; the SELECTION path only has
  ## to be comparable across m, is run at a deliberately small B_select, and
  ## would be entirely NA under the strict rule -- which is how this was found,
  ## as `m` came back empty and sample.int() refused it.
  need <- ceiling(2 / alpha)
  ends <- function(mm, reps, min_valid) {
    M <- do.call(cbind, .dea_replicate(one_rep(mm), reps, ncores, FALSE))
    Tst <- mm^kappa * (M - eff)                      ## n x reps
    nv <- rowSums(is.finite(Tst))
    qq <- function(pr) vapply(seq_len(n), function(i) {
      if (nv[i] < min_valid) return(NA_real_)
      stats::quantile(Tst[i, ], probs = pr, na.rm = TRUE, names = FALSE)
    }, numeric(1))
    list(lo = eff - n^(-kappa) * qq(1 - alpha / 2),
         hi = eff - n^(-kappa) * qq(alpha / 2),
         mean_T = ifelse(nv >= min_valid, rowMeans(Tst, na.rm = TRUE), NA_real_),
         n_valid = nv, boot = M)
  }

  grid <- NULL; vol <- NULL; path <- NULL
  if (is.null(m)) {
    grid <- .sub_grid(m_grid, n)
    if (progress) cat("  selecting m over", length(grid), "candidates\n")
    path <- matrix(NA_real_, length(grid), 2L,
                   dimnames = list(as.character(grid), c("lower", "upper")))
    floor_sel <- max(5L, ceiling(0.1 * as.integer(B_select)))
    for (i in seq_along(grid)) {
      e <- ends(grid[i], as.integer(B_select), floor_sel)
      path[i, ] <- c(mean(e$lo, na.rm = TRUE), mean(e$hi, na.rm = TRUE))
      if (progress) { cat("\r   m =", grid[i], "   "); utils::flush.console() }
    }
    if (progress) cat("\r", strrep(" ", 30), "\r", sep = "")
    vol <- .sub_volatility(path, window)
    el <- .sub_eligible(length(grid), window)
    if (el$window < as.integer(window)) {
      warning("The grid of ", length(grid), " candidates is too short for a ",
              "window of ", window, "; the volatility rule used a window of ",
              el$window, " instead. A narrower window is a weaker rule -- pass ",
              "a longer `m_grid`, or `m` directly.", call. = FALSE)
    }
    ## A candidate whose own endpoints could not be computed at all cannot be
    ## chosen, and neither can one whose window straddles such a gap.
    eligible <- el$idx[is.finite(vol[el$idx])]
    if (!length(eligible)) {
      stop("No candidate m produced a usable interval: at every m on the grid ",
           "too many DMUs fell outside the subsample's hull for the endpoints ",
           "to be computed. Raise `B_select`, give a coarser `m_grid` with ",
           "larger values, or pass `m` directly.", call. = FALSE)
    }
    m <- grid[eligible[which.min(vol[eligible])]]
    window_used <- el$window
  } else {
    if (!is.numeric(m) || length(m) != 1L || m < 2 || m >= n) {
      stop("`m` must be a single number with 2 <= m < n = ", n, ".", call. = FALSE)
    }
    m <- as.integer(m)
    window_used <- NA_integer_
  }

  fin <- ends(m, B, need)
  bc <- eff - n^(-kappa) * fin$mean_T

  tab <- data.frame(dmu = object$dmu, eff = eff,
                    bias_corrected = as.numeric(bc),
                    ci_lower = as.numeric(fin$lo),
                    ci_upper = as.numeric(fin$hi),
                    n_valid = as.integer(fin$n_valid),
                    row.names = NULL, stringsAsFactors = FALSE)
  if (any(fin$n_valid < need)) {
    bad <- which(fin$n_valid < need)
    warning(length(bad), " of ", n, " DMU(s) have fewer than ", need,
            " usable replications and are reported as NA (fewest: ",
            min(fin$n_valid), " of ", B, "). Scoring a point against a ",
            "SUBSAMPLE can be infeasible -- it may lie outside the smaller ",
            "hull -- and that falls hardest on the DMUs the fit calls ",
            "efficient, which lose 20-60% of their draws on ordinary designs. ",
            "A larger `m` reduces it, at the cost of the asymptotics the ",
            "method exists for. See ?dea_subsample.", call. = FALSE)
  }

  structure(list(
    table = tab, boot = fin$boot, m = m, m_grid = grid,
    endpoint_path = path, volatility = vol, window = as.integer(window),
    window_used = window_used,
    eligible = if (is.null(grid)) NULL else grid[eligible],
    kappa = kappa, B = B, B_select = as.integer(B_select), alpha = alpha,
    orientation = ori, rts = rts, n = n, p = p, q = q, dmu = object$dmu,
    fit = object, total_time = proc.time()[["elapsed"]] - t0, call = call
  ), class = "dea_subsample")
}

## The candidate grid, geometric between about n^0.4 and n^0.9. Geometric rather
## than arithmetic because the theory is about ORDERS of m, so equal ratios are
## the equal steps.
.sub_grid <- function(m_grid, n) {
  if (!is.null(m_grid)) {
    g <- sort(unique(as.integer(m_grid)))
    if (any(g < 2L) || any(g >= n)) {
      stop("`m_grid` must lie in 2:", n - 1L, ".", call. = FALSE)
    }
    if (length(g) < 3L) {
      stop("`m_grid` needs at least 3 candidates for the volatility rule to ",
           "mean anything; got ", length(g), ".", call. = FALSE)
    }
    return(g)
  }
  lo <- max(5, ceiling(n^0.4)); hi <- max(lo + 2, floor(n^0.9))
  hi <- min(hi, n - 1L)
  if (hi <= lo) {
    stop("n = ", n, " is too small to choose m by the volatility rule; pass ",
         "`m` directly.", call. = FALSE)
  }
  g <- unique(round(exp(seq(log(lo), log(hi), length.out = 12L))))
  g <- g[g >= 2L & g < n]
  if (length(g) < 3L) {
    stop("n = ", n, " leaves fewer than 3 usable candidates for m; pass `m` ",
         "directly.", call. = FALSE)
  }
  as.integer(g)
}

## Minimum volatility: the standard deviation of each endpoint over a window of
## neighbouring grid points, summed.
##
## ONLY CANDIDATES WITH A FULL WINDOW ARE ELIGIBLE, and that is not fussiness.
## A clipped window at the edge of the grid averages over fewer neighbours, and
## a standard deviation over fewer points is smaller on average for no reason
## connected to stability -- so clipping silently biases the rule toward the
## smallest and largest m. Observed before the fix: on a 100-DMU fit the rule
## picked the first grid point, m = 7, at volatility 0.0107 against 0.0109 for
## m = 16 in the middle of a genuinely flat stretch, purely because m = 7 was
## scored on 4 neighbours and m = 16 on 7. The eligible indices are therefore
## w+1 .. k-w; the volatility is still reported for every candidate so the whole
## path can be read, with the ineligible ends marked NA in the selection.
.sub_volatility <- function(path, window) {
  k <- nrow(path)
  w <- max(1L, as.integer(window))
  vapply(seq_len(k), function(i) {
    j <- seq(max(1L, i - w), min(k, i + w))
    if (length(j) < 2L) return(NA_real_)
    sum(apply(path[j, , drop = FALSE], 2L, stats::sd))
  }, numeric(1))
}

## Which grid points the rule may choose from: those whose window is complete.
## If the grid is too short for any, the window shrinks until something is
## eligible rather than the call failing -- a narrow window is a weaker rule,
## not a broken one, and the alternative is refusing a perfectly ordinary n.
.sub_eligible <- function(k, window) {
  for (w in seq(as.integer(window), 1L)) {
    e <- seq_len(k)[seq_len(k) > w & seq_len(k) <= k - w]
    if (length(e)) return(list(idx = e, window = w))
  }
  list(idx = seq_len(k), window = 0L)
}

print.dea_subsample <- function(x, ...) {
  cat("--- Subsampling intervals (m out of n) ---\n")
  cat("model: radial, ", x$rts, ", ",
      if (x$orientation == "in") "input" else "output", " orientation\n", sep = "")
  cat("n = ", x$n, "   m = ", x$m,
      if (is.null(x$m_grid)) "  (supplied)" else "  (minimum volatility)",
      "   kappa = ", format(x$kappa, digits = 4), "\n", sep = "")
  cat("B = ", x$B, if (!is.null(x$m_grid)) paste0("   B_select = ", x$B_select),
      "   ", format(round(x$total_time, 2)), " sec\n", sep = "")
  if (!is.null(x$m_grid)) {
    cat("\nm chosen from: ", paste(x$m_grid, collapse = " "), "\n", sep = "")
  }
  tb <- x$table
  cat("\nmean interval width: ",
      format(mean(tb$ci_upper - tb$ci_lower, na.rm = TRUE), digits = 4),
      "   mean correction: ",
      format(mean(tb$bias_corrected - tb$eff, na.rm = TRUE), digits = 4),
      "\n", sep = "")
  cat("usable replications: median ", stats::median(tb$n_valid), " of ", x$B,
      ", fewest ", min(tb$n_valid),
      if (anyNA(tb$ci_lower)) paste0("   (", sum(is.na(tb$ci_lower)),
                                     " DMU(s) have too few: NA)") else "",
      "\n", sep = "")
  cat("\nfirst rows\n")
  print(utils::head(data.frame(dmu = tb$dmu, eff = round(tb$eff, 4),
                               bias_corr = round(tb$bias_corrected, 4),
                               lower = round(tb$ci_lower, 4),
                               upper = round(tb$ci_upper, 4)), 10))
  invisible(x)
}

summary.dea_subsample <- function(object, ...) {
  print(object)
  if (!is.null(object$volatility)) {
    cat("\nthe m selection, as a path over the grid\n")
    print(data.frame(m = object$m_grid,
                     mean_lower = round(object$endpoint_path[, 1], 4),
                     mean_upper = round(object$endpoint_path[, 2], 4),
                     volatility = signif(object$volatility, 3),
                     eligible = object$m_grid %in% object$eligible,
                     chosen = object$m_grid == object$m, row.names = NULL))
    cat("\nA flat stretch is where the asymptotics have bitten; the rule takes\n",
        "the middle of the flattest one. Only candidates with a full window are\n",
        "eligible -- a clipped window has fewer neighbours and so a smaller\n",
        "standard deviation for no reason to do with stability. Inspect the path\n",
        "rather than trusting the choice.\n", sep = "")
  }
  invisible(object)
}
