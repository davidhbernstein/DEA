## ---------------------------------------------------------------------------
## dea_boot() -- the smoothed homogeneous bootstrap of Simar and Wilson (1998).
##
## WHY IT IS NEEDED.  The DEA frontier is spanned by the observed DMUs, so it
## lies inside the true one: every efficiency score is biased toward 1, and the
## bias is not small.  A DEA score reported on its own is a point estimate of
## unknown accuracy from an estimator with a known, one-signed bias.  This is
## the single largest gap between how DEA is used and what is known about it.
##
## WHAT THE BOOTSTRAP DOES.  Resample efficiencies, rebuild a pseudo-sample
## whose frontier is the ESTIMATED one, re-estimate against it, and read the
## sampling distribution of (theta_hat - theta) off the resulting spread.
##
## THE SMOOTHING AND THE REFLECTION.  A naive resample of the theta_hat values
## is inconsistent: it puts an atom of probability on each observed score, and
## in particular a large atom on 1, so the pseudo-frontier is a fixed set of
## points rather than a draw from a density.  Simar and Wilson smooth the
## resample with a Gaussian kernel, which fixes that but would then leak
## probability past the boundary at theta = 1 -- exactly where the mass is.
## The fix is reflection: smooth the 2n-point set {theta_hat, 2 - theta_hat},
## which is symmetric about 1 by construction, then fold draws back across 1.
## No mass crosses the boundary because the density is symmetric there.
##
## The variance correction that follows the smoothing matters for the same
## reason a kernel density is not a resample: adding h*epsilon inflates the
## variance by h^2, so the draws are rescaled to put it back.
##
## READING THE OUTPUT.  `bias_corrected` is 2*theta_hat - mean(theta_hat*).
## Simar and Wilson warn against using it unconditionally: the correction
## removes a bias but adds the variance of the estimate of that bias, and it is
## a net loss unless |bias| is large relative to the bootstrap standard error.
## Their rule of thumb -- correct only when |bias|/se > 1/sqrt(3) -- is
## evaluated here and reported per DMU as `correct_worthwhile`, rather than
## being applied silently.
## ---------------------------------------------------------------------------

dea_boot <- function(object, B = 2000, alpha = 0.05,
                     bw = c("silverman", "nrd0", "ucv", "sj"),
                     seed = NULL, ncores = 1L, progress = interactive()) {

  if (!inherits(object, "dea")) {
    stop("`object` must be a fit from dea(). The Simar-Wilson bootstrap is ",
         "defined for the radial efficiency estimator; dea_sbm() and ",
         "dea_ddf() have no equivalent published resampling scheme in this ",
         "package.", call. = FALSE)
  }
  if (!identical(object$model, "radial")) {
    stop("`object` is a ", object$model, " fit; dea_boot() needs a radial ",
         "fit from dea().", call. = FALSE)
  }
  if (isTRUE(object$super)) {
    stop("The bootstrap does not apply to a super-efficiency fit: the ",
         "Andersen-Petersen score is not an estimator of a distance to the ",
         "true frontier, so there is no bias to correct.", call. = FALSE)
  }
  if (identical(object$rts, "fdh")) {
    stop("The Simar-Wilson (1998) homogeneous bootstrap is built for the ",
         "convex DEA estimator. Under free disposal the relevant scheme is ",
         "the subsampling bootstrap of Jeong and Simar (2006), which is not ",
         "implemented here.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || B < 2) {
    stop("`B` must be a single number of bootstrap replications, at least 2.",
         call. = FALSE)
  }
  B <- as.integer(B)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be in (0, 1).", call. = FALSE)
  }

  t0 <- proc.time()[["elapsed"]]
  if (!isTRUE(object$self_ref)) {
    stop("The bootstrap resamples the technology, so it needs a fit whose ",
         "reference set IS the sample. This fit was scored against an ",
         "explicit `xref`/`yref`.", call. = FALSE)
  }
  X <- object$x; Y <- object$y
  n <- object$n; p <- object$p; q <- object$q
  ori <- object$orientation; rts <- object$rts
  eff <- object$eff
  if (any(!is.finite(eff))) {
    stop("The supplied fit has ", sum(!is.finite(eff)), " non-finite ",
         "efficiency score(s); the bootstrap has nothing to resample for ",
         "those DMUs.", call. = FALSE)
  }

  ## Reflect about 1. This holds in either orientation: theta lies in (0, 1]
  ## and phi in [1, Inf), and 1 is the boundary in both cases.
  refl <- c(eff, 2 - eff)
  h <- .dea_bandwidth(bw, refl, n)
  sd_refl <- stats::sd(refl)
  shrink <- sqrt(1 + h^2 / sd_refl^2)

  if (!is.null(seed)) {
    old <- .rng_snapshot(); on.exit(.rng_restore(old), add = TRUE)
    set.seed(seed)
  }

  ## Scaling is applied once, outside the replication loop: the pseudo-sample
  ## is a rescaling of the SAME columns, so the column means used here stay
  ## valid and the LP stays well conditioned for every replication.
  sc <- .dea_scale(X, Y, object$scaling)
  Xs <- sc$X; Ys <- sc$Y

  one_rep <- function(b) {
    beta <- sample(refl, n, replace = TRUE)
    tw   <- mean(beta) + (beta + h * stats::rnorm(n) - mean(beta)) / shrink
    ## Fold back across the boundary. The smoothed density is symmetric about
    ## 1, so folding is a measure-preserving map onto the correct half-line.
    st <- if (ori == "in") ifelse(tw > 1, 2 - tw, tw) else ifelse(tw < 1, 2 - tw, tw)
    st <- if (ori == "in") pmin(pmax(st, .DEA_CONSTANTS$MIN_POS), 1) else pmax(st, 1)

    ## The pseudo-sample: move each DMU along its own orientation so that its
    ## efficiency RELATIVE TO THE ESTIMATED FRONTIER is st instead of eff.
    if (ori == "in") {
      Xb <- Xs * (eff / st); Yb <- Ys
    } else {
      Xb <- Xs; Yb <- Ys * (eff / st)
    }
    Bl <- .lp_radial_build(Xb, Yb, rts, ori)
    ## Note which matrices go where: the TECHNOLOGY is the pseudo-sample, the
    ## points EVALUATED against it are the original observations. Bootstrapping
    ## the pseudo-points against themselves would just reproduce eff.
    ##
    ## The resulting scores routinely cross 1 -- theta* above it, phi* below --
    ## because the pseudo-technology is a shrunken copy of the estimated one and
    ## an original observation need not lie inside it. That is the whole point:
    ## the spread on the far side of 1 is what the interval is built from, and
    ## clamping these scores to the "valid" range would truncate exactly the
    ## tail being measured.
    vapply(seq_len(n), function(o) .lp_radial_at(Bl, Xs, Ys, o)$eff, numeric(1))
  }

  reps <- .dea_replicate(one_rep, B, ncores, progress)
  boot <- do.call(cbind, reps)            ## n x B
  rownames(boot) <- object$dmu

  mb   <- rowMeans(boot, na.rm = TRUE)
  bias <- mb - eff
  se   <- apply(boot, 1L, stats::sd, na.rm = TRUE)
  bc   <- eff - bias

  ## Percentile interval on (theta_hat* - theta_hat), inverted. Writing it as
  ## 2*theta_hat - quantile(theta_hat*) rather than quantile(theta_hat*)
  ## directly is what makes it an interval for theta and not for theta_hat.
  ql <- apply(boot, 1L, stats::quantile, probs = 1 - alpha/2, na.rm = TRUE)
  qu <- apply(boot, 1L, stats::quantile, probs = alpha/2,     na.rm = TRUE)
  ci_lo <- 2 * eff - ql
  ci_hi <- 2 * eff - qu

  tab <- data.frame(dmu = object$dmu, eff = as.numeric(eff),
                    bias = as.numeric(bias), se = as.numeric(se),
                    bias_corrected = as.numeric(bc),
                    ci_lower = as.numeric(ci_lo), ci_upper = as.numeric(ci_hi),
                    ## Simar and Wilson's rule: the correction pays only when
                    ## the bias it removes exceeds the noise it adds.
                    correct_worthwhile = is.finite(se) & se > 0 &
                                         abs(bias) / se > 1 / sqrt(3),
                    row.names = NULL, stringsAsFactors = FALSE)

  structure(list(
    table = tab, boot = boot, B = B, alpha = alpha, bw = h, bw_rule = bw[1],
    orientation = ori, rts = rts, n = n, p = p, q = q, dmu = object$dmu,
    fit = object, total_time = proc.time()[["elapsed"]] - t0,
    call = match.call()
  ), class = "dea_boot")
}

## Bandwidth for the reflected sample. "silverman" is the robust normal
## reference rule Simar and Wilson use, on 2n points; the others come straight
## from stats and are offered because the choice is genuinely open and the
## reflected density is bimodal enough that a normal reference oversmooths it.
.dea_bandwidth <- function(bw, refl, n) {
  if (is.numeric(bw)) {
    if (length(bw) != 1L || bw <= 0) stop("Numeric `bw` must be a single positive number.", call. = FALSE)
    return(bw)
  }
  bw <- .match_arg_ci(bw, c("silverman", "nrd0", "ucv", "sj"), "bw")
  h <- switch(bw,
    silverman = 0.9 * min(stats::sd(refl),
                          stats::IQR(refl) / 1.349) * (2 * n)^(-1/5),
    nrd0 = stats::bw.nrd0(refl),
    ucv  = suppressWarnings(stats::bw.ucv(refl)),
    sj   = suppressWarnings(stats::bw.SJ(refl)))
  if (!is.finite(h) || h <= 0) {
    ## Happens when every score is 1 -- a degenerate sample, but a real one.
    warning("Bandwidth rule \"", bw, "\" returned ", h, "; falling back to a ",
            "small positive value. This usually means the efficiency scores ",
            "have (almost) no spread.", call. = FALSE)
    h <- 0.01
  }
  h
}

## Replicate, on one core or several. parallel is in Suggests, so a request for
## more than one core degrades to a serial run with a message rather than an
## error if it is not installed -- and never on Windows, where mclapply is a
## serial stub anyway.
.dea_replicate <- function(fn, B, ncores, progress) {
  ncores <- as.integer(ncores)
  use_par <- ncores > 1L && .Platform$OS.type != "windows" &&
             requireNamespace("parallel", quietly = TRUE)
  if (ncores > 1L && !use_par) {
    message("ncores = ", ncores, " requested but parallel execution is not ",
            "available here; running serially.")
  }
  if (use_par) {
    return(parallel::mclapply(seq_len(B), fn, mc.cores = ncores))
  }
  out <- vector("list", B)
  for (b in seq_len(B)) {
    if (progress && (b %% 50 == 0 || b == 1L)) {
      cat(sprintf("\r  bootstrap %d/%d", b, B)); utils::flush.console()
    }
    out[[b]] <- fn(b)
  }
  if (progress) cat("\r", strrep(" ", 30), "\r", sep = "")
  out
}

print.dea_boot <- function(x, ...) {
  cat("--- Simar-Wilson bootstrap ---\n")
  cat("model: radial, ", x$rts, ", ",
      if (x$orientation == "in") "input" else "output", " orientation\n", sep = "")
  cat("B = ", x$B, "   bandwidth = ", format(x$bw, digits = 4),
      " (", x$bw_rule, ")   ", format(round(x$total_time, 2)), " sec\n", sep = "")
  tb <- x$table
  cat("\nmean bias: ", format(mean(tb$bias), digits = 4),
      "   mean se: ", format(mean(tb$se), digits = 4), "\n", sep = "")
  cat("bias correction worthwhile (|bias|/se > 1/sqrt(3)) for ",
      sum(tb$correct_worthwhile), " of ", x$n, " DMUs\n", sep = "")
  cat("\nfirst rows\n")
  print(utils::head(
    data.frame(dmu = tb$dmu, eff = round(tb$eff, 4),
               bias_corr = round(tb$bias_corrected, 4),
               lower = round(tb$ci_lower, 4), upper = round(tb$ci_upper, 4)), 10))
  invisible(x)
}

summary.dea_boot <- function(object, ...) {
  print(object)
  cat("\nefficiency, raw vs bias-corrected\n")
  print(rbind(raw = summary(object$table$eff),
              corrected = summary(object$table$bias_corrected)))
  invisible(object)
}
