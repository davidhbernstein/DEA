## ---------------------------------------------------------------------------
## dea_rts_test() -- the nonparametric test of returns to scale of Simar and
## Wilson (2002).
##
## WHAT IS MISSING WITHOUT IT.  dea_rts() classifies each DMU, and that is a
## comparison of three linear-programming optima with no standard error
## attached: a DMU whose crs and vrs scores differ by 0.001 is called "drs"
## with exactly the confidence of one whose scores differ by 0.4.  There is no
## statement anywhere in that output about the SAMPLE -- about whether the
## whole technology may be taken to be a cone.  That question has an answer,
## and the answer is a test.
##
## THE STATISTIC.  Under H0 the restricted technology is the true one, so the
## restricted and the unrestricted estimator are consistent for the same
## distance and their ratio tends to 1.  Under H1 the restricted technology is
## strictly larger than the true one, the restricted scores are strictly worse,
## and the ratio tends to something below 1.  So the test is one-sided and
## rejects for SMALL values.  Two aggregations of the per-DMU ratio are
## reported, both from Simar and Wilson (2002, section 4):
##
##   ratio of means    mean(s_H0) / mean(s_VRS)
##   mean of ratios    mean(s_H0 / s_VRS)
##
## They are not the same number and need not even move together -- the first
## weights a DMU by its own level, the second does not.
##
## MEASURED, THE MEAN OF RATIOS IS THE BETTER ONE, and by a long way.  On
## dea_sim(returns = 1), where H0 is exactly true, the two have size 0.050 and
## 0.025 at n = 60 against a nominal 0.05; at returns = 0.6 they have power 0.780
## and 0.155, a paired difference of -0.625 with t = -18.2.  Both are reported
## because Simar and Wilson give both and because that is one design, but nothing
## measured recommends the ratio of means.  horserace/README.md has the table,
## including the two cells where this package does not come out well: the size is
## conservative rather than nominal and becomes more so with n, and the power
## against a 20% departure from constant returns is 2.0%, which is BELOW the
## size.  rDEA sits at 3.5% on the same datasets, so that is the 2002 procedure
## and not this implementation of it.
##
## THE ORIENTATION, AND WHY `s` IS WRITTEN RATHER THAN `theta`.  `s` is the
## score on the package's common (0, 1] scale: theta input-oriented, 1/phi
## output-oriented, which is what efficiency(fit, "score") already returns.
## Writing it that way makes the statistic lie in (0, 1] and the test
## left-tailed in BOTH orientations, and makes `mean of ratios` the mean of
## dea_rts()'s own `scale_eff` column exactly.
##
## It also decides something that the two reference packages decide
## differently, and the difference is real rather than cosmetic.  Output
## oriented, `npsf::nptestrts` aggregates on the phi scale and this package and
## `rDEA::rts.test` aggregate on the 1/phi scale.  Measured on the output
## reference design (tools/make_reference_values.R, n = 40):
##
##                      this package / rDEA      npsf        1 / npsf
##   ratio of means           0.858356          1.221348     0.818767
##   mean of ratios           0.862094          1.244737     0.803382
##
## The reciprocal column is NOT the first one.  Taking a mean does not commute
## with taking a reciprocal, so neither aggregate transfers between the two
## scales -- and this holds for the mean of ratios as much as for the ratio of
## means, even though the PER-DMU ratios are exact reciprocals of each other.
## Both are legitimate statistics for the same hypothesis; this package follows
## rDEA, because (0, 1] is the scale the rest of the package already works on.
##
## INPUT oriented the question does not arise and all three agree exactly: on
## the input reference design, 0.856958226232 and 0.849487582832 from this
## package, from `rDEA::rts.test`'s `w_hat` and `w48_hat + 1`, and from
## `npsf::nptestrts`'s `sefficiencyMean` and `mean(sefficiency)`, to 1e-13.  The
## NIRS versions agree with `nrsOVERvrsMean` and `mean(nrsOVERvrs)` to the same
## tolerance.
##
## THE NULL DISTRIBUTION.  Simulate from a DGP that satisfies H0.  The smoothed
## homogeneous bootstrap of dea_boot() does exactly that if it is driven by the
## RESTRICTED scores rather than the vrs ones: resampling and smoothing
## {s_H0, 2 - s_H0} and moving each DMU so that its distance to the ESTIMATED
## RESTRICTED frontier is the resampled value produces a pseudo-sample whose
## own frontier satisfies H0 by construction.  Both technologies are then fitted
## to that pseudo-sample and the statistic recomputed; the p-value is the
## fraction of replications at or below the observed value.
##
## Two details of that are worth stating because they are choices:
##
## * the TECHNOLOGY is the pseudo-sample and the points EVALUATED against it are
##   the original observations, as in dea_boot() -- and as in both reference
##   packages, which is one convention this corner of the literature does agree
##   on.  Note what this buys here: in each orientation the untouched block of
##   data guarantees feasibility.  Input-oriented, outputs are unchanged, so
##   lambda = e_o satisfies every output row; output-oriented, inputs are.  So
##   the test never has to decide what an infeasible replication means.
##
## * one pseudo-sample per REPLICATION, shared by all n evaluated DMUs.  That
##   is what makes the replication a draw of the statistic: the statistic is an
##   average over DMUs, and averaging over n INDEPENDENT pseudo-samples instead
##   would shrink its spread by roughly sqrt(n) and leave a null distribution
##   far too narrow.  `npsf::nptestrts` draws a fresh pseudo-sample inside the
##   DMU loop; `rDEA::rts.test` draws one per replication, as here.  The
##   consequence is measurable and is measured in horserace/README.md.
##
## WHAT THE TEST CANNOT BE ASKED, AND THE GATE FOR IT.  Measured here
## (/tmp probes reproduced in horserace/rts_test_experiment.R): on
## dea_sim(returns = 1), where the technology IS a cone and H0 is therefore
## exactly true, the INPUT-oriented statistic does not converge to 1.  It sits
## at 0.887, 0.903, 0.894, 0.898 for n = 80, 200, 800, 3200 and stays there,
## while the output-oriented one goes 0.970, 0.985, 0.991, 0.997 as it should.
## So the input-oriented test rejects a true null with probability tending to
## one, and no null distribution can repair that -- the statistic itself is not
## consistent for 1 under H0.
##
## The cause is provable rather than guessed, and it is the `x_range` problem
## that convergence/README.md documents for the convergence rate, biting a test
## instead.  The vrs program cannot put sum(lambda_i x_ij) below min_i x_ij,
## because a convex combination of the sample's j-th input is at least its
## smallest value.  So a DMU whose H0 projection needs theta*x_oj < min_i x_ij
## in ANY input j has theta_vrs > theta_crs for a reason that has nothing to do
## with returns to scale, and no sample size fixes it.  Splitting the same
## designs on exactly that condition separates the two behaviours cleanly: the
## DMUs whose projection is reachable give 0.9815, 0.9914, 0.9991, 0.9995 at
## p = 1, and the 34% that is not reachable gives 0.755, 0.751, 0.771, 0.776.
## One consistent test statistic and one constant, mixed.
##
## So `reachable` is computed and reported per DMU -- the output-oriented
## mirror is phi*y_oj > max_i y_ij -- and a warning fires past a quarter of the
## sample.  The benign case is not "none flagged" but "a share that SHRINKS
## with n": output-oriented on the same cone it runs 0.188, 0.020, 0.050 while
## input-oriented holds at 0.362, 0.340, 0.338.  Neither `npsf::nptestrts` nor
## `rDEA::rts.test` reports anything of the kind, and both will happily run the
## input-oriented test on such a design.
##
## THE P-VALUE IS (count + 1)/(B + 1), not count/B.  The observed statistic is
## itself a draw under H0, so including it is the standard Monte Carlo
## correction, and it means the reported p-value is never exactly zero -- which
## is the honest thing for a quantity whose resolution is 1/(B + 1).
## `npsf::nptestrts` reports count/B and does print 0.0000 from 100 draws.
##
## THE PER-DMU TESTS COME FREE.  Each replication already yields the whole
## vector of pseudo-ratios, so a per-DMU p-value is one more comparison per
## row.  These are n simultaneous tests, so the table also carries the Sidak
## level 1 - (1 - alpha)^(1/n) at which they can be read together, rather than
## leaving the reader to apply a correction to a column of raw p-values.
## ---------------------------------------------------------------------------

dea_rts_test <- function(x, y, data = NULL, h0 = c("crs", "nirs", "ndrs"),
                         orientation = c("in", "out"), B = 2000, alpha = 0.05,
                         bw = c("silverman", "nrd0", "ucv", "sj"),
                         seed = NULL, ncores = 1L, progress = interactive(),
                         scaling = TRUE) {

  call <- match.call()
  .dea_check_data_arg(data)
  h0 <- .match_arg_ci(h0, c("crs", "nirs", "ndrs"), "h0")
  orientation <- .match_arg_ci(orientation, c("in", "out"), "orientation")
  if (!is.numeric(B) || length(B) != 1L || B < 2) {
    stop("`B` must be a single number of bootstrap replications, at least 2.",
         call. = FALSE)
  }
  B <- as.integer(B)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be in (0, 1).", call. = FALSE)
  }
  if (B < 200) {
    warning("B = ", B, " gives a p-value resolution of 1/(B+1) = ",
            format(1 / (B + 1), digits = 3), ", so the smallest value this ",
            "test can report is ", format(1 / (B + 1), digits = 3),
            ". Simar and Wilson use B = 1000 or more for a test; the ",
            "B = 2000 default is theirs.", call. = FALSE)
  }

  t0 <- proc.time()[["elapsed"]]
  X <- .dea_matrix(x, data, "x")
  Y <- .dea_matrix(y, data, "y")
  .dea_check(X, Y)
  dmu <- .dea_dmu_names(X, Y)
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  sc <- .dea_scale(X, Y, scaling)
  Xs <- sc$X; Ys <- sc$Y

  ## Both technologies, over the same reference set. peers and slacks are not
  ## wanted B times over, so they are off.
  fit <- function(XR, YR, rts) {
    e <- .dea_radial(Xs, Ys, XR, YR, rts, orientation, super = FALSE,
                     slack = FALSE, peers = FALSE, n, nrow(XR), p, q)$eff
    e[is.finite(e) & abs(e - 1) < .DEA_CONSTANTS$TOL_EFF] <- 1
    e
  }

  e0 <- fit(Xs, Ys, h0)
  ev <- fit(Xs, Ys, "vrs")
  if (any(!is.finite(e0)) || any(!is.finite(ev))) {
    stop(sum(!is.finite(e0) | !is.finite(ev)), " DMU(s) have a non-finite ",
         "score under ", h0, " or vrs, so the statistic is not defined. A ",
         "self-referenced radial fit should not produce these; inspect ",
         "dea(x, y, rts = \"", h0, "\")$status.", call. = FALSE)
  }

  stat <- .rtst_stat(e0, ev, orientation)
  ratio_obs <- .rtst_ratio(e0, ev, orientation)

  ## The support gate. See the long note above: this is not a heuristic, it is
  ## the exact condition under which the vrs program cannot follow the H0
  ## projection at all.
  reach <- .rtst_reachable(e0, Xs, Ys, orientation)
  if (mean(!reach) > 0.25) {
    warning(sum(!reach), " of ", n, " DMUs (",
            round(100 * mean(!reach)), "%) have an H0 projection that the vrs ",
            "program provably cannot reach: it would need ",
            if (orientation == "in") "an input below the smallest observed value"
            else "an output above the largest observed value",
            " in at least one dimension, which no convex combination of the ",
            "sample supplies. For those DMUs the two scores differ for a ",
            "reason unrelated to returns to scale, and the statistic is a ",
            "mixture of a consistent comparison and a constant. Check whether ",
            "the share shrinks with n -- if it does not, this orientation ",
            "cannot test this hypothesis on this support. See ?dea_rts_test.",
            call. = FALSE)
  }

  ## Reflect the RESTRICTED scores about 1: they are the consistent estimator
  ## under H0, and the pseudo-sample has to come from an H0 world.
  refl <- c(e0, 2 - e0)
  h <- .dea_bandwidth(bw, refl, n)
  shrink <- sqrt(1 + h^2 / stats::sd(refl)^2)

  if (!is.null(seed)) {
    old <- .rng_snapshot(); on.exit(.rng_restore(old), add = TRUE)
    set.seed(seed)
  }

  one_rep <- function(b) {
    beta <- sample(refl, n, replace = TRUE)
    tw   <- mean(beta) + (beta + h * stats::rnorm(n) - mean(beta)) / shrink
    st <- if (orientation == "in") ifelse(tw > 1, 2 - tw, tw) else ifelse(tw < 1, 2 - tw, tw)
    st <- if (orientation == "in") pmin(pmax(st, .DEA_CONSTANTS$MIN_POS), 1) else pmax(st, 1)
    if (orientation == "in") {
      Xb <- Xs * (e0 / st); Yb <- Ys
    } else {
      Xb <- Xs; Yb <- Ys * (e0 / st)
    }
    a <- fit(Xb, Yb, h0)
    v <- fit(Xb, Yb, "vrs")
    c(.rtst_stat(a, v, orientation), .rtst_ratio(a, v, orientation))
  }

  reps <- .dea_replicate(one_rep, B, ncores, progress)
  M <- do.call(cbind, reps)                     ## (2 + n) x B
  boot_stat  <- t(M[1:2, , drop = FALSE])       ## B x 2
  colnames(boot_stat) <- c("ratio_of_means", "mean_of_ratios")
  boot_ratio <- M[-(1:2), , drop = FALSE]       ## n x B
  rownames(boot_ratio) <- dmu

  ## One-sided, lower tail, with the observed value counted as a draw.
  pv <- function(bv, obs) (sum(bv <= obs, na.rm = TRUE) + 1) / (sum(is.finite(bv)) + 1)
  p_global <- c(ratio_of_means = pv(boot_stat[, 1], stat[["ratio_of_means"]]),
                mean_of_ratios = pv(boot_stat[, 2], stat[["mean_of_ratios"]]))

  p_dmu <- vapply(seq_len(n), function(i) pv(boot_ratio[i, ], ratio_obs[i]),
                  numeric(1))
  alpha_ind <- 1 - (1 - alpha)^(1 / n)

  tab <- data.frame(dmu = dmu,
                    h0 = as.numeric(e0), vrs = as.numeric(ev),
                    scale_eff = as.numeric(ratio_obs),
                    p_value = p_dmu,
                    scale_efficient = p_dmu > alpha_ind,
                    reachable = reach,
                    row.names = NULL, stringsAsFactors = FALSE)
  names(tab)[2] <- h0

  structure(list(
    statistic = stat, p_value = p_global, reject = p_global < alpha,
    table = tab, boot = boot_stat, boot_ratio = boot_ratio,
    h0 = h0, alternative = "vrs", orientation = orientation,
    B = B, alpha = alpha, alpha_individual = alpha_ind,
    bw = h, bw_rule = if (is.numeric(bw)) "user" else bw[1],
    eff_h0 = e0, eff_vrs = ev, reachable = reach,
    n = n, p = p, q = q, dmu = dmu,
    total_time = proc.time()[["elapsed"]] - t0, call = call
  ), class = "dea_rts_test")
}

## The common (0, 1] score: theta input-oriented, 1/phi output-oriented. See
## the note at the top on why the statistic is written on this scale.
.rtst_score <- function(e, orientation) if (identical(orientation, "out")) 1 / e else e

.rtst_stat <- function(e0, ev, orientation) {
  s0 <- .rtst_score(e0, orientation); sv <- .rtst_score(ev, orientation)
  c(ratio_of_means = mean(s0, na.rm = TRUE) / mean(sv, na.rm = TRUE),
    mean_of_ratios = mean(s0 / sv, na.rm = TRUE))
}

## The per-DMU ratio. Identical under either scale convention -- the phi's
## cancel -- which is why it is the one quantity every package agrees on.
.rtst_ratio <- function(e0, ev, orientation) {
  if (identical(orientation, "out")) ev / e0 else e0 / ev
}

## Can the vrs program reach the H0 projection at all?  Input-oriented it must
## put some input below min_i x_ij, output-oriented some output above
## max_i y_ij, and a convex combination can do neither.  Exact, not a rule of
## thumb -- which is what makes it worth reporting.
.rtst_reachable <- function(e0, Xs, Ys, orientation) {
  tol <- .DEA_CONSTANTS$TOL_SLACK
  if (identical(orientation, "in")) {
    !apply(sweep(Xs * e0, 2L, apply(Xs, 2L, min) - tol, "<"), 1L, any)
  } else {
    !apply(sweep(Ys * e0, 2L, apply(Ys, 2L, max) + tol, ">"), 1L, any)
  }
}

print.dea_rts_test <- function(x, ...) {
  cat("--- Simar-Wilson test of returns to scale ---\n")
  cat("H0: the technology is ", toupper(x$h0), "    H1: ", toupper(x$alternative),
      "    (", if (x$orientation == "in") "input" else "output",
      " orientation)\n", sep = "")
  cat("B = ", x$B, "   bandwidth = ", format(x$bw, digits = 4), " (", x$bw_rule,
      ")   n = ", x$n, "   ", format(round(x$total_time, 2)), " sec\n", sep = "")
  cat("\nstatistic lies in (0, 1]; the test rejects for SMALL values\n")
  st <- data.frame(statistic = round(as.numeric(x$statistic), 6),
                   p_value = round(as.numeric(x$p_value), 6),
                   reject = as.logical(x$reject),
                   row.names = names(x$statistic))
  print(st)
  cat("\nreject at alpha = ", x$alpha, ";  p-value resolution 1/(B+1) = ",
      format(1 / (x$B + 1), digits = 3), "\n", sep = "")
  ns <- sum(!x$table$scale_efficient)
  cat("\nindividually scale-INEFFICIENT (simultaneous level ",
      format(x$alpha_individual, digits = 3), "): ", ns, " of ", x$n, "\n", sep = "")
  cat("H0 projection unreachable inside the sample's support: ",
      sum(!x$reachable), " of ", x$n, " (see ?dea_rts_test)\n", sep = "")
  invisible(x)
}

summary.dea_rts_test <- function(object, ...) {
  print(object)
  cat("\nper-DMU scale efficiency\n")
  print(round(summary(object$table$scale_eff), 4))
  cat("\nfirst rows\n")
  print(utils::head(object$table, 10))
  invisible(object)
}

nobs.dea_rts_test <- function(object, ...) object$n
