## ---------------------------------------------------------------------------
## dea_ladder() -- the efficiency ladder: how much of one DMU's score is owed
## to which peer.
##
## The score of an inefficient DMU is a statement about the handful of
## efficient units that envelop it, and nothing else. Chapter 21 section 6 names
## sensitivity to a single observation -- an outlier, or a misplaced decimal
## point -- as one of DEA's two headline caveats, and there is no residual and
## no standard error in the estimator to warn anyone. The ladder answers the
## practical form of that caveat: drop the peer the score leans on hardest,
## re-score against what is left, and repeat. A DMU whose score barely moves is
## enveloped by a crowd; one that jumps to the frontier on the first rung was
## being measured against a single unit, and the analyst should want to know
## which.
##
## IT ONLY EVER MOVES ONE WAY, AND THAT IS THE CHECK, NOT A COINCIDENCE.
## Removing reference DMUs can only shrink the estimated technology, so the
## evaluated point can only become less enveloped: theta weakly INCREASES
## toward 1 under an input orientation and phi weakly DECREASES toward 1 under
## an output one. A ladder that moved the other way would mean the refit was
## not nested, which is the mistake this file is most able to make -- so
## test-ladder.R asserts monotonicity on every rung rather than spot-checking
## the endpoints.
##
## THE EVALUATED DMU STAYS IN THE REFERENCE SET. Dropping it would make the
## program unbounded in the usual case and meaningless in the rest; keeping it
## is also what makes the ladder terminate, since once its own peers are gone
## it references itself alone and scores exactly 1. The ladder therefore ends
## on the frontier by construction, and the interesting number is how many
## rungs that took.
##
## WHICH PEER GOES FIRST IS DEFINED, NOT LEFT TO THE SOLVER. The rung drops
## the peer with the largest lambda, ties broken by reference-set order. That
## matters because peer sets are not unique when the program has multiple
## optima -- which is common -- so an undefined rule would give a different
## ladder on a different machine. What the ladder is sensitive to is the peer
## set the fit actually returned; `peers(fit)` shows it, and the caveat printed
## by `summary()` applies here with full force.
##
## WHY IT TAKES A FIT RATHER THAN DATA. Every setting that defines the
## technology -- rts, orientation, scaling -- has to be identical to the one
## that produced the score being audited, and re-specifying them by hand is an
## invitation to audit a different model than the one fitted. The radial model
## only: the non-radial fits do not carry the arguments (a direction, a measure)
## that a refit would need, and guessing them is the same mistake one level
## down.
## ---------------------------------------------------------------------------

dea_ladder <- function(object, dmu, max_steps = NULL, threshold = 1e-8) {
  call <- match.call()
  if (!inherits(object, "dea")) {
    stop("`object` must be a fit from dea(); got an object of class \"",
         class(object)[1L], "\".", call. = FALSE)
  }
  if (!identical(object$model, "radial")) {
    stop("dea_ladder() supports the radial model only; this fit is \"",
         object$model, "\". The ladder itself is not radial-specific, but a ",
         "refit needs the arguments that define the model -- a direction, a ",
         "measure -- and the fit does not carry them.", call. = FALSE)
  }
  if (isTRUE(object$super)) {
    stop("This fit is a super-efficiency fit, which already excludes each DMU ",
         "from its own reference set. Refit with super = FALSE to audit the ",
         "ordinary score.", call. = FALSE)
  }
  if (!isTRUE(object$self_ref)) {
    stop("This fit was scored against an external reference set, so its peers ",
         "are not DMUs of the evaluated sample and dropping them is a ",
         "different question. Refit without `xref`/`yref`.", call. = FALSE)
  }
  if (is.null(object$lambda)) {
    stop("This fit was made with peers = FALSE, so the lambda matrix was not ",
         "kept and there is no peer to drop. Refit dea() with peers = TRUE.",
         call. = FALSE)
  }
  o <- .ladder_which(dmu, object$dmu)
  if (!is.finite(object$eff[o])) {
    stop("DMU \"", object$dmu[o], "\" was not solved in this fit (status ",
         object$status[o], "), so there is no score to audit.", call. = FALSE)
  }

  X <- object$x; Y <- object$y
  nref <- object$nref
  max_steps <- if (is.null(max_steps)) nref else {
    if (length(max_steps) != 1L || !is.finite(max_steps) || max_steps < 1)
      stop("`max_steps` must be a single positive number.", call. = FALSE)
    as.integer(max_steps)
  }

  ## Rung 0 is the fit as it stands, so the table reads as a ladder from the
  ## reported score rather than from a refit that might differ in the last
  ## digit.
  keep <- rep(TRUE, nref)
  lam  <- object$lambda[o, ]
  steps <- list(.ladder_row(0L, NA_character_, NA_real_, nref,
                            object$eff[o], object$status[o]))

  for (k in seq_len(max_steps)) {
    ## Peers other than the DMU itself. Its own column is what keeps the
    ## program feasible and is never a candidate.
    cand <- which(keep & lam > threshold)
    cand <- cand[cand != o]
    if (!length(cand)) break
    ## Largest lambda first; order() is stable, so a tie goes to the lower
    ## reference index and the ladder is reproducible.
    drop_at <- cand[order(-lam[cand])][1L]
    keep[drop_at] <- FALSE

    ## The weight recorded is the one the peer carried in the reference set it
    ## was dropped FROM, not in the original fit. Only rung 1 is dropping from
    ## the original fit; after that the peer set has been re-solved, and most
    ## units dropped later carried no weight at all to begin with.
    dropped_lambda <- lam[drop_at]

    re <- .ladder_refit(object, X, Y, o, keep)
    lam <- rep(0, nref)
    lam[keep] <- re$lambda
    steps[[length(steps) + 1L]] <-
      .ladder_row(k, object$ref[drop_at], dropped_lambda,
                  sum(keep), re$eff, re$status)
  }

  tab <- do.call(rbind, steps)
  rownames(tab) <- NULL
  bad <- which(!is.finite(tab$eff))
  if (length(bad)) {
    warning("The refit was not solved on ", length(bad), " rung(s) of the ",
            "ladder; their scores are NA. See the `status` column for the ",
            "lpSolve return codes.", call. = FALSE)
  }

  structure(list(
    table = tab, dmu = object$dmu[o], index = o,
    eff0 = object$eff[o], eff1 = tab$eff[nrow(tab)],
    n_step = nrow(tab) - 1L,
    exhausted = nrow(tab) - 1L < max_steps,
    rts = object$rts, orientation = object$orientation,
    threshold = threshold, fit_call = object$call, call = call
  ), class = "dea_ladder")
}

## One row of the ladder. Factored out so every rung is built the same way and
## rung 0 cannot drift from the rest.
.ladder_row <- function(step, dropped, lambda, n_ref, eff, status) {
  data.frame(step = as.integer(step), dropped = dropped, lambda = lambda,
             n_ref = as.integer(n_ref), eff = eff, status = as.integer(status),
             stringsAsFactors = FALSE)
}

## Re-score the one DMU against the surviving reference set.
##
## This is the external-reference form -- one program per rung instead of one
## per surviving DMU -- which is what `xref`/`yref` is for. The cost is that
## .lp_final() does not treat status 2 as retryable when the evaluated and
## reference matrices differ, and here they always do. But the evaluated DMU is
## itself in the reference set, so lambda = e_o is feasible and an
## "infeasible" is necessarily false. Rather than accept the NA, the rung is
## re-run self-referenced, where the guard applies; that path solves one
## program per surviving DMU and is only ever reached on a rung that already
## failed.
.ladder_refit <- function(object, X, Y, o, keep) {
  f <- suppressWarnings(dea(X[o, , drop = FALSE], Y[o, , drop = FALSE],
                            rts = object$rts, orientation = object$orientation,
                            slack = FALSE, peers = TRUE,
                            scaling = object$scaling,
                            xref = X[keep, , drop = FALSE],
                            yref = Y[keep, , drop = FALSE]))
  if (!is.finite(f$eff[1L]) && identical(as.integer(f$status[1L]), 2L)) {
    g <- suppressWarnings(dea(X[keep, , drop = FALSE], Y[keep, , drop = FALSE],
                              rts = object$rts,
                              orientation = object$orientation,
                              slack = FALSE, peers = TRUE,
                              scaling = object$scaling))
    i <- match(o, which(keep))
    return(list(eff = g$eff[i], status = g$status[i], lambda = g$lambda[i, ]))
  }
  list(eff = f$eff[1L], status = f$status[1L], lambda = f$lambda[1L, ])
}

## A DMU by name or by position, with the error naming what was asked for.
.ladder_which <- function(dmu, names) {
  if (length(dmu) != 1L)
    stop("`dmu` must name a single DMU; got ", length(dmu), ".", call. = FALSE)
  if (is.character(dmu)) {
    o <- match(dmu, names)
    if (is.na(o))
      stop("`dmu`: \"", dmu, "\" is not one of this fit's DMUs.", call. = FALSE)
    return(o)
  }
  o <- as.integer(dmu)
  if (!is.finite(o) || o < 1L || o > length(names))
    stop("`dmu`: index ", dmu, " is outside 1:", length(names), ".",
         call. = FALSE)
  o
}

print.dea_ladder <- function(x, ...) {
  cat("--- Efficiency ladder ---\n")
  cat("DMU: ", x$dmu, "   technology: ", toupper(x$rts), ", ",
      if (x$orientation == "in") "input" else "output", " oriented\n", sep = "")
  cat("score ", format(round(x$eff0, 4)), " -> ", format(round(x$eff1, 4)),
      " over ", x$n_step, " rung(s)\n", sep = "")
  tb <- x$table
  tb$eff <- round(tb$eff, 4)
  tb$lambda <- round(tb$lambda, 4)
  ## The ladder runs until the DMU references itself alone, so it is as long as
  ## the sample. The informative part is the top of it; $table has the rest.
  print(utils::head(tb, 10))
  if (nrow(tb) > 10L)
    cat("... ", nrow(tb) - 10L, " further rung(s) in $table\n", sep = "")
  cat("\nEach rung drops the remaining peer with the largest lambda and ",
      "re-scores\nagainst what is left. The DMU itself always stays in the ",
      "reference set, so\nthe ladder ends at 1. Peer sets are not unique: ",
      "this is the sensitivity of\nthe score to the peers THIS fit returned.\n",
      sep = "")
  invisible(x)
}

summary.dea_ladder <- function(object, ...) {
  print(object)
  d <- diff(object$table$eff)
  if (length(d) && any(is.finite(d))) {
    i <- which.max(abs(d))
    cat("\nlargest single move: rung ", object$table$step[i + 1L],
        ", dropping ", object$table$dropped[i + 1L], ", ",
        format(round(d[i], 4), nsmall = 4), "\n", sep = "")
    cat("share of the total move in that one rung: ",
        format(round(abs(d[i]) / abs(object$eff1 - object$eff0), 3)), "\n",
        sep = "")
  }
  invisible(object)
}

## The ladder as a step plot. The horizontal line is the frontier, which every
## ladder reaches; what the picture is for is the SHAPE of the climb.
plot.dea_ladder <- function(x, ...) {
  tb <- x$table
  plot(tb$step, tb$eff, type = "b", pch = 16,
       xlab = "peers dropped", ylab = "efficiency",
       main = paste0("Efficiency ladder, DMU ", x$dmu), ...)
  abline(h = 1, lty = 2)
  invisible(tb)
}

nobs.dea_ladder <- function(object, ...) nrow(object$table)
