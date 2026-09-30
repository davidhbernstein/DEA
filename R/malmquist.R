## ---------------------------------------------------------------------------
## dea_malmquist() -- the Malmquist productivity index of Fare, Grosskopf,
## Norris and Zhang (1994), on the panel layer in R/panel.R.
##
## Four distance functions per DMU per consecutive period pair, and the whole
## thing is already available: each one is a dea() fit with an `xref`/`yref`
## pointing at whichever period's technology is wanted. The two MIXED terms --
## one period's DMUs scored against the other period's frontier -- are the
## reason `xref`/`yref` exists at all.
##
## THE NAMING CONVENTION IS TWO LETTERS, OBSERVATION THEN TECHNOLOGY. `d_tf` is
## the TO-period observation scored against the FROM-period technology. Stating
## it that way is not fussiness: `Benchmarking::malmq()` returns the same four
## numbers as e00/e10/e11/e01 and the reader has to work out which index is
## which, and getting it backwards silently swaps efficiency change with its
## reciprocal.
##
##   efficiency change  EC = d_tt / d_ff
##   technical change   TC = sqrt( (d_tf / d_tt) * (d_ff / d_ft) )
##   Malmquist          M  = EC * TC = sqrt( (d_tf / d_ff) * (d_tt / d_ft) )
##
## with M > 1 productivity GROWTH. Every ratio above is on a common (0, 1]
## efficiency scale, so the output-oriented case is handled by converting phi
## to 1/phi once, at the source, rather than by carrying a second set of
## formulas. The mixed terms routinely leave (0, 1] -- an observation scored
## against another period's frontier need not be inside it -- and that is the
## measurement, not an error.
##
## WHY THE DEFAULT IS CONSTANT RETURNS, and it is not laziness. The index is
## defined on a cone by Fare et al., for a reason that shows up immediately in
## practice: the mixed distance functions can be INFEASIBLE under variable
## returns, because a convex hull cannot extrapolate to a point outside it,
## while a cone can always be scaled to reach one. Software that quietly
## defaults to `vrs` produces a table with holes in it and rarely says so.
## Here `vrs` is allowed, the infeasible terms are NA, and the count is
## reported per period pair.
##
## THE SCALE DECOMPOSITION IS AN IDENTITY, AND IS TESTED AS MORE THAN ONE.
## Fare et al. split efficiency change into a pure part and a scale part,
##
##   PEC = v_tt / v_ff                      (vrs own-period scores)
##   SEC = (d_tt / v_tt) / (d_ff / v_ff)    (scale efficiency, then and now)
##
## and EC = PEC * SEC follows by cancellation. Checking that product is
## therefore worth nothing on its own -- it cannot fail. What the test checks
## is that PEC and SEC equal the same quantities computed from independent
## vrs and crs fits.
## ---------------------------------------------------------------------------

dea_malmquist <- function(x, y = NULL, id = NULL, period = NULL, data = NULL,
                          rts = c("crs", "vrs", "nirs", "ndrs"),
                          orientation = c("in", "out"),
                          scale_decomposition = TRUE,
                          scaling = TRUE) {
  call <- match.call()
  rts <- .match_arg_ci(rts, c("crs", "vrs", "nirs", "ndrs"), "rts")
  orientation <- .match_arg_ci(orientation, c("in", "out"), "orientation")

  pn <- if (inherits(x, "dea_panel")) x else
        dea_panel(x, y, id, period, data = data)

  if (!identical(rts, "crs") && isTRUE(scale_decomposition)) {
    ## Fare et al.'s pure/scale split is defined for the CRS-BASED index and
    ## only for it: the scale part is what the crs frontier says and the vrs
    ## frontier does not. Under vrs it is identically 1, a column of ones
    ## dressed as a result; under nirs or ndrs it would be some other
    ## decomposition of some other index, and calling it theirs would be wrong.
    if (!missing(scale_decomposition))
      warning("The pure/scale split is defined for the crs-based index. Under ",
              "rts = \"", rts, "\" there is no such decomposition to report: ",
              if (identical(rts, "vrs"))
                "the scale part is the gap between the crs and vrs frontiers, and this index already sits on the vrs one, so it would be 1 for every DMU."
              else
                "splitting a nirs or ndrs index against the vrs frontier gives some other decomposition of some other index, not Fare et al.'s.",
              " Ignoring scale_decomposition.", call. = FALSE)
    scale_decomposition <- FALSE
  }

  fit <- function(X, Y, XR, YR) {
    f <- suppressWarnings(dea(X, Y, rts = rts, orientation = orientation,
                              slack = FALSE, peers = FALSE, scaling = scaling,
                              xref = XR, yref = YR))
    e <- as.numeric(f$eff)
    ## One conversion, at the source: every formula below is written on the
    ## (0, 1] input scale.
    if (identical(orientation, "out")) e <- 1 / e
    list(e = e, status = f$status)
  }

  per <- pn$periods
  rows <- list(); dist <- list(); infeas <- integer(0)

  for (k in seq_len(length(per) - 1L)) {
    t0 <- per[k]; t1 <- per[k + 1L]
    pr <- .panel_pair(pn, t0, t1)
    if (!length(pr$id)) {
      warning("No DMU is observed in both ", t0, " and ", t1,
              "; that pair has no index.", call. = FALSE)
      next
    }
    A <- pr$a; B <- pr$b

    ff <- fit(A$X, A$Y, A$X, A$Y)   # from-period obs, from-period technology
    tt <- fit(B$X, B$Y, B$X, B$Y)   # to   on to
    tf <- fit(B$X, B$Y, A$X, A$Y)   # to   on from   -- mixed
    ft <- fit(A$X, A$Y, B$X, B$Y)   # from on to     -- mixed

    d_ff <- ff$e[pr$ia]; d_tt <- tt$e[pr$ib]
    d_tf <- tf$e[pr$ib]; d_ft <- ft$e[pr$ia]

    ec <- d_tt / d_ff
    tc <- sqrt((d_tf / d_tt) * (d_ff / d_ft))
    m  <- ec * tc

    pec <- sec <- rep(NA_real_, length(pr$id))
    if (scale_decomposition) {
      vf <- .malm_fit_vrs(A$X, A$Y, orientation, scaling)[pr$ia]
      vt <- .malm_fit_vrs(B$X, B$Y, orientation, scaling)[pr$ib]
      pec <- vt / vf
      sec <- (d_tt / vt) / (d_ff / vf)
    }

    bad <- sum(!is.finite(d_tf)) + sum(!is.finite(d_ft))
    if (bad) infeas[paste(t0, "->", t1)] <- bad

    rows[[length(rows) + 1L]] <- data.frame(
      id = pr$id, from = t0, to = t1,
      malmquist = m, effch = ec, techch = tc,
      purech = pec, scalech = sec, stringsAsFactors = FALSE)
    dist[[length(dist) + 1L]] <- data.frame(
      id = pr$id, from = t0, to = t1,
      d_ff = d_ff, d_tt = d_tt, d_tf = d_tf, d_ft = d_ft,
      stringsAsFactors = FALSE)
  }

  if (!length(rows))
    stop("No period pair had a DMU observed in both periods.", call. = FALSE)

  tab <- do.call(rbind, rows)
  if (!scale_decomposition) { tab$purech <- NULL; tab$scalech <- NULL }
  rownames(tab) <- NULL
  dst <- do.call(rbind, dist); rownames(dst) <- NULL

  if (length(infeas)) {
    warning("The mixed distance functions were infeasible for ", sum(infeas),
            " DMU-period(s) and are reported as NA (",
            paste(sprintf("%s: %d", names(infeas), infeas), collapse = "; "),
            "). A convex hull cannot extrapolate to a point outside it, which ",
            "is why the index is usually computed under rts = \"crs\".",
            call. = FALSE)
  }

  structure(list(
    table = tab, distances = dst, panel = pn,
    rts = rts, orientation = orientation,
    scale_decomposition = scale_decomposition,
    infeasible = infeas,
    periods = per, n_pair = length(rows),
    call = call
  ), class = "dea_malmquist")
}

## The vrs own-period scores the scale split needs. Separate so that the four
## crs distances above stay readable, and on the same (0, 1] scale as they are.
.malm_fit_vrs <- function(X, Y, orientation, scaling) {
  f <- suppressWarnings(dea(X, Y, rts = "vrs", orientation = orientation,
                            slack = FALSE, peers = FALSE, scaling = scaling))
  e <- as.numeric(f$eff)
  if (identical(orientation, "out")) e <- 1 / e
  e
}

print.dea_malmquist <- function(x, ...) {
  cat("--- Malmquist productivity index (",
      if (x$orientation == "in") "input" else "output", " oriented, rts = ",
      x$rts, ") ---\n", sep = "")
  cat("DMUs: ", x$panel$n_id, "   periods: ", x$panel$n_period,
      "   period pairs: ", x$n_pair, "\n", sep = "")
  cat("period order: ", paste(x$periods, collapse = " < "), "\n", sep = "")
  cat("\ngeometric mean by period pair (> 1 is growth)\n")
  print(round(.malm_gm(x$table), 4))
  if (length(x$infeasible))
    cat("\ninfeasible mixed distances: ", sum(x$infeasible), "\n", sep = "")
  invisible(x)
}

summary.dea_malmquist <- function(object, ...) {
  print(object)
  cat("\ndistribution across all DMU-pairs\n")
  cols <- intersect(c("malmquist", "effch", "techch", "purech", "scalech"),
                    names(object$table))
  for (cl in cols) {
    v <- object$table[[cl]]; v <- v[is.finite(v)]
    cat(sprintf("  %-10s n %4d   gm %7.4f   min %7.4f   max %7.4f\n",
                cl, length(v), exp(mean(log(v))), min(v), max(v)))
  }
  invisible(object)
}

## The geometric mean is the right average for a ratio index: the mean of the
## logs is what makes "no change" sit at 1 and makes a doubling and a halving
## cancel. An arithmetic mean of ratios does neither.
.malm_gm <- function(tab) {
  cols <- intersect(c("malmquist", "effch", "techch", "purech", "scalech"),
                    names(tab))
  pair <- paste(tab$from, "->", tab$to)
  out <- lapply(cols, function(cl)
    tapply(tab[[cl]], factor(pair, levels = unique(pair)),
           function(v) { v <- v[is.finite(v) & v > 0]
                         if (length(v)) exp(mean(log(v))) else NA_real_ }))
  m <- do.call(cbind, out); colnames(m) <- cols
  m
}

nobs.dea_malmquist <- function(object, ...) nrow(object$table)
