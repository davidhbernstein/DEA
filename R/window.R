## ---------------------------------------------------------------------------
## dea_window() -- window analysis (chapter 3 section 5), on the panel layer in
## R/panel.R.
##
## A DMU observed in each of T periods is treated as T DISTINCT DMUs inside a
## rolling window of w consecutive periods, so each of the T - w + 1 windows
## estimates one technology from n_id * w points. That is the whole purpose:
## DEA discriminates badly when the sample is small relative to p + q -- the
## standard heuristic is n >= 3(p + q), which .dea_check() warns about -- and a
## single period of a short panel often fails it while a window of three
## periods does not.
##
## THE OUTPUT IS A SPREAD, NOT A SCORE. A DMU-period that sits in the interior
## of the panel appears in up to w windows and therefore gets up to w different
## scores, because each window's technology is estimated from a different set of
## points. Reading any one of them as "the" score throws away what the method
## produces. The stability of a DMU across its windows is the result, and
## summary() reports exactly that: per DMU, the mean score and the range of it
## across windows.
##
## WHAT POOLING ASSUMES, AND IT IS THE OPPOSITE OF MALMQUIST. Putting several
## periods into one technology assumes the frontier did not move within the
## window; a DMU that looks worse in a later window may simply be facing a
## frontier that advanced. Window analysis measures stability under a frontier
## held fixed, dea_malmquist() measures the movement of the frontier itself, and
## neither substitutes for the other. The documentation says this because the
## two are routinely confused.
##
## THE WIDTH IS AD HOC AND THE BOOK SAYS SO. There is no theory that picks w;
## chapter 3 section 5 advises running two or three widths and seeing whether
## the answer moves. `width` therefore takes a VECTOR, and the table carries a
## `width` column, so the advice costs one argument rather than a loop the user
## has to write and then reconcile by hand.
##
## PERIOD ORDER IS INHERITED FROM dea_panel(), DELIBERATELY. A window is a
## contiguous block of the panel's period order, so a wrong order does not
## error here -- it silently groups periods that are not neighbours. The panel
## object prints the order it settled on, and that is where to check it.
## ---------------------------------------------------------------------------

dea_window <- function(x, y = NULL, id = NULL, period = NULL, data = NULL,
                       width = 3L,
                       rts = c("vrs", "crs", "nirs", "ndrs", "fdh"),
                       orientation = c("in", "out"),
                       scaling = TRUE) {
  call <- match.call()
  rts <- .match_arg_ci(rts, c("vrs", "crs", "nirs", "ndrs", "fdh"), "rts")
  orientation <- .match_arg_ci(orientation, c("in", "out"), "orientation")

  pn <- if (inherits(x, "dea_panel")) x else
        dea_panel(x, y, id, period, data = data)

  width <- .window_widths(width, pn$n_period)

  rows <- list()
  for (w in width) {
    nwin <- pn$n_period - w + 1L
    for (k in seq_len(nwin)) {
      per <- pn$periods[k + seq_len(w) - 1L]
      i <- unlist(pn$rows[per], use.names = FALSE)
      ## Each DMU-period is its own DMU here, so the fit is an ordinary
      ## self-referenced one on the stacked window; nothing about it knows
      ## that some of its rows are the same unit at different dates.
      f <- suppressWarnings(dea(pn$X[i, , drop = FALSE],
                                pn$Y[i, , drop = FALSE],
                                rts = rts, orientation = orientation,
                                slack = FALSE, peers = FALSE,
                                scaling = scaling))
      rows[[length(rows) + 1L]] <- data.frame(
        width = as.integer(w), window = as.integer(k),
        win_from = per[1L], win_to = per[w],
        id = pn$id[i], period = as.character(pn$period)[i],
        eff = as.numeric(f$eff), status = as.integer(f$status),
        stringsAsFactors = FALSE)
    }
  }

  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  bad <- sum(!is.finite(tab$eff))
  if (bad) {
    warning(bad, " DMU-period(s) across all windows were not solved and are ",
            "NA. See the `status` column for the lpSolve return codes.",
            call. = FALSE)
  }

  structure(list(
    table = tab, panel = pn,
    width = width, rts = rts, orientation = orientation,
    n_window = vapply(width, function(w) pn$n_period - w + 1L, integer(1)),
    unsolved = bad, call = call
  ), class = "dea_window")
}

## Widths that make a window: at least two periods, at most the panel's length,
## whole numbers only. A width equal to the panel's length is one window -- a
## pooled fit, not an error, and the only case the method degenerates in.
.window_widths <- function(width, n_period) {
  if (!length(width) || anyNA(width) || !is.numeric(width))
    stop("`width` must be one or more whole numbers.", call. = FALSE)
  if (any(width != as.integer(width)))
    stop("`width` must be whole: a window is a count of periods.", call. = FALSE)
  width <- sort(unique(as.integer(width)))
  if (any(width < 2L)) {
    stop("`width` must be at least 2; a window of one period is just a ",
         "period-by-period fit, which dea() already does.", call. = FALSE)
  }
  if (any(width > n_period)) {
    stop("`width` of ", paste(width[width > n_period], collapse = ", "),
         " exceeds the panel's ", n_period, " period(s).", call. = FALSE)
  }
  width
}

## Per-DMU stability across windows, which is what the method is for. Computed
## from the long table rather than carried alongside it, so there is one
## definition of it and print(), summary() and the user all see the same one.
##
## TWO RANGES, AND ONLY ONE OF THEM IS ABOUT THE WINDOWS. `range` is the
## spread of every score the DMU received, over all its periods, windows and
## widths -- the window literature's "total range". It includes the DMU's
## genuine change from one period to the next, so a unit that improved steadily
## has a large one under a perfectly stable technology. `ldr`, the largest
## difference between scores of the SAME DMU-period in different windows of the
## SAME width, holds the observation fixed and moves only the window, so it is
## the part that measures the technology's instability. The first version of
## summary() reported `range` alone as "range across windows", which read the
## first as the second; a paper illustration built on it reported a number
## that meant nothing.
.window_stability <- function(tab) {
  key <- factor(tab$id, levels = unique(tab$id))
  f <- function(fun) vapply(split(tab$eff, key), function(v) {
    v <- v[is.finite(v)]
    if (length(v)) fun(v) else NA_real_
  }, numeric(1))
  ## Range of each (width, DMU, period) cell, then the largest per DMU.
  ok <- is.finite(tab$eff)
  cell <- interaction(tab$width[ok], tab$id[ok], tab$period[ok], drop = TRUE,
                      lex.order = TRUE)
  cr <- vapply(split(tab$eff[ok], cell), function(v) max(v) - min(v),
               numeric(1))
  cid <- vapply(split(as.character(tab$id[ok]), cell), `[`, character(1), 1L)
  ldr <- vapply(levels(key), function(k) {
    v <- cr[cid == k]
    if (length(v)) max(v) else NA_real_
  }, numeric(1))
  data.frame(id = levels(key), n = as.integer(table(key)),
             mean = f(mean), min = f(min), max = f(max),
             range = f(max) - f(min), ldr = unname(ldr), row.names = NULL,
             stringsAsFactors = FALSE)
}

print.dea_window <- function(x, ...) {
  cat("--- Window analysis (",
      if (x$orientation == "in") "input" else "output", " oriented, rts = ",
      x$rts, ") ---\n", sep = "")
  cat("DMUs: ", x$panel$n_id, "   periods: ", x$panel$n_period,
      "   DMU-periods: ", x$panel$n, "\n", sep = "")
  cat("period order: ", paste(x$panel$periods, collapse = " < "), "\n", sep = "")
  cat("width(s): ", paste(x$width, collapse = ", "), "   windows: ",
      paste(x$n_window, collapse = ", "), "   fits scored: ", nrow(x$table),
      "\n", sep = "")
  cat("\nmean efficiency by window\n")
  print(round(.window_means(x$table), 4))
  if (x$unsolved) cat("\nnot solved: ", x$unsolved, " DMU-period(s)\n", sep = "")
  if (length(x$width) == 1L) {
    cat("\nThe width is not pinned down by theory. Try two or three widths -- ",
        "`width` takes\na vector -- and read the answer as stable only where ",
        "it does not move.\n", sep = "")
  }
  invisible(x)
}

## One mean per (width, window), labelled by the periods the window spans.
.window_means <- function(tab) {
  key <- paste0("w", tab$width, ": ", tab$win_from, "-", tab$win_to)
  key <- factor(key, levels = unique(key))
  vapply(split(tab$eff, key), function(v) {
    v <- v[is.finite(v)]
    if (length(v)) mean(v) else NA_real_
  }, numeric(1))
}

summary.dea_window <- function(object, ...) {
  print(object)
  st <- .window_stability(object$table)
  cat("\nstability across windows, by DMU (first ", min(10L, nrow(st)), " of ",
      nrow(st), ")\n", sep = "")
  pr <- st
  for (cl in c("mean", "min", "max", "range", "ldr")) pr[[cl]] <- round(pr[[cl]], 4)
  print(utils::head(pr, 10))
  v <- st$ldr[is.finite(st$ldr)]
  if (length(v)) {
    cat("\nsame DMU-period, different windows (ldr): median ",
        format(round(stats::median(v), 4)), ", max ", format(round(max(v), 4)),
        " (DMU ", st$id[which.max(st$ldr)], ")\n", sep = "")
    cat("A DMU whose ldr is large is not being measured by a stable ",
        "technology;\nits score depends on which periods the window happens ",
        "to contain. `range`\nalso includes its change between periods, ",
        "which is not instability.\n", sep = "")
  }
  invisible(object)
}

## Mean efficiency by window, one line per width, against the window's last
## period -- which is where a window is conventionally dated.
plot.dea_window <- function(x, ...) {
  tab <- x$table
  pos <- match(tab$win_to, x$panel$periods)
  key <- paste(tab$width, pos)
  agg <- tapply(tab$eff, factor(key, levels = unique(key)),
                function(v) { v <- v[is.finite(v)]
                              if (length(v)) mean(v) else NA_real_ })
  wd <- as.integer(sub(" .*", "", names(agg)))
  px <- as.integer(sub(".* ", "", names(agg)))
  plot(range(px), range(agg, na.rm = TRUE), type = "n",
       xlab = "window (last period)", ylab = "mean efficiency",
       main = "Window analysis", xaxt = "n", ...)
  graphics::axis(1, at = sort(unique(px)), labels = x$panel$periods[sort(unique(px))])
  for (i in seq_along(x$width)) {
    k <- wd == x$width[i]
    graphics::lines(px[k], agg[k], type = "b", pch = 16, lty = i)
  }
  if (length(x$width) > 1L) {
    legend("topright", legend = paste("width", x$width),
           lty = seq_along(x$width), pch = 16, bty = "n")
  }
  invisible(data.frame(width = wd, period = x$panel$periods[px],
                       mean_eff = as.numeric(agg), row.names = NULL))
}

nobs.dea_window <- function(object, ...) nrow(object$table)
