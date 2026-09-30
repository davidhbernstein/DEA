## ---------------------------------------------------------------------------
## dea_panel() -- the panel layer.
##
## This exists as its own object rather than as arguments to dea_malmquist()
## because it is the prerequisite for four different things, not one: the
## Malmquist index, window analysis, the Malmquist-Luenberger index and the
## Hicks-Moorsteen family. Designing it around whichever arrives first is how
## an interface ends up with a `period` argument that means something slightly
## different in each of them.
##
## What a panel IS, here: a long-format stack of DMU-periods. One row per
## (id, period), the same input and output columns throughout, and nothing
## assumed about balance. Every index built on top of it works period-pair by
## period-pair, and each pair uses the DMUs observed in BOTH -- so an
## unbalanced panel is an ordinary case rather than an error, and the object
## reports which pairs lost DMUs and how many.
##
## THE ORDER OF THE PERIODS IS DATA, NOT A GUESS. A productivity index is a
## statement about consecutive periods, so getting the order wrong silently
## reverses every index it produces. A factor's levels are honoured as given;
## anything else is sorted, which is right for numbers and for ISO-like dates
## and is WRONG for labels such as "Q1", "Q10", "Q2". The object prints the
## order it settled on, and a factor is the way to override it.
##
## THE TECHNOLOGY IS ESTIMATED PER PERIOD, from that period's DMUs alone. That
## is what makes a frontier shift measurable; pooling periods would define it
## away.
## ---------------------------------------------------------------------------

dea_panel <- function(x, y, id, period, data = NULL) {
  call <- match.call()
  .dea_check_data_arg(data)
  X <- .dea_matrix(x, data, "x")
  Y <- .dea_matrix(y, data, "y")
  if (nrow(X) != nrow(Y)) {
    stop("`x` has ", nrow(X), " rows and `y` has ", nrow(Y),
         "; a panel needs one row of each per DMU-period.", call. = FALSE)
  }
  id <- .panel_key(id, data, nrow(X), "id")
  pd <- .panel_key(period, data, nrow(X), "period")

  ## Period order: a factor's levels are the user's statement of order, and
  ## sort() is the fallback. See the header -- this is the one place where a
  ## silent mistake reverses every index downstream.
  lev <- if (is.factor(pd)) levels(droplevels(pd)) else sort(unique(pd))
  pd  <- factor(as.character(pd), levels = as.character(lev))
  idc <- as.character(id)

  if (nlevels(pd) < 2L) {
    stop("`period` has only ", nlevels(pd), " distinct value(s). A panel index ",
         "compares periods, so it needs at least two.", call. = FALSE)
  }
  dup <- duplicated(data.frame(idc, as.character(pd)))
  if (any(dup)) {
    w <- which(dup)
    stop(sum(dup), " duplicated (id, period) pair(s), first at row ", w[1],
         " (id ", idc[w[1]], ", period ", as.character(pd)[w[1]], "). A panel ",
         "has one observation per DMU per period; aggregate or subset first.",
         call. = FALSE)
  }

  rows <- split(seq_len(nrow(X)), pd)
  ids  <- lapply(rows, function(i) idc[i])
  all_ids <- unique(idc)
  balanced <- all(vapply(ids, function(v) setequal(v, all_ids), logical(1)))

  structure(list(
    X = X, Y = Y, id = idc, period = pd,
    periods = levels(pd), ids = all_ids,
    rows = rows, period_ids = ids,
    n = nrow(X), n_id = length(all_ids), n_period = nlevels(pd),
    p = ncol(X), q = ncol(Y),
    balanced = balanced,
    call = call
  ), class = "dea_panel")
}

## Resolve `id`/`period` the same way x and y are resolved: a bare vector, or a
## column name when `data` is given. Deliberately NOT .dea_matrix(), which
## coerces to double -- an id is a label and a period may be one too.
.panel_key <- function(arg, data, n, what) {
  if (inherits(arg, "formula")) {
    if (length(arg) != 2L)
      stop("`", what, "` must be a ONE-sided formula such as ~ year.", call. = FALSE)
    arg <- all.vars(arg)
  }
  if (is.character(arg) && length(arg) == 1L && !is.null(data)) {
    d <- as.data.frame(data)
    if (!arg %in% names(d))
      stop("`", what, "`: column \"", arg, "\" not found in `data`.", call. = FALSE)
    arg <- d[[arg]]
  }
  if (length(arg) != n) {
    stop("`", what, "` has length ", length(arg), " but the data has ", n,
         " rows.", call. = FALSE)
  }
  if (anyNA(arg))
    stop("`", what, "` contains NA. Every row must say which DMU and which ",
         "period it is.", call. = FALSE)
  arg
}

## The (x, y) of one period, and the ids in the order those rows appear.
.panel_slice <- function(pn, t) {
  i <- pn$rows[[t]]
  list(X = pn$X[i, , drop = FALSE], Y = pn$Y[i, , drop = FALSE], id = pn$id[i])
}

## The DMUs observed in BOTH periods, in the order of the first, with the row
## positions needed to line the two slices up. An index is a comparison of a
## DMU with itself, so anything seen in only one period has no index for that
## pair -- reported, never silently dropped.
.panel_pair <- function(pn, t0, t1) {
  a <- .panel_slice(pn, t0); b <- .panel_slice(pn, t1)
  keep <- a$id[a$id %in% b$id]
  list(a = a, b = b, id = keep,
       ia = match(keep, a$id), ib = match(keep, b$id),
       dropped = length(a$id) + length(b$id) - 2L * length(keep))
}

print.dea_panel <- function(x, ...) {
  cat("--- Panel ---\n")
  cat("DMU-periods: ", x$n, "   DMUs: ", x$n_id, "   periods: ", x$n_period,
      "\n", sep = "")
  cat("inputs: ", x$p, "   outputs: ", x$q, "\n", sep = "")
  cat("period order: ", paste(x$periods, collapse = " < "), "\n", sep = "")
  cat("balanced: ", if (x$balanced) "yes" else "no", "\n", sep = "")
  if (!x$balanced) {
    k <- vapply(x$rows, length, integer(1))
    cat("  DMUs per period: ", paste(sprintf("%s=%d", names(k), k),
                                     collapse = "  "), "\n", sep = "")
  }
  invisible(x)
}

summary.dea_panel <- function(object, ...) {
  print(object)
  cat("\nconsecutive pairs\n")
  p <- object$periods
  for (k in seq_len(length(p) - 1L)) {
    pr <- .panel_pair(object, p[k], p[k + 1L])
    cat(sprintf("  %s -> %s : %d DMU(s) in both", p[k], p[k + 1L],
                length(pr$id)))
    if (pr$dropped) cat(sprintf(", %d observation(s) not matched", pr$dropped))
    cat("\n")
  }
  invisible(object)
}

nobs.dea_panel <- function(object, ...) object$n
