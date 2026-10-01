## ---------------------------------------------------------------------------
## dea_categorical() -- categorical environmental variables (Banker and Morey
## 1986b), as NESTED REFERENCE SETS rather than as a mixed-integer program.
##
## THE PROBLEM.  A school in a deprived catchment, a hospital in a rural
## district, a branch in a small town: the DMU operates in an environment it did
## not choose and cannot change, and scoring it against units in easier
## environments measures the environment rather than the management.  The
## environment is ORDERED but not a quantity -- "more deprived" is a rank, not a
## number of units of anything -- so it cannot simply be added as an input.
##
## THE SOLUTION IS A RESTRICTION ON WHO MAY BE A PEER, AND NOTHING ELSE.  Order
## the categories from least to most favourable.  A DMU in category k is scored
## against every DMU in categories 1..k -- everyone whose environment was no
## better than its own.  That is `xref`/`yref`, which this package has had from
## the start, applied once per category.  No new program, no integer variables.
##
## WHY NOT BANKER AND MOREY'S OWN FORMULATION.  Theirs is a mixed-integer
## program, and Charnes, Cooper, Lewin and Seiford (1994, ch. 3) are explicit
## that it "incorrectly prescribes categories that are undefined or
## meaningless" -- the MIP can hand a DMU a target in a category that does not
## exist.  Adapting the SOLUTION PROCEDURE rather than the model is both correct
## and allows several categorical variables at once, which the MIP does not.  So
## the MIP is deliberately not built.
##
## WHICH END IS WHICH, AND WHY IT IS ASKED FOR EXPLICITLY.  Higher category =
## MORE FAVOURABLE environment.  Category 1 is scored against itself alone;
## the top category is scored against everybody.  Getting this backwards does
## not error -- it silently answers the opposite question -- so `category` must
## be an ordered factor or a numeric rank, an unordered factor is refused, and
## the reference-set size per category is reported so that the direction is
## visible in the output rather than taken on trust.
##
## THE COST OF NESTING IS SAMPLE SIZE, AND IT FALLS ENTIRELY ON CATEGORY 1.
## That category is scored against only its own members, so a design that is
## comfortable pooled can be far below the n >= 3(p+q) rule of thumb once split.
## Warned about per category, because the scores stay finite and plausible and
## nothing else would reveal it.
##
## IF THE CATEGORIES ARE NOT NESTED, THIS IS THE WRONG FUNCTION.  Public versus
## private universities are not ordered -- neither environment is a subset of
## the other's advantages -- and inventing an order produces a number that
## answers nothing.  The right analysis there is separate ones, and `print()`
## says so rather than leaving the reader to wonder.
## ---------------------------------------------------------------------------

dea_categorical <- function(x, y, category, data = NULL,
                            rts = c("vrs", "crs", "nirs", "ndrs"),
                            orientation = c("in", "out"),
                            slack = TRUE, peers = TRUE, scaling = TRUE) {

  call <- match.call()
  .dea_check_data_arg(data)
  rts <- .match_arg_ci(rts, c("vrs", "crs", "nirs", "ndrs"), "rts")
  orientation <- .match_arg_ci(orientation, c("in", "out"), "orientation")
  t0 <- proc.time()[["elapsed"]]

  X <- .dea_matrix(x, data, "x")
  Y <- .dea_matrix(y, data, "y")
  .dea_check(X, Y)
  dmu <- .dea_dmu_names(X, Y)
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)

  cat_in <- .cat_resolve(category, data, n)
  k <- cat_in$index                      ## 1..K, 1 = least favourable
  lev <- cat_in$levels
  K <- length(lev)

  ## The pooled fit, for the comparison that says whether the categorical
  ## structure did anything at all. One extra sweep, and without it the output
  ## cannot distinguish "the environment matters" from "the environment was
  ## recorded but is irrelevant here".
  pooled <- suppressWarnings(
    dea(X, Y, rts = rts, orientation = orientation, slack = slack,
        peers = FALSE, scaling = scaling))

  eff <- rep(NA_real_, n); st <- integer(n)
  L <- if (peers) matrix(0, n, n) else NULL
  sx <- if (slack) matrix(NA_real_, n, p) else NULL
  sy <- if (slack) matrix(NA_real_, n, q) else NULL
  refsize <- integer(n)
  small <- integer(0)

  for (kk in seq_len(K)) {
    ev <- which(k == kk)
    rf <- which(k <= kk)                 ## everyone no better off
    refsize[ev] <- length(rf)
    if (length(rf) < 3 * (p + q)) small <- c(small, kk)
    ## The dimensionality warning is raised once below, per category and with
    ## the category named, rather than K times with no indication of which.
    f <- suppressWarnings(
      dea(X[ev, , drop = FALSE], Y[ev, , drop = FALSE],
          rts = rts, orientation = orientation, slack = slack,
          peers = peers, scaling = scaling,
          xref = X[rf, , drop = FALSE], yref = Y[rf, , drop = FALSE]))
    eff[ev] <- as.numeric(f$eff)
    st[ev] <- f$status
    if (peers) L[ev, rf] <- f$lambda
    if (slack) { sx[ev, ] <- f$slack_x; sy[ev, ] <- f$slack_y }
  }

  if (length(small)) {
    warning(length(small), " categor", if (length(small) == 1L) "y" else "ies",
            " (", paste(lev[small], collapse = ", "), ") ",
            if (length(small) == 1L) "has" else "have",
            " a reference set below the n >= 3(p+q) = ", 3 * (p + q),
            " rule of thumb; the smallest is ", min(refsize),
            " DMUs. Nesting costs sample size and the cost falls on the ",
            "lowest categories, where almost everything can be efficient by ",
            "dimension alone. The scores stay finite and plausible, so ",
            "nothing else would show this.", call. = FALSE)
  }
  .dea_report_unsolved(st, n, FALSE, "dea_categorical")

  eff[is.finite(eff) & abs(eff - 1) < .DEA_CONSTANTS$TOL_EFF] <- 1
  names(eff) <- dmu
  if (!is.null(L)) { L[abs(L) < .DEA_CONSTANTS$TOL_LAMBDA] <- 0; dimnames(L) <- list(dmu, dmu) }
  if (slack) { dimnames(sx) <- list(dmu, colnames(X)); dimnames(sy) <- list(dmu, colnames(Y)) }

  ## Both scores on the common (0, 1] scale, so that the effect has one sign
  ## whichever orientation was used. A smaller reference set is a smaller
  ## technology, so the categorical score is never the WORSE of the two.
  sc_cat <- if (orientation == "in") eff else 1 / eff
  sc_pool <- if (orientation == "in") as.numeric(pooled$eff) else 1 / as.numeric(pooled$eff)

  structure(list(
    eff = eff, model = "radial", rts = rts, orientation = orientation,
    super = FALSE,
    lambda = L, sum_lambda = if (is.null(L)) NULL else rowSums(L),
    slack_x = sx, slack_y = sy,
    efficient = if (is.null(sx)) eff == 1 else
      eff == 1 & rowSums(sx) <= .DEA_CONSTANTS$TOL_SLACK &
      rowSums(sy) <= .DEA_CONSTANTS$TOL_SLACK,
    status = st,
    ## The categorical layer's own output.
    category = factor(lev[k], levels = lev, ordered = TRUE),
    category_levels = lev, n_category = K,
    reference_size = refsize,
    pooled_eff = stats::setNames(as.numeric(pooled$eff), dmu),
    category_effect = stats::setNames(sc_cat - sc_pool, dmu),
    x = X, y = Y, xref = X, yref = Y, self_ref = FALSE,
    dmu = dmu, ref = dmu, n = n, nref = n, p = p, q = q,
    scaling = scaling, slack = slack,
    total_time = proc.time()[["elapsed"]] - t0,
    call = call
  ), class = c("dea_categorical", "dea"))
}

## Resolve `category` into an ordered index. An unordered factor is refused
## rather than coerced: alphabetical order is not an ordering of environments,
## and silently inventing one is exactly the failure this function exists to
## avoid.
.cat_resolve <- function(category, data, n) {
  if (is.character(category) && length(category) == 1L && !is.null(data)) {
    data <- as.data.frame(data)
    if (!category %in% names(data)) {
      stop("`category`: column \"", category, "\" not found in `data`.", call. = FALSE)
    }
    category <- data[[category]]
  }
  if (length(category) != n) {
    stop("`category` has length ", length(category), " and there are ", n,
         " DMUs.", call. = FALSE)
  }
  if (anyNA(category)) {
    stop("`category` has ", sum(is.na(category)), " missing value(s). A DMU ",
         "with no recorded environment cannot be placed in a nest; drop it or ",
         "assign it deliberately.", call. = FALSE)
  }
  if (is.factor(category) && !is.ordered(category)) {
    stop("`category` is an unordered factor. This model nests the reference ",
         "sets, so it needs to know which environment is more favourable, and ",
         "the level order of a plain factor is alphabetical rather than ",
         "meaningful. Pass an ordered factor -- factor(z, levels = ..., ",
         "ordered = TRUE) -- or a numeric rank.\n  If the groups are not ",
         "ordered at all (public vs private, urban vs coastal), no nesting is ",
         "correct and the right analysis is a separate one per group, or a ",
         "metafrontier via `xref`/`yref`.", call. = FALSE)
  }
  if (is.ordered(category)) {
    category <- droplevels(category)
    lev <- levels(category)
    idx <- as.integer(category)
  } else {
    if (!is.numeric(category)) {
      stop("`category` must be an ordered factor or a numeric rank; got ",
           class(category)[1], ".", call. = FALSE)
    }
    if (any(!is.finite(category))) {
      stop("`category` contains non-finite values.", call. = FALSE)
    }
    u <- sort(unique(category))
    lev <- as.character(u)
    idx <- match(category, u)
  }
  if (length(lev) < 2L) {
    stop("`category` has only one distinct value, so every DMU is scored ",
         "against every other and this is an ordinary dea() fit.", call. = FALSE)
  }
  list(index = idx, levels = lev)
}

print.dea_categorical <- function(x, ...) {
  cat("--- Categorical environment: nested reference sets ---\n")
  cat("technology:  ", toupper(x$rts), ", ",
      if (x$orientation == "in") "input" else "output", " orientation\n", sep = "")
  cat("DMUs: ", x$n, "   inputs: ", x$p, "   outputs: ", x$q,
      "   categories: ", x$n_category, "   (", format(round(x$total_time, 3)),
      " sec)\n", sep = "")
  tb <- data.frame(
    category = x$category_levels,
    n = as.integer(table(x$category)),
    reference = vapply(seq_len(x$n_category),
                       function(i) x$reference_size[which(as.integer(x$category) == i)[1]],
                       integer(1)),
    mean_eff = round(tapply(x$eff, x$category, mean, na.rm = TRUE), 4),
    mean_pooled = round(tapply(x$pooled_eff, x$category, mean, na.rm = TRUE), 4),
    effect = round(tapply(x$category_effect, x$category, mean, na.rm = TRUE), 4),
    row.names = NULL)
  cat("\ncategory 1 is the LEAST favourable environment and is scored against\n",
      "itself alone; each later one is scored against every category up to it\n", sep = "")
  print(tb)
  cat("\n`effect` is how much a DMU's efficiency improves on the (0, 1] scale\n",
      "by being compared only with units no better placed than itself.\n", sep = "")
  if (max(abs(x$category_effect), na.rm = TRUE) < 1e-8) {
    cat("\nIt is zero everywhere: the nesting changed nothing on this data, so\n",
        "the environment, as recorded, is doing no work here.\n", sep = "")
  }
  invisible(x)
}

summary.dea_categorical <- function(object, ...) {
  print(object)
  cat("\nIf the categories are not genuinely nested -- if no environment's\n",
      "advantages are a subset of another's -- this model is the wrong one and\n",
      "separate analyses are the right one.\n", sep = "")
  cat("\nfirst rows\n")
  print(utils::head(data.frame(
    dmu = object$dmu, category = as.character(object$category),
    reference = object$reference_size,
    eff = round(object$eff, 4), pooled = round(object$pooled_eff, 4),
    effect = round(object$category_effect, 4), row.names = NULL), 10))
  invisible(object)
}
