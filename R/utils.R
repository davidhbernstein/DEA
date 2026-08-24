## ---------------------------------------------------------------------------
## Shared constants, input coercion and validation.
## ---------------------------------------------------------------------------

## Tolerances. DEA is a linear program, so the answer is exact up to the
## solver's tolerance -- these are about deciding when a number that "should"
## be 0 or 1 actually is, not about statistical precision.
.DEA_CONSTANTS <- list(
  TOL_EFF    = 1e-9,   ## snap an efficiency this close to 1 onto 1
  TOL_LAMBDA = 1e-8,   ## a lambda below this is not a peer
  TOL_SLACK  = 1e-8,   ## a slack below this is zero
  LP_EPSEL   = 1e-12,  ## lpSolveAPI rounding epsilon
  MIN_POS    = 1e-12   ## floor for a denominator that must be positive
)

`%||%` <- function(a, b) if (is.null(a)) b else a

## Resolve an x/y argument into a numeric matrix with one row per DMU.
##
## Accepts a matrix, a data frame, a bare numeric vector (one variable), or --
## when `data` is supplied -- a character vector of column names or a one-sided
## formula.  The formula/character forms exist so that a call reads like the
## rest of R's modelling code; the matrix form exists because that is what the
## DEA literature and every other DEA package use.
.dea_matrix <- function(arg, data, what) {
  if (inherits(arg, "formula")) {
    if (is.null(data)) {
      stop("`", what, "` was given as a formula, which needs `data`.", call. = FALSE)
    }
    if (length(arg) != 2L) {
      stop("`", what, "` must be a ONE-sided formula such as ~ x1 + x2.", call. = FALSE)
    }
    tt <- stats::terms(arg, data = data)
    m  <- stats::model.matrix(tt, data = as.data.frame(data))
    m  <- m[, colnames(m) != "(Intercept)", drop = FALSE]
    if (!ncol(m)) stop("`", what, "` selected no variables.", call. = FALSE)
    return(.dea_num(m, what))
  }
  if (is.character(arg)) {
    if (is.null(data)) {
      stop("`", what, "` was given as column names, which needs `data`.", call. = FALSE)
    }
    data <- as.data.frame(data)
    miss <- setdiff(arg, names(data))
    if (length(miss)) {
      stop("`", what, "`: column(s) not found in `data`: ",
           paste(miss, collapse = ", "), call. = FALSE)
    }
    return(.dea_num(as.matrix(data[, arg, drop = FALSE]), what))
  }
  if (is.data.frame(arg)) return(.dea_num(as.matrix(arg), what))
  if (is.vector(arg) && is.numeric(arg)) {
    return(.dea_num(matrix(arg, ncol = 1L, dimnames = list(names(arg), what)), what))
  }
  if (is.matrix(arg)) return(.dea_num(arg, what))
  stop("`", what, "` must be a matrix, data frame, numeric vector, or -- with ",
       "`data` -- column names or a one-sided formula.", call. = FALSE)
}

.dea_num <- function(m, what) {
  storage.mode(m) <- "double"
  if (!is.numeric(m) || any(!is.finite(m))) {
    stop("`", what, "` contains missing or non-finite values. DEA has no ",
         "missing-data theory: drop or impute those rows before calling.",
         call. = FALSE)
  }
  if (is.null(colnames(m))) {
    colnames(m) <- paste0(what, seq_len(ncol(m)))
  }
  m
}

## Shared validation for every entry point.  `require_positive` is model
## specific: radial DEA tolerates zeros, the slacks-based measure divides by
## the DMU's own levels and does not.
.dea_check <- function(X, Y, require_positive = FALSE, model = "dea") {
  model <- if (is.null(model)) "dea" else model
  if (nrow(X) != nrow(Y)) {
    stop("`x` and `y` describe different numbers of DMUs (", nrow(X), " vs ",
         nrow(Y), ").", call. = FALSE)
  }
  n <- nrow(X)
  if (n < 1L) stop("No DMUs supplied.", call. = FALSE)
  if (any(X < 0) || any(Y < 0)) {
    stop("Negative values in `x` or `y`. The radial and slacks-based models ",
         "assume a non-negative technology; use dea_ddf(), whose translation ",
         "property handles data of either sign.", call. = FALSE)
  }
  if (require_positive && (any(X <= 0) || any(Y <= 0))) {
    stop(model, "() divides by each DMU's own inputs and outputs, so every ",
         "value must be strictly positive. Found ", sum(X <= 0), " zero input ",
         "and ", sum(Y <= 0), " zero output entries. Use dea() (radial ",
         "efficiency tolerates zeros) or dea_ddf() with an explicit direction.",
         call. = FALSE)
  }
  if (any(rowSums(Y) <= 0)) {
    bad <- which(rowSums(Y) <= 0)
    stop("DMU(s) ", paste(utils::head(bad, 5), collapse = ", "),
         if (length(bad) > 5) ", ..." else "",
         " produce no output at all. Such a DMU has no well-defined output ",
         "efficiency and distorts the input-oriented frontier; drop it.",
         call. = FALSE)
  }
  ## Not an error -- DEA is defined for any n -- but with few DMUs relative to
  ## the number of dimensions almost everything is efficient by construction
  ## and the scores carry no information.  The rule of thumb is n >= 3(p+q).
  d <- ncol(X) + ncol(Y)
  if (n < 3 * d) {
    warning("Only ", n, " DMUs for ", ncol(X), " inputs and ", ncol(Y),
            " outputs. With n < 3(p+q) = ", 3 * d, " most DMUs are efficient ",
            "by dimension alone, not by performance -- see ?dea, 'the curse ",
            "of dimensionality'.", call. = FALSE)
  }
  invisible(TRUE)
}

## Column scaling.  Radial efficiency, the SBM and a proportional-direction DDF
## are all invariant to the units of each input and output, so rescaling every
## column to mean one leaves the answer alone and leaves the LP far better
## conditioned.  Slacks and peers are mapped back to the original units on the
## way out.  A DDF with a FIXED direction is NOT units invariant, so its
## direction is scaled by the same factors -- which restores invariance rather
## than breaking it.
.dea_scale <- function(X, Y, on = TRUE) {
  if (!on) {
    return(list(X = X, Y = Y, sx = rep(1, ncol(X)), sy = rep(1, ncol(Y))))
  }
  sx <- colMeans(X); sy <- colMeans(Y)
  sx[!is.finite(sx) | sx <= 0] <- 1
  sy[!is.finite(sy) | sy <= 0] <- 1
  list(X = sweep(X, 2L, sx, "/"), Y = sweep(Y, 2L, sy, "/"), sx = sx, sy = sy)
}

.dea_dmu_names <- function(X, Y) {
  nm <- rownames(X) %||% rownames(Y)
  if (is.null(nm)) as.character(seq_len(nrow(X))) else nm
}


## The evaluated set faces a weaker standard than the reference set: it does
## not span anything, so the dimensionality rule of thumb does not apply to it
## and a single evaluated DMU is perfectly well defined.
.dea_check_eval <- function(X, Y, require_positive = FALSE, model = "dea") {
  if (any(X < 0) || any(Y < 0)) {
    stop("Negative values in the evaluated `x`/`y`. Use dea_ddf(), whose ",
         "translation property handles data of either sign.", call. = FALSE)
  }
  if (require_positive && (any(X <= 0) || any(Y <= 0))) {
    stop(model, "() divides by each evaluated DMU's own inputs and outputs, ",
         "so every value must be strictly positive.", call. = FALSE)
  }
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
## Case-insensitive match.arg, so that rts = "VRS" and rts = "vrs" are the same
## choice. Applied work writes these acronyms in capitals and it is not worth
## an error.
##
## Aliases for the two one-sided technologies: Benchmarking writes
## sum(lambda) <= 1 as "drs" and sum(lambda) >= 1 as "irs". This package writes
## them as "nirs"/"ndrs" because those name the RESTRICTION rather than the
## conclusion -- a DMU evaluated against a non-increasing technology has not
## thereby been found to have decreasing returns. Both spellings are accepted
## so that switching packages does not silently change the model.
## ---------------------------------------------------------------------------
.RTS_ALIAS <- c(drs = "nirs", irs = "ndrs", nirts = "nirs", ndrts = "ndrs")

.match_arg_ci <- function(arg, choices, what) {
  if (missing(arg) || is.null(arg)) return(choices[1L])
  if (length(arg) > 1L) {
    if (identical(as.character(arg), as.character(choices))) return(choices[1L])
    stop("`", what, "` must be a single value; got ", length(arg), ".", call. = FALSE)
  }
  a <- tolower(as.character(arg))
  if (identical(what, "rts") && a %in% names(.RTS_ALIAS)) a <- .RTS_ALIAS[[a]]
  hit <- match(a, tolower(choices))
  if (is.na(hit)) {
    stop("`", what, "` must be one of ", paste0("\"", choices, "\"", collapse = ", "),
         "; got \"", arg, "\".", call. = FALSE)
  }
  choices[hit]
}
