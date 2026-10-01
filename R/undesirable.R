## ---------------------------------------------------------------------------
## dea_undesirable() -- efficiency when some outputs are BAD.
##
## Pollution, defects, non-performing loans, patient mortality.  A bad output is
## not an output: free disposal runs the wrong way for it, so listing it among
## `y` and solving the ordinary program asks the technology to let a unit
## produce MORE of it, and rewards a unit for doing so.
##
## There are three routes in the literature and this function offers the two
## that need no new solver.  Which one is in use changes what the answer MEANS,
## so `method` has no default that is right for every question and both are
## documented rather than ranked.
##
## ROUTE 1, `method = "translate"` -- Seiford and Zhu (2002).  Replace the bad
## b by w - b for a constant w large enough to keep it positive, so that more of
## the transformed variable is better, and run the ordinary radial program.  The
## transformed variable is a genuine output: a unit with less pollution has more
## of w - b, and free disposal of it says a unit could always pollute more,
## which is true.
##
## THE PRECONDITION IS THE WHOLE STORY, AND IT IS NARROWER THAN IT LOOKS.  w is
## arbitrary, so the answer is only meaningful where the model is invariant to
## translating that output.  Measured on a 40-DMU design, running exactly this
## transformation at w = max(b) + 1, 20, 100 and 1000, the largest change in ANY
## DMU's score:
##
##     vrs  in    7.0e-11        crs  in    1.8e-02
##     vrs  out   2.5e-01        crs  out   2.6e-02
##     nirs in    4.0e-03        ndrs in    1.8e-02
##     nirs out   2.5e-01        ndrs out   2.6e-02
##
## So it holds in EXACTLY ONE CELL, vrs and input-oriented, and the gap is eight
## orders of magnitude -- not a tolerance question.  Not "a convex technology"
## either: nirs is convex and moves by 4e-3.  `deaR::undesirable_basic`'s help says
## rts must be vrs and says nothing about the orientation, which is the half of
## the precondition that bites hardest: studying outputs is exactly when a user
## reaches for an output orientation.  So the other combinations are REFUSED
## here rather than served with a number, the way dea_profit() refuses a cone.
##
## ROUTE 2, `method = "input"` -- Hailu and Veeman (2001).  Put the bad among
## the inputs.  Free disposal then runs the right way (a unit could always use
## more of it), it needs no constant and no invariance, and it works under every
## technology and orientation.  What it does not encode is NULL-JOINTNESS: the
## program is free to recommend less bad at no cost in good output, which the
## physics of most production processes forbids.  That is the objection, it is a
## real one, and it is the reason route 3 exists.
##
## ROUTE 3, NOT IMPLEMENTED -- the weak-disposability technology of Fare and
## Grosskopf, which puts the bads on EQUALITY rows with an abatement factor and
## is the only one of the three that encodes null-jointness.  It needs a new
## program rather than a rearrangement of this one.  ROADMAP.md carries it.
##
## THE SECOND PRECONDITION IS NUMERICAL, AND THE FAILURE IS SILENT.  Invariance
## to w is exact in algebra and not in double precision: a large w makes the
## transformed column nearly constant, and after the column is scaled to mean
## one its values differ in the last few digits.  Measured on the same design,
## against the default w = max(b) + 1, writing `rel` for the column's relative
## spread range(b) / mean(w - b):
##
##     rel 1.4e+00  (default)   exact          0 unsolved
##     rel 8.1e-02              3e-12          0
##     rel 7.7e-04              1.7e-08        0
##     rel 7.7e-05              5.1e-07        0
##     rel 7.7e-06              2.5e-06        1 unsolved
##     rel 7.7e-07              0.188          3 unsolved, 7 failures
##     rel 7.7e-09              0.188          0 unsolved
##
## The last row is why this is gated rather than merely documented: at w = 1e9
## every program reports a clean optimum, nothing is NA, no warning fires
## anywhere, and the scores are wrong by 0.188.  Same shape as the set.column()
## trap in .lp_radial_at() -- the solver is content and the answer is not.  So a
## relative spread below 1e-3 warns, which is three orders of magnitude inside
## the first observed failure and keeps the error near 1e-8.
##
## WHAT THE SLACK ON A BAD MEANS, which is the one piece of good luck here.
## Under route 1 the output slack s on w - b satisfies (w - b_peer) =
## (w - b_o) + s, i.e. b_peer = b_o - s: the slack IS the reducible bad, in the
## caller's own units, with no sign flip.  Under route 2 it is an input slack on
## b, which is the same quantity.  So `slack_bad` means one thing under both
## routes even though it is stored in different places by each.
## ---------------------------------------------------------------------------

dea_undesirable <- function(x, y, bad, data = NULL,
                            method = c("translate", "input"),
                            rts = c("vrs", "crs", "nirs", "ndrs"),
                            orientation = c("in", "out"),
                            trans = NULL, ...) {

  call <- match.call()
  .dea_check_data_arg(data)
  method <- .match_arg_ci(method, c("translate", "input"), "method")
  rts <- .match_arg_ci(rts, c("vrs", "crs", "nirs", "ndrs"), "rts")
  orientation <- .match_arg_ci(orientation, c("in", "out"), "orientation")

  X <- .dea_matrix(x, data, "x")
  Y <- .dea_matrix(y, data, "y")
  B <- .dea_matrix(bad, data, "bad")
  if (nrow(B) != nrow(X)) {
    stop("`bad` has ", nrow(B), " rows and `x` has ", nrow(X),
         "; they describe different numbers of DMUs.", call. = FALSE)
  }
  if (ncol(B) < 1L) stop("`bad` selected no variables.", call. = FALSE)
  if (identical(colnames(B), paste0("bad", seq_len(ncol(B))))) {
    ## .dea_matrix() names an unnamed matrix after its argument; keep that, but
    ## a single column reads better without the index.
    if (ncol(B) == 1L) colnames(B) <- "bad"
  }
  nb <- ncol(B)

  if (identical(method, "translate")) {
    if (!identical(rts, "vrs") || !identical(orientation, "in")) {
      stop("method = \"translate\" is the Seiford-Zhu (2002) transformation, ",
           "which replaces each bad output by `trans - bad`. `trans` is an ",
           "ARBITRARY constant, so the result is only meaningful where the ",
           "model does not depend on it -- that is, where the model is ",
           "invariant to translating an output. Measured, that is rts = ",
           "\"vrs\" with orientation = \"in\" and nothing else -- where the ",
           "score moves by 7e-11 across four constants, against 4e-3 to 0.25 ",
           "in every other combination, including the other convex ones. At ",
           "rts = \"", rts, "\", orientation = \"", orientation, "\" the ",
           "answer would be a function of the constant you pick. Use ",
           "method = \"input\", which needs no constant and works under every ",
           "technology, or fit rts = \"vrs\" input-oriented. See ",
           "?dea_undesirable.", call. = FALSE)
    }
    tr <- .undes_trans(trans, B)
    Bt <- sweep(-B, 2L, tr, "+")            ## trans - bad, column by column
    if (any(Bt <= 0)) {
      bad_j <- which(apply(Bt, 2L, function(v) any(v <= 0)))
      stop("`trans` leaves the transformed bad output non-positive in column(s) ",
           paste(colnames(B)[bad_j], collapse = ", "),
           ". Every transformed value must be strictly positive -- the ",
           "default, max(bad) + 1 per column, guarantees it.", call. = FALSE)
    }
    ## The numerical half of the precondition. See the ladder above: the error
    ## is invisible until it is large, and then it is large with every program
    ## reporting success.
    rel <- apply(B, 2L, function(v) diff(range(v))) / colMeans(Bt)
    rel[!is.finite(rel)] <- Inf
    if (any(rel < 1e-3)) {
      j <- which.min(rel)
      warning("`trans` is large relative to the spread of \"", colnames(B)[j],
              "\": the transformed column's relative spread is ",
              format(rel[j], digits = 3), ". Invariance to `trans` is exact in ",
              "algebra and not in double precision -- measured on a 40-DMU ",
              "design, a relative spread of 7.7e-4 moves scores by 1.7e-8, ",
              "7.7e-6 loses a DMU to a numerical failure, and 7.7e-9 moves ",
              "them by 0.188 while every program still reports a clean ",
              "optimum. Use the default, max(bad) + 1 per column, unless there ",
              "is a reason not to.", call. = FALSE)
    }
    colnames(Bt) <- paste0(colnames(B), "_t")
    fit <- dea(X, cbind(Y, Bt), rts = rts, orientation = orientation, ...)
    bad_cols <- ncol(Y) + seq_len(nb)
    sb <- if (is.null(fit$slack_y)) NULL else fit$slack_y[, bad_cols, drop = FALSE]
  } else {
    tr <- NULL
    fit <- dea(cbind(X, B), Y, rts = rts, orientation = orientation, ...)
    bad_cols <- ncol(X) + seq_len(nb)
    sb <- if (is.null(fit$slack_x)) NULL else fit$slack_x[, bad_cols, drop = FALSE]
  }
  if (!is.null(sb)) colnames(sb) <- colnames(B)

  fit$bad <- B
  fit$bad_method <- method
  fit$bad_trans <- tr
  fit$slack_bad <- sb
  fit$nbad <- nb
  fit$call <- call
  class(fit) <- c("dea_undesirable", class(fit))
  fit
}

## deaR's default, per column, and it is a good one: max + 1 is the smallest
## round choice that keeps every transformed value strictly positive whatever
## the data. The answer does not depend on it -- that is the precondition above
## -- so the default is a convenience and not an assumption.
.undes_trans <- function(trans, B) {
  nb <- ncol(B)
  if (is.null(trans)) return(apply(B, 2L, max) + 1)
  if (!is.numeric(trans) || !all(is.finite(trans))) {
    stop("`trans` must be a finite numeric vector.", call. = FALSE)
  }
  if (length(trans) == 1L) return(rep(trans, nb))
  if (length(trans) != nb) {
    stop("`trans` must have length 1 or one entry per undesirable output (",
         nb, "); got ", length(trans), ".", call. = FALSE)
  }
  as.numeric(trans)
}

print.dea_undesirable <- function(x, ...) {
  cat("--- Undesirable outputs: ",
      if (x$bad_method == "translate") "Seiford-Zhu translation"
      else "bad outputs as inputs (Hailu-Veeman)", " ---\n", sep = "")
  cat(x$nbad, " undesirable output(s): ",
      paste(colnames(x$bad), collapse = ", "), "\n", sep = "")
  if (!is.null(x$bad_trans)) {
    cat("translation constant(s): ",
        paste(format(x$bad_trans, digits = 6), collapse = ", "),
        "   (the score does not depend on these; see ?dea_undesirable)\n", sep = "")
  }
  if (!is.null(x$slack_bad)) {
    cat("reducible bad output (mean slack, caller's units): ",
        paste(format(colMeans(x$slack_bad), digits = 4), collapse = ", "),
        "\n", sep = "")
  }
  cat("\n")
  NextMethod()
}
