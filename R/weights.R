## ---------------------------------------------------------------------------
## Weight restrictions: assurance regions and the cone-ratio model.
##
## An unrestricted DEA score lets every DMU choose the prices that flatter it
## most, including prices that value an input at zero. Chapter 3 s4 calls
## fixing that "perhaps the most significant of the proposed extensions to
## DEA", and the machinery is already here: a restriction is extra ROWS on the
## multiplier program, which R/lp.R has built since 1.0.0.
##
## WHAT CHANGES, AND IT IS NOT SMALL: once restrictions are imposed the
## envelopment and multiplier forms STOP BEING INTERCHANGEABLE. Strong duality
## still holds, but the dual of the restricted multiplier program is not the
## ordinary envelopment program -- it acquires extra terms, one per
## restriction, which enter the input and output constraints exactly where
## slacks do:
##
##   X'lambda + A_x z <= theta x_o,     Y'lambda - A_y z >= y_o,    z >= 0
##
## Those `A z` terms are the RESIDUE. They are not slacks, they are not peers'
## contributions, and a package that reports them as either is reporting a
## projection onto a point that is not in the technology. So dea_wr() returns
## the score and the weights and does NOT return lambda, peers or slacks, and
## says so rather than returning something plausible. That is the whole reason
## this is a separate entry point instead of an argument to dea().
##
## THE SCORE IS THE MULTIPLIER OPTIMUM. With no restrictions it equals dea()'s
## to solver precision, and the tests assert that -- it is the one case where
## the two forms can still be compared, and it is what makes the rest
## believable.
##
## UNITS. The solver works on columns divided by their means, so a weight the
## solver reports applies to x_i / sx_i and the caller's weight is v_i / sx_i
## (see .dea_multipliers()). A restriction is written in the CALLER's units, so
## its coefficients are divided by the same factors on the way in. Getting this
## wrong would not error: it would quietly impose a different restriction, on
## data whose scale nobody stated. test-weights.R fixes it by requiring the
## restricted score to be identical with scaling on and off.
## ---------------------------------------------------------------------------

## A restriction set is a list of rows, each `a'v + b'u  (type)  rhs`, with the
## coefficients in the caller's units. The constructors below all produce that
## shape; they exist so a caller writes what they mean rather than assembling
## a matrix and hoping the column order is the one the package used.
.wr_new <- function(v, u, type, rhs, label) {
  structure(list(rows = list(list(v = v, u = u, type = type, rhs = rhs,
                                  label = label))),
            class = "dea_restriction")
}

wr_bound <- function(on = c("v", "u"), index, lower = NA, upper = NA) {
  on <- .match_arg_ci(on, c("v", "u"), "on")
  index <- .wr_index(index, "index")
  if (is.na(lower) && is.na(upper))
    stop("wr_bound() needs at least one of `lower` and `upper`.", call. = FALSE)
  if (!is.na(lower) && !is.na(upper) && lower > upper)
    stop("wr_bound(): lower = ", lower, " is above upper = ", upper, ".",
         call. = FALSE)
  out <- list()
  side <- function(coef, type, rhs, lab)
    .wr_new(v = if (on == "v") coef else NULL,
            u = if (on == "u") coef else NULL, type = type, rhs = rhs,
            label = lab)
  if (!is.na(lower))
    out <- c(out, list(side(c(index, 1), ">=", lower,
                            sprintf("%s[%d] >= %g", on, index, lower))))
  if (!is.na(upper))
    out <- c(out, list(side(c(index, 1), "<=", upper,
                            sprintf("%s[%d] <= %g", on, index, upper))))
  do.call(c, out)
}

## L <= w[num] / w[den] <= U, written without dividing. A ratio constraint
## divided out is nonlinear; multiplied through it is a plain row, which is why
## assurance regions are cheap and why they are stated this way in Thompson et
## al. (1986).
wr_ratio <- function(on = c("v", "u"), numerator, denominator,
                     lower = NA, upper = NA) {
  on <- .match_arg_ci(on, c("v", "u"), "on")
  numerator <- .wr_index(numerator, "numerator")
  denominator <- .wr_index(denominator, "denominator")
  if (numerator == denominator)
    stop("wr_ratio(): numerator and denominator are the same index (",
         numerator, "), so the ratio is 1 whatever the weights are.",
         call. = FALSE)
  if (is.na(lower) && is.na(upper))
    stop("wr_ratio() needs at least one of `lower` and `upper`.", call. = FALSE)
  if (!is.na(lower) && !is.na(upper) && lower > upper)
    stop("wr_ratio(): lower = ", lower, " is above upper = ", upper, ".",
         call. = FALSE)
  side <- function(coef, lab)
    .wr_new(v = if (on == "v") coef else NULL,
            u = if (on == "u") coef else NULL, type = "<=", rhs = 0, label = lab)
  out <- list()
  ## w[num] - U w[den] <= 0  and  L w[den] - w[num] <= 0
  if (!is.na(upper))
    out <- c(out, list(side(rbind(c(numerator, 1), c(denominator, -upper)),
                            sprintf("%s[%d]/%s[%d] <= %g", on, numerator, on,
                                    denominator, upper))))
  if (!is.na(lower))
    out <- c(out, list(side(rbind(c(numerator, -1), c(denominator, lower)),
                            sprintf("%s[%d]/%s[%d] >= %g", on, numerator, on,
                                    denominator, lower))))
  do.call(c, out)
}

## The general case of Wong and Beasley (1990): any linear inequality in the
## weights. `v` and `u` are full-length coefficient vectors (or NULL for none).
wr_linear <- function(v = NULL, u = NULL, type = c("<=", ">=", "="), rhs = 0) {
  type <- match.arg(type)
  if (is.null(v) && is.null(u))
    stop("wr_linear() needs coefficients on `v`, on `u`, or on both.",
         call. = FALSE)
  if (!is.numeric(rhs) || length(rhs) != 1L || !is.finite(rhs))
    stop("wr_linear(): `rhs` must be one finite number.", call. = FALSE)
  chk <- function(z, what) {
    if (is.null(z)) return(NULL)
    if (!is.numeric(z) || anyNA(z) || any(!is.finite(z)))
      stop("wr_linear(): `", what, "` must be finite numbers.", call. = FALSE)
    as.numeric(z)
  }
  .wr_new(v = chk(v, "v"), u = chk(u, "u"), type = type, rhs = rhs,
          label = "linear")
}

.wr_index <- function(i, what) {
  if (!is.numeric(i) || length(i) != 1L || is.na(i) || i != as.integer(i) || i < 1)
    stop("`", what, "` must be a single positive whole number (a column ",
         "position), got ", paste(utils::head(i, 3), collapse = ", "), ".",
         call. = FALSE)
  as.integer(i)
}

c.dea_restriction <- function(...) {
  parts <- list(...)
  ok <- vapply(parts, inherits, logical(1), "dea_restriction")
  if (!all(ok))
    stop("Weight restrictions combine only with other weight restrictions; ",
         "element ", which(!ok)[1L], " is a ", class(parts[!ok][[1L]])[1L], ".",
         call. = FALSE)
  structure(list(rows = do.call(c, lapply(parts, `[[`, "rows"))),
            class = "dea_restriction")
}

print.dea_restriction <- function(x, ...) {
  cat("--- Weight restrictions (", length(x$rows), ") ---\n", sep = "")
  for (r in x$rows) cat("  ", r$label, "\n", sep = "")
  invisible(x)
}

## Expand the restriction set into a k x (p+q) coefficient matrix, then divide
## each column by its scale factor so that the rows mean in the SOLVER's
## coordinates what the caller wrote in theirs.
.wr_rows <- function(rs, p, q, sx, sy) {
  if (is.null(rs)) return(NULL)
  if (!inherits(rs, "dea_restriction"))
    stop("`restrictions` must be built by wr_bound(), wr_ratio() or ",
         "wr_linear(), combined with c(); got a ", class(rs)[1L], ".",
         call. = FALSE)
  k <- length(rs$rows)
  mat <- matrix(0, k, p + q)
  type <- character(k); rhs <- numeric(k)
  for (j in seq_len(k)) {
    r <- rs$rows[[j]]
    mat[j, ] <- .wr_row_one(r, p, q)
    type[j] <- r$type; rhs[j] <- r$rhs
  }
  if (any(rowSums(abs(mat)) == 0))
    stop("Restriction ", which(rowSums(abs(mat)) == 0)[1L],
         " has no non-zero coefficient on any weight.", call. = FALSE)
  mat <- sweep(mat, 2L, c(sx, sy), "/")
  list(mat = mat, type = type, rhs = rhs,
       label = vapply(rs$rows, `[[`, character(1), "label"))
}

.wr_row_one <- function(r, p, q) {
  row <- numeric(p + q)
  put <- function(z, off, n, what) {
    if (is.null(z)) return(invisible(NULL))
    if (is.matrix(z)) {
      for (m in seq_len(nrow(z))) {
        i <- z[m, 1L]
        if (i > n) stop("A restriction refers to ", what, "[", i, "] but there ",
                        "are only ", n, ".", call. = FALSE)
        row[off + i] <<- row[off + i] + z[m, 2L]
      }
    } else if (length(z) == 2L && z[1L] == as.integer(z[1L]) && z[1L] >= 1 &&
               z[1L] <= n && !is.matrix(z)) {
      ## the (index, coefficient) shorthand the bound/ratio constructors use
      row[off + as.integer(z[1L])] <<- row[off + as.integer(z[1L])] + z[2L]
    } else {
      if (length(z) != n)
        stop("A restriction gives ", length(z), " coefficient(s) on ", what,
             " but there are ", n, ".", call. = FALSE)
      row[off + seq_len(n)] <<- row[off + seq_len(n)] + z
    }
    invisible(NULL)
  }
  put(r$v, 0L, p, "v")
  put(r$u, p, q, "u")
  row
}

## ---------------------------------------------------------------------------
## dea_wr() -- the restricted radial model.
## ---------------------------------------------------------------------------

dea_wr <- function(x, y, restrictions, data = NULL,
                   rts = c("vrs", "crs", "nirs", "ndrs"),
                   orientation = c("in", "out"),
                   scaling = TRUE,
                   xref = NULL, yref = NULL, dataref = NULL) {
  call <- match.call()
  rts <- .match_arg_ci(rts, c("vrs", "crs", "nirs", "ndrs"), "rts")
  orientation <- .match_arg_ci(orientation, c("in", "out"), "orientation")
  t0 <- proc.time()[["elapsed"]]

  d <- .dea_data(x, y, data, xref, yref, dataref, FALSE, "dea_wr")
  X <- d$X; Y <- d$Y; XR <- d$XR; YR <- d$YR
  n <- nrow(X); nr <- nrow(XR); p <- ncol(X); q <- ncol(Y)

  sc  <- .dea_scale(XR, YR, scaling)
  XRs <- sc$X; YRs <- sc$Y
  Xs  <- sweep(X, 2L, sc$sx, "/"); Ys <- sweep(Y, 2L, sc$sy, "/")
  R <- .wr_rows(restrictions, p, q, sc$sx, sc$sy)

  M <- .lp_mult_build(XRs, YRs, rts, orientation, restrict = R)
  eff <- numeric(n); st <- integer(n)
  V <- matrix(NA_real_, n, p); U <- matrix(NA_real_, n, q); u0 <- numeric(n)
  for (o in seq_len(n)) {
    r <- .lp_mult_at(M, Xs, Ys, o)
    if (!r$status %in% c(0L, 1L)) {
      ## A fresh object, for the same reason the envelopment sweep needs one:
      ## inherited basis state, not a statement about the data. An infeasible
      ## program (2) is a statement about the data and is not retried.
      if (!identical(r$status, 2L)) {
        M2 <- .lp_mult_build(XRs, YRs, rts, orientation, restrict = R)
        r <- .lp_mult_at(M2, Xs, Ys, o)
      }
    }
    st[o] <- r$status
    eff[o] <- if (r$status %in% c(0L, 1L)) r$eff else NA_real_
    V[o, ] <- r$v; U[o, ] <- r$u; u0[o] <- r$u0
  }

  .wr_report_unsolved(st, n, length(R$label))
  V <- sweep(V, 2L, sc$sx, "/")
  U <- sweep(U, 2L, sc$sy, "/")
  eff[is.finite(eff) & abs(eff - 1) < .DEA_CONSTANTS$TOL_EFF] <- 1
  names(eff) <- names(u0) <- d$dmu
  dimnames(V) <- list(d$dmu, colnames(X)); dimnames(U) <- list(d$dmu, colnames(Y))

  structure(list(
    eff = eff, model = "radial", rts = rts, orientation = orientation,
    super = FALSE, restricted = TRUE, restrictions = restrictions,
    restriction_labels = R$label,
    ## Deliberately absent: lambda, peers and slacks. The dual of a restricted
    ## multiplier program is not the ordinary envelopment program -- it carries
    ## one extra non-negative term per restriction, in exactly the place a
    ## slack sits -- so a lambda read off it does not project onto a point of
    ## the technology and a "slack" read off it is not one. See the header.
    lambda = NULL, sum_lambda = NULL, slack_x = NULL, slack_y = NULL,
    v = V, u = U, u0 = u0, multipliers = TRUE, mult_status = st,
    efficient = is.finite(eff) & eff == 1, status = st,
    x = X, y = Y, xref = XR, yref = YR, self_ref = d$self,
    dmu = d$dmu, ref = d$ref, n = n, nref = nr, p = p, q = q,
    scaling = scaling, slack = FALSE,
    total_time = proc.time()[["elapsed"]] - t0,
    call = call
  ), class = c("dea_wr", "dea"))
}

## Restrictions are the one thing here that can make a DEA program genuinely
## infeasible -- ask for weights no price vector can satisfy and there is no
## answer, which is information rather than a failure. So the message names the
## restrictions as the likely cause instead of the generic conditioning advice.
.wr_report_unsolved <- function(st, n, k) {
  bad <- !st %in% c(0L, 1L)
  if (!any(bad)) return(invisible(NULL))
  infeas <- sum(st == 2L)
  if (infeas) {
    warning(infeas, " of ", n, " DMU(s) had an INFEASIBLE restricted program ",
            "and are reported as NA. With ", k, " restriction(s) this usually ",
            "means no price vector satisfies all of them at once -- check them ",
            "against each other before checking the data. An unrestricted fit ",
            "of the same data cannot be infeasible, so the restrictions are ",
            "the place to look.", call. = FALSE)
  }
  num <- sum(bad & st != 2L)
  if (num) {
    warning(num, " of ", n, " DMU(s) FAILED NUMERICALLY (lpSolve status other ",
            "than 0, 1 or 2) and are reported as NA.", call. = FALSE)
  }
  invisible(NULL)
}

print.dea_wr <- function(x, ...) {
  cat("--- DEA with weight restrictions (radial, ",
      if (x$orientation == "in") "input" else "output", " oriented, rts = ",
      x$rts, ") ---\n", sep = "")
  cat("DMUs: ", x$n, "   inputs: ", x$p, "   outputs: ", x$q,
      "   restrictions: ", length(x$restriction_labels), "\n", sep = "")
  for (l in x$restriction_labels) cat("  ", l, "\n", sep = "")
  e <- x$eff[is.finite(x$eff)]
  cat("\nefficiency\n"); print(round(summary(e), 4))
  cat("\non the restricted frontier: ", sum(x$efficient, na.rm = TRUE), " of ",
      x$n, "\n", sep = "")
  if (anyNA(x$eff))
    cat("infeasible or unsolved: ", sum(is.na(x$eff)), "\n", sep = "")
  cat("\nNo peers or slacks: with restrictions the envelopment dual carries a\n",
      "residue term per restriction, so neither is the quantity it would be\n",
      "in an unrestricted fit. ?dea_wr explains.\n", sep = "")
  invisible(x)
}

peers.dea_wr <- function(object, ...) {
  stop("A restricted fit has no peers to report. The dual of the restricted ",
       "multiplier program is not the ordinary envelopment program: it carries ",
       "one non-negative residue term per restriction, sitting where a slack ",
       "sits, so the lambda it yields does not project the DMU onto a point of ",
       "the technology. Fit the same data with dea() if unrestricted peers are ",
       "what is wanted, and do not mix the two in one table.", call. = FALSE)
}

slacks.dea_wr <- function(object, ...) {
  stop("A restricted fit has no slacks to report -- see peers() on the same ",
       "object for why. The residue terms the restricted dual carries occupy ",
       "the same place in the constraints as slacks and are not slacks.",
       call. = FALSE)
}

## ---------------------------------------------------------------------------
## dea_cone() -- the cone-ratio model of Charnes, Cooper, Wei and Huang (1989),
## which needs no new linear program at all.
##
## When the admissible weight sets are polyhedral cones in SUM form,
## V = {A'alpha : alpha >= 0} and U = {B'gamma : gamma >= 0}, the cone-ratio
## model is exactly a CCR model on transformed data: Xt = X A', Yt = Y B'
## (chapter 3 s4). Substituting v = A'alpha into v'x gives alpha'(A x), so the
## program in alpha on the transformed data IS the program in v on the
## original -- and alpha >= 0 is what an ordinary DEA program already imposes.
##
## So this is a preprocessing step and a call to dea(). What it is NOT is
## interchangeable with dea_wr(): a cone in sum form becomes a set of
## INEQUALITIES only through its polar, which is a separate computation and is
## not generally as small. Both are provided because both appear in the
## literature and neither is a wrapper for the other.
## ---------------------------------------------------------------------------

dea_cone <- function(x, y, A = NULL, B = NULL, data = NULL,
                     rts = c("crs", "vrs", "nirs", "ndrs"),
                     orientation = c("in", "out"), ...) {
  call <- match.call()
  rts <- .match_arg_ci(rts, c("crs", "vrs", "nirs", "ndrs"), "rts")
  d <- .dea_data(x, y, data, NULL, NULL, NULL, FALSE, "dea_cone")
  A <- .cone_matrix(A, ncol(d$X), "A", "input")
  B <- .cone_matrix(B, ncol(d$Y), "B", "output")
  fit <- dea(d$X %*% t(A), d$Y %*% t(B), rts = rts, orientation = orientation,
             ...)
  fit$cone <- list(A = A, B = B)
  fit$call <- call
  fit
}

.cone_matrix <- function(M, k, what, side) {
  if (is.null(M)) return(diag(k))
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  if (ncol(M) != k)
    stop("`", what, "` must have one column per ", side, " (", k, "), got ",
         ncol(M), ".", call. = FALSE)
  if (any(!is.finite(M)))
    stop("`", what, "` contains missing or non-finite values.", call. = FALSE)
  if (nrow(M) < k) {
    ## Fewer generators than dimensions is a cone with empty interior: every
    ## weight vector is forced onto a lower-dimensional face, which is legal
    ## and is almost never what someone meant to type.
    warning("`", what, "` has ", nrow(M), " generator(s) for ", k, " ", side,
            "(s), so the admissible weights lie in a ", nrow(M),
            "-dimensional subspace. That is a much stronger restriction than ",
            "it looks; check it is intended.", call. = FALSE)
  }
  M
}
