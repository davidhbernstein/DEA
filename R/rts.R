## ---------------------------------------------------------------------------
## dea_rts() -- scale efficiency and the returns-to-scale classification.
##
## Three technologies, one sweep, because the classification is a COMPARISON
## and reading it off one fit is not possible:
##
##   scale efficiency   SE = theta_crs / theta_vrs   (input orientation)
##                      SE = phi_vrs   / phi_crs     (output orientation)
##
## SE is 1 when the DMU operates at the most productive scale size, and the
## shortfall is the part of its distance to the CRS frontier that comes from
## being the wrong SIZE rather than from being badly run.
##
## THE CLASSIFICATION.  If theta_crs = theta_vrs the DMU is already at the most
## productive scale size.  Otherwise the non-increasing technology decides.
##
## Getting this the right way round takes one careful sentence.  NIRS allows
## sum(lambda) <= 1, i.e. it permits scaling a reference point DOWN but not up.
## A DMU in the INCREASING-returns region is smaller than the most productive
## scale size, so the constant-returns frontier above it is reached by scaling
## the MPSS peer DOWN -- which NIRS permits.  There, NIRS coincides with CRS.
## A DMU in the DECREASING-returns region is larger than the MPSS, so reaching
## the CRS frontier would mean scaling a peer UP, which NIRS forbids; there it
## coincides with VRS.  Hence
##
##   theta_crs  = theta_vrs                    constant returns
##   theta_nirs = theta_crs  != theta_vrs      increasing returns
##   theta_nirs = theta_vrs != theta_crs       decreasing returns
##
## The intuition runs the opposite way to the label -- "non-increasing returns"
## agreeing with the constant-returns frontier is what marks the INCREASING
## region -- and an earlier version of this file had it inverted.  What caught
## it was the check below: the independent sum(lambda) criterion agreed with
## the inverted rule on 1.7% of DMUs, which is the signature of a sign error
## rather than of a tolerance being wrong.  See Fare, Grosskopf and Lovell
## (1994) and Coelli et al. (1998, section 6.5).
##
## That alternative -- read sum(lambda) from the CRS program and call it IRS
## below 1, DRS above -- gives the same answer when the CRS optimum is unique
## and an arbitrary one when it is not, which is common.  sum(lambda) is
## reported here so the two can be compared, but the verdict uses the NIRS
## rule.  On the reference design the two agree on 100% of DMUs under either
## orientation, which is what makes the rule above checkable at all.
##
## ORIENTATION CHANGES THE ANSWER, AND IS MEANT TO.  Returns to scale are a
## property of the frontier POINT a DMU is benchmarked against, not of the
## DMU's own coordinates.  An output-oriented projection holds inputs fixed, so
## the verdict is about the DMU's own scale.  An input-oriented projection
## moves inputs down, and for a badly inefficient DMU it can land below the
## most productive scale size even though the DMU itself sits above it -- so
## the same DMU is classified DRS by the output sweep and IRS by the input one.
## On a 60-DMU design, 34 of them flip this way.  Both classifications are
## correct about their own projection: the input one matches the projected
## input level against the MPSS on 100% of DMUs, and the output one matches the
## observed input level on 100%.  Pick the orientation that matches the
## decision being made, and do not compare the two tables row by row.
## ---------------------------------------------------------------------------

dea_rts <- function(x, y, data = NULL,
                    orientation = c("in", "out"),
                    tol = 1e-6,
                    scaling = TRUE) {

  call <- match.call()
  orientation <- .match_arg_ci(orientation, c("in", "out"), "orientation")
  t0 <- proc.time()[["elapsed"]]

  X <- .dea_matrix(x, data, "x")
  Y <- .dea_matrix(y, data, "y")
  .dea_check(X, Y)
  dmu <- .dea_dmu_names(X, Y)
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  sc <- .dea_scale(X, Y, scaling)

  one <- function(r) .dea_radial(sc$X, sc$Y, sc$X, sc$Y, r, orientation,
                                 super = FALSE, slack = FALSE, peers = TRUE,
                                 n, n, p, q)
  f_crs  <- one("crs"); f_vrs <- one("vrs"); f_nirs <- one("nirs")
  snap <- function(v) { v[is.finite(v) & abs(v - 1) < .DEA_CONSTANTS$TOL_EFF] <- 1; v }
  e_crs <- snap(f_crs$eff); e_vrs <- snap(f_vrs$eff); e_nirs <- snap(f_nirs$eff)

  se <- if (orientation == "in") e_crs / e_vrs else e_vrs / e_crs
  se[!is.finite(se)] <- NA_real_

  cls <- ifelse(abs(e_crs - e_vrs) < tol, "crs",
         ifelse(abs(e_nirs - e_crs) < tol, "irs", "drs"))
  cls[!is.finite(e_crs) | !is.finite(e_vrs) | !is.finite(e_nirs)] <- NA_character_

  out <- data.frame(dmu = dmu, crs = e_crs, vrs = e_vrs, nirs = e_nirs,
                    scale_eff = se, rts = cls,
                    sum_lambda_crs = f_crs$sum_lambda,
                    row.names = NULL, stringsAsFactors = FALSE)
  structure(list(table = out, orientation = orientation, tol = tol,
                 n = n, p = p, q = q, dmu = dmu,
                 total_time = proc.time()[["elapsed"]] - t0,
                 call = call),
            class = "dea_rts")
}

print.dea_rts <- function(x, ...) {
  cat("--- Returns to scale (", if (x$orientation == "in") "input" else "output",
      " orientation) ---\n", sep = "")
  cat("DMUs: ", x$n, "   inputs: ", x$p, "   outputs: ", x$q, "\n", sep = "")
  tb <- table(factor(x$table$rts, levels = c("irs", "crs", "drs")), useNA = "ifany")
  cat("\nclassification\n"); print(tb)
  se <- x$table$scale_eff
  cat("\nscale efficiency\n")
  print(round(summary(se[is.finite(se)]), 4))
  cat("\nat most productive scale size (SE = 1): ",
      sum(is.finite(se) & abs(se - 1) < x$tol), " of ", x$n, "\n", sep = "")
  invisible(x)
}

summary.dea_rts <- function(object, ...) {
  print(object)
  cat("\nfirst rows\n")
  print(utils::head(object$table, 10))
  invisible(object)
}
