## ---------------------------------------------------------------------------
## dea_mult() -- the multiplicative, or piecewise Cobb-Douglas, model of
## Charnes, Cooper, Seiford and Stutz (1982, 1983).
##
## EVERY OTHER MODEL IN THIS PACKAGE ENVELOPS THE DATA WITH FLAT PIECES.  The
## radial, additive, slacks-based and directional models all build a technology
## from CONVEX combinations, so the frontier is piecewise LINEAR and a peer
## group is an arithmetic average of its members.  That is an assumption about
## the production function and it is usually made silently.  This one makes the
## other assumption: the frontier is piecewise LOG-linear, so each facet is a
## Cobb-Douglas function y = A prod x_i^b_i, and a peer group is a weighted
## GEOMETRIC average.
##
##   target x_i = prod_j x_ij ^ lambda_j      rather than   sum_j lambda_j x_ij
##
## For an economist that is the natural shape: Cobb-Douglas is the functional
## form most production studies assume, and this estimates its envelope without
## assuming the exponents.  No package on CRAN implements it.
##
## IT IS THE UNWEIGHTED ADDITIVE MODEL ON LOGS, which is why it needs no new
## solver.  Taking logs turns a product into a sum and a power into a
## coefficient, so
##
##   max sum s^-  +  sum s^+
##   s.t. sum_j lambda_j log x_ij + s_i^- = log x_io
##        sum_j lambda_j log y_rj - s_r^+ = log y_ro,  sum lambda_j = 1
##
## is exactly .dea_additive() on log data.  The slacks are then LOG RATIOS --
## s_i^- = log(x_io / x_i*) -- so exp(-s_i^-) is the factor input i can be
## multiplied by, and that is reported as `ratio_x` rather than left for the
## reader to exponentiate.  `eff` is exp(-objective), in (0, 1] with 1 exactly
## on the frontier, which equals prod_i (x_i*/x_io) prod_r (y_ro/y_r*).
##
## IT CANNOT BE BUILT ON THE PUBLIC dea_add(), and the reason is worth knowing:
## .dea_check() refuses negative data, and log x is negative for every x below
## 1 -- which is most data, depending only on the unit it is recorded in.  So
## this calls .dea_additive() directly, after validating the ORIGINAL data,
## where positivity is a real requirement rather than an artefact.
##
## UNITS INVARIANCE IS THE WHOLE POINT OF THE 1983 PAPER, AND IT NEEDS vrs.
## Multiplying an input column by c adds log c to that column in log space, so
## units invariance HERE is the additive model's translation invariance THERE
## (Ali and Seiford 1990) -- and that needs the convexity row.  Measured on a
## 40-DMU design, the mean objective with one input column multiplied by
## c = 1, 10, 1000:
##
##   vrs    0.370673  0.370673  0.370673     (max change 1.3e-11)
##   crs    0.441869  0.950590  2.187686     (max change 5.02)
##   nirs   0.424713  0.949094  2.187686     (max change 5.02)
##   ndrs   0.387829  0.372169  0.370673     (max change 0.247)
##
## Under vrs the answer does not move; under crs it moves by a factor of five
## because nothing was done to the data except record an input in different
## units.  The reason is geometric: without the convexity row the log-space
## technology is a cone through the log-space origin, which is the point
## x = y = 1 -- a point with no meaning that moves whenever a unit changes.
##
## SO WHY IS crs OFFERED AT ALL?  Because it is the published 1982 model, whose
## defect is exactly what the 1983 paper's title ("Invariant multiplicative
## efficiency...") says it fixed, and reproducing a published result is a real
## use.  It WARNS rather than refuses, which is the opposite of the choice made
## in dea_undesirable(), and the distinction is deliberate: there,
## non-vrs was never a model anyone proposed, so serving it would have invented
## a quantity; here it is somebody's model and the caller is told what it costs.
##
## NOTE WHAT rts DOES NOT MEAN HERE.  These are restrictions on sum(lambda) in
## LOG space.  "vrs" does give a variable-returns Cobb-Douglas envelope, but
## "crs" does NOT give a constant-returns technology in the original space -- it
## gives the cone anchored at x = y = 1 described above.  The argument keeps its
## name because it is the same row of the same program, not because the
## returns-to-scale reading survives the transformation.
## ---------------------------------------------------------------------------

dea_mult <- function(x, y, data = NULL,
                     rts = c("vrs", "crs", "nirs", "ndrs"),
                     peers = TRUE,
                     xref = NULL, yref = NULL, dataref = NULL) {

  call <- match.call()
  rts <- .match_arg_ci(rts, c("vrs", "crs", "nirs", "ndrs"), "rts")
  t0 <- proc.time()[["elapsed"]]

  ## Validate the data the caller actually supplied. Positivity is checked here
  ## because the model takes logs, not because the additive program needs it.
  d <- .dea_data(x, y, data, xref, yref, dataref, FALSE, "dea_mult",
                 require_positive = TRUE)
  X <- d$X; Y <- d$Y; XR <- d$XR; YR <- d$YR
  dmu <- d$dmu
  n <- nrow(X); nr <- nrow(XR); p <- ncol(X); q <- ncol(Y)

  if (!identical(rts, "vrs")) {
    warning("rts = \"", rts, "\" drops the convexity row, which is what makes ",
            "this model units invariant. Multiplying one input column by 1000 ",
            "moved the mean objective from 0.442 to 2.188 under \"crs\" on a ",
            "40-DMU design, against no change at all under \"vrs\": without ",
            "that row the log-space technology is a cone anchored at ",
            "x = y = 1, a point that moves whenever a unit changes. This is ",
            "the 1982 model, and the 1983 paper exists to fix exactly this. ",
            "See ?dea_mult.", call. = FALSE)
  }

  lX <- log(X); lY <- log(Y); lXR <- log(XR); lYR <- log(YR)

  ## Straight to the internal sweep: the public dea_add() would refuse this
  ## data for having negative entries, which every value below 1 produces.
  ## Solved in log space with no column scaling -- rescaling a log column is a
  ## POWER transformation of the original variable, which is not a change of
  ## units and would change the model.
  res <- .dea_additive(lX, lY, lXR, lYR, rts, "unweighted", peers, n, nr, p, q)

  z <- res$eff                              ## the total log-slack, >= 0
  z[is.finite(z) & abs(z) < .DEA_CONSTANTS$TOL_EFF] <- 0
  eff <- exp(-z)
  names(eff) <- dmu

  lsx <- res$sx; lsy <- res$sy
  lsx[abs(lsx) < .DEA_CONSTANTS$TOL_SLACK] <- 0
  lsy[abs(lsy) < .DEA_CONSTANTS$TOL_SLACK] <- 0
  dimnames(lsx) <- list(dmu, colnames(X)); dimnames(lsy) <- list(dmu, colnames(Y))

  ## The back-transformation. These are the model's own quantities: a factor,
  ## not a difference. ratio_x is in (0, 1] and ratio_y in [1, Inf).
  rx <- exp(-lsx); ry <- exp(lsy)
  ## ... and the same thing as a difference, so that the inherited methods and
  ## anyone reading `slack_x` out of habit get the caller's units.
  sx <- X * (1 - rx); sy <- Y * (ry - 1)
  dimnames(sx) <- dimnames(rx) <- list(dmu, colnames(X))
  dimnames(sy) <- dimnames(ry) <- list(dmu, colnames(Y))

  L <- res$lambda
  if (!is.null(L)) {
    L[abs(L) < .DEA_CONSTANTS$TOL_LAMBDA] <- 0
    dimnames(L) <- list(dmu, d$ref)
  }

  .dea_report_unsolved(res$status, n, d$self, "The multiplicative program")

  structure(list(
    eff = eff, model = "multiplicative", measure = "multiplicative",
    rts = rts, orientation = "none", super = FALSE,
    lambda = L, sum_lambda = if (is.null(L)) NULL else rowSums(L),
    ## The model's own quantities first, then the derived ones.
    ratio_x = rx, ratio_y = ry,
    log_slack_x = lsx, log_slack_y = lsy,
    slack_x = sx, slack_y = sy,
    efficient = is.finite(eff) & rowSums(lsx) == 0 & rowSums(lsy) == 0,
    status = res$status,
    x = X, y = Y, xref = XR, yref = YR, self_ref = d$self,
    dmu = dmu, ref = d$ref, n = n, nref = nr, p = p, q = q,
    scaling = FALSE, slack = TRUE,
    total_time = proc.time()[["elapsed"]] - t0,
    call = call
  ), class = c("dea_mult", "dea"))
}

print.dea_mult <- function(x, ...) {
  cat("--- Multiplicative (piecewise Cobb-Douglas) model ---\n")
  cat("technology:  ", toupper(x$rts),
      if (identical(x$rts, "vrs")) "  (units invariant)" else
        "  (NOT units invariant -- see ?dea_mult)", "\n", sep = "")
  cat("DMUs: ", x$n, "   inputs: ", x$p, "   outputs: ", x$q,
      "   (", format(round(x$total_time, 3)), " sec)\n", sep = "")
  cat("\nefficiency, exp(-total log slack); 1 = on the frontier\n")
  print(round(summary(x$eff[is.finite(x$eff)]), 4))
  cat("\nefficient DMUs: ", sum(x$efficient, na.rm = TRUE), " of ", x$n, "\n", sep = "")
  cat("\nmean factor each input could be multiplied by\n")
  print(round(colMeans(x$ratio_x, na.rm = TRUE), 4))
  cat("\nmean factor each output could be multiplied by\n")
  print(round(colMeans(x$ratio_y, na.rm = TRUE), 4))
  invisible(x)
}

summary.dea_mult <- function(object, ...) {
  print(object)
  cat("\nThe peer target is a weighted GEOMETRIC mean of the peers, not an\n",
      "arithmetic one: target x_i = prod_j x_ij^lambda_j.\n", sep = "")
  cat("\nfirst rows: efficiency and the input factors\n")
  print(utils::head(round(cbind(eff = object$eff, object$ratio_x), 4), 10))
  invisible(object)
}
