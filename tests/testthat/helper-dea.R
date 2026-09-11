## Shared fixtures. Deliberately small: these tests are about correctness and
## contracts, not about statistical power, and R CMD check has a time budget.

toy <- function(n = 40, p = 2, q = 1, seed = 1, returns = 1) {
  s <- dea_sim(n, p = p, q = q, returns = returns, seed = seed)
  list(x = s$x, y = s$y, theta = s$theta, phi = s$phi, data = s$data)
}

## A hand-worked one-input one-output example. The constant-returns frontier is
## the ray through the DMU with the largest y/x, so every score is computable
## with a pencil -- which is the point: it does not depend on any package
## agreeing with any other.
tiny <- list(
  x = matrix(c(2, 4, 3, 5, 6), ncol = 1, dimnames = list(NULL, "x")),
  y = matrix(c(1, 3, 2, 4, 3), ncol = 1, dimnames = list(NULL, "y")))
## max y/x = 4/5 = 0.8, attained by DMU 4.
tiny_crs_in <- (tiny$y[, 1] / tiny$x[, 1]) / 0.8

## ---------------------------------------------------------------------------
## Reference values from independent implementations.
##
## `reference/*.csv` holds what Benchmarking 0.33 and DJL 3.9 computed on these
## same fixtures, recorded once by ../../tools/make_reference_values.R. The
## comparison is against those recorded numbers rather than against a live call,
## for three reasons set out in full in that script:
##
##   * a skip_if_not_installed() check is skipped on most of CRAN's machines,
##     so the old form of this test usually did not run at all;
##   * a DEA package should not make other DEA packages a dependency, even a
##     suggested one;
##   * "agrees with what Benchmarking computed on a stated date from a stated
##     version" stays true, where "agrees with whatever it does now" drifts.
##
## The live, current comparison against these and five other packages still
## happens in ../../../horserace/, which is not shipped.
##
## These are not hand-checked constants and are not a substitute for the
## hand-worked `tiny` fixture above or for recovery of dea_sim()'s closed-form
## truth. They catch coding errors. Two implementations sharing a
## misunderstanding would agree perfectly, which is exactly why all three kinds
## of check exist.
## ---------------------------------------------------------------------------
ref_values <- function(name) {
  utils::read.csv(test_path("reference", name), stringsAsFactors = FALSE)
}
