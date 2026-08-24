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
