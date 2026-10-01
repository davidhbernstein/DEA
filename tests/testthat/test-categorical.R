## ---------------------------------------------------------------------------
## Categorical environmental variables as nested reference sets.
##
## No reference implementation to pin against -- deaR has no categorical model
## and the Banker-Morey MIP is deliberately not built -- so what is checked is
## the structure the model IS: who may be a peer, and the two identities that
## follow from it. These are provable from the construction, which is what makes
## them worth asserting rather than a cross-check.
## ---------------------------------------------------------------------------

cat_fixture <- function(n = 90, seed = 21) {
  set.seed(seed)
  x <- matrix(runif(n * 2, 2, 10), n, 2,
              dimnames = list(NULL, c("staff", "budget")))
  env <- rep(1:3, length.out = n)
  boost <- c(1.0, 1.15, 1.32)[env]
  y <- matrix((x[, 1] * x[, 2])^0.4 * boost * exp(-abs(rnorm(n, 0, 0.2))), n, 1,
              dimnames = list(NULL, "output"))
  list(x = x, y = y, env = env,
       cf = factor(c("deprived", "average", "affluent")[env],
                   levels = c("deprived", "average", "affluent"), ordered = TRUE))
}

test_that("a DMU is never given a peer from a more favourable category", {
  ## This is the entire model. If it fails, nothing else matters.
  f <- cat_fixture()
  r <- dea_categorical(f$x, f$y, f$cf)
  ki <- as.integer(r$category)
  forbidden <- outer(ki, ki, function(a, b) b > a)   ## peer better placed
  expect_true(all(r$lambda[forbidden] == 0))
  expect_true(any(r$lambda[!forbidden] > 0))         ## and the rest is used
  ## The reference set is the cumulative count, by construction.
  expect_equal(r$reference_size, cumsum(table(f$cf))[ki], ignore_attr = TRUE)
})

test_that("the two ends of the nest match the fits they must equal", {
  f <- cat_fixture()
  for (ori in c("in", "out")) {
    r <- dea_categorical(f$x, f$y, f$cf, orientation = ori)
    ki <- as.integer(r$category)

    ## The TOP category is scored against everybody, so its scores are the
    ## pooled ones exactly.
    top <- which(ki == 3L)
    expect_equal(unname(r$eff[top]), unname(r$pooled_eff[top]), tolerance = 1e-9)
    expect_equal(unname(r$category_effect[top]), rep(0, length(top)),
                 tolerance = 1e-9)

    ## The BOTTOM category is scored against itself alone, so its scores are an
    ## ordinary fit on that subset.
    bot <- which(ki == 1L)
    solo <- suppressWarnings(dea(f$x[bot, ], f$y[bot, ], orientation = ori))
    expect_equal(unname(r$eff[bot]), unname(solo$eff), tolerance = 1e-9)
  }
})

test_that("nesting can only help, never hurt", {
  ## A smaller reference set is a smaller technology, so the categorical score
  ## is never the worse of the two on the common (0, 1] scale -- in either
  ## orientation, which is why the effect is defined on that scale.
  f <- cat_fixture()
  for (ori in c("in", "out")) {
    r <- dea_categorical(f$x, f$y, f$cf, orientation = ori)
    expect_true(all(r$category_effect >= -1e-9))
    expect_gt(mean(r$category_effect), 0)
    ## and the gain is larger the worse the environment
    m <- tapply(r$category_effect, r$category, mean)
    expect_true(m[[1]] > m[[2]] && m[[2]] >= m[[3]])
  }
})

test_that("a numeric rank and an ordered factor are the same thing", {
  f <- cat_fixture()
  a <- dea_categorical(f$x, f$y, f$cf)
  b <- dea_categorical(f$x, f$y, f$env)
  expect_equal(unname(a$eff), unname(b$eff))
  expect_equal(a$reference_size, b$reference_size)
  ## The labels differ; the nesting does not.
  expect_equal(a$category_levels, c("deprived", "average", "affluent"))
  expect_equal(b$category_levels, c("1", "2", "3"))
})

test_that("an unordered grouping is refused rather than ordered alphabetically", {
  ## The failure this function exists to prevent: alphabetical order is not an
  ## ordering of environments, and inventing one answers a question nobody asked.
  f <- cat_fixture()
  plain <- factor(as.character(f$cf))
  expect_error(dea_categorical(f$x, f$y, plain), "unordered factor")
  expect_error(dea_categorical(f$x, f$y, plain), "separate analysis|separate one")
  expect_error(dea_categorical(f$x, f$y, as.character(f$cf)), "ordered factor")
})

test_that("the other refusals fire", {
  f <- cat_fixture(60)
  expect_error(dea_categorical(f$x, f$y, f$cf[-1]), "length")
  bad <- f$cf; bad[3] <- NA
  expect_error(dea_categorical(f$x, f$y, bad), "missing value")
  expect_error(dea_categorical(f$x, f$y, rep(1L, 60)), "only one distinct")
  expect_error(dea_categorical(f$x, f$y, "nope", data = as.data.frame(f$x)),
               "not found")
})

test_that("a category whose reference set is too small is reported", {
  ## Nesting costs sample size and the cost falls on the lowest category. The
  ## scores stay finite and plausible, so the warning is the only notice there
  ## is.
  ## 3(p+q) = 9 here, and the lowest category has 6 members to be scored
  ## against -- comfortable pooled at 30 DMUs, well under the rule once split.
  set.seed(4)
  x <- matrix(runif(30 * 2, 2, 10), 30, 2)
  y <- matrix(runif(30, 1, 5), 30, 1)
  cf <- factor(rep(c("a", "b"), times = c(6, 24)), levels = c("a", "b"),
               ordered = TRUE)
  expect_warning(r <- dea_categorical(x, y, cf), "rule of thumb")
  expect_warning(dea_categorical(x, y, cf), "falls on the lowest")
  expect_equal(min(r$reference_size), 6L)
  ## and more of that category is efficient than would be pooled -- which is
  ## the risk the warning is about. A fixed threshold here would be a guess;
  ## the comparison against the pooled fit on the same DMUs is the property.
  lo <- as.integer(r$category) == 1L
  expect_gte(mean(r$eff[lo] == 1), mean(r$pooled_eff[lo] == 1))
  expect_gt(mean(r$eff[lo] == 1), 0)
})

test_that("the methods run and the no-effect case is called out", {
  f <- cat_fixture(60)
  r <- dea_categorical(f$x, f$y, f$cf)
  expect_output(print(r), "nested reference sets")
  expect_output(print(r), "LEAST favourable")
  expect_output(summary(r), "not genuinely nested")
  expect_equal(nobs(r), 60L)
  expect_s3_class(r, "dea")

  ## When the environment does no work, the print method says so rather than
  ## leaving a column of zeros to be noticed.
  set.seed(8)
  n <- 60
  xx <- matrix(runif(n * 2, 2, 10), n, 2)
  yy <- matrix((xx[, 1] * xx[, 2])^0.4, n, 1)        ## everyone on one frontier
  cc <- factor(rep(c("a", "b"), length.out = n), levels = c("a", "b"),
               ordered = TRUE)
  z <- dea_categorical(xx, yy, cc)
  expect_lt(max(abs(z$category_effect)), 1e-8)
  expect_output(print(z), "doing no work here")
})
