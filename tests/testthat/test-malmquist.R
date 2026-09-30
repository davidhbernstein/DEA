## ---------------------------------------------------------------------------
## The panel layer and the Malmquist index.
## ---------------------------------------------------------------------------

panel_fixture <- function(n = 20, seed0 = 21, seed1 = 22) {
  d0 <- toy(n, 2, 1, seed0)
  d1 <- toy(n, 2, 1, seed1)
  list(x = rbind(d0$x, d1$x * 0.92), y = rbind(d0$y, d1$y * 1.12),
       id = rep(seq_len(n), 2L), period = rep(1:2, each = n))
}

test_that("a panel validates its keys and reports its own shape", {
  f <- panel_fixture(8)
  pn <- dea_panel(f$x, f$y, f$id, f$period)
  expect_s3_class(pn, "dea_panel")
  expect_equal(pn$n_id, 8L)
  expect_equal(pn$n_period, 2L)
  expect_true(pn$balanced)
  expect_equal(nobs(pn), 16L)

  ## An unbalanced panel is an ordinary case, not an error.
  drop <- !(f$id == 1 & f$period == 2)
  un <- dea_panel(f$x[drop, ], f$y[drop, ], f$id[drop], f$period[drop])
  expect_false(un$balanced)

  ## The things that ARE errors.
  expect_error(dea_panel(f$x, f$y, f$id, rep(1, length(f$id))), "at least two")
  expect_error(dea_panel(f$x, f$y, rep(1L, length(f$id)), f$period), "duplicated")
  expect_error(dea_panel(f$x, f$y, f$id[-1], f$period), "length")
  bad <- f$id; bad[3] <- NA
  expect_error(dea_panel(f$x, f$y, bad, f$period), "contains NA")
})

test_that("period ORDER comes from the data, and a factor is how to set it", {
  ## Sorting is right for numbers and wrong for labels such as Q1, Q10, Q2.
  ## Getting it wrong silently reverses every index built on the panel, so the
  ## object is required to say which order it settled on.
  f <- panel_fixture(6)
  lab <- c("Q9", "Q10")[f$period]
  expect_equal(dea_panel(f$x, f$y, f$id, lab)$periods, c("Q10", "Q9"))
  expect_equal(
    dea_panel(f$x, f$y, f$id, factor(lab, levels = c("Q9", "Q10")))$periods,
    c("Q9", "Q10"))

  ## And reversing the order really does invert the index, rather than being
  ## cosmetic.
  a <- dea_malmquist(f$x, f$y, f$id, factor(lab, levels = c("Q9", "Q10")))
  b <- dea_malmquist(f$x, f$y, f$id, factor(lab, levels = c("Q10", "Q9")))
  expect_equal(a$table$malmquist, 1 / b$table$malmquist, tolerance = 1e-8)
})

test_that("the Malmquist index matches productivity::malm in both orientations", {
  ## The pinned values are productivity 1.1.0's. Benchmarking 0.33 agrees too
  ## on input orientation and gives the RECIPROCAL on output; see
  ## tools/make_reference_values.R and R/malmquist.R for why this package
  ## follows the Shephard-distance convention, under which M > 1 is growth
  ## whichever side it is measured from.
  f <- panel_fixture(20)
  ref <- ref_values("malmquist_productivity.csv")
  for (ori in c("in", "out")) {
    r <- ref[ref$orientation == ori, ]
    o <- dea_malmquist(f$x, f$y, f$id, f$period, rts = "crs", orientation = ori)
    expect_equal(o$table$id, as.character(r$id), info = ori)
    expect_equal(o$table$malmquist, r$malmquist, tolerance = 1e-7, info = ori)
    expect_equal(o$table$effch,     r$effch,     tolerance = 1e-7, info = ori)
    expect_equal(o$table$techch,    r$techch,    tolerance = 1e-7, info = ori)
  }
})

test_that("M is EC times TC, and the four distances are the ones claimed", {
  f <- panel_fixture(15, 3, 4)
  o <- dea_malmquist(f$x, f$y, f$id, f$period, rts = "crs")
  expect_equal(o$table$malmquist, o$table$effch * o$table$techch,
               tolerance = 1e-10)

  ## The naming convention is observation-then-technology, and that is checked
  ## against fits built by hand rather than taken on trust -- reading the two
  ## mixed terms the other way round silently swaps efficiency change for its
  ## reciprocal.
  i0 <- seq_len(15); i1 <- 16:30
  X0 <- f$x[i0, ]; Y0 <- f$y[i0, ]; X1 <- f$x[i1, ]; Y1 <- f$y[i1, ]
  d <- o$distances
  expect_equal(d$d_ff, as.numeric(dea(X0, Y0, rts = "crs", slack = FALSE)$eff))
  expect_equal(d$d_tt, as.numeric(dea(X1, Y1, rts = "crs", slack = FALSE)$eff))
  expect_equal(d$d_tf, as.numeric(dea(X1, Y1, rts = "crs", slack = FALSE,
                                      xref = X0, yref = Y0)$eff))
  expect_equal(d$d_ft, as.numeric(dea(X0, Y0, rts = "crs", slack = FALSE,
                                      xref = X1, yref = Y1)$eff))
  expect_equal(o$table$effch, d$d_tt / d$d_ff)
})

test_that("the pure/scale split equals independently computed fits", {
  ## EC = PEC * SEC holds by cancellation, so checking the product proves
  ## nothing -- it cannot fail. What can fail is either piece being the wrong
  ## quantity, so each is compared against its own fit.
  f <- panel_fixture(20, 5, 6)
  o <- dea_malmquist(f$x, f$y, f$id, f$period, rts = "crs")
  i0 <- seq_len(20); i1 <- 21:40
  c0 <- dea(f$x[i0, ], f$y[i0, ], rts = "crs", slack = FALSE)$eff
  v0 <- dea(f$x[i0, ], f$y[i0, ], rts = "vrs", slack = FALSE)$eff
  c1 <- dea(f$x[i1, ], f$y[i1, ], rts = "crs", slack = FALSE)$eff
  v1 <- dea(f$x[i1, ], f$y[i1, ], rts = "vrs", slack = FALSE)$eff
  expect_equal(o$table$purech,  as.numeric(v1 / v0))
  expect_equal(o$table$scalech, as.numeric((c1 / v1) / (c0 / v0)))
  expect_equal(o$table$effch, o$table$purech * o$table$scalech,
               tolerance = 1e-10)

  ## and it is not offered where it would be meaningless
  ## Two warnings fire here -- the split being dropped, and the vrs mixed
  ## distances being infeasible -- so both are collected rather than letting
  ## expect_warning() match one and leak the other.
  w <- character()
  g <- withCallingHandlers(
    dea_malmquist(f$x, f$y, f$id, f$period, rts = "vrs",
                  scale_decomposition = TRUE),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("defined for the crs-based index", w)))
  expect_false("scalech" %in% names(g$table))
})

test_that("the mixed distances can be infeasible under vrs, and say so", {
  ## This is the reason the default is crs: a convex hull cannot extrapolate
  ## to a point outside it, while a cone can always be scaled to reach one.
  f <- panel_fixture(20, 5, 6)
  expect_warning(g <- dea_malmquist(f$x, f$y, f$id, f$period, rts = "vrs"),
                 "infeasible")
  expect_true(anyNA(g$table$malmquist))
  expect_true(sum(g$infeasible) > 0)

  expect_false(anyNA(dea_malmquist(f$x, f$y, f$id, f$period,
                                   rts = "crs")$table$malmquist))
})

test_that("an unbalanced panel indexes the DMUs present in both periods", {
  f <- panel_fixture(12)
  keep <- !(f$id %in% c(1, 2) & f$period == 2)
  o <- dea_malmquist(f$x[keep, ], f$y[keep, ], f$id[keep], f$period[keep],
                     rts = "crs")
  expect_equal(nrow(o$table), 10L)
  expect_false(any(c("1", "2") %in% o$table$id))
})

test_that("column names and formulas select the same panel as matrices", {
  f <- panel_fixture(10)
  df <- data.frame(id = f$id, t = f$period, x1 = f$x[, 1], x2 = f$x[, 2],
                   y1 = f$y[, 1])
  a <- dea_malmquist(f$x, f$y, f$id, f$period, rts = "crs")
  b <- dea_malmquist(c("x1", "x2"), "y1", "id", "t", data = df, rts = "crs")
  c_ <- dea_malmquist(~ x1 + x2, ~ y1, ~ id, ~ t, data = df, rts = "crs")
  expect_equal(a$table$malmquist, b$table$malmquist)
  expect_equal(a$table$malmquist, c_$table$malmquist)

  ## and a prebuilt panel is the same thing again
  expect_equal(a$table$malmquist,
               dea_malmquist(dea_panel(f$x, f$y, f$id, f$period),
                             rts = "crs")$table$malmquist)
})

test_that("three periods give two consecutive pairs, chained in order", {
  d <- lapply(1:3, function(k) toy(12, 2, 1, 30 + k))
  X <- do.call(rbind, lapply(seq_along(d), function(k) d[[k]]$x * 0.95^(k - 1)))
  Y <- do.call(rbind, lapply(seq_along(d), function(k) d[[k]]$y * 1.05^(k - 1)))
  o <- dea_malmquist(X, Y, rep(1:12, 3), rep(1:3, each = 12), rts = "crs")
  expect_equal(o$n_pair, 2L)
  expect_equal(unique(paste(o$table$from, o$table$to)), c("1 2", "2 3"))
  expect_equal(nrow(o$table), 24L)
})

test_that("print and summary work on both classes", {
  f <- panel_fixture(10)
  pn <- dea_panel(f$x, f$y, f$id, f$period)
  expect_output(print(pn), "Panel")
  expect_output(summary(pn), "consecutive pairs")
  o <- dea_malmquist(pn, rts = "crs")
  expect_output(print(o), "Malmquist")
  expect_output(summary(o), "techch")
  expect_equal(nobs(o), nrow(o$table))
})
