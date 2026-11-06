## ---------------------------------------------------------------------------
## The efficiency ladder, and the isolated-DMU flag that is its other half.
## ---------------------------------------------------------------------------

test_that("the ladder only ever moves away from the frontier, on every rung", {
  ## This is the property that would break if a refit were not nested in the
  ## one before it, which is the mistake dea_ladder() is most able to make.
  ## Asserted rung by rung rather than at the endpoints, because a ladder that
  ## wandered and came back would pass an endpoint check.
  d <- toy(30, 2, 1, seed = 7)
  for (rts in c("vrs", "crs", "nirs")) {
    fin  <- dea(d$x, d$y, rts = rts, orientation = "in")
    fout <- dea(d$x, d$y, rts = rts, orientation = "out")

    lin <- dea_ladder(fin, dmu = which.min(fin$eff))
    expect_true(all(diff(lin$table$eff) >= -1e-9),
                info = paste("input,", rts))

    lout <- dea_ladder(fout, dmu = which.max(fout$eff))
    expect_true(all(diff(lout$table$eff) <= 1e-9),
                info = paste("output,", rts))
  }
})

test_that("rung 0 is the fit as reported, and the ladder ends at 1", {
  d <- toy(25, 2, 1, seed = 3)
  f <- dea(d$x, d$y, rts = "vrs")
  o <- which.min(f$eff)
  lad <- dea_ladder(f, dmu = o)

  expect_equal(lad$table$eff[1], as.numeric(f$eff[o]))
  expect_equal(lad$table$step[1], 0L)
  expect_true(is.na(lad$table$dropped[1]))
  expect_equal(lad$table$n_ref[1], f$nref)

  ## The DMU keeps its own column, so once every other peer is gone it
  ## references itself alone and the score is exactly 1. That is what makes
  ## the ladder terminate.
  expect_equal(lad$table$eff[nrow(lad$table)], 1, tolerance = 1e-9)
  expect_equal(diff(lad$table$n_ref), rep(-1L, nrow(lad$table) - 1L))
})

test_that("a rung's score matches an independent refit from the definition", {
  ## The ladder re-scores one DMU against a reduced reference set using
  ## xref/yref. Here the same number is recomputed the other way -- fit the
  ## reduced set against ITSELF and read the row -- so the check does not
  ## depend on the external-reference path being right.
  d <- toy(20, 2, 1, seed = 11)
  f <- dea(d$x, d$y, rts = "vrs")
  o <- which.min(f$eff)
  lad <- dea_ladder(f, dmu = o, max_steps = 3)

  keep <- rep(TRUE, f$nref)
  for (k in seq_len(nrow(lad$table) - 1L)) {
    keep[match(lad$table$dropped[k + 1L], f$ref)] <- FALSE
    g <- dea(d$x[keep, , drop = FALSE], d$y[keep, , drop = FALSE], rts = "vrs")
    i <- match(o, which(keep))
    expect_equal(lad$table$eff[k + 1L], as.numeric(g$eff[i]),
                 tolerance = 1e-8,
                 info = paste("rung", k))
  }
})

test_that("which peer is dropped is defined and reproducible", {
  d <- toy(25, 2, 1, seed = 5)
  f <- dea(d$x, d$y, rts = "vrs")
  o <- which.min(f$eff)

  ## Rung 1 drops the largest lambda in the ORIGINAL fit, and records that
  ## lambda -- not the one the peer had in some later reference set.
  p <- peers(f)
  mine <- p[p$dmu == f$dmu[o] & p$peer != f$dmu[o], ]
  top <- mine$peer[which.max(mine$lambda)]
  lad <- dea_ladder(f, dmu = o)
  expect_equal(lad$table$dropped[2], top)
  expect_equal(lad$table$lambda[2], max(mine$lambda))

  ## Same call, same table: nothing here is left to the solver's choice of
  ## vertex between runs.
  expect_equal(dea_ladder(f, dmu = o)$table, lad$table)
  ## By name and by index are the same DMU.
  expect_equal(dea_ladder(f, dmu = f$dmu[o])$table, lad$table)
})

test_that("max_steps truncates and the object reports the stopping reason", {
  d <- toy(20, 2, 1, seed = 9)
  f <- dea(d$x, d$y, rts = "vrs")
  o <- which.min(f$eff)
  short <- dea_ladder(f, dmu = o, max_steps = 2)
  expect_equal(nrow(short$table), 3L)
  expect_equal(short$n_step, 2L)
  expect_false(short$exhausted)

  full <- dea_ladder(f, dmu = o)
  expect_true(full$exhausted)
  expect_equal(short$table, head(full$table, 3))
})

test_that("the ladder refuses the fits it cannot audit", {
  d <- toy(15, 2, 1, seed = 2)
  f <- dea(d$x, d$y, rts = "vrs")

  expect_error(dea_ladder(d$x, dmu = 1), "must be a fit")
  expect_error(dea_ladder(dea_sbm(d$x, d$y), dmu = 1), "radial model only")
  expect_error(dea_ladder(dea(d$x, d$y, super = TRUE), dmu = 1),
               "super-efficiency")
  expect_error(dea_ladder(dea(d$x, d$y, peers = FALSE), dmu = 1),
               "peers = FALSE")
  ## An external reference set makes "drop a peer" a different question.
  ext <- dea(d$x[1:5, ], d$y[1:5, ], xref = d$x, yref = d$y)
  expect_error(dea_ladder(ext, dmu = 1), "external reference set")

  expect_error(dea_ladder(f, dmu = "nope"), "not one of this fit")
  expect_error(dea_ladder(f, dmu = 99), "outside 1:15")
  expect_error(dea_ladder(f, dmu = c(1, 2)), "a single DMU")
  expect_error(dea_ladder(f, dmu = 1, max_steps = 0), "positive")
})

test_that("print, summary and plot say what the object is", {
  d <- toy(18, 2, 1, seed = 4)
  f <- dea(d$x, d$y, rts = "vrs")
  lad <- dea_ladder(f, dmu = which.min(f$eff))
  expect_output(print(lad), "Efficiency ladder")
  expect_output(print(lad), "ends at 1")
  expect_output(summary(lad), "largest single move")
  expect_equal(nobs(lad), nrow(lad$table))
  pdf(NULL); on.exit(dev.off())
  expect_s3_class(plot(lad), "data.frame")
})

## ---------------------------------------------------------------------------
## The isolated-DMU flag in summary.dea().
## ---------------------------------------------------------------------------

test_that("the isolated-DMU flag is not a second name for `efficient`", {
  ## This is the bug the flag was born with. "The DMU's own lambda is e_o" is
  ## true of EVERY efficient DMU, so the first version of this flagged exactly
  ## charnes1981's 27 efficient units and reported it as news. The flag is
  ## about the OTHER direction -- efficient and in nobody else's reference set
  ## -- so it must be a STRICT subset of the efficient ones on a design where
  ## some efficient DMUs do envelop others.
  data("charnes1981", package = "DEA")
  X <- as.matrix(charnes1981[, paste0("x", 1:5)])
  Y <- as.matrix(charnes1981[, paste0("y", 1:3)])
  f <- dea(X, Y, rts = "vrs")
  iso <- DEA:::.dea_self_referencing(f)

  expect_true(all(f$efficient[iso]))
  expect_gt(sum(iso), 0L)
  expect_lt(sum(iso), sum(f$efficient))

  ## And the same set computed the other way, off the peer table.
  p <- peers(f)
  by_peers <- setdiff(f$dmu[f$efficient], unique(p$peer[p$dmu != p$peer]))
  expect_setequal(f$dmu[iso], by_peers)
})

test_that("the flag is empty for the fits where it has no meaning", {
  d <- toy(20, 2, 1, seed = 6)
  ## Super-efficiency excludes each DMU from its own reference set.
  expect_false(any(DEA:::.dea_self_referencing(
    dea(d$x, d$y, super = TRUE))))
  ## With an external reference set the DMU is not in the set at all.
  expect_false(any(DEA:::.dea_self_referencing(
    dea(d$x[1:6, ], d$y[1:6, ], xref = d$x, yref = d$y))))
  ## And with no lambda there is nothing to count.
  expect_false(any(DEA:::.dea_self_referencing(
    dea(d$x, d$y, peers = FALSE))))
})

test_that("summary() reports the isolated DMUs when there are any", {
  data("charnes1981", package = "DEA")
  X <- as.matrix(charnes1981[, paste0("x", 1:5)])
  Y <- as.matrix(charnes1981[, paste0("y", 1:3)])
  f <- dea(X, Y, rts = "vrs")
  expect_output(summary(f), "efficient but enveloping no one")
  expect_output(summary(f), "dea_ladder")
})
