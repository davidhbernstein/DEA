## ---------------------------------------------------------------------------
## Window analysis.
## ---------------------------------------------------------------------------

window_fixture <- function(n = 12, t = 4) {
  d <- lapply(seq_len(t), function(k) toy(n, 2, 1, seed = 300 + k))
  list(x = do.call(rbind, lapply(d, `[[`, "x")),
       y = do.call(rbind, lapply(d, `[[`, "y")),
       id = rep(seq_len(n), times = t),
       period = rep(seq_len(t), each = n),
       n = n, t = t)
}

test_that("a window's scores are an ordinary dea() fit on its stacked rows", {
  ## The definition: each DMU-period is a distinct DMU inside the window, so
  ## every window is a plain self-referenced fit on those rows. Checked
  ## window by window against dea() called directly, which is the only thing
  ## that makes the rest of this file worth anything.
  f <- window_fixture()
  pn <- dea_panel(f$x, f$y, f$id, f$period)
  w <- dea_window(pn, width = 3, rts = "vrs")

  for (k in 1:2) {
    per <- as.character(k:(k + 2))
    i <- unlist(pn$rows[per], use.names = FALSE)
    direct <- dea(f$x[i, , drop = FALSE], f$y[i, , drop = FALSE], rts = "vrs",
                  slack = FALSE, peers = FALSE)
    got <- w$table[w$table$window == k, ]
    expect_equal(got$eff, as.numeric(direct$eff), tolerance = 1e-10,
                 info = paste("window", k))
    expect_equal(got$id, as.character(f$id[i]))
    expect_equal(got$period, as.character(f$period[i]))
  }
})

test_that("the number of windows and of rows is T - w + 1 and n*w each", {
  f <- window_fixture(n = 10, t = 5)
  for (wd in 2:5) {
    w <- dea_window(f$x, f$y, f$id, f$period, width = wd)
    expect_equal(as.integer(w$n_window), 5L - wd + 1L)
    expect_equal(nrow(w$table), (5L - wd + 1L) * 10L * wd)
    expect_setequal(unique(w$table$width), wd)
  }
})

test_that("a width equal to the panel length is the pooled fit", {
  ## The one case the method degenerates in: one window, every DMU-period in
  ## it. Allowed, and it had better equal the pooled fit exactly.
  f <- window_fixture(n = 8, t = 3)
  w <- dea_window(f$x, f$y, f$id, f$period, width = 3)
  expect_equal(as.integer(w$n_window), 1L)
  pooled <- dea(f$x, f$y, rts = "vrs", slack = FALSE, peers = FALSE)
  expect_equal(w$table$eff, as.numeric(pooled$eff), tolerance = 1e-10)
})

test_that("a vector width runs every width and labels which is which", {
  ## The width is not pinned down by theory, so running two or three of them
  ## is the documented advice; `width` takes a vector so that costs an
  ## argument rather than a loop.
  f <- window_fixture(n = 9, t = 4)
  w <- dea_window(f$x, f$y, f$id, f$period, width = c(2, 4))
  expect_equal(sort(unique(w$table$width)), c(2L, 4L))
  expect_equal(as.integer(w$n_window), c(3L, 1L))

  ## Each width's rows are exactly what that width alone would have produced.
  for (wd in c(2L, 4L)) {
    one <- dea_window(f$x, f$y, f$id, f$period, width = wd)
    expect_equal(w$table[w$table$width == wd, ], one$table,
                 ignore_attr = TRUE, info = paste("width", wd))
  }
  ## Widths are sorted and de-duplicated rather than taken as given.
  expect_equal(dea_window(f$x, f$y, f$id, f$period,
                          width = c(3, 2, 3))$width, c(2L, 3L))
})

test_that("a DMU-period in the interior gets one score per window it is in", {
  ## The output is a spread, not a score: with width 3 over 4 periods, period
  ## 3 sits in both windows and period 1 in only the first.
  f <- window_fixture(n = 6, t = 4)
  w <- dea_window(f$x, f$y, f$id, f$period, width = 3)
  tb <- w$table[w$table$id == 1, ]
  expect_equal(as.integer(table(tb$period)[c("1", "3", "4")]), c(1L, 2L, 1L))

  ## And the two scores for the same DMU-period generally differ, because the
  ## two windows are different technologies. If they never did, the method
  ## would have nothing to report.
  both <- tb$eff[tb$period == "3"]
  expect_length(both, 2L)
  st <- DEA:::.window_stability(w$table)
  expect_equal(nrow(st), 6L)
  expect_true(all(st$range >= 0))
  expect_gt(max(st$range), 0)
})

test_that("ldr moves only the window, and range also moves the period", {
  ## `ldr` holds the DMU-period fixed and compares windows of one width;
  ## `range` pools over periods too, so it includes genuine change over time.
  ## Pinned by hand for one DMU, because the two were once conflated.
  f <- window_fixture(n = 10, t = 4)
  w <- dea_window(f$x, f$y, f$id, f$period, width = 3)
  st <- DEA:::.window_stability(w$table)
  tb <- w$table[w$table$id == 1, ]
  by_period <- tapply(tb$eff, tb$period, function(v) max(v) - min(v))
  expect_equal(st$ldr[st$id == "1"], max(by_period))
  expect_equal(st$range[st$id == "1"], max(tb$eff) - min(tb$eff))
  expect_true(all(st$ldr <= st$range + 1e-12))
  expect_gt(max(st$ldr), 0)

  ## One window per width means every DMU-period is scored once: no window
  ## variation exists, so ldr is exactly zero while range need not be.
  w1 <- dea_window(f$x, f$y, f$id, f$period, width = 4)
  st1 <- DEA:::.window_stability(w1$table)
  expect_equal(st1$ldr, rep(0, nrow(st1)))
  expect_gt(max(st1$range), 0)

  ## Two widths: a DMU-period's width-2 and width-3 scores are not compared
  ## with each other. With width 4 added to width 2, ldr equals width 2's.
  w24 <- dea_window(f$x, f$y, f$id, f$period, width = c(2, 4))
  w2 <- dea_window(f$x, f$y, f$id, f$period, width = 2)
  expect_equal(DEA:::.window_stability(w24$table)$ldr,
               DEA:::.window_stability(w2$table)$ldr)
})

test_that("an unbalanced panel is an ordinary case", {
  f <- window_fixture(n = 8, t = 3)
  keep <- !(f$id == 1 & f$period == 2)
  w <- dea_window(f$x[keep, ], f$y[keep, ], f$id[keep], f$period[keep],
                  width = 2)
  ## DMU 1 is missing from period 2, so it appears in fewer window rows than
  ## a DMU present throughout.
  st <- DEA:::.window_stability(w$table)
  expect_lt(st$n[st$id == "1"], st$n[st$id == "2"])
})

test_that("the width has to be a count of periods that fits the panel", {
  f <- window_fixture(n = 6, t = 3)
  expect_error(dea_window(f$x, f$y, f$id, f$period, width = 1),
               "at least 2")
  expect_error(dea_window(f$x, f$y, f$id, f$period, width = 4),
               "exceeds the panel")
  expect_error(dea_window(f$x, f$y, f$id, f$period, width = 2.5),
               "must be whole")
  expect_error(dea_window(f$x, f$y, f$id, f$period, width = NA),
               "whole numbers")
})

test_that("the period order comes from the panel, not from dea_window()", {
  ## A window is a contiguous block of the panel's order, so the order is the
  ## panel's business -- and getting it wrong groups periods that are not
  ## neighbours instead of failing.
  f <- window_fixture(n = 6, t = 3)
  lab <- c("Q9", "Q10", "Q11")[f$period]
  sorted <- dea_window(f$x, f$y, f$id, lab, width = 2)
  expect_equal(sorted$panel$periods, c("Q10", "Q11", "Q9"))
  stated <- dea_window(f$x, f$y, f$id,
                       factor(lab, levels = c("Q9", "Q10", "Q11")), width = 2)
  expect_equal(stated$panel$periods, c("Q9", "Q10", "Q11"))
  ## Different grouping, therefore different scores.
  expect_false(isTRUE(all.equal(sorted$table$eff, stated$table$eff)))
})

test_that("print, summary and plot say what the object is", {
  f <- window_fixture(n = 8, t = 4)
  w <- dea_window(f$x, f$y, f$id, f$period, width = 3)
  expect_output(print(w), "Window analysis")
  expect_output(print(w), "mean efficiency by window")
  expect_output(print(w), "two or three widths")
  expect_output(summary(w), "stability across windows")
  expect_equal(nobs(w), nrow(w$table))

  w2 <- dea_window(f$x, f$y, f$id, f$period, width = c(2, 3))
  pdf(NULL); on.exit(dev.off())
  out <- plot(w2)
  expect_s3_class(out, "data.frame")
  expect_setequal(unique(out$width), c(2L, 3L))
})

test_that("a dea_panel and raw vectors are the same call", {
  f <- window_fixture(n = 7, t = 3)
  pn <- dea_panel(f$x, f$y, f$id, f$period)
  a <- dea_window(pn, width = 2)
  b <- dea_window(f$x, f$y, f$id, f$period, width = 2)
  expect_equal(a$table, b$table)

  ## Column names through `data`, as everywhere else in the package.
  df <- data.frame(firm = f$id, yr = f$period, k = f$x[, 1], l = f$x[, 2],
                   out = f$y[, 1])
  c3 <- dea_window(c("k", "l"), "out", "firm", "yr", data = df, width = 2)
  expect_equal(c3$table$eff, b$table$eff, tolerance = 1e-10)
})
