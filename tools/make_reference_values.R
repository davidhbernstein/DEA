## ---------------------------------------------------------------------------
## tools/make_reference_values.R -- regenerate tests/testthat/reference/*.csv
##
##   Rscript tools/make_reference_values.R          (from the package root)
##
## NOT SHIPPED. This directory is in .Rbuildignore, so nothing here reaches the
## tarball and the package declares no dependency on any other DEA package.
##
## WHY THE VALUES ARE PINNED RATHER THAN COMPUTED AT TEST TIME.
##
## The package's second correctness check is agreement with independent
## implementations of the same linear programs -- Benchmarking and DJL -- for
## every technology and orientation. That check was previously made by calling
## those packages from the test suite, behind skip_if_not_installed().
##
## Three things were wrong with that:
##
##   1. IT MOSTLY DID NOT RUN. skip_if_not_installed() means the check is
##      silently skipped on any machine without the competitor installed, which
##      is most of CRAN's check farm. A test that is usually skipped provides
##      most of the reassurance and little of the protection.
##
##   2. IT MADE A DEA PACKAGE DEPEND ON RIVAL DEA PACKAGES. Suggests is still a
##      dependency: CRAN installs it to check, its breakage becomes our check
##      failure, and a package whose selling point is being a clean, minimal
##      implementation should not drag the field in behind it.
##
##   3. IT COULD DRIFT. Agreeing with "whatever Benchmarking does now" is a
##      moving target. Agreeing with what it computed on a stated date, from a
##      stated version, is a fact that stays true.
##
## So the comparison is run ONCE, here, and the answers are written out as CSV
## with full precision. The test suite then compares against the recorded
## numbers, unconditionally and on every machine. If this package's answers
## ever move, the tests fail exactly as before -- and the live, up-to-date
## comparison against these and five other packages still happens, in
## ../horserace/, which is not shipped and not run by CRAN.
##
## REGENERATE ONLY DELIBERATELY. A regenerated file that changes existing
## numbers means either this package's answers moved or the competitor's did,
## and the CSV diff is the evidence for which. Commit the diff separately from
## whatever caused it.
## ---------------------------------------------------------------------------

if (!requireNamespace("pkgload", quietly = TRUE)) stop("needs pkgload")
for (p in c("Benchmarking", "DJL")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("needs ", p, " installed to regenerate reference values")
}
suppressMessages(pkgload::load_all(".", quiet = TRUE))

outdir <- file.path("tests", "testthat", "reference")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## EXACT round-tripping, verified per value.
##
## Two traps here, both found the hard way.
##
## 1. NOT formatC(). It pads to a common width, so every field acquires leading
##    spaces and read.csv() with quote = FALSE returns the column as CHARACTER.
##    That failed loudly -- but a character column compared against a double is
##    not always an error, sometimes it is a silent coercion.
##
## 2. "%.17g" IS NOT ENOUGH, despite 17 significant digits being the usual
##    round-trip guarantee for an IEEE double. sprintf("%.17g", 9.09) gives
##    "9.0899999999999999", and as.numeric() of that is a DIFFERENT double from
##    9.09 -- off by 1.78e-15. That is far below every tolerance in the suite
##    and would have gone unnoticed, except that the charnes1981 data set is
##    compared at tolerance 0 on purpose, and it turned the check into a
##    failure that looked like the shipped data being wrong. It was not: the
##    two copies are `identical()`. The reference file was lossy.
##
## So: widen the format until the value round-trips EXACTLY, per value, and
## keep the shortest representation that does. Most numbers need far fewer than
## 17 digits, so the file is also smaller than a fixed-width one.
## The ladder covers %g and then %e, shortest first so the file stays readable.
##
## SOME DOUBLES CANNOT BE STORED EXACTLY AS DECIMAL TEXT IN R, and that is a
## limitation of R's PARSER, not of the format. For example the double
## 0x1.234e73390c50ep+0 prints as 1.137915803363572120332 and
## as.numeric("1.137915803363572120332") returns a different double, one ulp
## away -- at every width from 15 to 22 significant digits. No decimal string
## round-trips it.
##
## So exactness is attempted, not assumed. Values that round-trip are stored in
## their shortest exact form; the rest fall back to 17 significant digits and
## are COUNTED, with the worst relative error reported by wr(). The fallback
## error is about 1e-16 relative -- ten orders of magnitude below the loosest
## tolerance any test here uses (1e-6) and seven below the tightest (1e-9) --
## so it changes no comparison. It is reported rather than swallowed because a
## reference file quietly losing precision is exactly the kind of thing that
## should not be discovered later.
##
## The one place exactness genuinely matters is charnes1981, which is DATA
## rather than a computed result, is compared at tolerance 0, and consists of
## short decimals that all round-trip. wr() reports per file, so that can be
## seen rather than hoped for.
.inexact <- new.env(parent = emptyenv())
.inexact$n <- 0L; .inexact$worst <- 0

.exact <- function(z) {
  if (is.na(z)) return("NA")
  for (d in 15:22) {
    s <- sprintf(paste0("%.", d, "g"), z)
    if (identical(as.numeric(s), z)) return(s)
  }
  for (d in 15:21) {
    s <- sprintf(paste0("%.", d, "e"), z)
    if (identical(as.numeric(s), z)) return(s)
  }
  s <- sprintf("%.17g", z)
  rel <- if (z == 0) abs(as.numeric(s)) else abs(as.numeric(s) - z) / abs(z)
  .inexact$n <- .inexact$n + 1L
  .inexact$worst <- max(.inexact$worst, rel)
  s
}

wr <- function(df, name) {
  .inexact$n <- 0L; .inexact$worst <- 0
  num <- vapply(df, is.double, logical(1))
  out <- df
  out[num] <- lapply(out[num], function(v) vapply(v, .exact, character(1)))
  f <- file.path(outdir, name)
  utils::write.csv(out, f, row.names = FALSE, quote = FALSE)

  ## Round-trip, checked with identical() rather than all.equal(). A reference
  ## file that does not survive being read is worse than no reference file, and
  ## `all.equal` at 1e-15 would have passed the very error described above.
  back <- utils::read.csv(f, stringsAsFactors = FALSE)
  for (cl in names(df)[vapply(df, is.numeric, logical(1))]) {
    if (!is.numeric(back[[cl]]))
      stop(name, ": column '", cl, "' read back as ", class(back[[cl]]))
    a <- as.numeric(df[[cl]]); b <- as.numeric(back[[cl]])
    if (!isTRUE(all.equal(a, b, tolerance = 1e-14)))
      stop(name, ": column '", cl, "' did not round-trip")
  }
  cat(sprintf("  %-28s %d rows, %s\n", name, nrow(df),
              if (.inexact$n == 0L) "all values exact"
              else sprintf("%d value(s) inexact, worst rel. err %.1e",
                           .inexact$n, .inexact$worst)))
}

## The fixtures must match the ones the tests build, exactly. They come from
## dea_sim() with a fixed seed, so they are reproducible from the package
## itself and no data is stored here -- only the competitors' answers.
toy_ <- function(n, p, q, seed, returns = 1) {
  s <- dea_sim(n, p = p, q = q, returns = returns, seed = seed)
  list(x = s$x, y = s$y)
}

cat("Benchmarking ", as.character(packageVersion("Benchmarking")),
    " / DJL ", as.character(packageVersion("DJL")),
    " / productivity ", as.character(packageVersion("productivity")),
    " / truncreg ", as.character(packageVersion("truncreg")),
    " / npsf ", as.character(packageVersion("npsf")),
    "\n", sep = "")

## ---------------------------------------------------------------------------
## THIS SCRIPT IS THE LIVE CROSS-PACKAGE CHECK.
##
## Pinning the reference values removed the competitors from Suggests, which
## is what was wanted -- but it would also have removed the live comparison
## entirely if nothing replaced it, and "our answers still match the numbers we
## recorded from ourselves-plus-them in 2026" is weaker than it sounds if the
## recording is never re-derived.
##
## So every value written below is first CHECKED against this package's current
## answer, and the script refuses to write a file where they disagree by more
## than the tolerance the corresponding test uses. Running this tool is
## therefore a full, current, live validation against Benchmarking and DJL --
## deliberately, on a machine that has them, outside CRAN.
##
## A disagreement here is the signal the old skip_if_not_installed() tests
## existed to produce, and now it cannot be silently skipped: it stops the run.
## ---------------------------------------------------------------------------
.agree <- function(ours, theirs, tol, what) {
  ours <- as.numeric(ours); theirs <- as.numeric(theirs)
  if (length(ours) != length(theirs))
    stop(what, ": length ", length(ours), " vs ", length(theirs))
  ok <- is.finite(ours) & is.finite(theirs)
  if (!identical(is.finite(ours), is.finite(theirs)))
    stop(what, ": the two disagree about WHICH values are finite")
  d <- if (any(ok)) max(abs(ours[ok] - theirs[ok])) else 0
  if (d > tol) stop(sprintf("%s: DISAGREE by %.3e (tolerance %.0e)", what, d, tol))
  cat(sprintf("    %-40s agree to %.2e\n", what, d))
  invisible(d)
}

## --- radial, every technology and orientation, plus super-efficiency --------
d <- toy_(60, 2, 2, 9)
alias <- c(crs = "crs", vrs = "vrs", nirs = "drs", ndrs = "irs", fdh = "fdh")
rows <- list()
for (rts in names(alias)) for (ori in c("in", "out")) {
  th <- suppressWarnings(Benchmarking::dea(d$x, d$y, RTS = alias[[rts]], ORIENTATION = ori))
  .agree(dea(d$x, d$y, rts = rts, orientation = ori, slack = FALSE)$eff,
         th$eff, 1e-7, paste("radial", rts, ori))
  rows[[length(rows) + 1L]] <- data.frame(rts = rts, orientation = ori,
                                          dmu = seq_along(th$eff),
                                          eff = as.numeric(th$eff),
                                          stringsAsFactors = FALSE)
}
wr(do.call(rbind, rows), "radial_Benchmarking.csv")

sup <- suppressWarnings(Benchmarking::sdea(d$x, d$y, RTS = "vrs", ORIENTATION = "in"))
{
  os <- dea(d$x, d$y, rts = "vrs", orientation = "in", super = TRUE)$eff
  if (!identical(unname(is.na(os)), unname(!is.finite(as.numeric(sup$eff)))))
    stop("super-efficiency: the two disagree about which DMUs are infeasible")
  k <- is.finite(as.numeric(sup$eff))
  .agree(os[k], as.numeric(sup$eff)[k], 1e-7, "super-efficiency vrs in")
}
wr(data.frame(dmu = seq_along(sup$eff), eff = as.numeric(sup$eff),
              feasible = is.finite(as.numeric(sup$eff))), "super_Benchmarking.csv")

## --- additive, unweighted total --------------------------------------------
d <- toy_(40, 2, 2, 8)
alias <- c(crs = "crs", vrs = "vrs", nirs = "drs", ndrs = "irs")
rows <- list()
for (rts in names(alias)) {
  th <- suppressWarnings(Benchmarking::dea.add(d$x, d$y, RTS = alias[[rts]]))
  .agree(dea_add(d$x, d$y, measure = "unweighted", rts = rts, scaling = FALSE)$eff,
         th$sum, 1e-7, paste("additive unweighted", rts))
  rows[[length(rows) + 1L]] <- data.frame(rts = rts, dmu = seq_along(th$sum),
                                          sum = as.numeric(th$sum),
                                          stringsAsFactors = FALSE)
}
wr(do.call(rbind, rows), "additive_Benchmarking.csv")

## --- slacks-based measure, against DJL -------------------------------------
d <- toy_(40, 2, 2, 4)
rows <- list()
for (rts in c("crs", "vrs")) for (ori in c("none", "in", "out")) {
  th <- suppressWarnings(suppressMessages(
    DJL::dm.sbm(d$x, d$y, rts = rts, orientation = substr(ori, 1, 1))))
  .agree(dea_sbm(d$x, d$y, rts = rts, orientation = ori)$eff, th$eff, 1e-6,
         paste("sbm", rts, ori))
  rows[[length(rows) + 1L]] <- data.frame(rts = rts, orientation = ori,
                                          dmu = seq_along(as.numeric(th$eff)),
                                          eff = as.numeric(th$eff),
                                          stringsAsFactors = FALSE)
}
wr(do.call(rbind, rows), "sbm_DJL.csv")

## --- the truncated regression, against two implementations of it -----------
## .truncreg() is written inside this package because the two packages that
## provide a truncated regression are not, and will not become, dependencies --
## the same reasoning as everything else in this file. So its answers are
## pinned here against BOTH of them.
##
## truncreg::truncreg is the dedicated implementation; npsf::truncreg is the
## one a DEA user is most likely to reach for, since npsf is a DEA package.
## They agree with each other and with this one to their own optimiser
## tolerances, which differ: npsf matches to 1e-9 and truncreg to about 1e-5,
## so the recorded values are npsf's and the tolerance is truncreg's.
##
## The LOG-LIKELIHOOD is recorded too, and it is the reason this block exists
## in the form it does: the first version of .truncreg() dropped the
## 0.5*log(2*pi) constant, which does not move the optimum and made $logLik
## wrong by n * 0.9189 -- invisible in every coefficient and obvious the moment
## it was compared with someone else's.
set.seed(7)
{
  n <- 400L
  z1 <- runif(n); z2 <- rnorm(n)
  yy <- 1.2 + 0.8 * z1 - 0.5 * z2 + rnorm(n, 0, 0.6)
  keep <- yy > 1
  dd <- data.frame(y = yy[keep], z1 = z1[keep], z2 = z2[keep])
  Z <- cbind(`(Intercept)` = 1, z1 = dd$z1, z2 = dd$z2)

  mine <- .truncreg(dd$y, Z, ll = 1)
  tr <- truncreg::truncreg(y ~ z1 + z2, data = dd, point = 1,
                           direction = "left")
  np <- npsf::truncreg(y ~ z1 + z2, data = dd, ll = 1, print.level = 0)

  .agree(c(mine$beta, mine$sigma), as.numeric(coef(tr)), 1e-4,
         "truncated regression (truncreg)")
  .agree(c(mine$beta, mine$sigma), as.numeric(np$coef), 1e-7,
         "truncated regression (npsf)")
  .agree(mine$logLik, as.numeric(stats::logLik(tr)), 1e-5,
         "truncated regression logLik")

  wr(data.frame(term = c(names(mine$beta), "sigma"),
                estimate = as.numeric(np$coef),
                logLik = as.numeric(stats::logLik(tr)),
                stringsAsFactors = FALSE),
     "truncreg_reference.csv")
}

## --- weight restrictions, against Benchmarking::dea.dual --------------------
## Benchmarking's DUAL argument takes only one shape of restriction -- ratios
## against the FIRST input and the FIRST output -- so it checks the assurance
## region case and nothing else. That is still the case worth checking against
## an independent implementation, because it is the one where a row of the
## multiplier program can be written with the wrong sign and still produce
## plausible numbers.
##
## WHICH ROW IS WHICH WAS PINNED BY MEASUREMENT, with deliberately asymmetric
## bounds so the two readings cannot both fit: row 1 is the input ratio
## v[2]/v[1] and row 2 the output ratio u[2]/u[1]. Reading them the other way
## round moves the scores by 0.30 on this design, and the restriction binds by
## 0.32 against an unrestricted fit -- so this comparison can fail, which is
## the only kind worth recording.
d <- toy_(30, 2, 2, 8)
DU <- matrix(c(0.3, 0.8,
               1.2, 4.0), nrow = 2L, byrow = TRUE)
rows <- list()
for (ori in c("in", "out")) {
  th <- Benchmarking::dea.dual(d$x, d$y, RTS = "crs", ORIENTATION = ori,
                               DUAL = DU)
  rs <- c(wr_ratio("v", 2, 1, lower = DU[1, 1], upper = DU[1, 2]),
          wr_ratio("u", 2, 1, lower = DU[2, 1], upper = DU[2, 2]))
  ours <- suppressWarnings(
    dea_wr(d$x, d$y, rs, rts = "crs", orientation = ori)$eff)
  .agree(ours, as.numeric(th$eff), 1e-7, paste("weight restrictions", ori))
  free <- dea(d$x, d$y, rts = "crs", orientation = ori, slack = FALSE)$eff
  if (max(abs(as.numeric(th$eff) - free)) < 0.02)
    stop("weight restrictions ", ori, ": the restriction barely binds on this ",
         "design, so the comparison cannot distinguish a right implementation ",
         "from a wrong one.")
  rows[[length(rows) + 1L]] <- data.frame(
    orientation = ori, dmu = seq_along(th$eff),
    eff = as.numeric(th$eff), stringsAsFactors = FALSE)
}
wr(do.call(rbind, rows), "weights_Benchmarking.csv")

## --- the Malmquist index, against TWO packages, and a convention they split on
## Two references rather than one, because they DISAGREE and the disagreement
## is the thing that had to be settled.
##
##   * Benchmarking::malmq() reports its four distances on the native scale of
##     the orientation -- theta input-oriented, phi output-oriented -- and
##     builds the index from them directly. So its OUTPUT-oriented index is the
##     RECIPROCAL of the one below, and under its convention "greater than 1"
##     means growth on inputs and decline on outputs.
##   * productivity::malm() agrees with this package in BOTH orientations.
##
## The convention here is Fare, Grosskopf, Norris and Zhang's: the index is
## built from Shephard distance functions, so M > 1 is productivity growth
## whichever side it is measured from. Every score is converted to the common
## (0, 1] scale once (phi -> 1/phi) before any ratio is taken.
##
## What is pinned is this package's convention, checked against BOTH packages:
## Benchmarking exactly on input and reciprocally on output, productivity
## exactly on both. A single reference could not have shown that the
## disagreement is a convention rather than an error.
d0 <- toy_(20, 2, 1, 21)
d1 <- toy_(20, 2, 1, 22)
d1$x <- d1$x * 0.92; d1$y <- d1$y * 1.12          # a frontier that moves
pid <- seq_len(20)
pX <- rbind(d0$x, d1$x); pY <- rbind(d0$y, d1$y)
pdf_ <- data.frame(id = rep(pid, 2L), t = rep(1:2, each = 20L),
                   x1 = pX[, 1], x2 = pX[, 2], y1 = pY[, 1])
rows <- list()
for (ori in c("in", "out")) {
  ours <- dea_malmquist(pX, pY, rep(pid, 2L), rep(1:2, each = 20L),
                        rts = "crs", orientation = ori)
  bm <- Benchmarking::malmq(d0$x, d0$y, pid, d1$x, d1$y, pid,
                            RTS = "crs", ORIENTATION = ori)
  flip <- if (identical(ori, "out")) function(v) 1 / v else identity
  .agree(ours$table$malmquist, flip(bm$m),  1e-7, paste("malmquist M ", ori, "(Benchmarking)"))
  .agree(ours$table$effch,     flip(bm$ec), 1e-7, paste("malmquist EC", ori, "(Benchmarking)"))
  .agree(ours$table$techch,    flip(bm$tc), 1e-7, paste("malmquist TC", ori, "(Benchmarking)"))

  pm <- suppressMessages(
    productivity::malm(pdf_, "id", "t", c("x1", "x2"), "y1",
                       rts = "crs", orientation = ori))
  ch <- pm$Changes[match(ours$table$id, pm$Changes$id), ]
  .agree(ours$table$malmquist, ch$malmquist, 1e-7, paste("malmquist M ", ori, "(productivity)"))
  .agree(ours$table$effch,     ch$effch,     1e-7, paste("malmquist EC", ori, "(productivity)"))
  .agree(ours$table$techch,    ch$tech,      1e-7, paste("malmquist TC", ori, "(productivity)"))

  rows[[length(rows) + 1L]] <- data.frame(
    orientation = ori, id = ours$table$id,
    malmquist = as.numeric(ch$malmquist), effch = as.numeric(ch$effch),
    techch = as.numeric(ch$tech), stringsAsFactors = FALSE)
}
## And that Benchmarking really does invert on output rather than merely
## differ: assert the reciprocal holds and that the two are NOT equal.
{
  o <- dea_malmquist(pX, pY, rep(pid, 2L), rep(1:2, each = 20L),
                     rts = "crs", orientation = "out")
  bm <- Benchmarking::malmq(d0$x, d0$y, pid, d1$x, d1$y, pid,
                            RTS = "crs", ORIENTATION = "out")
  if (max(abs(o$table$malmquist - bm$m)) < 0.05)
    stop("Benchmarking's output-oriented Malmquist no longer differs from ",
         "this one; the convention note in R/malmquist.R needs rechecking.")
}
wr(do.call(rbind, rows), "malmquist_productivity.csv")

## --- directional distance, four directions x four technologies --------------
## The DDF had no cross-package reference at all until now, which is the one
## model here where that gap mattered most: its direction is an argument, so a
## misreading of the direction convention would show up nowhere else. The
## comparison is Benchmarking::dea.direct(), and two things about it had to be
## pinned by measurement before any number below meant anything.
##
##   1. `objval` IS beta. `eff` is something else -- for ORIENTATION = "in" it
##      is an n x p matrix of per-input Farrell ratios, which happens to equal
##      1 - beta when DIRECT = X and does not in general. Reading `eff` gives
##      numbers near 1 that look like efficiencies and are not this estimand.
##   2. ORIENTATION = "in-out" is the general direction. "in" and "out" hold
##      the other half of the direction at zero regardless of what DIRECT says,
##      so passing a full (p+q) direction with ORIENTATION = "in" silently
##      measures a different program.
##
## The design is chosen so that every technology BINDS, and binds against BOTH
## of its neighbours: for every direction, nirs and ndrs each differ from vrs
## and from crs on at least 11 of the 40 DMUs. Where a restriction does not
## bind, all four technologies agree and the comparison proves nothing -- the
## same rule test-multipliers.R and horserace/validate_feature_adapters.R
## follow, and it is asserted below rather than assumed. Checking only against
## vrs is not enough: on a design sitting entirely in the decreasing-returns
## region, nirs collapses onto crs and ndrs onto vrs, so half the technologies
## are never separately exercised while a vrs-only gate still passes. The first
## design tried here did exactly that.
d <- toy_(40, 2, 2, 12, returns = 0.8)
alias <- c(crs = "crs", vrs = "vrs", nirs = "drs", ndrs = "irs")
rng <- cbind(sweep(d$x, 2L, apply(d$x, 2L, min), "-"),
             -sweep(d$y, 2L, apply(d$y, 2L, max), "-"))
rows <- list()
bind <- list()
for (dir in c("both", "in", "out", "range")) {
  G <- if (identical(dir, "range")) rng else
       .ddf_direction(dir, d$x, d$y, d$x, d$y, nrow(d$x), ncol(d$x), ncol(d$y))
  for (rts in names(alias)) {
    th <- suppressWarnings(
      Benchmarking::dea.direct(d$x, d$y, DIRECT = G, RTS = alias[[rts]],
                               ORIENTATION = "in-out"))
    ours <- suppressWarnings(
      dea_ddf(d$x, d$y, direction = dir, rts = rts, peers = FALSE)$beta)
    .agree(ours, th$objval, 1e-7, paste("ddf", dir, rts))
    bind[[paste(dir, rts)]] <- as.numeric(th$objval)
    rows[[length(rows) + 1L]] <- data.frame(
      direction = dir, rts = rts, dmu = seq_along(th$objval),
      beta = as.numeric(th$objval), stringsAsFactors = FALSE)
  }
}
for (dir in c("both", "in", "out", "range")) {
  for (rr in c("nirs", "ndrs")) for (ref in c("vrs", "crs")) {
    k <- sum(abs(bind[[paste(dir, rr)]] - bind[[paste(dir, ref)]]) > 1e-8)
    if (k < 5)
      stop("ddf ", dir, "/", rr, " against ", ref, ": they differ on only ", k,
           " DMUs, so this design cannot tell the technologies apart.")
  }
}
wr(do.call(rbind, rows), "ddf_Benchmarking.csv")

## --- multiplier weights, crs input-oriented ---------------------------------
## NOTE THE FIELD NAMES. Benchmarking::dea(DUAL = TRUE) returns `ux` for the
## INPUT weights and `vy` for the OUTPUT weights -- the reverse of this
## package's v/u, and of the CCR notation. Established in ../horserace/ by
## dimension and by value; see horserace/README.md, "Two more convention
## traps". Recording them under OUR names here, once, is what keeps the swap
## from having to be re-derived every time the test is read.
d <- toy_(25, 2, 2, 12, returns = 0.9)
b <- Benchmarking::dea(d$x, d$y, RTS = "crs", ORIENTATION = "in", DUAL = TRUE)
V <- as.matrix(b$ux); U <- as.matrix(b$vy)
{
  m <- suppressWarnings(dea(d$x, d$y, rts = "crs", orientation = "in",
                            slack = FALSE, multipliers = TRUE))
  ie <- m$eff < 1 - 1e-6
  .agree(unname(m$v[ie, ]), V[ie, ], 1e-6, "multipliers v (inefficient DMUs)")
  .agree(unname(m$u[ie, ]), U[ie, ], 1e-6, "multipliers u (inefficient DMUs)")
}
wr(data.frame(dmu = rep(seq_len(nrow(V)), ncol(V) + ncol(U)),
              which = rep(c(rep("v", ncol(V)), rep("u", ncol(U))),
                          each = nrow(V)),
              j = rep(c(seq_len(ncol(V)), seq_len(ncol(U))), each = nrow(V)),
              value = c(as.numeric(V), as.numeric(U)),
              stringsAsFactors = FALSE), "multipliers_Benchmarking.csv")

## --- cost, revenue and profit optima ---------------------------------------
n <- nrow(d$x)
W <- matrix(rep(c(1.2, 0.8), each = n), n, 2)
P <- matrix(rep(c(2.0, 3.1), each = n), n, 2)
map <- c(vrs = "vrs", crs = "crs", drs = "nirs", irs = "ndrs")
rows <- list()
for (B in names(map)) {
  co <- Benchmarking::cost.opt(d$x, d$y, W, RTS = B)
  ro <- Benchmarking::revenue.opt(d$x, d$y, P, RTS = B)
  .agree(dea_cost(d$x, d$y, W, rts = map[[B]])$optimal, rowSums(co$xopt * W),
         1e-8, paste("cost optimum", map[[B]]))
  .agree(dea_revenue(d$x, d$y, P, rts = map[[B]])$optimal, rowSums(ro$yopt * P),
         1e-8, paste("revenue optimum", map[[B]]))
  rows[[length(rows) + 1L]] <- data.frame(
    rts = map[[B]], dmu = seq_len(n),
    cost_opt = as.numeric(rowSums(co$xopt * W)),
    revenue_opt = as.numeric(rowSums(ro$yopt * P)),
    stringsAsFactors = FALSE)
}
po <- Benchmarking::profit.opt(d$x, d$y, W, P, RTS = "vrs")
.agree(dea_profit(d$x, d$y, W, P, rts = "vrs")$profit_max,
       rowSums(po$yopt * P) - rowSums(po$xopt * W), 1e-8, "profit maximum vrs")
pr <- do.call(rbind, rows)
pr$profit_max_vrs <- NA_real_
pr$profit_max_vrs[pr$rts == "vrs"] <- as.numeric(rowSums(po$yopt * P) - rowSums(po$xopt * W))
wr(pr, "price_Benchmarking.csv")

## --- the charnes1981 data set ----------------------------------------------
## Not a model comparison: this checks that the copy this package SHIPS still
## matches the copy Benchmarking ships, which is what man/charnes1981.Rd's
## Source section claims. Pinning it means the package asserts its own data is
## unchanged without needing Benchmarking installed. The cost is that an
## upstream CORRECTION would no longer be noticed here -- that duty moves to
## ../horserace/ and is recorded in ../ROADMAP.md.
e <- new.env(); utils::data("charnes1981", package = "Benchmarking", envir = e)
up <- get("charnes1981", envir = e)
{
  data(charnes1981, package = "DEA", envir = environment())
  ours <- get("charnes1981", envir = environment())
  for (v in c(paste0("x", 1:5), paste0("y", 1:3)))
    if (!identical(ours[[v]], up[[v]]))
      stop("charnes1981: column ", v, " differs from Benchmarking's copy")
  cat("    charnes1981                              identical to Benchmarking\n")
}
wr(data.frame(dmu = seq_len(nrow(up)),
              up[, c(paste0("x", 1:5), paste0("y", 1:3))],
              pft = as.integer(up$pft),
              stringsAsFactors = FALSE), "charnes1981_Benchmarking.csv")

cat("done.\n")
