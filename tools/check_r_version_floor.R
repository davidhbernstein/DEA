## ---------------------------------------------------------------------------
## Does the package actually run on the oldest R it CLAIMS to support?
##
## `DESCRIPTION` says `Depends: R (>= 4.0.0)`.  Nothing in the CI matrix tests
## that: the oldest cell is `oldrel-1`, which is two releases behind current and
## therefore years NEWER than the declared floor.  So the floor is an
## undefended promise, and the failure mode is silent on every machine we own.
##
## It has already been broken once.  `tests/testthat/test-interface.R` used the
## native pipe `|>`, which is R 4.1.0 syntax, so on R 4.0.x that file did not
## fail a test -- it failed to PARSE, taking all of its assertions with it
## before one of them ran.  Found by reading, not by running, which is the
## point of this script.
##
## WHY THE PARSER AND NOT grep.  `|>` and `\(x)` vanish at parse time: the
## pipe becomes an ordinary call and the lambda an ordinary `function`, so
## nothing survives into the fitted object to test for.  They have to be caught
## in the SOURCE.  But a grep over this package's source is hopeless -- the
## comments are prose, and a sentence mentioning the pipe, or a warning string
## containing a backslash, reads identically to the real thing.  So the tokens
## come from `getParseData()`, which classifies a `|>` inside a string as
## STR_CONST and one inside a comment as COMMENT, and only the genuine operator
## as PIPE.
##
## WHAT THIS DOES NOT CATCH.  Syntax it knows about, and a hand-kept list of
## base functions, below -- and R ships no database of when each base function
## appeared, so that list is a floor on what is checked, not a ceiling.  A
## function added to base R in 4.3 and used here would pass this script and
## still break a 4.0 user.  The only complete check is an actual old R; this is
## the cheap check that runs today.
##
## Not shipped: `tools/` is in `.Rbuildignore`.  Run from the package root:
##     Rscript tools/check_r_version_floor.R
## Exits non-zero if anything in the package needs a newer R than it declares.
## ---------------------------------------------------------------------------

## Syntax, by the R version that introduced it, keyed on getParseData()'s token
## name.  The lambda token really is the 4-character string '\\'.
SYNTAX <- list(
  list(token = "PIPE",   version = "4.1.0", what = "native pipe `|>`"),
  list(token = "'\\\\'", version = "4.1.0", what = "lambda shorthand `\\(x)`")
)

## Base functions newer than 4.0.0 that are plausible to reach for here.  See
## the caveat above: this list is hand-kept and deliberately short.
NEWER_FUNS <- c(
  "...names"       = "4.1.0",
  "Sys.setLanguage"= "4.2.0",
  "array2DF"       = "4.3.0",
  "sort_by"        = "4.4.0",
  "%||%"           = "4.4.0",
  "grepv"          = "4.5.0"
)

## --- the declared floor ------------------------------------------------------
dcf <- read.dcf("DESCRIPTION")
dep <- if ("Depends" %in% colnames(dcf)) dcf[1, "Depends"] else ""
m   <- regmatches(dep, regexpr("R \\(>=\\s*([0-9.-]+)\\)", dep))
if (!length(m)) stop("DESCRIPTION declares no `R (>= ...)` in Depends.")
floor_str <- gsub("[^0-9.-]", "", sub("R \\(>=", "", m))
floor_ver <- as.numeric_version(floor_str)
cat("Declared floor: R >=", format(floor_ver), "\n")

## --- everything with R code in it -------------------------------------------
## Rd examples are included because R CMD check RUNS them, so 4.1 syntax in an
## example breaks a 4.0 check exactly as it would in R/.  Rd2ex writes each
## file's examples out as plain R, which is then parsed like any other source.
files <- c(Sys.glob("R/*.R"), Sys.glob("tests/*.R"), Sys.glob("tests/testthat/*.R"))
rd <- Sys.glob("man/*.Rd")
ex_of <- character(0)
for (f in rd) {
  out <- tempfile(fileext = ".R")
  ok <- tryCatch({ tools::Rd2ex(f, out); TRUE }, error = function(e) FALSE)
  if (ok && file.exists(out) && file.size(out) > 0) { files <- c(files, out); ex_of[out] <- f }
}

## Vignettes: pull the R out of fenced chunks.  Crude on purpose -- it only has
## to be good enough to parse, and a chunk it mangles shows up as a parse error
## here rather than passing silently.
for (v in c(Sys.glob("vignettes/*.Rmd"), Sys.glob("vignettes/*.Rnw"))) {
  ln <- readLines(v, warn = FALSE)
  open <- grepl("^\\s*```\\{r", ln); close <- grepl("^\\s*```\\s*$", ln)
  keep <- logical(length(ln)); inside <- FALSE
  for (i in seq_along(ln)) {
    if (open[i]) { inside <- TRUE; next }
    if (inside && close[i]) { inside <- FALSE; next }
    keep[i] <- inside
  }
  if (any(keep)) {
    out <- tempfile(fileext = ".R"); writeLines(ln[keep], out)
    files <- c(files, out); ex_of[out] <- v
  }
}

label <- function(f) if (!is.na(ex_of[f]) && nzchar(ex_of[f])) ex_of[f] else f

## Names the package DEFINES is not the same as names it borrows from base, and
## the difference decides the answer.  `%||%` entered base R in 4.4.0, and this
## package uses it -- but it also defines its own at the top of R/utils.R, which
## shadows base's, so it costs nothing on R 4.0.  Vendoring a newer base
## function is the ordinary way to hold an old floor, so a checker blind to the
## vendored copy reports exactly the wrong thing.  That was this script's first
## finding on its first run, and it was a false positive.
defined <- character(0)
for (f in Sys.glob("R/*.R")) {
  for (e in as.list(parse(f))) {
    if (is.call(e) && length(e) >= 3L &&
        as.character(e[[1L]]) %in% c("<-", "=", "<<-") &&
        (is.symbol(e[[2L]]) || is.character(e[[2L]]))) {
      defined <- c(defined, as.character(e[[2L]]))
    }
  }
}
vendored <- intersect(names(NEWER_FUNS), defined)
if (length(vendored)) {
  cat("Vendored in R/ (so the floor holds):", paste(vendored, collapse = ", "), "\n")
  NEWER_FUNS <- NEWER_FUNS[setdiff(names(NEWER_FUNS), vendored)]
}

## --- scan --------------------------------------------------------------------
bad <- data.frame(file = character(0), line = integer(0),
                  needs = character(0), what = character(0))
for (f in files) {
  pd <- tryCatch(getParseData(parse(f, keep.source = TRUE)),
                 error = function(e) { cat("  ! could not parse", label(f), ":",
                                           conditionMessage(e), "\n"); NULL })
  if (is.null(pd) || !nrow(pd)) next
  tok <- pd[pd$terminal, , drop = FALSE]

  for (s in SYNTAX) {
    hit <- tok$token == s$token
    if (any(hit)) bad <- rbind(bad, data.frame(
      file = label(f), line = tok$line1[hit], needs = s$version, what = s$what))
  }
  ## Function calls and infix operators, by name.  SYMBOL_FUNCTION_CALL covers
  ## f(), SPECIAL covers %||% and friends.
  nm <- tok$text[tok$token %in% c("SYMBOL_FUNCTION_CALL", "SPECIAL")]
  ln <- tok$line1[tok$token %in% c("SYMBOL_FUNCTION_CALL", "SPECIAL")]
  hit <- nm %in% names(NEWER_FUNS)
  if (any(hit)) bad <- rbind(bad, data.frame(
    file = label(f), line = ln[hit], needs = unname(NEWER_FUNS[nm[hit]]),
    what = paste0("`", nm[hit], ifelse(grepl("^%", nm[hit]), "", "()"), "`")))
}

## --- verdict -----------------------------------------------------------------
if (nrow(bad)) bad <- bad[as.numeric_version(bad$needs) > floor_ver, , drop = FALSE]
cat("Scanned", length(files), "file(s) of R code.\n")
if (!nrow(bad)) {
  cat("OK: nothing used needs more than R", format(floor_ver), ".\n")
  quit(status = 0)
}
bad <- bad[order(bad$file, bad$line), ]
cat("\nFAIL:", nrow(bad), "use(s) of features newer than the declared floor:\n\n")
for (i in seq_len(nrow(bad))) {
  cat(sprintf("  %s:%d  %s  needs R >= %s\n",
              bad$file[i], bad$line[i], bad$what[i], bad$needs[i]))
}
cat("\nEither rewrite these, or raise Depends in DESCRIPTION to R (>= ",
    format(max(as.numeric_version(bad$needs))), ").\n", sep = "")
quit(status = 1)
