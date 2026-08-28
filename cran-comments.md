## Resubmission

This is a resubmission of `DEA` 1.0.0, rejected on 2026-08-25 with:

> Is this related to the archived package?
> If yes: Rejected as copyright and authors are not correctly cited.
> If no: Rejected: Please choose a unique name.

**The answer is yes**, and the citation has been corrected.

`DEA` 0.1-2 was published in February 2008 by Zuleyka Diaz-Martinez and Jose
Fernandez-Menendez and archived thereafter. This package is its successor: it
covers the same models -- CCR, BCC, additive and slacks-based efficiency -- and
takes over the name with the original authors' agreement. I contacted them and
both have given permission to be added to the package. I can forward that
correspondence on request.

What has changed since the rejected version:

* **Both original authors are now listed as authors** in `Authors@R`, each with
  a `comment` recording precisely what they authored:

      Zuleyka Diaz-Martinez [aut] (author of the original DEA package (2008),
        which this package succeeds)
      Jose Fernandez-Menendez [aut] (author of the original DEA package (2008),
        which this package succeeds)

* **The Description field states the relationship**, so a user arriving from
  the archive is told what this is without having to open NEWS.

* **`NEWS.md` maps the old API onto the new one** function by function, and
  all sixteen entry points of 0.1-2 are now covered, envelopment and multiplier
  form alike. The multiplier form was the one genuine gap when 1.0.0 was first
  submitted: the score was there, since the two programs are duals and attain
  the same value, but the optimal weights `v` and `u` that
  `dea.ccr.io.mul()` and its siblings returned were not. `dea(..., multipliers
  = TRUE)` and `multipliers()` now return them. The mapping was checked
  against the archived tarball function by function, not from memory.

### On copyright

No copyright holder is declared beyond the authors, and deliberately so. **No
code from 0.1-2 is included** -- not adapted, not translated, not referenced.
The old package bundled a C implementation of the simplex method taken from
GLPK; this one has no compiled code at all and solves through **lpSolveAPI**.
The implementation is new throughout, so there is no third-party copyright to
declare. The relationship is one of succession in name, subject and
maintainership, not of derived code, and I would rather state that plainly than
have it inferred.

If you would prefer the original authors recorded as `ctb` rather than `aut`,
or a `cph` entry added, I am happy to make either change.

## Test environments

* local macOS 26.5 (aarch64), R 4.5.2
* GitHub Actions, all passing:
  - macOS-latest, R release
  - windows-latest, R release
  - ubuntu-latest, R devel
  - ubuntu-latest, R release
  - ubuntu-latest, R oldrel-1

## R CMD check results

0 errors | 0 warnings | 1 note.

The note is from the incoming check and reads:

    New submission
    Package was archived on CRAN

Both lines are expected. There are no notes arising from the package itself.

A local run also emits notes for `checking HTML version of manual` (this
machine's `tidy` is too old and `V8` is not installed) and occasionally
`unable to verify current time`. Both are properties of the machine, not of the
package.

## Other notes

* The version is 1.0.0 rather than 0.1.0 because R orders `0.1.0` *below* the
  archived `0.1-2`, so a first release numbered 0.1.0 would have been a
  downgrade.

* The package solves n linear programs per call through **lpSolveAPI** and has
  no compiled code of its own.

* Examples and the vignette use small samples and low bootstrap replication
  counts deliberately, so the whole check runs in well under a minute.
  `dea_boot()`'s example uses `B = 50` and the vignette `B = 100`, each noting
  that 2000 is the usual recommendation. The one slower example is marked
  `\donttest{}`.

* The bundled data set `charnes1981` is the Program Follow Through data of
  Charnes, Cooper and Rhodes (1981), the original DEA application. The same
  figures are already distributed on CRAN in `Benchmarking` (as
  `charnes1981`) and `npsf` (as `ccr81`), both GPL-2; the two copies were
  checked against each other and agree on every value before this one was made.
  Columns are renamed and the programme indicator is logical rather than
  integer, but no measurement has been altered, and the source is credited in
  the help page.

* `Benchmarking` and `DJL` appear in Suggests only. They are used in the test
  suite to cross-check that this package's linear programs agree with
  independent implementations of the same models, and every such test is
  guarded by `skip_if_not_installed()`.
