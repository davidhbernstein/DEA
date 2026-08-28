## Resubmission

This is a resubmission of `DEA` 1.0.0, rejected on 2026-08-25 with:

> Is this related to the archived package?
> If yes: Rejected as copyright and authors are not correctly cited.
> If no: Rejected: Please choose a unique name.

**The answer is yes**, and the authorship has been corrected.

`DEA` 0.1-2 was published in February 2008 by Zuleyka Diaz-Martinez and Jose
Fernandez-Menendez and archived thereafter. **Both are authors of this package.**
I contacted them after the rejection, and they have joined the project rather
than merely consenting to the name being reused; `Authors@R` now reads

    David Bernstein [aut, cre] (ORCID: 0000-0002-2267-5741)
    Zuleyka Diaz-Martinez [aut]
    Jose Fernandez-Menendez [aut]

with no qualifying comment, because none is warranted: they are authors of this
package on the same footing as I am, not merely of the one it succeeds. I can
forward the correspondence on request.

The rest of what has changed since the rejected version:

* **The Description field states the succession**, so a user arriving from the
  archive is told what this is without having to open NEWS.

* **`NEWS.md` maps the old API onto the new one** function by function, and
  all sixteen entry points of 0.1-2 are now covered, envelopment and multiplier
  form alike. The multiplier form was the one genuine gap when 1.0.0 was first
  submitted: the score was there, since the two programs are duals and attain
  the same value, but the optimal weights `v` and `u` that `dea.ccr.io.mul()`
  and its siblings returned were not. `dea(..., multipliers = TRUE)` and
  `multipliers()` now return them. The mapping was checked against the archived
  tarball function by function, not from memory.

### On copyright

No copyright holder is declared beyond the authors, and deliberately so. **No
code from 0.1-2 is carried over** -- not adapted, not translated, not
referenced. The 2008 package bundled a C implementation of the simplex method
taken from GLPK; this one has no compiled code at all and solves through
**lpSolveAPI**. The implementation is new throughout, so there is no
third-party copyright to declare, and I would rather state that plainly than
have it inferred from the shared authorship.

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
