## Test environments

* local macOS 26.5 (aarch64), R 4.5.2
* GitHub Actions: macOS-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes

## This submission reuses an archived package name

`DEA` 0.1-2 was published in February 2008 by Zuleyka Diaz-Martinez and Jose
Fernandez-Menendez and archived thereafter; there has been no release in
eighteen years. This submission reuses that name for an entirely new package.

I want to be explicit about what that means, so nothing is discovered later:

* **The code is unrelated.** Nothing from the 2008 package is reused, adapted
  or referenced. There is no shared authorship and no shared source. The old
  package also carried a bundled C implementation of the simplex algorithm
  taken from GLPK; this one has no compiled code at all and solves through
  **lpSolveAPI**.

* **The subject matter is close, which is the reason for wanting the name.**
  `DEA` 0.1-2 provided CCR, BCC, additive and slacks-based models in both
  envelopment and multiplier form. This package provides radial CCR/BCC
  efficiency (plus non-increasing, non-decreasing and free-disposal
  technologies), the slacks-based measure, the directional distance function,
  and the Simar-Wilson bootstrap. It is a successor in topic, not in code.

* **The API is not backwards compatible**, and cannot be. The old package
  exported sixteen functions of the form `dea.ccr.io.env()`,
  `dea.sbm.bcc.oo()` and so on; none of them exist here, where the same choices
  are arguments to `dea()` and `dea_sbm()`. Anyone with code written against
  0.1-2 will find that it does not run. Given that the package has been off
  CRAN since 2008, I judge the population of such code to be effectively empty,
  but the incompatibility is real and I would rather state it plainly.
  `NEWS.md` carries a table mapping the old names onto the new ones.

* **The version is 1.0.0**, which is greater than the archived 0.1-2. I noticed
  in checking that an apparently natural 0.1.0 for a first release would in
  fact have been a *downgrade* under R's version ordering.

* `NEWS.md` opens by saying all of the above, so a user arriving from the
  archive is not misled.

If you would prefer this not reuse the archived name, I am happy to rename and
resubmit — `deafit`, `deatools`, `npdea` and `deainfer` are all free, and the
change is mechanical. I would also gladly make contact with the original
maintainers first if that is the process you would rather I follow.

## Other notes

* The package solves n linear programs per call through **lpSolveAPI** and has
  no compiled code of its own.

* Examples and the vignette use small samples and low bootstrap replication
  counts deliberately, so the whole check runs in well under a minute.
  `dea_boot()`'s example uses `B = 50` and the vignette `B = 100`, each with a
  note that 2000 is the usual recommendation in practice. The one slower
  example is marked `\donttest{}`.

* The bundled data set `charnes1981` is the Program Follow Through data of
  Charnes, Cooper and Rhodes (1981), the original DEA application. The same
  figures are already distributed on CRAN in `Benchmarking` (as
  `charnes1981`) and `npsf` (as `ccr81`), both GPL-2; the two copies were
  checked against each other and agree on every value before this one was
  made. Columns are renamed and the programme indicator is logical rather than
  integer, but no measurement has been altered, and the source is credited in
  the help page.

* `Benchmarking` and `DJL` appear in Suggests only. They are used in the test
  suite to cross-check that this package's linear programs agree with
  independent implementations of the same models, and every such test is
  guarded by `skip_if_not_installed()`.
