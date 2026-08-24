# DEA 1.0.0

First release.

**This package shares only its name with `DEA` 0.1-2**, which was published in
2008 by Zuleyka Diaz-Martinez and Jose Fernandez-Menendez and archived shortly
afterwards. The name has been reused; no code has. The two cover much the same
ground — the 2008 package had CCR, BCC, additive and slacks-based models in
both envelopment and multiplier form — but nothing carries over, and the API is
entirely different:

| `DEA` 0.1-2 | here |
|---|---|
| `dea.ccr.io.env()`, `dea.ccr.oo.env()`, `dea.bcc.io.env()`, `dea.bcc.oo.env()` | `dea(rts = "crs" / "vrs", orientation = "in" / "out")` |
| `dea.ccr.io.mul()`, `dea.bcc.oo.mul()`, ... | not exposed separately; the multiplier form is the dual of the same program |
| `dea.sbm.ccr()`, `dea.sbm.bcc.io()`, ... | `dea_sbm(rts =, orientation =)` |
| `dea.add.env()`, `dea.add.mul()` | no direct equivalent; see `dea_sbm()` |

Code written against 0.1-2 will not run. The version number starts at 1.0.0
rather than 0.1.0 because 0.1.0 would rank *below* the archived 0.1-2 under R's
version ordering.

## Estimators

* `dea()` — radial Debreu-Farrell technical efficiency under constant
  (`"crs"`), variable (`"vrs"`), non-increasing (`"nirs"`) and non-decreasing
  (`"ndrs"`) returns to scale, and under free disposal (`"fdh"`), in either
  orientation. Optional second-stage slack maximization, Andersen-Petersen
  super-efficiency, and an explicit reference technology through
  `xref`/`yref`.

* `dea_sbm()` — the slacks-based measure of Tone (2001), non-oriented,
  input-oriented or output-oriented. The non-oriented measure is a fractional
  program and is solved exactly by the Charnes-Cooper linearization.

* `dea_ddf()` — the directional distance function of Chambers, Chung and Fare
  (1996), with the direction given as one of five shorthands, a fixed vector,
  or one row per DMU. Accepts zero and negative data, which the radial and
  slacks-based measures cannot.

* `dea_rts()` — scale efficiency and the returns-to-scale classification, from
  a single sweep of the constant-, variable- and non-increasing-returns
  technologies.

## Inference

* `dea_boot()` — the smoothed homogeneous bootstrap of Simar and Wilson (1998):
  bias, bootstrap standard error, bias-corrected efficiency and percentile
  confidence intervals. Reports Simar and Wilson's own `|bias|/se > 1/sqrt(3)`
  rule per DMU as `correct_worthwhile` rather than applying the correction
  silently. Parallel over replications through `parallel::mclapply`.

* `dea_rate()` — the rate a DEA estimator can attain, as the slope of log mean
  squared error on log sample size, given the number of inputs and outputs and
  the returns-to-scale assumption.

* `dea_sim()` — simulate a convex, freely disposable technology whose input-
  and output-oriented Farrell efficiencies have a closed form, so that an
  estimator can be scored against a truth rather than against another
  implementation.

## Data

* `charnes1981` — the Program Follow Through data of Charnes, Cooper and Rhodes
  (1981), 70 school sites with five inputs and three outputs. The original DEA
  application, and the standard worked example in the field. Identical to the
  copies distributed as `Benchmarking::charnes1981` and `npsf::ccr81`, which
  were checked against each other before this one was made.

## Interface

* `x` and `y` accept matrices, data frames, bare numeric vectors, or -- with
  `data` -- character column names or one-sided formulas.

* One result class, `"dea"`, shared by all three estimators, with
  `print()`, `summary()`, `plot()`, `fitted()`, `nobs()`, `efficiency()`,
  `peers()` and `slacks()` methods. `efficiency(fit, "score")` maps every
  model's measure onto a common `(0, 1]` scale, and deliberately refuses to do
  so for a directional model whose direction is neither purely input nor
  purely output, because no such mapping exists.

* `"drs"` and `"irs"` are accepted as aliases for `"nirs"` and `"ndrs"`, the
  spellings **Benchmarking** uses, so switching packages does not silently
  change the model.

* Input and output columns are rescaled to mean one before solving. Radial
  efficiency, the slacks-based measure and a proportional-direction
  directional distance are all invariant to this, so no answer changes and the
  linear program is much better conditioned; slacks come back in the caller's
  own units.
