# DEA 1.0.1

Documentation and testing only. **No user-visible behaviour changes**: every
function returns exactly what it returned in 1.0.0.

## Exact property tests

A new test file asserts properties rather than values — statements that must
hold for every estimator, on any data, with no reference answer to compare
against. A value test catches a wrong formula; a property test catches a wrong
*index*, which produces plausible numbers on the fixture it was written against
and survives comparison with another package that makes the same slip.

Every estimator is registered in one table and every applicable property
applied to each: invariance to rescaling each column separately, invariance to
the order of the rows, invariance to duplicating a DMU, monotonicity in the
reference set, agreement between the envelopment and multiplier solutions,
and the cross-model identities tying `dea_ddf()` to `dea()` (`beta = 1 - theta`
for a direction of `x`, `beta = phi - 1` for a direction of `y`). The
estimators that are *not* units invariant — the unweighted additive model and a
fixed-direction directional model, both by design — are asserted to remain so,
since silently gaining an invariance is as much a defect as losing one.

One property is asserted to **fail**, and it is the reason `dea_cross()`
defaults to `secondary = "benevolent"`. Without a secondary goal the
cross-efficiency weights are not determined, and reordering the rows of the
data moves the scores by about 1e-2 — first order relative to the scores
themselves. With either Doyle-Green secondary goal the answer is invariant to
1e-8.

## Fewer dependencies

**`Benchmarking` and `DJL` have been removed from `Suggests`.** A DEA package
should not make rival DEA packages a dependency of any kind: CRAN installs
Suggests in order to check, so their breakage becomes this package's check
failure, and a reader assessing what `DEA` costs to install should not find the
field listed underneath it.

Installing `DEA` now pulls in exactly one package that is not part of R itself
— `lpSolveAPI`, which has no dependencies of its own.

Nothing was lost from the test suite. The two packages were used only to
confirm that this one solves the same linear programs, behind
`skip_if_not_installed()` — which meant the comparison was **silently skipped
on any machine without them**, including most of CRAN's. Their answers are now
recorded in `tests/testthat/reference/`, and the suite compares against those
unconditionally, on every machine. The live comparison still happens, and still
covers every technology and orientation; it just happens in a maintainer script
that is not shipped, and which refuses to update the recorded values if this
package and theirs ever disagree.

## `?dea_boot` now maps onto the published algorithm

Two new sections. **The algorithm, step by step** ties each of the eight steps
of Simar and Wilson's (1998) Algorithm #1 to what the code does, including the
bandwidth (computed on the reflected `2n` points, not on the `n` scores), the
variance correction after smoothing, the folding rule at the boundary, and the
exact form of the interval.

**What is assumed, and when it fails** states, for the first time, that this is
the *homogeneous* bootstrap: it assumes the inefficiency distribution does not
depend on the input-output mix, and where that fails the procedure is not
consistent. The remedy is the heterogeneous bootstrap of Simar and Wilson
(2000), which this package does not implement. The existing refusals for
free-disposal-hull and super-efficiency fits are given their reasons in the
same place.

# DEA 1.0.0

First release.

**This package succeeds `DEA` 0.1-2**, published in 2008 by Zuleyka
Diaz-Martinez and Jose Fernandez-Menendez and archived shortly afterwards. Both
are authors of this package.

The two cover the same ground — the 2008 package had CCR, BCC, additive and
slacks-based models in envelopment and multiplier form — but **no code is
carried over**. The 2008 package bundled a C implementation of the simplex
method from GLPK; this one has no compiled code at all and solves through
**lpSolveAPI**. The implementation is new throughout and the API is entirely
different:

| `DEA` 0.1-2 | here |
|---|---|
| `dea.ccr.io.env()`, `dea.ccr.oo.env()`, `dea.bcc.io.env()`, `dea.bcc.oo.env()` | `dea(rts = "crs" / "vrs", orientation = "in" / "out")` |
| `dea.ccr.io.mul()`, `dea.ccr.oo.mul()`, `dea.bcc.io.mul()`, `dea.bcc.oo.mul()` | `dea(..., multipliers = TRUE)`, then `multipliers(fit)` for `v` and `u` |
| `dea.sbm.ccr()`, `dea.sbm.ccr.io()`, `dea.sbm.ccr.oo()`, `dea.sbm.bcc()`, `dea.sbm.bcc.io()`, `dea.sbm.bcc.oo()` | `dea_sbm(rts = "crs" / "vrs", orientation = "none" / "in" / "out")` |
| `dea.add.env()` | `dea_add(measure = "unweighted")`, plus the normalized RAM and MIP measures |
| `dea.add.mul()` | `dea_add()` for the slacks; `dea(..., multipliers = TRUE)` for the supporting weights |

All sixteen entry points of 0.1-2 are covered, envelopment and multiplier form
alike. Code written against them will not run: the names, the arguments and the
return values are all different. The version number starts at 1.0.0 rather than
0.1.0 because 0.1.0 would rank *below* the archived 0.1-2 under R's version
ordering.

## Prices, and weights

* `dea(..., multipliers = TRUE)` solves the **multiplier (dual) program** at
  each DMU and returns the optimal weights `v` on the inputs, `u` on the
  outputs and the returns-to-scale intercept `u0`, in the caller's own units.
  `multipliers(fit)` puts them in one matrix. For an *efficient* DMU these are
  not unique -- there is a whole face of optimal weight vectors -- and the
  documentation says so rather than leaving it to be discovered.

* `dea_cost()` and `dea_revenue()` -- cost and revenue efficiency, each
  factoring **exactly** into technical times allocative, so a DMU that sits on
  the frontier while buying the wrong mix for its prices is visible as such.

* `dea_profit()` -- Nerlovian profit inefficiency, the profit gap normalized by
  the value of a direction. It *adds* into technical plus allocative rather
  than multiplying, because observed profit is routinely zero or negative and a
  ratio is undefined exactly where the question matters. Constant and
  non-decreasing returns are refused: maximum profit over a cone is unbounded
  as soon as one reference DMU is profitable.

* `dea_cross()` -- cross-efficiency, with the Doyle and Green (1994)
  benevolent and aggressive secondary goals. Plain cross-efficiency depends on
  which vertex of an optimal face the solver stopped at, and the spread is
  first-order rather than numerical; running both goals brackets it. On
  `charnes1981` the constant-returns model calls 19 of 70 sites efficient and
  cannot rank them, and cross-efficiency ranks all 70 with no ties.

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

* `dea_add()` — the additive model of Charnes, Cooper, Golany, Seiford and
  Stutz (1985), with three objective weightings: the Range Adjusted Measure of
  Cooper, Park and Pastor (1999), the Measure of Inefficiency Proportions, and
  the original unweighted total. One program on three scales; the projection
  and the efficient set are identical across them, and that set is exactly the
  Pareto-Koopmans set `dea()` and `dea_sbm()` find.

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

* A stray positional argument is refused rather than absorbed. `data` is the
  third positional argument of every entry point, so `dea_add(x, y, "ram")`
  used to put `"ram"` there and silently take the default measure. It now
  errors.

* Input and output columns are rescaled to mean one before solving. Radial
  efficiency, the slacks-based measure and a proportional-direction
  directional distance are all invariant to this, so no answer changes and the
  linear program is much better conditioned; slacks come back in the caller's
  own units.
