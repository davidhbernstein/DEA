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
