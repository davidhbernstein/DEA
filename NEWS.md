# DEA 1.0.3

## Free disposal now reports an unsolved DMU as unsolved

Reported by GitHub Copilot on issue #1, and confirmed: `rts = "fdh"` returned
`integer(n)` as its status vector, so **every** DMU looked like an optimal
solve -- including those with no dominating peer, whose efficiency is `NA`
because the score does not exist.

The effect was that `$status` could not distinguish a valid free-disposal score
from one that could not be computed, and `dea()` had nothing to warn about. On
identical data both `rts = "vrs"` and `dea_sbm()` returned status 2 and warned;
free disposal was alone in staying silent. It now returns 2 and warns like
everything else, with wording that says free disposal needs one dominating
reference DMU rather than a combination of them.

This can only fire against an external `xref`/`yref`. A self-referenced fit is
unaffected, because every DMU dominates itself.

`?dea` now documents what the three status codes mean rather than only what 0
means.

## A startup message, and a citation that keeps up with the version

Attaching the package now prints its version and points at `citation('DEA')`,
matching `sfa`.

`inst/CITATION` used to write the version out by hand, so it went on claiming
1.0.0 after the package moved past it. It now reads `meta$Version` from the
DESCRIPTION, which is what makes the startup message's advice true. Author
names are also split into given and family parts, so BibTeX no longer braces
each one as a literal and author-year styles cite "Bernstein et al." rather
than "David H. Bernstein et al.".

## Why you cannot test the bootstrap's homogeneity assumption by regression

`?dea_boot` has said since 1.0.1 that the Simar-Wilson (1998) bootstrap assumes
inefficiency does not depend on the input-output mix. It now also says what
happens if you try to check that by regressing the estimated scores on the
inputs and outputs: the test rejects on data that satisfy the assumption, in
100 per cent of 200 samples drawn from a DGP that satisfies it by construction.

The estimator manufactures the association. Scores are all measured against one
frontier estimated from the whole sample, and the DEA bias is larger in sparse
regions than in dense ones, so the correlation between an input and the
*estimated* score is -0.36 where the correlation with the *true* score is +0.04.

A null simulated from the fitted homogeneous model removes that objection and
still does not give a usable test. The measurement is kept, reproducible, in
`horserace/homogeneity_experiment.R` of the development repository, along with
why neither route works.

## Complementary slackness now ties the peers to the multipliers

The `u0` sign table is the part of this package most likely to be silently
wrong: getting a sign backwards does not error, it solves a different
technology's program and returns a plausible number. It had two guards. Neither
was enough on its own.

Feasibility of the weights has **no** teeth against a sign error at all --- with
the table deliberately inverted the weights stay feasible to 1e-12, because
they are then a valid price vector for a different technology. Value agreement
against the envelopment form does catch it, but only where the returns-to-scale
restriction binds, which is why that test carries a `skip_if`.

The new check has teeth everywhere. For a linear program, complementary
slackness holds between *any* optimal primal solution and *any* optimal dual
solution, so a reference DMU carrying positive `lambda` must have its multiplier
constraint exactly tight --- even though the two programs are solved
independently here, and even though the weights are non-unique at an efficient
DMU. Measured against an inverted table, the residual at a positive `lambda`
goes from 1e-12 to between 0.4 and 1.6 in all four of nirs/ndrs by in/out.

The feasibility test also now runs in both orientations rather than the input
one alone; the output form's constraint has the opposite sense, and it was
untested. And the same slackness check now runs against a reference set of a
different size from the evaluated set --- 12 scored against 50 --- because with
`n == nref` every place that indexes a reference column by `n` still works,
which is how one such bug once survived every self-referenced test.

## A fourth silent failure, this time in the multiplier program

Found while adding the complementary slackness check above, which is the point
of adding it.

1.0.2 fixed two places where a linear program failed and said nothing, and
1.0.3 fixed a third. The multiplier (dual) sweep was the fourth and last: it
had neither the retry nor the reporting. A dual that failed returned `NA`
weights while `status` --- which describes the *envelopment* program --- went
on saying 0. The dual's own code was recorded in `mult_status`, and nothing
ever read it.

**Only super-efficiency fits were affected.** Ordinary fits show no dual
failures at all, including at n = 800 under either technology.

**The retry.** The same lpSolveAPI basis carryover as in 1.0.2: on a 40-DMU
constant-returns super-efficiency fit, four DMUs came back infeasible on the
reused program and solved to optimality on a fresh one. Note this is retried
where the envelopment sweep deliberately does *not* retry an infeasible
program. In the dual that reasoning inverts --- if the envelopment program is
feasible and bounded, strong duality makes its dual feasible and bounded too,
so an infeasible dual is a contradiction rather than an answer.

**The reporting** masks the documented case: under super-efficiency an
infeasible envelopment program has an unbounded dual, and warning about it
would be noise about an answer already given as `NA`. The mask is on the
*pairing* --- a dual failure whose envelopment program was infeasible --- and
not on the status code, because the solver does not name that case
consistently: one DMU reported it as 5 on five successive fresh programs. A
dual that fails where the envelopment program solved is still reported.

`?dea` now documents `mult_status`, which was returned but never described.

## Three more bugs, from auditing every solve site

The multiplier fix above was the fourth program that could fail and say
nothing. That made the pattern worth chasing to the end rather than one more
time, so every one of the package's ten linear-program solve sites was checked
for the three things the earlier four were missing: the status recorded, a
retry on a fresh program, and the failure reported.

**`dea_ddf()` never reported an unsolved DMU.** The fifth instance, and
reachable from ordinary code rather than a contrived program:
`direction = "in"` holds the output direction at zero, so a DMU whose outputs
exceed everything the reference set can produce has no feasible point at all.
It returned `beta = NA` in silence. Established the same way as the free
disposal bug in 1.0.3 --- by contrast on identical data, where `dea()` warned
and `dea_ddf()` did not, with both reporting the same status vector.

**The additive sweep lost a DMU it could have solved.** It reported failures
but never retried them, and the failure is the same lpSolveAPI basis carryover
fixed elsewhere in 1.0.2: one DMU of 1200 under variable returns gave status 5
on the reused program and status 0 on a fresh one, five times out of five.

This changes a documented behaviour. The unweighted measure at a 1000:1 column
spread used to lose a DMU outright; it is now recovered, and the recovered
answer is right --- solved independently under four `lpSolve` scaling modes the
objective is 8e-9 to 2.5e-8, zero to within the noise that spread produces, and
the DMU is efficient. The conditioning warning still fires. Users simply no
longer lose a DMU on top of it.

**`dea_cross()` returned `NaN` and `-Inf` for a DMU nobody could appraise.**
Every rater's weights value an all-zero-input DMU at zero, so the
cross-efficiency ratio has no denominator, and averaging with `na.rm = TRUE`
turned that into `NaN` with `spread` of `-Inf` --- explained to the user only
by base R's "no non-missing arguments to max". Both are now `NA` with a warning
that says what happened. A DMU appraised by only *part* of the panel is now
warned about too: that case produced not a missing answer but a plausible wrong
one, a mean over fewer raters reported as though the whole panel had spoken.

Suite 839 -> 857 assertions, none skipped.

## The retry guard now covers every sweep, and is tested directly

1.0.3 left `dea_sbm()` and the price models without the rebuild-on-failure
guard, on the grounds that no failure could be produced in either. **That was
wrong, and wrong because the search was incomplete:** it exercised the
non-oriented slacks-based program and never the oriented one. A wider search
finds the oriented measure failing readily --- nine configurations, including
n = 150 at seed 12, where DMU 108 returns lpSolve status 5 on the reused
program.

As everywhere else, the failure is recoverable and the recovered answer is
right: solved independently under four `lpSolve` scaling modes, that DMU's rho
is exactly 1 in all four, and the package now reports exactly 1 where it
previously reported `NA`.

**The guard is now one function.** `.lp_solve_retry()` holds the rule that was
being re-derived at each sweep, including which statuses are answers rather
than failures --- infeasible is an answer for an envelopment program and is not
retried, while for the multiplier program an infeasible dual contradicts strong
duality and is. It is used by the slacks-based, additive, directional, cost,
revenue and profit sweeps.

**And it is tested directly**, not only through a sweep. A rare branch whose
first real execution is the day it matters has never been run, so the retry is
driven with one-variable programs of known status: a solved program is not
rebuilt, an infeasible one is not rebuilt, anything else is rebuilt once, and
the caller receives the *fresh* program rather than the failed one.

Nothing else moved. On 800 DMUs the cost, revenue, profit and non-oriented
slacks-based scores are bit-identical before and after; the only changed number
in the package is the one DMU that used to be lost.

## plot() now draws the object, on every class

`plot.dea` was the only plot method, and what happened to the other five
classes was worse than having none.

`dea_price`, `dea_cross` and `dea_sim` all carry components named `x` and `y`
--- the input and output matrices. R's generic therefore fell through to
`plot.default`, which **found** them and silently drew the raw inputs against
the raw outputs: a plausible-looking scatter with nothing to do with a price
decomposition or a cross-efficiency ranking. `dea_rts` and `dea_boot` have no
such components, so the same fall-through failed instead with *"'x' is a list,
but does not have components 'x' and 'y'"* --- a true statement about
`plot.default` and a baffling one about a returns-to-scale object.

Five methods now draw what each object is for, and each returns the plotted
data invisibly so it can be redrawn another way:

* `plot.dea_boot` --- the standard caterpillar: DMUs sorted, interval as a
  segment, bias-corrected estimate on it, and the raw score marked too, because
  the distance between the two *is* the bias being corrected.
* `plot.dea_rts` --- scale efficiency against `sum(lambda)` under constant
  returns, coloured by class. That puts the classification and its own evidence
  in one picture: the rule is that `sum(lambda)` below 1 is the
  increasing-returns region, and an early version of this package had that
  inverted.
* `plot.dea_cross` --- the range of appraisals each DMU receives, with its
  self-appraisal marked. The mean is what `eff` already reports; the spread is
  what the matrix adds.
* `plot.dea_price` --- technical against allocative efficiency, with contours
  of constant overall efficiency taken from the data rather than a fixed
  ladder, which would otherwise fall off the panel whenever efficiency is
  confined near 1.
* `plot.dea_sim` --- the sample under its true frontier for one input and one
  output, and the distribution of true efficiency otherwise.

All five handle the `NA`s these objects are documented to produce: an
infeasible super-efficiency DMU, or a rated DMU that no rater could appraise.

Found by measuring test coverage rather than by guesswork --- the package sits
at 88%, and the print, summary and plot methods were most of what the suite
never reached.

# DEA 1.0.2

Two silent defects in the radial path, both in how solver failures were
handled. No interface changes.

## Failed DMUs were not reported

`dea()` was the only entry point that never called the internal reporter that
`dea_sbm()`, `dea_add()` and the price models have always used. A DMU whose
linear program failed came back as `NA` with nothing said about it.

Worse, `.dea_radial()` never read the **second stage's** status at all. The
radial score is found first and the slacks by a second program from that
projection; if the second one failed, its `NA` slacks were stored while
`status` kept the first stage's `0`. The fit then reported a DMU as solved,
with missing slacks, and `efficient` quietly `NA`.

Both are fixed: stage two's status is recorded, and `dea()` reports. Under
`super = TRUE` an infeasible program is the documented result rather than a
fault, so those stay silent.

## The failures themselves are fixed

Those two bugs were hiding **56 of 1200 DMUs failing numerically** on a
variable-returns fit — about 5%, growing with `n`, and absent under constant
returns.

The cause is not tolerance. The package builds one linear program per
(technology, orientation) and rewrites only the right-hand side for each DMU,
which is where its speed comes from — but `lpSolveAPI` carries basis and
factorisation state on that object, and for some right-hand sides the inherited
state is bad enough that the solve gives up. Of those 56 failures, loosening
`epsel` from 1e-12 to 1e-9 fixed one, every scaling mode fixed at most a
quarter, and `guess.basis()` fixed one — while rebuilding the object fixed all
56.

So a failed solve is now retried once on a fresh program. Failures are rare
enough that the rebuild costs little, and none of the tested designs — up to
n = 1200 across five technologies and both orientations — now leaves any DMU
unsolved.

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
