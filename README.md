# DEA

<!-- badges: start -->
[![R-CMD-check](https://github.com/davidhbernstein/DEA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/davidhbernstein/DEA/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/DEA)](https://CRAN.R-project.org/package=DEA)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://cran.r-project.org/web/licenses/GPL-2)
<!-- badges: end -->

**Data envelopment analysis in R.** Five nonparametric efficiency estimators
behind one interface and one result object; cost, revenue and profit efficiency
where prices are known, and cross-efficiency where they are not; and the
sampling theory that DEA practice usually leaves out: bias correction,
confidence intervals, and an honest statement of how fast any of it can
converge.

By David H. Bernstein, Zuleyka Díaz-Martínez and José Fernández-Menéndez.
A companion to [`sfa`](https://github.com/davidhbernstein/sfa), which does the
parametric half of the same problem.

> **On the name.** This package succeeds `DEA` 0.1-2, published in 2008 by
> Zuleyka Díaz-Martínez and José Fernández-Menéndez and archived shortly
> afterwards; both are authors here. It covers the same models, but no code is
> carried over — the 2008 package bundled a GLPK simplex in C, this one has no
> compiled code at all, and the API is entirely different. See
> [`NEWS.md`](NEWS.md) for a table mapping the old function names onto the new
> arguments.

## Installation

Not yet on CRAN — a first submission is pending. From GitHub:

```r
# install.packages("remotes")
remotes::install_github("davidhbernstein/DEA")
```

Once it is accepted, `install.packages("DEA")`.

## Quick start

```r
library(DEA)

## The original DEA application: 70 school sites, 5 inputs, 3 outputs
## (Charnes, Cooper and Rhodes 1981), bundled with the package.
data(charnes1981)
x <- charnes1981[, paste0("x", 1:5)]
y <- charnes1981[, paste0("y", 1:3)]

fit <- dea(x, y, rts = "crs", orientation = "in")
fit
#> --- Data envelopment analysis ---
#> model:       radial (Debreu-Farrell)
#> technology:  CRS, input orientation
#> DMUs: 70   inputs: 5   outputs: 3
#>
#> theta (input-oriented, 1 = on the frontier)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#>  0.7883  0.9002  0.9404  0.9378  1.0000  1.0000
#>
#> efficient DMUs: 19 of 70

head(efficiency(fit))   # theta, in (0, 1]
head(peers(fit))        # who each DMU is benchmarked against
head(slacks(fit))       # what the radial cut leaves behind
dea_rts(x, y)           # scale efficiency and returns to scale
```

`x` and `y` can be matrices, data frames, bare vectors, or — with `data` —
column names or one-sided formulas:

```r
dea(~ x1 + x2 + x3 + x4 + x5, ~ y1 + y2 + y3,
    data = charnes1981, rts = "crs")
```

For checking an estimator rather than using one, `dea_sim()` draws a technology
whose true efficiency is known in closed form:

```r
sim <- dea_sim(n = 200, p = 2, q = 1, returns = 0.9, seed = 1)
fit <- dea(sim$x, sim$y, rts = "vrs", orientation = "out")
mean(fit$eff - sim$phi)   # the bias, against a truth written down first
plot(fit)
```

## What is here

| function | what it estimates |
|---|---|
| `dea()` | radial Debreu–Farrell efficiency: CCR, BCC, non-increasing and non-decreasing returns, and the free disposal hull; either orientation; second-stage slacks; Andersen–Petersen super-efficiency |
| `dea_sbm()` | the slacks-based measure of Tone (2001), non-oriented or oriented |
| `dea_ddf()` | the directional distance function of Chambers, Chung and Färe (1996), which handles zero and negative data |
| `dea_add()` | the additive model of Charnes et al. (1985), unweighted or as the Range Adjusted Measure |
| `dea_rts()` | scale efficiency and the returns-to-scale classification |
| `dea_cost()`, `dea_revenue()` | cost and revenue efficiency, each factoring exactly into technical × allocative |
| `dea_profit()` | Nerlovian profit inefficiency, which decomposes by addition rather than multiplication |
| `dea_cross()` | cross-efficiency, with the Doyle–Green benevolent and aggressive secondary goals |
| `multipliers()` | the optimal weights `v`, `u`, `u0` from the multiplier form of any radial fit |
| `dea_boot()` | bias correction and confidence intervals, by the Simar–Wilson (1998) smoothed homogeneous bootstrap |
| `dea_sim()` | a technology with a closed-form answer, for testing an estimator against a truth |
| `dea_rate()` | the rate at which any of this can converge |
| `charnes1981` | the Program Follow Through data the CCR model was introduced on |

All the estimators return an object of class `"dea"`, so `print()`,
`summary()`, `plot()`, `efficiency()`, `peers()`, `slacks()`, `multipliers()`
and `fitted()` work the same way across them. The price models and
cross-efficiency return their own classes, with the same methods where the same
question makes sense.

```r
## On the frontier, but buying the wrong mix for the prices you face?
ce <- dea_cost(x, y, w = c(1.4, 0.9, 2.1, 1.2, 1.0), rts = "crs")
sum(ce$technical > 1 - 1e-9)   # 19 sites are technically efficient
sum(ce$eff == 1)               # 2 of them are also allocatively efficient

## 19 sites tie at 1 and cannot be ranked. Cross-efficiency ranks all 70.
ben <- dea_cross(x, y, secondary = "benevolent")
agg <- dea_cross(x, y, secondary = "aggressive")
range(ben$eff - agg$eff)       # how much the answer depends on whose weights
```

Cross-efficiency is reported through a *pair* of secondary goals on purpose.
An efficient DMU has a whole face of optimal weight vectors, all giving it the
same score of 1 but giving everyone else different ones, so a single
cross-efficiency number is an artefact of the vertex a solver stopped at. The
benevolent and aggressive answers bracket it: a ranking that survives from one
to the other is in the data, and one that does not was in the solver.

## Four things this package insists on

**A DEA score is an estimate.** The frontier is spanned by the observed DMUs,
so it lies inside the true one and every score is biased toward 1. On a
200-unit simulated sample the mean bias is around 0.05 on a (0, 1] scale — not
a rounding error. `dea_boot()` estimates it:

```r
b <- dea_boot(fit, B = 2000, seed = 1)
head(b$table)   # eff, bias, se, bias_corrected, ci_lower, ci_upper
```

It also reports `correct_worthwhile`, which applies Simar and Wilson's own
warning: the correction removes a bias and adds the variance of the estimate of
that bias, so it is a net loss where the bias is small relative to the noise.
The package reports that judgement per DMU rather than making it silently.

**Convergence is slow, and slower with every variable added.** DEA is not
root-*n* consistent. `dea_rate()` gives the slope that log mean squared error
can attain against log sample size:

```r
dea_rate(1, 1, "vrs")   # -1.333
dea_rate(4, 4, "vrs")   # -0.444
```

At four inputs and four outputs, a hundredfold increase in sample size buys
about a factor of 8 in mean squared error where a parametric estimator would
buy 100. `dea()` warns when `n < 3(p+q)`, where most DMUs are efficient by
dimension alone rather than by performance.

**Radial efficiency is not Pareto–Koopmans efficiency.** A DMU can score 1 and
still be dominated, because no *common* factor cuts every input even though one
input alone could come down. The second-stage slack program is on by default,
`print()` reports both counts when they differ, and `$efficient` is the
stricter judgement. `dea_sbm()` folds the distinction into a single number.

**Peer sets are not unique.** These linear programs frequently have multiple
optima, so `peers()` returns one valid benchmark set rather than *the*
benchmark set. Scores are unique; the weights achieving them generally are not,
and no DEA package can make them so. `summary()` says this out loud.

## Correctness

Efficiency scores are checked three ways, in the package's own test suite:

* against **hand-worked answers** on a five-DMU example, so at least one check
  depends on no other package agreeing;
* against **`Benchmarking`** for every radial technology and orientation, and
  against **`DJL`** for all three slacks-based orientations — agreement to
  around 1e-12, which is linear-programming tolerance;
* against **`dea_sim()`'s closed-form truth**, which is the only check that
  catches a misunderstanding rather than a coding error. Two packages making
  the same mistake agree perfectly.

Two further pipelines live in the development repository rather than the
package: a convergence study that measures each estimator's rate against the
theoretical target from `dea_rate()`, and a horserace benchmarking this package
against the other DEA packages on CRAN for speed, agreement and edge-case
behaviour.

## References

Banker, R. D., Charnes, A. and Cooper, W. W. (1984). Some models for estimating
technical and scale inefficiencies in data envelopment analysis.
*Management Science* 30, 1078–1092.

Chambers, R. G., Chung, Y. and Färe, R. (1996). Benefit and distance functions.
*Journal of Economic Theory* 70, 407–419.

Charnes, A., Cooper, W. W. and Rhodes, E. (1978). Measuring the efficiency of
decision making units. *European Journal of Operational Research* 2, 429–444.

Kneip, A., Park, B. U. and Simar, L. (1998). A note on the convergence of
nonparametric DEA estimators for production efficiency scores.
*Econometric Theory* 14, 783–793.

Simar, L. and Wilson, P. W. (1998). Sensitivity analysis of efficiency scores:
how to bootstrap in nonparametric frontier models. *Management Science* 44,
49–61.

Tone, K. (2001). A slacks-based measure of efficiency in data envelopment
analysis. *European Journal of Operational Research* 130, 498–509.

## License

GPL (>= 2).
