\name{dea_reg}
\alias{dea_reg}
\alias{print.dea_reg}
\alias{summary.dea_reg}
\title{Two-stage efficiency regression (Simar and Wilson 2007)}
\description{
Regress efficiency on covariates, correctly. \code{dea_reg()} implements
Algorithms #1 and #2 of Simar and Wilson (2007): a truncated regression of the
estimated distance functions on \code{z}, with a bootstrap for the standard
errors, and --- under Algorithm #2 --- a bias correction of the scores before
the regression is run at all.
}
\usage{
dea_reg(object, z, data = NULL,
        algorithm = c("sw2", "sw1"),
        L1 = 100L, L2 = 2000L, alpha = 0.05,
        seed = NULL, progress = interactive())
}
\arguments{
\item{object}{A radial fit from \code{\link{dea}}.}
\item{z}{The covariates, in any of the forms \code{\link{dea}} accepts for
\code{x}: a matrix, a data frame, column names or a one-sided formula. An
intercept is added unless a constant column is already present.}
\item{data}{Optional data frame the names or formula refer to.}
\item{algorithm}{\code{"sw2"} for the double bootstrap of Algorithm #2, the
default and the one the paper recommends; \code{"sw1"} for the single
bootstrap of Algorithm #1.}
\item{L1}{Replications in the bias-correction loop, Algorithm #2 only. Each
one is a full DEA sweep, so this is where the time goes.}
\item{L2}{Replications in the bootstrap of the coefficients.}
\item{alpha}{One minus the nominal coverage of the intervals.}
\item{seed}{Passed to \code{set.seed()} for reproducibility.}
\item{progress}{Print a counter.}
}
\value{
An object of class \code{"dea_reg"}. \code{table} has one row per coefficient
with the estimate, the bootstrap standard error, the interval and the
bootstrap's own estimate of the bias in that coefficient. \code{delta} is the
response the regression used, \code{delta_bc} and \code{bias} the corrected
scores and the correction under Algorithm #2, and \code{boot} the full matrix
of bootstrap replications.
}
\details{
\strong{Why not OLS, and why not Tobit.} Regressing DEA scores on covariates is
one of the most common things done with DEA and one of the most common things
done wrong. The scores share an estimated frontier, so they are correlated by
construction in a way that depends on the sample, and conventional standard
errors do not describe their sampling distribution. Tobit does not fix this.
Censoring is a statement about how a variable is \emph{observed}; the mass at
the frontier here is an artefact of a boundary estimator, not a censored
observation of something else. The efficient DMUs are therefore
\emph{dropped}, not censored --- once it is known the frontier was hit, that
observation says nothing further about \code{z}.

\strong{The response is Shephard's distance, so read the signs carefully.} The
regression is in \eqn{\delta \ge 1}, with \eqn{\delta = 1/\theta}
input-oriented and \eqn{\delta = \phi} output-oriented, truncated from the left
at 1. A \strong{positive} coefficient therefore means the covariate is
associated with \emph{less} efficiency. That is the opposite of what a reader
expects from a regression on \eqn{\theta}, it is what the paper reports, and
\code{print()} says so on every fit.

\strong{What Algorithm #2 adds, and why it is the default.} \eqn{\hat\delta} is
a biased estimate of \eqn{\delta}, biased toward 1, so Algorithm #1 regresses
on something systematically shifted. The bias is roughly common across DMUs, so
it lands mostly on the \emph{intercept}. Algorithm #2's inner loop draws a
\eqn{\delta^*} from the fitted model, builds a pseudo-sample whose true
distance is that \eqn{\delta^*} --- scaling a DMU's inputs by a factor scales
its input distance by the same factor, so the construction is exact ---
re-estimates the original points against the pseudo-technology, and corrects
\eqn{\hat\delta} before the regression is fitted.

As in \code{\link{dea_boot}}, the pseudo-sample is the \emph{technology} and
the points \emph{evaluated} against it are the originals. Getting that
backwards produces a number for every DMU and measures nothing.

\strong{On a design whose \eqn{\beta} is known}, 200 replications at
\eqn{n = 150}, against a nominal 95\% interval:

\tabular{lrrr}{
                        \tab Algorithm #1 \tab Algorithm #2 \tab OLS on \eqn{\hat\delta} \cr
  intercept, bias       \tab  -0.070 \tab  -0.039 \tab  +0.052 \cr
  intercept, coverage   \tab   0.650 \tab   0.785 \tab   0.440 \cr
  slope 1, bias         \tab  -0.015 \tab  +0.005 \tab  -0.160 \cr
  slope 1, coverage     \tab   0.900 \tab   0.905 \tab   0.065 \cr
  slope 2, bias         \tab  +0.009 \tab  -0.001 \tab  +0.083 \cr
  slope 2, coverage     \tab   0.885 \tab   0.920 \tab   0.000
}

Read the slope rows first. \strong{The naive route's intervals essentially
never contain the truth} --- 0.065 and 0.000 against a nominal 0.95 --- and its
slope bias is ten times Algorithm #2's. Two things are wrong with it at once
and the table cannot separate them: the likelihood is misspecified, because the
data generating process is truncated, and the standard errors describe a
sampling distribution the scores do not have.

\strong{One cell goes the other way, and it is worth saying so.} On the
\emph{intercept point estimate} OLS is no worse: its bias is +0.052 against
Algorithm #2's -0.039, in the opposite direction and of much the same size, and
paired across the 200 datasets the difference in absolute error does not
resolve. The case against the naive route rests on the slopes and on the
intervals, not on every number in the table.

Between the two algorithms: \strong{#2 is nearer the truth for every
coefficient}, and comparing them on the same 200 datasets rather than as two
samples, decisively so (paired \eqn{|t|} of 13.4, 3.8 and 4.4). On interval
\emph{coverage} it beats #1 on the intercept, +0.135 with a paired standard
error of 0.027; for the slopes both sit near 0.90 and the difference is not
resolvable at this many replications. \strong{Neither reaches nominal coverage
on the intercept}, which is where the frontier bias concentrates.

The script that measures this is \code{horserace/twostage_experiment.R} in the
package's development repository; a smaller version of it is a test.

\strong{What is not bootstrapped.} The covariates. \code{z} is fixed, exactly
as in the paper; the resampling is of the error term and, under Algorithm #2,
of the frontier.

\strong{Cost.} Algorithm #1 is \code{L2} truncated regressions and is fast.
Algorithm #2 adds \code{L1} full DEA sweeps, which is where essentially all of
its time goes; the paper's suggested \code{L1 = 100} is the default.
}
\references{
Simar, L. and Wilson, P. W. (2007). Estimation and inference in two-stage,
semi-parametric models of production processes.
\emph{Journal of Econometrics} 136, 31--64.

Simar, L. and Wilson, P. W. (2011). Two-stage DEA: caveat emptor.
\emph{Journal of Productivity Analysis} 36, 205--218.
}
\seealso{\code{\link{dea}}, \code{\link{dea_boot}}}
\examples{
## A technology whose true distances depend on z.
set.seed(1)
n <- 80
X <- matrix(runif(2 * n, 1, 10), n, 2)
z <- cbind(z1 = runif(n), z2 = rnorm(n))
delta <- 1 + rexp(n, 1 / (0.1 + 0.25 * z[, 1]))
y <- matrix(apply(X, 1, function(r) prod(r)^0.45) / delta, ncol = 1)

fit <- dea(X, y, rts = "vrs", orientation = "out", slack = FALSE)

## Algorithm #1 is quick; Algorithm #2 costs L1 DEA sweeps on top.
dea_reg(fit, z, algorithm = "sw1", L2 = 200, seed = 1, progress = FALSE)

\donttest{
r <- dea_reg(fit, z, algorithm = "sw2", L1 = 50, L2 = 500, seed = 1,
             progress = FALSE)
summary(r)
}
}
