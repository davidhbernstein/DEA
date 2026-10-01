\name{dea_subsample}
\alias{dea_subsample}
\alias{print.dea_subsample}
\alias{summary.dea_subsample}
\title{Confidence intervals by the m out of n bootstrap}
\description{
Subsampling intervals for radial efficiency scores: draw \eqn{m < n}
observations without replacement, rescale by the estimator's convergence rate,
and invert. Simar and Wilson (2011), with the data-driven choice of \eqn{m}
following Politis, Romano and Wolf.
}
\usage{
dea_subsample(object, m = NULL, B = 2000, alpha = 0.05,
              m_grid = NULL, B_select = 200L, window = 3L,
              seed = NULL, ncores = 1L, progress = interactive())

\method{print}{dea_subsample}(x, \dots)
\method{summary}{dea_subsample}(object, \dots)
}
\arguments{
\item{object}{A self-referenced radial fit from \code{\link{dea}}. Unlike
\code{\link{dea_boot}}, \code{rts = "fdh"} is accepted.}
\item{m}{Subsample size. \code{NULL} chooses it by minimum volatility.}
\item{B}{Replications at the chosen \code{m}.}
\item{alpha}{Interval level.}
\item{m_grid}{Candidate values of \code{m}. The default is twelve points
spaced geometrically between about \eqn{n^{0.4}} and \eqn{n^{0.9}}.}
\item{B_select}{Replications used at each candidate while choosing \code{m}.
Smaller than \code{B} on purpose: the selection only needs the endpoints to be
comparable across \code{m}.}
\item{window}{Half-width of the neighbourhood the volatility is computed over.}
\item{seed, ncores, progress}{As in \code{\link{dea_boot}}.}
\item{x}{An object of class \code{"dea_subsample"}.}
\item{\dots}{Ignored.}
}
\value{
An object of class \code{"dea_subsample"} whose \code{table} has one row per DMU
with \code{eff}, \code{bias_corrected}, \code{ci_lower}, \code{ci_upper} and
\code{n_valid}. Also \code{m}, \code{m_grid}, \code{endpoint_path},
\code{volatility}, \code{eligible} and \code{kappa}.
}
\details{
\code{\link{dea_boot}} resamples all \eqn{n} observations. For a boundary
estimator that is not known to be consistent: the score's limiting distribution
is driven by how many observations land near the frontier, and an
\eqn{n}-out-of-\eqn{n} resample reproduces the sample's own boundary rather than
the population's. Simar and Wilson (2011) propose subsampling instead.

Writing \eqn{\kappa} for the convergence exponent that \code{\link{dea_rate}}
already knows --- \eqn{2/(p+q+1)} under variable returns, \eqn{2/(p+q)} under
constant, \eqn{1/(p+q)} under free disposal ---
\eqn{n^{\kappa}(\hat\theta_n - \theta)} has a limiting distribution, and
\deqn{T^*_m = m^{\kappa}(\hat\theta^*_m - \hat\theta_n)}
estimates it, where \eqn{\hat\theta^*_m} scores the same observation against a
subsample of size \eqn{m}. Inverting gives the reported interval. No orientation
branch is needed: input-oriented \eqn{\theta} is biased up and \eqn{\theta^*_m}
is larger still, output-oriented \eqn{\phi} is biased down and \eqn{\phi^*_m}
smaller still, and the signs carry themselves.

It also covers the case \code{\link{dea_boot}} refuses. Free disposal has no
smoothed homogeneous bootstrap --- the hull is not convex --- and subsampling is
the scheme Jeong and Simar (2006) propose for it.

\strong{Choosing \eqn{m} is the whole difficulty.} The theory asks for
\eqn{m \to \infty} with \eqn{m/n \to 0} and says nothing about which \eqn{m} to
use at the \eqn{n} in front of you. The minimum-volatility rule computes the
interval across a grid and takes the \eqn{m} whose neighbourhood moves the
endpoints least, on the argument that a flat stretch is where the asymptotics
have bitten. Two things about the implementation:

\itemize{
\item Only candidates with a \emph{full} window are eligible. A clipped window
averages over fewer neighbours, and a standard deviation over fewer points is
smaller for no reason connected to stability --- so clipping biases the rule
toward the extreme \eqn{m}. Observed before this was fixed: on a 100-DMU fit the
rule picked the first grid point at volatility 0.0107 against 0.0109 for a
candidate in the middle of a genuinely flat stretch, purely because one was
scored on four neighbours and the other on seven.
\item Politis, Romano and Wolf choose \eqn{m} for \emph{one} parameter; here
there are \eqn{n} intervals at once. The volatility is therefore computed on the
\emph{mean} endpoints across DMUs and one \eqn{m} serves all of them. That
aggregation is a choice the method does not make, so \code{m_grid},
\code{endpoint_path} and \code{volatility} are returned for inspection and
\code{summary()} prints the path.
}
}
\section{What subsampling costs, and who pays}{
Scoring an original observation against a subsample can have no solution at all:
under variable returns the point may lie outside the smaller hull, and then there
is no answer rather than a large one. That is not an edge case, and the burden
falls on exactly the DMUs an analyst cares about. Share of replications with no
solution, on \code{dea_sim(returns = 0.9)}, output-oriented, variable returns:
\tabular{lrrr}{
 \tab \eqn{m=n^{0.4}} \tab \eqn{m=n^{0.6}} \tab \eqn{m=n^{0.8}} \cr
\eqn{n=60}, all DMUs \tab 31.3\% \tab 15.7\% \tab 6.3\% \cr
\eqn{n=60}, called efficient \tab \strong{58.1\%} \tab 39.2\% \tab 22.1\% \cr
\eqn{n=60}, interior \tab 26.5\% \tab 11.5\% \tab 3.5\% \cr
\eqn{n=150}, all DMUs \tab 27.3\% \tab 13.3\% \tab 5.0\% \cr
\eqn{n=150}, called efficient \tab \strong{59.6\%} \tab 42.7\% \tab 24.9\%
}
The worst single DMU lost between 63\% and 97\% of its draws. So the interval for
a frontier DMU is built from the fewest observations, and a smaller \eqn{m} makes
that worse at exactly the rate that makes the asymptotics better.

\code{n_valid} is therefore reported per DMU, an interval needs at least
\eqn{\lceil 2/\alpha \rceil} valid draws to be reported at all --- below that the
\eqn{\alpha/2} quantile is extrapolation into an empty tail rather than an order
statistic --- and anything short comes back \code{NA} with a warning naming the
count.
}
\section{Does it cover better than dea_boot}{
Measured against \code{dea_sim}'s closed-form truth, 100 replications,
output-oriented variable returns, nominal 95\%, both methods on the same fit in
every replication:
\tabular{lrrrr}{
 \tab coverage \tab frontier \tab interior \tab width \cr
\eqn{n=50}, \code{dea_boot} \tab 0.814 \tab 0.794 \tab 0.819 \tab 0.170 \cr
\eqn{n=50}, \code{dea_subsample} \tab 0.811 \tab \strong{0.602} \tab 0.871 \tab 0.157 \cr
\eqn{n=200}, \code{dea_boot} \tab 0.879 \tab 0.710 \tab 0.897 \tab 0.089 \cr
\eqn{n=200}, \code{dea_subsample} \tab \strong{0.925} \tab 0.695 \tab \strong{0.950} \tab 0.097
}
Three things to take from it, and only one of them is a recommendation.

\strong{At \eqn{n=200} it is a real improvement, and it is the interior that
improves.} Overall coverage rises 0.046 with a paired standard error of 0.0047
(\eqn{t = 9.8}), and interior coverage reaches 0.950 against a nominal 0.95 ---
which is the package's largest documented inference weakness closed, not
narrowed.

\strong{At \eqn{n=50} it is not an improvement at all.} The paired difference is
-0.003 and does not resolve. Subsampling needs sample size, which is what
\eqn{m \to \infty} with \eqn{m/n \to 0} asks for and what \eqn{n=50} cannot give.

\strong{It never helps the frontier DMUs} --- 0.794 to 0.602 at \eqn{n=50} and
0.710 to 0.695 at \eqn{n=200} --- because that is where the infeasibility above
concentrates. If the units of interest are the ones the fit calls efficient,
this method has nothing to offer them and the smaller \eqn{m} the selection
prefers makes it worse.
}
\references{
Jeong, S.-O. and Simar, L. (2006). Linearly interpolated FDH estimators of
multivariate efficiency scores. \emph{Journal of Multivariate Analysis} 97,
2141--2161.

Politis, D. N., Romano, J. P. and Wolf, M. (1999). \emph{Subsampling}. Springer.

Simar, L. and Wilson, P. W. (2011). Inference by the \eqn{m} out of \eqn{n}
bootstrap in nonparametric frontier models. \emph{Journal of Productivity
Analysis} 36, 33--53.
}
\seealso{\code{\link{dea_boot}} for the full-sample bootstrap this is an
alternative to, \code{\link{dea_rate}} for the exponent, \code{\link{dea}}.}
\examples{
sim <- dea_sim(80, p = 2, q = 1, returns = 0.9, seed = 5)
fit <- dea(sim$x, sim$y, rts = "vrs", orientation = "out", slack = FALSE)

r <- dea_subsample(fit, m = 20, B = 200, seed = 1, progress = FALSE)
r

## The m selection, as a path to be inspected rather than trusted.
\donttest{
sel <- dea_subsample(fit, B = 200, B_select = 50, seed = 1, progress = FALSE)
summary(sel)
}
}
