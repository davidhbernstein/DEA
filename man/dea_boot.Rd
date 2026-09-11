\name{dea_boot}
\alias{dea_boot}
\alias{print.dea_boot}
\alias{summary.dea_boot}
\title{Bias correction and confidence intervals by the Simar-Wilson bootstrap}
\description{
The smoothed homogeneous bootstrap of Simar and Wilson (1998): bias-corrected
efficiency estimates and confidence intervals for a radial DEA fit.
}
\usage{
dea_boot(object, B = 2000, alpha = 0.05,
         bw = c("silverman", "nrd0", "ucv", "sj"),
         seed = NULL, ncores = 1L, progress = interactive())

\method{print}{dea_boot}(x, \dots)
\method{summary}{dea_boot}(object, \dots)
}
\arguments{
\item{object}{A fit from \code{\link{dea}}, scored against its own sample
(no \code{xref}), not a super-efficiency or free-disposal-hull fit. For the
\code{print} method, an object of class \code{"dea_boot"}.}
\item{B}{Bootstrap replications. 2000 is the usual recommendation; each one
costs a full DEA sweep, so the run is \code{B} times the cost of the fit.}
\item{alpha}{One minus the nominal coverage of the interval.}
\item{bw}{Bandwidth for the kernel smoothing, as a rule name or a positive
number. \code{"silverman"} is the robust normal-reference rule Simar and Wilson
use, applied to the reflected sample; the others are \code{\link[stats]{bw.nrd0}},
\code{\link[stats]{bw.ucv}} and \code{\link[stats]{bw.SJ}}.}
\item{seed}{Integer seed. The caller's random-number stream is restored
afterwards.}
\item{ncores}{Replications run in parallel through \code{parallel::mclapply}
when this exceeds 1. Not available on Windows; a request there runs serially
with a message.}
\item{progress}{Print a counter while running serially.}
\item{x}{An object of class \code{"dea_boot"}.}
\item{\dots}{Ignored.}
}
\value{
An object of class \code{"dea_boot"}. Its \code{table} component is a data
frame with one row per DMU and columns \code{dmu}, \code{eff}, \code{bias},
\code{se}, \code{bias_corrected}, \code{ci_lower}, \code{ci_upper} and
\code{correct_worthwhile}. The full \code{n} by \code{B} matrix of bootstrap
scores is returned as \code{boot}, and the bandwidth actually used as
\code{bw}.
}
\details{
\strong{Why it is needed.} The DEA frontier is spanned by the observed DMUs, so
it lies inside the true one: every efficiency score is biased toward 1, and the
bias is not small. A DEA score reported on its own is a point estimate of
unknown accuracy from an estimator with a known, one-signed bias. This is the
largest gap between how DEA is used in practice and what is known about it.

\strong{Smoothing and reflection.} Resampling the estimated efficiencies
directly is inconsistent: it puts an atom of probability on each observed
score, in particular a large one on 1, so the pseudo-frontier is a fixed set of
points rather than a draw from a density. Smoothing with a Gaussian kernel
fixes that but would leak probability past the boundary at 1, exactly where the
mass is. The remedy is reflection: smooth the \eqn{2n}-point set
\eqn{\{\hat\theta_i\} \cup \{2 - \hat\theta_i\}}, which is symmetric about 1 by
construction, then fold draws back across it.

\strong{Whether to use the correction.} \code{bias_corrected} is
\eqn{2\hat\theta - \mathrm{mean}(\hat\theta^*)}. Simar and Wilson warn against
applying it unconditionally: it removes a bias but adds the variance of the
estimate of that bias, and is a net loss unless the bias is large relative to
the bootstrap standard error. Their rule of thumb --- correct only when
\eqn{|bias|/se > 1/\sqrt{3}} --- is evaluated per DMU and reported as
\code{correct_worthwhile}, rather than being applied silently.

\strong{Intervals.} The interval is the percentile interval for
\eqn{\hat\theta - \theta}, inverted, and so is centred near
\code{bias_corrected} rather than near \code{eff}.
}
\section{The algorithm, step by step}{
Each step of Algorithm #1 of Simar and Wilson (1998), and what this function
does for it. Stated explicitly because the details are where bootstrap
procedures for frontier models differ from one another, and two implementations
that describe themselves the same way can be doing different things.

\enumerate{
\item \emph{Fit.} \eqn{\hat\theta_i} for \eqn{i = 1, \dots, n} comes from the
supplied \code{\link{dea}} object; it is not recomputed here. The fit must be
scored against its own sample, because the procedure resamples the same units
that span the technology.

\item \emph{Reflect.} Form the \eqn{2n}-point set
\eqn{\{\hat\theta_i\} \cup \{2 - \hat\theta_i\}}. This is exact in either
orientation: \eqn{\theta \in (0, 1]} and \eqn{\phi \in [1, \infty)} both have
their boundary at 1.

\item \emph{Bandwidth.} \eqn{h} is computed on the reflected set of \eqn{2n}
points, not on the \eqn{n} original scores. The default \code{"silverman"} is
the robust normal-reference rule \eqn{0.9 \min(s, IQR/1.349)(2n)^{-1/5}}.
The alternatives are offered because the choice is genuinely open --- the
reflected density is bimodal enough that a normal reference oversmooths it ---
and a positive number may be given directly.

\item \emph{Resample and smooth.} Draw \eqn{\beta_i^*} with replacement from
the reflected set and set
\eqn{\tilde\theta_i^* = \bar\beta^* + (\beta_i^* + h\varepsilon_i - \bar\beta^*)/\sqrt{1 + h^2/s^2_{refl}}},
with \eqn{\varepsilon_i} standard normal. The denominator is the variance
correction: adding \eqn{h\varepsilon} inflates the variance by \eqn{h^2}, and
without rescaling the pseudo-sample is overdispersed.

\item \emph{Fold.} Map draws back across the boundary ---
\eqn{\theta^* = 2 - \tilde\theta^*} where \eqn{\tilde\theta^* > 1} for an input
orientation, and the mirror for an output orientation. The smoothed density is
symmetric about 1 by construction, so folding is measure preserving and no mass
crosses the boundary.

\item \emph{Build the pseudo-sample.} Move each DMU along its own orientation
so that its efficiency \emph{relative to the estimated frontier} becomes
\eqn{\theta^*_i} instead of \eqn{\hat\theta_i}: \eqn{x^*_i = x_i \hat\theta_i / \theta^*_i}
for an input orientation, \eqn{y^*_i = y_i \hat\theta_i / \theta^*_i} for an
output one.

\item \emph{Re-estimate.} The \strong{technology} is the pseudo-sample; the
points \strong{evaluated} against it are the original observations. Scoring the
pseudo-points against themselves would simply reproduce \eqn{\hat\theta}.

\item \emph{Read off.} Over \code{B} replications:
\eqn{bias = \mathrm{mean}(\hat\theta^*) - \hat\theta},
\eqn{se = \mathrm{sd}(\hat\theta^*)}, and
\code{bias_corrected} \eqn{= \hat\theta - bias = 2\hat\theta - \mathrm{mean}(\hat\theta^*)}.
The interval is the percentile interval on \eqn{\hat\theta^* - \hat\theta}
inverted, i.e. \eqn{2\hat\theta - q_{1-\alpha/2}} to
\eqn{2\hat\theta - q_{\alpha/2}}, which is an interval for \eqn{\theta} rather
than for \eqn{\hat\theta}.
}

Bootstrap scores routinely cross 1 --- \eqn{\theta^*} above it, \eqn{\phi^*}
below --- because the pseudo-technology is a shrunken copy of the estimated one
and an original observation need not lie inside it. That tail \emph{is} the
interval, and it is not clamped.
}

\section{What is assumed, and when it fails}{
\strong{Homogeneity.} This is the \emph{homogeneous} bootstrap: it assumes the
distribution of inefficiency does not depend on the input-output mix, so that
one common density of \eqn{\theta} may be resampled for every DMU. Where
inefficiency varies systematically with the inputs --- large units
systematically closer to the frontier than small ones, say --- the assumption
fails and the procedure is not consistent. The remedy is the heterogeneous
bootstrap of Simar and Wilson (2000), which conditions the resampled density on
\eqn{(x, y)} and is \strong{not implemented here}.

\strong{Convexity.} The resampling scheme is built for the convex DEA
estimator. Under \code{rts = "fdh"} the relevant procedure is the subsampling
bootstrap of Jeong and Simar (2006), not this one, and \code{dea_boot} refuses
a free-disposal-hull fit rather than returning something inapplicable.

\strong{Super-efficiency.} Refused. The Andersen-Petersen score is not an
estimator of a distance to the true frontier, so there is no bias to correct.

\strong{Boundary DMUs.} The sampling distribution of a DMU estimated as
efficient is one-sided, and intervals are hardest exactly there. Coverage
should be judged separately for frontier and interior units rather than
averaged; the harness in \code{horserace/} of the package's development
repository reports the split for that reason.

\strong{Choosing \code{B}.} 2000 is the usual recommendation. Interval
endpoints are sample quantiles, so a small \code{B} interpolates the tails from
few order statistics. Cost is linear in \code{B} and the whole run is \code{B}
times the cost of one fit.
}

\references{
Simar, L. and Wilson, P. W. (1998). Sensitivity analysis of efficiency scores:
how to bootstrap in nonparametric frontier models. \emph{Management Science}
44, 49--61.

Simar, L. and Wilson, P. W. (2000). A general methodology for bootstrapping in
nonparametric frontier models. \emph{Journal of Applied Statistics} 27,
779--802.

Jeong, S.-O. and Simar, L. (2006). Linearly interpolated FDH efficiency score
for nonconvex frontiers. \emph{Journal of Multivariate Analysis} 97,
2141--2161.
}
\seealso{\code{\link{dea}}, \code{\link{dea_rate}}}
\examples{
sim <- dea_sim(40, p = 1, q = 1, seed = 3)
fit <- dea(sim$x, sim$y, rts = "vrs", orientation = "in")

## B is tiny here so the example runs quickly; use 2000 in practice.
b <- dea_boot(fit, B = 50, seed = 1, progress = FALSE)
b
head(b$table)

## The correction moves the estimates toward the truth.
c(raw = mean(fit$eff - sim$theta),
  corrected = mean(b$table$bias_corrected - sim$theta))
}
