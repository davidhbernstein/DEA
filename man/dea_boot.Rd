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
\references{
Simar, L. and Wilson, P. W. (1998). Sensitivity analysis of efficiency scores:
how to bootstrap in nonparametric frontier models. \emph{Management Science}
44, 49--61.

Simar, L. and Wilson, P. W. (2000). A general methodology for bootstrapping in
nonparametric frontier models. \emph{Journal of Applied Statistics} 27,
779--802.
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
