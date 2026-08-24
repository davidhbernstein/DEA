\name{dea_sim}
\alias{dea_sim}
\alias{print.dea_sim}
\alias{dea_rate}
\title{Simulate a technology with known efficiency, and the rate it can be
estimated at}
\description{
\code{dea_sim} draws a sample from a production technology whose distance
functions have a closed form, so that an estimator can be scored against a
truth written down before it ran. \code{dea_rate} gives the rate at which the
DEA estimator can converge on such a design, expressed as the slope of log mean
squared error on log sample size.
}
\usage{
dea_sim(n, p = 1L, q = 1L, returns = 1, ineff = c("exp", "hnorm"),
        mean_ineff = 0.3, x_range = c(1, 2), seed = NULL)

dea_rate(p, q, rts = c("vrs", "crs", "nirs", "ndrs", "fdh"),
         what = c("mse", "estimator"))

\method{print}{dea_sim}(x, \dots)
}
\arguments{
\item{n}{Number of DMUs.}
\item{p, q}{Numbers of inputs and outputs.}
\item{returns}{Elasticity of scale of the true frontier, in \eqn{(0, 1]}. At 1
the technology is a cone, so a constant-returns fit is consistent; below 1 it
is strictly concave and a constant-returns fit is inconsistent by
construction, which makes it a useful negative control.}
\item{ineff}{Distribution of the inefficiency term \eqn{u \ge 0}: exponential
or half-normal. Both are parameterized so that \eqn{E[u]} is
\code{mean_ineff}, so the two designs are comparable.}
\item{mean_ineff}{\eqn{E[u]}.}
\item{x_range}{Support of each input, as two increasing positive numbers.
Widening it downward matters for any estimand whose projection moves inputs
down --- see \sQuote{Details}.}
\item{seed}{Integer seed. The caller's random-number stream is restored
afterwards.}
\item{rts}{Returns-to-scale assumption of the estimator whose rate is wanted.}
\item{what}{\code{"mse"} returns the slope for mean squared error;
\code{"estimator"} returns it for the estimator itself, which is half as
steep.}
\item{x}{An object of class \code{"dea_sim"}.}
\item{\dots}{Ignored.}
}
\value{
\code{dea_sim} returns an object of class \code{"dea_sim"}: a list with
\code{x} and \code{y} (matrices), \code{data} (both as one data frame),
\code{theta} and \code{phi} (the true input- and output-oriented Farrell
efficiencies), \code{u}, \code{frontier} and \code{design}.

\code{dea_rate} returns a single negative number.
}
\details{
\strong{The design.} Inputs are uniform on \code{x_range}\eqn{^p}, by default
\eqn{[1, 2]^p}. The frontier is the
Cobb-Douglas aggregate \eqn{f(x) = \prod_j x_j^{r/p}} with \eqn{r =}
\code{returns}, and the output set at \eqn{x} is the non-negative part of the
ball of radius \eqn{f(x)}:
\deqn{T = \{(x, y) : \|y\| \le f(x)\}.}
This is a genuine production technology in the sense DEA assumes -- \eqn{f} is
concave for \eqn{r \le 1} so \eqn{T} is convex, both \eqn{f} and the norm are
monotone so \eqn{T} is freely disposable, and at \eqn{r = 1} it is a cone.
Setting \eqn{\|y\| = f(x)e^{-u}} then gives
\deqn{\phi = f(x)/\|y\| = e^{u}, \qquad \theta = e^{-u/r},}
the second because scaling every input by \eqn{\theta} scales \eqn{f} by
\eqn{\theta^{r}}. With \eqn{q > 1} the direction of \eqn{y} within the ball is
drawn uniformly and independently of \eqn{u}, so the output mix carries no
information about efficiency; fixing it instead would put every DMU on one ray,
and a \eqn{q}-output problem on a single ray is a one-output problem wearing
\eqn{q} columns.

\strong{\code{theta} needs care, and \code{x_range} is why it is an argument.}
\code{theta} is the distance to the technology as defined for \emph{all}
\eqn{x > 0}, and reaching it means scaling inputs down --- for an inefficient
DMU, below \eqn{\min(x)}. A variable-returns hull cannot extrapolate there, so
on the default support it does not estimate \code{theta} consistently however
large the sample: the mean squared error hits a floor and a log-log sweep
produces a slope near zero. That is a property of the estimand, not a bug.
Three ways out: use \code{phi}, since outputs expand at fixed \eqn{x} and stay
inside the support; fit \code{rts = "crs"} against \code{theta} on a
\code{returns = 1} design, where the technology is a cone and can be scaled
without bound; or widen \code{x_range} downward so the sample has room beneath
the points being scored. The same applies to the directional distance
function's two-sided direction, which also moves inputs down.

\strong{The rate.} DEA is consistent but not root-\eqn{n} consistent, and the
rate worsens with every variable added. For the radial estimator with \eqn{p}
inputs and \eqn{q} outputs,
\deqn{|\hat\theta - \theta| = O_p(n^{-2/(p+q+1)})}
under variable returns (Kneip, Park and Simar 1998),
\eqn{O_p(n^{-2/(p+q)})} under constant returns (Park, Simar and Weiner 2000),
and \eqn{O_p(n^{-1/(p+q)})} for the free disposal hull. Squaring doubles the
exponent, so regressing log mean squared error on log \eqn{n} over a grid of
sample sizes has slope \eqn{-4/(p+q+1)}, \eqn{-4/(p+q)} and \eqn{-2/(p+q)}
respectively. That number, not \eqn{-1}, is the benchmark: the familiar
root-\eqn{n} slope is reached only at \eqn{p + q = 3} under variable returns,
and at \eqn{p + q = 8} the slope is \eqn{-4/9}, so a hundredfold increase in
sample size buys about a factor of 8 in mean squared error where a parametric
estimator would buy 100.

These are worst-case rates over a class of technologies. On a frontier that
happens to be linear the estimator does better, so a slope \emph{steeper} than
the target is not a defect; only a slope materially flatter is.
}
\references{
Kneip, A., Park, B. U. and Simar, L. (1998). A note on the convergence of
nonparametric DEA estimators for production efficiency scores.
\emph{Econometric Theory} 14, 783--793.

Park, B. U., Simar, L. and Weiner, C. (2000). The FDH estimator for
productivity efficiency scores. \emph{Econometric Theory} 16, 855--877.

Simar, L. and Wilson, P. W. (2015). Statistical approaches for non-parametric
frontier models: a guided tour. \emph{International Statistical Review} 83,
77--110.
}
\seealso{\code{\link{dea}}, \code{\link{dea_boot}}}
\examples{
sim <- dea_sim(200, p = 1, q = 1, returns = 1, seed = 1)
sim

## A constant-returns fit recovers the truth on a constant-returns design.
fit <- dea(sim$x, sim$y, rts = "crs", orientation = "in", slack = FALSE)
c(bias = mean(fit$eff - sim$theta), rmse = sqrt(mean((fit$eff - sim$theta)^2)))

## What the rate ought to be, here and in a bigger problem.
dea_rate(1, 1, "crs")
dea_rate(4, 4, "vrs")
}
