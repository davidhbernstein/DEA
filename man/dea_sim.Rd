\name{dea_sim}
\alias{dea_sim}
\alias{print.dea_sim}
\alias{plot.dea_sim}
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
        mean_ineff = 0.3, x_range = c(1, 2), seed = NULL,
        frontier = c("cobb", "ces", "facet"), frontier_par = NULL)

dea_rate(p, q, rts = c("vrs", "crs", "nirs", "ndrs", "fdh"),
         what = c("mse", "estimator"))

\method{print}{dea_sim}(x, \dots)
\method{plot}{dea_sim}(x, \dots)
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
\item{frontier}{Shape of the true frontier: \code{"cobb"} the Cobb-Douglas
aggregate (smooth and strictly concave; the default and the original),
\code{"ces"} a CES aggregate whose curvature \code{frontier_par} sets, or
\code{"facet"} a piecewise-linear frontier with \code{frontier_par} facets.
All three are homogeneous of degree \code{returns} and all agree on the
diagonal, \eqn{f(c, \ldots, c) = c^r}, so \code{returns} means the same thing
across them. See \strong{Frontier geometry} below.}
\item{frontier_par}{The shape's one parameter, or \code{NULL} for its default.
For \code{"ces"} it is \eqn{\rho \le 1} (default \code{0.5}); for
\code{"facet"} the number of facets, at least 2 (default \code{3}). Refused
for \code{"cobb"}, which has no free parameter.}
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
efficiencies), \code{u}, \code{frontier} and \code{design}. A
\code{frontier = "facet"} design also carries \code{facet}, the index of the
facet each DMU sits under, and \code{facet_normals}; for the other shapes both
are \code{NULL}.

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

\strong{Frontier geometry.} A convergence rate is a claim about a class of
technologies, and one member of that class cannot confirm it, so the shape of
\eqn{f} is an argument. Changing it costs nothing, which is not obvious: both
closed forms above survive \emph{any} change of shape provided \eqn{f} stays
homogeneous of degree \eqn{r}, because
\deqn{\theta = \min\{t : f(tx) \ge \|y\|\}
            = \min\{t : t^{r} f(x) \ge f(x)e^{-u}\} = e^{-u/r},}
which uses homogeneity and nothing else --- not smoothness, not strict
concavity, and nothing Cobb-Douglas --- while \eqn{\phi = f(x)/\|y\|} holds for
any \eqn{f} whatever. So each shape below is built homogeneous of degree
\code{returns} and its truth needs no re-derivation:
\describe{
\item{\code{"cobb"}}{\eqn{f(x) = (\prod_j x_j)^{r/p}}. Smooth and strictly
concave.}
\item{\code{"ces"}}{\eqn{f(x) = ((1/p)\sum_j x_j^{\rho})^{r/\rho}}. Nests the
others: \eqn{\rho \to 0} \emph{is} Cobb-Douglas, \eqn{\rho = 1} is an arithmetic
mean and so a flat frontier, and \eqn{\rho \to -\infty} is Leontief
\eqn{\min_j x_j}. So \eqn{\rho} moves curvature continuously. \eqn{\rho > 1} is
refused rather than clamped: it makes \eqn{f} convex, so \eqn{T} is not convex
and a variable-returns estimator converges to the convex hull of \eqn{T} instead
of to \eqn{T}. The \code{theta} returned would still be the true distance and
would no longer be the estimand, so the design would quietly score the
estimator against the wrong number.}
\item{\code{"facet"}}{\eqn{f(x) = (\min_k a_k'x / p)^{r}} over
\code{frontier_par} facets. Piecewise linear, which is the shape a DEA hull can
represent exactly, and the shape with flat faces, so the optimal face need not
be a point and the projection need not be unique. The normals are
\emph{deterministic}: drawn from the random stream they would move with
\code{seed}, and \dQuote{the same geometry at two sample sizes} --- which any
rate study depends on --- would be false. Each is scaled to sum to \eqn{p}, so
every facet passes through the diagonal ray and the frontier is a fan of flat
faces meeting along it.}
}
At \eqn{p = 2} the scaled normals are collinear, and a minimum of \eqn{a_k'x}
over collinear \eqn{a_k} is always attained at an endpoint, so exactly two
facets can ever bind however many are asked for; more than two at \eqn{p = 2} is
refused for that reason. At \eqn{p \ge 3} interior facets do bind, and the
returned \code{facet} says which, so a design whose geometry is simpler than it
claims can be seen rather than assumed.

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
