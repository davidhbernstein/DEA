\name{dea-methods}
\alias{efficiency}
\alias{efficiency.dea}
\alias{peers}
\alias{peers.dea}
\alias{slacks}
\alias{slacks.dea}
\alias{print.dea}
\alias{summary.dea}
\alias{plot.dea}
\alias{fitted.dea}
\alias{nobs.dea}
\title{Methods for objects of class "dea"}
\description{
Extractors and display methods shared by \code{\link{dea}},
\code{\link{dea_sbm}} and \code{\link{dea_ddf}}, which all return the same
class.
}
\usage{
efficiency(object, \dots)
\method{efficiency}{dea}(object, type = c("natural", "score"), \dots)

peers(object, \dots)
\method{peers}{dea}(object, threshold = 1e-08, \dots)

slacks(object, \dots)
\method{slacks}{dea}(object, \dots)

\method{print}{dea}(x, \dots)
\method{summary}{dea}(object, \dots)
\method{plot}{dea}(x, ngrid = 400, \dots)
\method{fitted}{dea}(object, \dots)
\method{nobs}{dea}(object, \dots)
}
\arguments{
\item{object, x}{An object of class \code{"dea"}.}
\item{type}{\code{"natural"} returns the model's own measure --- \eqn{\theta},
\eqn{\phi}, \eqn{\rho} or \eqn{\beta}. \code{"score"} maps all of them onto
\eqn{(0, 1]} with 1 meaning efficient, which is what a table comparing models
needs.}
\item{threshold}{Intensity weights at or below this are not reported as peers.}
\item{ngrid}{Number of grid points used to trace the frontier in
\code{plot()}.}
\item{\dots}{Passed to the underlying plot, or ignored.}
}
\value{
\code{efficiency} and \code{nobs} return numeric vectors; \code{peers} returns
a data frame with columns \code{dmu}, \code{peer} and \code{lambda};
\code{slacks} and \code{fitted} return matrices with one row per DMU.
\code{print} and \code{summary} return their argument invisibly.
}
\details{
\code{efficiency(object, "score")} deliberately \emph{fails} for a directional
fit whose direction is neither purely input nor purely output. \eqn{\beta} is
an additive distance measured in units of \eqn{g} and there is no ratio it
corresponds to; a formula such as \eqn{(1-\beta)/(1+\beta)} would produce a
number that looks comparable to a Farrell score and is not.

\code{fitted()} returns the frontier projection: where each DMU would sit if it
were efficient. For a radial fit with slacks that is the radial projection with
the remaining slack removed, so the result is a Pareto-Koopmans point.

\code{plot()} draws the estimated frontier when there is one input and one
output, and the distribution of scores otherwise. The frontier is traced by
\emph{evaluating} the fitted technology on a grid of input levels rather than
by joining up the efficient DMUs, so the same code draws it correctly under
every returns-to-scale assumption --- including the free disposal hull's
staircase --- and puts the vertical and horizontal faces where they actually
are.

\strong{Peer sets are not unique.} These linear programs frequently have
multiple optima, so \code{peers()} reports one valid benchmark set rather than
the benchmark set. Efficiency scores are unique; the weights achieving them
generally are not.
}
\seealso{\code{\link{dea}}, \code{\link{dea_sbm}}, \code{\link{dea_ddf}}}
\examples{
sim <- dea_sim(60, p = 1, q = 1, seed = 1)
fit <- dea(sim$x, sim$y, rts = "vrs")

summary(fit)
head(efficiency(fit, "score"))
head(peers(fit))
head(slacks(fit))
head(fitted(fit))

plot(fit)
}
