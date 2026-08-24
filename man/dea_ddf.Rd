\name{dea_ddf}
\alias{dea_ddf}
\title{Directional distance function}
\description{
The directional technology distance function of Chambers, Chung and Fare
(1996): the number of units of a chosen direction vector by which a DMU could
move toward the frontier, contracting inputs and expanding outputs at the same
time.
}
\usage{
dea_ddf(x, y, data = NULL,
        direction = c("both", "in", "out", "unit", "mean"),
        rts = c("vrs", "crs", "nirs", "ndrs"),
        peers = TRUE, scaling = TRUE,
        xref = NULL, yref = NULL, dataref = NULL)
}
\arguments{
\item{x, y, data}{Inputs and outputs, in any of the forms accepted by
\code{\link{dea}}. Unlike the other models here, values may be zero or
negative.}
\item{direction}{The direction \eqn{g = (g_x, g_y)} along which distance is
measured. One of
\describe{
  \item{\code{"both"}}{\eqn{g = (x_o, y_o)}, each DMU's own data. The default.}
  \item{\code{"in"}}{\eqn{g = (x_o, 0)}; then \eqn{\beta = 1 - \theta}.}
  \item{\code{"out"}}{\eqn{g = (0, y_o)}; then \eqn{\beta = \phi - 1}.}
  \item{\code{"unit"}}{\eqn{g = (1, \ldots, 1)}, one unit of every variable.}
  \item{\code{"mean"}}{the sample mean of each variable, the same for all DMUs.}
}
Or a numeric vector of length \code{p + q} giving one fixed direction, inputs
first; or a \code{n} by \code{p + q} matrix giving one per DMU.}
\item{rts}{Returns to scale, as in \code{\link{dea}}.}
\item{peers, scaling, xref, yref, dataref}{As in \code{\link{dea}}.}
}
\value{
An object of class \code{"dea"} with \code{model = "ddf"}. Note that
\code{eff} (also available as \code{beta}) measures \emph{inefficiency}:
\strong{0 is on the frontier and larger is worse}, the reverse of every other
model in this package. The direction actually used is returned in \code{g}.
There are no second-stage slacks -- the projection is along \eqn{g}.
}
\details{
\strong{Why not just use a radial measure.} Three reasons, in order of how
often they bite.

Inputs and outputs move together. A radial measure has to pick a side;
\code{"both"} contracts inputs and expands outputs at once, which is what most
efficiency questions actually ask.

Zeros and negatives are fine. \eqn{\beta} is an additive contraction, so
nothing is divided by a DMU's own level. A zero output leaves the radial output
score undefined and the slacks-based measure undefined; here it is just a
number. Net income, net exports and abatement credits can be negative, and the
directional distance function is the standard way to handle them.

The direction is visible. A radial model also chooses a direction; it simply
does not say which.

\strong{Units.} \eqn{\beta} is scale invariant when \eqn{g} is proportional to
the DMU's own data and is \emph{not} when \eqn{g} is fixed, since doubling an
input's units halves the \eqn{\beta} needed to move one unit of \eqn{g}. That
is a property of the estimand rather than a defect, but a fixed direction must
be stated in the same units as the data. Internal rescaling is applied to the
direction as well, so the reported \eqn{\beta} is the one the caller's units
imply.

\strong{No common ratio scale.} \code{\link{efficiency}(fit, "score")} refuses
to convert \eqn{\beta} to a \eqn{(0, 1]} score for a direction that is not
purely input or purely output, because no such conversion exists.
}
\references{
Chambers, R. G., Chung, Y. and Fare, R. (1996). Benefit and distance functions.
\emph{Journal of Economic Theory} 70, 407--419.

Chung, Y. H., Fare, R. and Grosskopf, S. (1997). Productivity and undesirable
outputs: a directional distance function approach.
\emph{Journal of Environmental Management} 51, 229--240.
}
\seealso{\code{\link{dea}}, \code{\link{dea_sbm}}}
\examples{
sim <- dea_sim(60, p = 2, q = 1, seed = 1)

dea_ddf(sim$x, sim$y, direction = "both", rts = "vrs")

## The radial models are the two one-sided directions.
b <- dea_ddf(sim$x, sim$y, direction = "in", rts = "vrs")
t <- dea(sim$x, sim$y, orientation = "in", rts = "vrs", slack = FALSE)
max(abs(b$beta - (1 - t$eff)))

## A fixed direction: one unit of every variable.
dea_ddf(sim$x, sim$y, direction = "unit", rts = "vrs")
}
