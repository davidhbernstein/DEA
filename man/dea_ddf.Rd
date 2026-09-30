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
        direction = c("both", "in", "out", "unit", "mean", "range"),
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
  \item{\code{"range"}}{\eqn{g_o = (x_o - \min x, \max y - y_o)}, the vector
    from the DMU to the technology's ideal point. This is the \strong{range
    directional model} of Portela, Thanassoulis and Simpson (2004); see
    \sQuote{The range directional model} below. The minimum and maximum are
    taken over the \emph{reference} set, unlike \code{"mean"}, which uses the
    evaluated set.}
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

A DMU whose program does not solve is reported as \code{NA} with a warning,
and its \code{status} says which kind of failure it was: 2 for an infeasible
program, anything other than 0 or 1 for a numerical one. Infeasibility is
reachable in ordinary use --- \code{direction = "in"} holds \eqn{g_y} at zero,
so a DMU whose outputs exceed everything the reference set can produce has no
feasible point at all.
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

\strong{The range directional model.} \code{direction = "range"} is the model
of Portela, Thanassoulis and Simpson (2004). The direction is the vector from
each DMU to the \emph{ideal point} of the reference set --- the corner of its
bounding box that is best in every coordinate at once --- so \eqn{\beta} reads
as the fraction of the distance to that point the DMU could travel. Two
properties follow, and both are \code{rts = "vrs"} properties.

\eqn{\beta} lies in \eqn{[0, 1]}. Under variable returns \eqn{X'\lambda} is a
convex combination of the reference inputs and so cannot fall below their
minimum, which forces \eqn{\beta \le 1}; the same argument applies on the
output side. The bound also held under every cone tested here, but the
elementary argument does not cover them, so it is stated for what it is proved
for.

\eqn{\beta} is translation invariant, which is the reason to use this model at
all. Adding a constant to a column leaves \eqn{g} unchanged, since it is a
difference of two values in that column, and shifts both sides of the
constraint by \eqn{t \sum \lambda} --- equal only when \eqn{\sum \lambda = 1}.
So negative inputs or outputs need no translation before fitting, and the
answer does not depend on one. Under a cone this fails, and by a wide margin
rather than subtly.

The range is taken over the \emph{reference} set, because it is the reference
DMUs that span the technology and the \eqn{[0, 1]} bound is a statement about
that technology. With a separate \code{xref}/\code{yref} an evaluated DMU can
lie outside that set's box; a direction component that would point backwards is
clamped to zero, with a warning, since a coordinate already better than
anything in the technology has no room to improve. A DMU that attains the best
value in \emph{every} input and output is the ideal point itself: its direction
is zero, its \eqn{\beta} is 0, no program is solved --- it would be unbounded,
not infeasible --- and its peers are reported as \code{NA}, with a warning.

\code{\link{dea_profit}} takes the same argument, but has no answer at such a
DMU: the Nerlovian measure divides the profit gap by \eqn{w'g_x + r'g_y}, which
is zero there.
}
\references{
Chambers, R. G., Chung, Y. and Fare, R. (1996). Benefit and distance functions.
\emph{Journal of Economic Theory} 70, 407--419.

Chung, Y. H., Fare, R. and Grosskopf, S. (1997). Productivity and undesirable
outputs: a directional distance function approach.
\emph{Journal of Environmental Management} 51, 229--240.

Portela, M. C. A. S., Thanassoulis, E. and Simpson, G. (2004). Negative data in
DEA: a directional distance approach applied to bank branches.
\emph{Journal of the Operational Research Society} 55, 1111--1121.
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

## The range directional model: beta is the fraction of the distance to the
## ideal point, so it lies in [0, 1], and it survives data of either sign.
r <- dea_ddf(sim$x, sim$y, direction = "range", rts = "vrs")
range(r$beta)
shifted <- dea_ddf(sim$x - 3, sim$y - 1, direction = "range", rts = "vrs")
max(abs(shifted$beta - r$beta))    # translation invariant under vrs
}
