\name{dea_sbm}
\alias{dea_sbm}
\title{Slacks-based measure of efficiency}
\description{
Tone's (2001) slacks-based measure: a non-radial efficiency score that is 1 if
and only if the DMU is Pareto-Koopmans efficient, units invariant, and monotone
in every input and output slack.
}
\usage{
dea_sbm(x, y, data = NULL,
        rts = c("vrs", "crs", "nirs", "ndrs"),
        orientation = c("none", "in", "out"),
        peers = TRUE, scaling = TRUE,
        xref = NULL, yref = NULL, dataref = NULL)
}
\arguments{
\item{x, y, data}{Inputs and outputs, in any of the forms accepted by
\code{\link{dea}}. Every value must be strictly positive: the measure divides
each slack by the DMU's own level.}
\item{rts}{Returns to scale, as in \code{\link{dea}}. There is no \code{"fdh"}
option -- the measure is defined on a convex technology.}
\item{orientation}{\code{"none"} (the default) scores input and output slack
together, as a ratio; \code{"in"} scores input slack only and \code{"out"}
output slack only. All three return a score in \eqn{(0, 1]}.}
\item{peers, scaling, xref, yref, dataref}{As in \code{\link{dea}}.}
}
\value{
An object of class \code{"dea"}, with the same components as
\code{\link{dea}}'s return value; \code{eff} holds \eqn{\rho} and
\code{model} is \code{"sbm"}. Slacks are always computed, so
\code{efficient} is \code{eff == 1}.
}
\details{
The non-oriented measure is
\deqn{\rho = \frac{1 - (1/p)\sum_i s_i^-/x_{io}}{1 + (1/q)\sum_r s_r^+/y_{ro}}}
minimized over the technology. It answers a different question from the radial
score: not \dQuote{by what common factor can every input be cut} but
\dQuote{what fraction of each input is being wasted}. The two coincide only
when no slack remains at the radial projection, and \eqn{\rho \le \theta}
always.

\strong{The fractional program is solved exactly.} A ratio of linear forms is
not a linear program, and the Charnes-Cooper change of variable used here --
introduce \eqn{t > 0}, substitute \eqn{\Lambda = t\lambda} and \eqn{S = ts},
and normalize the denominator to 1 -- turns it into one without approximation.

\strong{The oriented measures are one-sided by construction.} \eqn{\rho_I} sees
only input waste and \eqn{\rho_O} only forgone output, so a DMU can score 1 on
one of them while carrying slack on the other side. Use \code{"none"} unless
the orientation is a substantive assumption about what the DMU controls.

\strong{An external reference set can make the program infeasible.} The
balance constraints are equalities with non-negative slacks, so they require
some convex combination of the reference DMUs to weakly dominate the evaluated
point in every input \emph{and} every output. A DMU scored against its own
sample dominates itself, so this is free; a point scored against someone
else's frontier may lie outside the estimated technology, and then there is no
solution. \code{dea_sbm()} reports \code{NA} and warns, naming the cause.
\code{\link{dea}} has no such problem --- its radial score simply falls below
1 for a point outside the technology --- so use it when scoring out of sample.

\strong{Zeros.} Because every slack is divided by the DMU's own level, a zero
input or output has no defined contribution. \code{dea_sbm()} stops rather than
substituting a small number, which would make the score depend on the number
chosen. Use \code{\link{dea}}, whose radial score tolerates zeros, or
\code{\link{dea_ddf}}, whose additive measure never divides by the data.
}
\references{
Tone, K. (2001). A slacks-based measure of efficiency in data envelopment
analysis. \emph{European Journal of Operational Research} 130, 498--509.
}
\seealso{\code{\link{dea}}, \code{\link{dea_ddf}}}
\examples{
sim <- dea_sim(60, p = 2, q = 1, seed = 1)

rho <- dea_sbm(sim$x, sim$y, rts = "vrs")
rho

## rho <= theta, with equality exactly where the radial projection leaves no
## slack behind.
th <- dea(sim$x, sim$y, rts = "vrs", orientation = "in")
all(rho$eff <= th$eff + 1e-8)
identical(rho$eff == 1, unname(th$efficient))
}
