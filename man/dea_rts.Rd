\name{dea_rts}
\alias{dea_rts}
\alias{print.dea_rts}
\alias{summary.dea_rts}
\title{Scale efficiency and returns to scale}
\description{
Fit the constant-, variable- and non-increasing-returns technologies to the
same data in one sweep, and use the three scores to decompose efficiency into
its technical and scale parts and to classify each DMU's local returns to
scale.
}
\usage{
dea_rts(x, y, data = NULL, orientation = c("in", "out"), tol = 1e-6,
        scaling = TRUE)

\method{print}{dea_rts}(x, \dots)
\method{summary}{dea_rts}(object, \dots)
}
\arguments{
\item{x, y, data}{Inputs and outputs, in any of the forms accepted by
\code{\link{dea}}. For the \code{print} method, an object of class
\code{"dea_rts"}.}
\item{orientation}{\code{"in"} or \code{"out"}, as in \code{\link{dea}}.}
\item{tol}{How close two efficiency scores must be to count as equal. The
classification is a comparison of three linear-programming optima, so it needs
a tolerance; \code{1e-6} is loose relative to the solver and tight relative to
any economically meaningful difference.}
\item{scaling}{As in \code{\link{dea}}.}
\item{object}{An object of class \code{"dea_rts"}.}
\item{\dots}{Ignored.}
}
\value{
An object of class \code{"dea_rts"} whose \code{table} component is a data
frame with one row per DMU and columns \code{dmu}, \code{crs}, \code{vrs},
\code{nirs}, \code{scale_eff}, \code{rts} (one of \code{"irs"}, \code{"crs"},
\code{"drs"}) and \code{sum_lambda_crs}.
}
\details{
Scale efficiency is \eqn{SE = \theta_{CRS}/\theta_{VRS}} under an input
orientation and \eqn{\phi_{VRS}/\phi_{CRS}} under an output orientation, so it
lies in \eqn{(0, 1]} either way. It is 1 when the DMU operates at the most
productive scale size, and the shortfall is the part of its distance to the
constant-returns frontier that comes from being the wrong \emph{size} rather
than from being badly run.

The classification uses the non-increasing technology, which permits scaling a
reference point down but not up. A DMU in the increasing-returns region is
smaller than the most productive scale size, so the constant-returns frontier
above it is reached by scaling its peer \emph{down} -- which the non-increasing
technology permits, and there it coincides with the constant-returns frontier.
A DMU in the decreasing-returns region would need its peer scaled \emph{up},
which it forbids, and there it coincides with the variable-returns frontier:
\deqn{\theta_{CRS} = \theta_{VRS} \Rightarrow \mbox{constant returns}}
\deqn{\theta_{NIRS} = \theta_{CRS} \ne \theta_{VRS} \Rightarrow \mbox{increasing returns}}
\deqn{\theta_{NIRS} = \theta_{VRS} \ne \theta_{CRS} \Rightarrow \mbox{decreasing returns}}
Note that the rule runs opposite to what the label suggests: it is agreement
with the \emph{constant}-returns frontier that marks the increasing-returns
region.

The often-quoted alternative -- read \eqn{\sum\lambda} off the
constant-returns program and call it increasing below 1, decreasing above --
gives the same answer when that program has a unique optimum and an arbitrary
one when it does not, which is common. \code{sum_lambda_crs} is reported so
the two can be compared, but the verdict uses the rule above.

\strong{The orientation changes the answer, and is meant to.} Returns to scale
are a property of the frontier \emph{point} a DMU is benchmarked against, not
of the DMU's own coordinates. An output-oriented projection holds inputs fixed,
so the verdict describes the DMU's own scale. An input-oriented projection
moves inputs down, and for a badly inefficient DMU it can land below the most
productive scale size even though the DMU itself sits above it --- so the same
DMU can be classified decreasing by the output sweep and increasing by the
input one. Both are correct about their own projection. Choose the orientation
that matches the decision being made, and do not read the two tables against
each other row by row.

These are point estimates from an estimator with a known bias, and the
classification inherits that. A formal test of constant returns needs a
bootstrap; see Simar and Wilson (2002).
}
\references{
Coelli, T. J., Rao, D. S. P., O'Donnell, C. J. and Battese, G. E. (1998).
\emph{An Introduction to Efficiency and Productivity Analysis}. Springer.

Fare, R., Grosskopf, S. and Lovell, C. A. K. (1994).
\emph{Production Frontiers}. Cambridge University Press.

Simar, L. and Wilson, P. W. (2002). Non-parametric tests of returns to scale.
\emph{European Journal of Operational Research} 139, 115--132.
}
\seealso{\code{\link{dea}}, \code{\link{dea_boot}}}
\examples{
sim <- dea_sim(80, p = 1, q = 1, returns = 0.7, seed = 2)
r <- dea_rts(sim$x, sim$y)
r
head(r$table)
}
