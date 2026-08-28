\name{dea_cross}
\alias{dea_cross}
\alias{print.dea_cross}
\alias{summary.dea_cross}
\alias{efficiency.dea_cross}
\alias{nobs.dea_cross}
\title{Cross-efficiency}
\description{
Peer appraisal instead of self-appraisal. Ordinary DEA lets every DMU pick the
weights that flatter it most, which is why so many of them score 1 and why the
ranking is so often empty. Cross-efficiency scores each DMU under \emph{every
other} DMU's optimal weights and averages the result.
}
\usage{
dea_cross(x, y, data = NULL,
          secondary = c("benevolent", "aggressive", "none"),
          rts = c("crs", "vrs"),
          self = TRUE, scaling = TRUE,
          xref = NULL, yref = NULL, dataref = NULL)
}
\arguments{
\item{x, y, data}{Inputs and outputs, in any of the forms accepted by
\code{\link{dea}}.}
\item{secondary}{How to choose among a rater's alternate optimal weights.
\code{"benevolent"} maximizes the average cross-efficiency it awards everyone
else, \code{"aggressive"} minimizes it, and \code{"none"} takes whatever the
solver returned. See \strong{Details}: \code{"none"} is not reproducible and is
not the default.}
\item{rts}{Returns to scale. \code{"vrs"} is accepted only with
\code{secondary = "none"}.}
\item{self}{Whether a DMU's own appraisal counts in its mean. \code{FALSE}
drops the diagonal, which is Doyle and Green's convention.}
\item{scaling, xref, yref, dataref}{As in \code{\link{dea}}. The reference set
is the panel of \emph{raters}; separating it from \code{x}/\code{y} lets a
fixed panel score newcomers.}
}
\value{
An object of class \code{"dea_cross"}, a list whose components include:
\item{eff}{The cross-efficiency of each DMU: the mean of the appraisals it
receives.}
\item{cross_matrix}{The \eqn{n_{ref} \times n} matrix \eqn{E}, one row per
rater and one column per rated DMU. Its diagonal is the ordinary DEA score.}
\item{own}{Each DMU's ordinary DEA score, or \code{NULL} when the rated DMUs
are not among the raters.}
\item{maverick}{Doyle and Green's index, \code{(own - eff) / eff}: how much a
DMU's self-appraisal exceeds what its peers award it. Large means a DMU that
looks efficient only under weights nobody else would choose.}
\item{spread}{The range of appraisals each DMU receives -- a direct measure of
how much the answer depends on whose weights are used.}
\item{v, u, u0}{The weights each rater used.}
}
\details{
\deqn{E_{do} = \frac{u_d'y_o}{v_d'x_o}}
is DMU \eqn{o}'s efficiency priced at DMU \eqn{d}'s own optimal weights. The
column mean is \eqn{o}'s cross-efficiency. It very rarely ties, which is the
practical reason to use it: on the 70 sites of \code{\link{charnes1981}} the
constant-returns model calls 19 of them efficient and cannot rank them, while
cross-efficiency ranks all 70.

\strong{The weights are not unique, and that is the whole difficulty.} An
efficient DMU generally has an entire face of optimal weight vectors, all
giving it the same score of 1 but giving \emph{other} DMUs quite different
cross-efficiencies. A cross-efficiency computed from whichever vertex the
solver happened to stop at is therefore not reproducible across solvers, let
alone across packages, and the spread is first-order rather than numerical.

Doyle and Green (1994) resolve this by holding the rater's own score at
\eqn{\theta_d} and choosing among the remaining optima with a stated secondary
goal. Reporting both the benevolent and the aggressive answer brackets the
result: a ranking that survives from one to the other is a property of the
data, and one that does not was a property of the solver. That is why
\code{"none"} is available but is not the default.

\strong{Constant returns is the default for a reason.} Under variable returns
the ratio carries the intercept \eqn{u_0}, is no longer bounded by 1 and can be
negative, so the average of a column is not an efficiency. \code{dea_cross()}
allows \code{rts = "vrs"} for description but refuses to attach a secondary
goal to it, because \dQuote{maximize everyone else's average} would then be
optimizing something that is not an efficiency.
}
\references{
Sexton, T. R., Silkman, R. H. and Hogan, A. J. (1986). Data envelopment
analysis: critique and extensions. \emph{New Directions for Program
Evaluation} 32, 73--105.

Doyle, J. and Green, R. (1994). Efficiency and cross-efficiency in DEA:
derivations, meanings and uses. \emph{Journal of the Operational Research
Society} 45, 567--578.
}
\seealso{\code{\link{dea}}, \code{\link{multipliers}}}
\examples{
X <- as.matrix(charnes1981[, paste0("x", 1:5)])
Y <- as.matrix(charnes1981[, paste0("y", 1:3)])
rownames(X) <- rownames(Y) <- charnes1981$site

ben <- dea_cross(X, Y, secondary = "benevolent")
agg <- dea_cross(X, Y, secondary = "aggressive")
summary(ben)

## The bracket. Where it is wide, the ranking is an artefact of the weights.
range(ben$eff - agg$eff)

## 19 sites tie at 1 under ordinary DEA; none tie under cross-efficiency.
sum(ben$own == 1)
length(unique(round(ben$eff, 8)))
}
