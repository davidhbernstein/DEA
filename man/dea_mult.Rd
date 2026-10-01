\name{dea_mult}
\alias{dea_mult}
\alias{print.dea_mult}
\alias{summary.dea_mult}
\title{The multiplicative, or piecewise Cobb-Douglas, model}
\description{
Envelop the data with piecewise \emph{log}-linear facets rather than linear
ones, so that each facet is a Cobb-Douglas function and a peer group is a
weighted geometric mean. The multiplicative model of Charnes, Cooper, Seiford
and Stutz (1982, 1983).
}
\usage{
dea_mult(x, y, data = NULL, rts = c("vrs", "crs", "nirs", "ndrs"),
         peers = TRUE, xref = NULL, yref = NULL, dataref = NULL)

\method{print}{dea_mult}(x, \dots)
\method{summary}{dea_mult}(object, \dots)
}
\arguments{
\item{x, y, data}{Inputs and outputs, in any of the forms accepted by
\code{\link{dea}}. Every value must be \strong{strictly positive}: the model
takes logs. For the \code{print} method, an object of class \code{"dea_mult"}.}
\item{rts}{Restriction on \eqn{\sum\lambda}, \emph{in log space}. Only
\code{"vrs"} gives a units-invariant model; the others warn. See below.}
\item{peers, xref, yref, dataref}{As in \code{\link{dea}}.}
\item{object}{An object of class \code{"dea_mult"}.}
\item{\dots}{Ignored.}
}
\value{
An object of class \code{c("dea_mult", "dea")}, so the methods for \code{"dea"}
apply. Its own components:
\item{eff}{\eqn{\exp(-\mbox{total log slack})}, in \eqn{(0, 1]} with 1 exactly
on the frontier. It equals \eqn{\prod_i (x_i^*/x_{io}) \prod_r (y_{ro}/y_r^*)}.}
\item{ratio_x, ratio_y}{The factor each input could be multiplied by (in
\eqn{(0,1]}) and each output could be multiplied by (in \eqn{[1,\infty)}).
These are the model's native quantities.}
\item{log_slack_x, log_slack_y}{The underlying additive slacks in log space.}
\item{slack_x, slack_y}{The same improvements as differences in the caller's
units, so that code written against \code{"dea"} still reads.}
}
\details{
Every other estimator here builds its technology from \emph{convex}
combinations, so the frontier is piecewise linear and a peer target is an
arithmetic average of its peers. That is an assumption about the production
function, and it is usually made silently. This model makes the other one:
\deqn{\mbox{target } x_i = \prod_j x_{ij}^{\lambda_j} \quad\mbox{rather than}\quad \sum_j \lambda_j x_{ij}}
so each facet of the frontier is a Cobb-Douglas function
\eqn{y = A \prod_i x_i^{b_i}} whose exponents are estimated rather than assumed.
For most production studies that is the natural shape. No other package on CRAN
implements it.

It needs no new solver: taking logs turns a product into a sum, so the program
is exactly the unweighted additive model of \code{\link{dea_add}} on
\code{log(x)} and \code{log(y)}. The slacks are then log ratios, and
\code{ratio_x} and \code{ratio_y} are their exponentials, reported so that
nobody has to remember to exponentiate.

\strong{Units invariance needs \code{rts = "vrs"}.} Multiplying an input column
by \eqn{c} adds \eqn{\log c} to that column in log space, so units invariance
here \emph{is} the additive model's translation invariance (Ali and Seiford
1990), which needs the convexity row. The mean objective on a 40-DMU design with
one input column multiplied by \eqn{c}:
\tabular{lrrr}{
rts \tab c = 1 \tab c = 10 \tab c = 1000 \cr
vrs \tab 0.370673 \tab 0.370673 \tab 0.370673 \cr
crs \tab 0.441869 \tab 0.950590 \tab 2.187686 \cr
nirs \tab 0.424713 \tab 0.949094 \tab 2.187686 \cr
ndrs \tab 0.387829 \tab 0.372169 \tab 0.370673
}
Under \code{"vrs"} nothing moves; under \code{"crs"} the answer changes by a
factor of five although nothing was done but record an input in different units.
Geometrically, without the convexity row the log-space technology is a cone
through the log-space origin --- the point \eqn{x = y = 1}, which has no meaning
and moves whenever a unit changes.

The other technologies are still offered, with a warning, because
\code{"crs"} is the published 1982 model and the 1983 paper's title
("Invariant multiplicative efficiency\ldots") says what it was written to fix;
reproducing a published result is a real use. That is the opposite of the choice
made in \code{\link{dea_undesirable}}, where the invalid combinations are
refused outright --- there, no one had proposed them, so serving them would have
invented a quantity.

\strong{\code{rts} does not mean returns to scale here.} These are restrictions
on \eqn{\sum\lambda} in \emph{log} space. \code{"vrs"} does give a
variable-returns Cobb-Douglas envelope, but \code{"crs"} does not give a
constant-returns technology in the original space; it gives the cone described
above. The argument keeps its name because it is the same row of the same
program, not because the returns-to-scale reading survives the transformation.
}
\section{The bias runs the other way when the technology is not log-convex}{
This is the one thing to understand before using the model in place of a
linear-envelope one.

The whole of DEA's inferential story rests on the estimated frontier lying
\emph{inside} the true one, so that every score is biased toward 1 in a known
direction --- that is what \code{\link{dea_boot}} corrects. For this model that
holds only if the true technology is \strong{log-convex}, i.e. if
\eqn{\log f(e^z)} is concave. It is not implied by ordinary convexity.

For a linear frontier \eqn{f(x) = 0.5x_1 + 0.5x_2}, \eqn{\log f(e^z)} is
log-sum-exp, which is \emph{convex}, so the log chord runs above the surface:
between \eqn{x = (2, 9)} and \eqn{(9, 2)} the chord gives 5.50 where the truth
is 4.24 at the midpoint, 30\% high. No estimator is involved in that, it is a
property of the hull.

Measured in \code{horserace/multiplicative_experiment.R}, which asks both
envelopes for \eqn{\max\{y : (x_0,y) \in T\}} at the same fixed interior grid, so
that only the envelope differs:
100 replications at \eqn{n = 150}:
\tabular{lrrr}{
true shape \tab envelope \tab mean signed shortfall \tab \% of points above the truth \cr
cobb \tab log-convex \tab 0.0089 \tab 0.0 \cr
cobb \tab convex \tab 0.0211 \tab 0.0 \cr
linear \tab log-convex \tab \strong{-0.0804} \tab \strong{97.1} \cr
linear \tab convex \tab 0.0091 \tab 0.0 \cr
ces \tab log-convex \tab 0.0174 \tab 0.0 \cr
ces \tab convex \tab 0.0239 \tab 0.0
}
On a technology that really is Cobb-Douglas the log-convex envelope is 2.4 times
tighter and still inside; on a CES technology, which is neither shape, it is 1.4
times tighter and still inside. On a \emph{linear} technology it sits above the
truth at 97\% of the evaluation points. That is not a larger error but a
different kind of error: the estimated technology contains points the true one
does not.

And it does not improve with more data --- it cannot, because the estimate is
converging to the log-convex hull of the technology, which strictly contains it.
Between \eqn{n = 80} and \eqn{n = 150} the overshoot grew from -0.046 to -0.080
and the share of violated points from 87.5\% to 97.1\%. More data makes this
model more confidently wrong about the wrong technology.

So this is a maintained hypothesis about the production function, not a drop-in
alternative. Use it when a Cobb-Douglas or log-convex technology is what you
would have assumed anyway, and do not read a score from it as conservative.
}
\references{
Ali, A. I. and Seiford, L. M. (1990). Translation invariance in data
envelopment analysis. \emph{Operations Research Letters} 9, 403--405.

Charnes, A., Cooper, W. W., Seiford, L. and Stutz, J. (1982). A multiplicative
model for efficiency analysis. \emph{Socio-Economic Planning Sciences} 16,
223--224.

Charnes, A., Cooper, W. W., Seiford, L. and Stutz, J. (1983). Invariant
multiplicative efficiency and piecewise Cobb-Douglas envelopments.
\emph{Operations Research Letters} 2, 101--103.
}
\seealso{\code{\link{dea_add}}, which is the same program on untransformed data,
\code{\link{dea}}, \code{\link{dea_boot}}.}
\examples{
set.seed(3)
n <- 40
x <- matrix(runif(n * 2, 2, 10), n, 2,
            dimnames = list(NULL, c("labour", "capital")))
y <- matrix((x[, 1] * x[, 2])^0.4 * exp(-abs(rnorm(n, 0, 0.25))), n, 1,
            dimnames = list(NULL, "output"))

m <- dea_mult(x, y)
m
head(m$ratio_x)        # the factor each input could be multiplied by

## The target is a weighted GEOMETRIC mean of the peers.
all.equal(unname(x * m$ratio_x), unname(exp(m$lambda \%*\% log(x))),
          tolerance = 1e-8)

## And the score does not change when an input is recorded in other units.
x2 <- x; x2[, 1] <- x2[, 1] * 1000
all.equal(dea_mult(x2, y)$eff, m$eff, tolerance = 1e-9)
}
