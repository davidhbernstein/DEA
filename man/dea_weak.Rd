\name{dea_weak}
\alias{dea_weak}
\alias{print.dea_weak}
\title{Weak disposability of undesirable outputs}
\description{
The directional distance function on the weak-disposability technology, where
bad outputs can be reduced only by scaling the activity back rather than thrown
away freely. Chung, Fare and Grosskopf (1997), with the variable-returns
formulation of Kuosmanen (2005).
}
\usage{
dea_weak(x, y, bad, data = NULL,
         direction = c("goods_bads", "all", "bads"),
         rts = c("vrs", "crs"),
         technology = c("kuosmanen", "fare"),
         peers = TRUE, scaling = TRUE)

\method{print}{dea_weak}(x, \dots)
}
\arguments{
\item{x, y, data}{Inputs and the \emph{good} outputs, in any of the forms
accepted by \code{\link{dea}}. For the \code{print} method, an object of class
\code{"dea_weak"}.}
\item{bad}{The undesirable outputs, one column each.}
\item{direction}{\code{"goods_bads"} raises good output and lowers bad output
together, holding inputs --- the Chung-Fare-Grosskopf default. \code{"all"} also
contracts inputs; \code{"bads"} moves only the bads. A numeric vector of length
\eqn{p+q+b} or a matrix with one row per DMU gives an explicit direction. Every
entry is a \strong{magnitude}: the signs are already in the program, so a
negative entry is refused.}
\item{rts}{\code{"vrs"} or \code{"crs"}. The one-sided technologies have no
standard weak-disposability form and are not offered.}
\item{technology}{Which variable-returns formulation. Ignored under
\code{"crs"}, where the two are the same program.}
\item{peers, scaling}{As in \code{\link{dea_ddf}}.}
\item{\dots}{Ignored.}
}
\value{
An object of class \code{c("dea_weak", "dea")}. \code{eff} is \eqn{\beta}: the
proportion by which good output could rise and bad output fall at the same time,
0 on the frontier and larger worse, as in \code{\link{dea_ddf}}. Also
\code{lambda} and, for the Kuosmanen form, \code{mu}; \code{bad}; and
\code{null_joint} with \code{n_clean}.
}
\details{
\code{\link{dea_undesirable}} offers two routes that are honest rearrangements
of an ordinary program, and neither encodes the one thing everybody believes
about pollution: that you cannot simply stop producing it. Under free disposal a
DMU may always cut its bad output at no cost, so the model will cheerfully
prescribe less pollution and the same product.

Weak disposability is an equality row:
\deqn{\sum_j \lambda_j x_{ij} \le x_{io} - \beta g_{xi}}
\deqn{\sum_j \lambda_j y_{rj} \ge y_{ro} + \beta g_{yr}}
\deqn{\sum_j \lambda_j b_{fj} = b_{fo} - \beta g_{bf}}
so the bads can be reduced only by moving along the technology rather than by
assumption.

\strong{Null-jointness is a property of the data, not a constraint.} The
technology means what it is supposed to mean only if no good output can be had
without some bad; a DMU with positive goods and zero bads says the opposite, and
one such row lets the frontier claim that clean production is already feasible.
It cannot be imposed by a constraint, so it is checked, reported by
\code{print()}, and warned about.
}
\section{Which variable-returns formulation, and how much it matters}{
Fare and Grosskopf's variable-returns version applies one abatement factor to the
whole sample, which in linear form is the single-intensity program above with
\eqn{\sum\lambda = 1}. Kuosmanen (2005) showed this is not the minimal technology
satisfying weak disposability and variable returns, and gave a two-vector form in
which each DMU carries its own abatement factor:
\deqn{\sum_j (\lambda_j + \mu_j) x_{ij} \le x_{io} - \beta g_{xi}, \quad \sum_j (\lambda_j + \mu_j) = 1}
with the good and bad rows still in \eqn{\lambda} alone, so \eqn{\sum\lambda} may
fall below 1 while the input side still averages to one. That is strictly more
freedom --- \eqn{\mu = 0} is always available --- so Kuosmanen's technology
contains Fare and Grosskopf's and \eqn{\beta_{K} \ge \beta_{FG}} is an identity.

How much larger is not an identity, and the answer depends almost entirely on the
\emph{direction}, which neither paper discusses. From
\code{horserace/weak_disposability_experiment.R}, 50 replications at \eqn{n=60}:
\tabular{lrrrr}{
direction \tab coupling \tab \eqn{\mu} used \tab mean gap \tab max gap \cr
goods_bads \tab prop \tab 6.6\% \tab 0.00091 \tab 0.113 \cr
goods_bads \tab steep \tab 0.7\% \tab 0.00017 \tab 0.107 \cr
bads \tab prop \tab 16.6\% \tab 0.00435 \tab 0.203 \cr
bads \tab flat \tab \strong{30.0\%} \tab 0.01347 \tab \strong{0.393} \cr
all \tab flat \tab 2.1\% \tab 0.00039 \tab 0.097
}
On the standard direction the second vector is used for a few per cent of DMUs
and the mean gap is a thousandth --- but the largest gap is 0.1 even there, so
the choice bites a few units hard rather than everyone a little. It matters most
when the direction asks for bads to be cut and nothing else, which is exactly
when scaling the activity back is the attractive move and \eqn{\mu} is what
permits it.

What it does not much change is \emph{who} is on the frontier: the efficient set
differs between the formulations for 0.2\% to 0.7\% of DMUs. The formulation is a
question about the size of measured inefficiency far more than about its
classification.

Kuosmanen's is the default because it is the one satisfying minimal
extrapolation. Fare and Grosskopf's is kept for reproducing published results,
the same reason \code{\link{dea_mult}} keeps the 1982 model. Under \code{"crs"}
there is no convexity row, \eqn{\mu} only adds slack that free disposal already
allows, and the two are the same program --- \code{technology} is then ignored
and says so.
}
\section{Does weak disposability change the answer}{
If it did not, the honest thing would be to say so and point at the one-liner in
\code{\link{dea_undesirable}}. Against the same direction and the same data, with
the bad treated as a freely disposable input instead:
\tabular{lrrrr}{
coupling \tab mean weak \tab mean strong \tab difference \tab rank correlation \cr
prop \tab 0.0508 \tab 0.0519 \tab -0.0011 \tab 0.955 \cr
steep \tab 0.0663 \tab 0.0697 \tab -0.0034 \tab 0.951 \cr
flat \tab 0.0866 \tab 0.0950 \tab -0.0084 \tab 0.917
}
Weak disposability gives the smaller \eqn{\beta} --- it is the smaller technology,
so less improvement is available --- by between 2\% and 10\% of the level, and
the rankings agree at 0.92 to 0.96 rather than exactly.

The gap is largest where the good-bad coupling is \emph{weakest}, which is the
opposite of the natural guess. When bads move almost independently of output, the
free-disposal model finds plenty of apparently costless abatement --- peers with
similar output and much lower emissions --- and weak disposability forbids
exactly that. When bads are tightly tied to output there is little free abatement
to find and the two models nearly agree.
}
\references{
Chung, Y. H., Fare, R. and Grosskopf, S. (1997). Productivity and undesirable
outputs: a directional distance function approach.
\emph{Journal of Environmental Management} 51, 229--240.

Fare, R., Grosskopf, S., Lovell, C. A. K. and Pasurka, C. (1989). Multilateral
productivity comparisons when some outputs are undesirable.
\emph{Review of Economics and Statistics} 71, 90--98.

Kuosmanen, T. (2005). Weak disposability in nonparametric production analysis
with undesirable outputs. \emph{American Journal of Agricultural Economics} 87,
1077--1082.
}
\seealso{\code{\link{dea_undesirable}} for the two cheaper routes,
\code{\link{dea_ddf}} for the directional model on an ordinary technology.}
\examples{
set.seed(7)
n <- 50
x <- matrix(runif(n * 2, 5, 20), n, 2,
            dimnames = list(NULL, c("capital", "energy")))
product <- matrix((x[, 1] * x[, 2])^0.4 * exp(-abs(rnorm(n, 0, 0.2))), n, 1,
                  dimnames = list(NULL, "product"))
## Null-joint by construction: emissions rise with output, so nobody is clean.
co2 <- matrix(product[, 1]^1.1 * exp(rnorm(n, 0, 0.15)), n, 1,
              dimnames = list(NULL, "CO2"))

w <- dea_weak(x, product, co2)
w

## The bad-output row is an equality, which is the whole model.
peer <- as.numeric(w$lambda \%*\% co2)
max(abs(peer - (co2[, 1] - w$eff * co2[, 1])))
}
