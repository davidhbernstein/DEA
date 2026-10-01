\name{dea_undesirable}
\alias{dea_undesirable}
\alias{print.dea_undesirable}
\title{Efficiency with undesirable outputs}
\description{
Radial efficiency when some outputs are bad --- pollution, defects,
non-performing loans, mortality --- by either the Seiford and Zhu (2002)
transformation or by treating the bad outputs as inputs. Both reduce to
\code{\link{dea}} on rearranged data; what this function adds is the refusals,
the diagnostics and a slack that reads directly as reducible bad output.
}
\usage{
dea_undesirable(x, y, bad, data = NULL,
                method = c("translate", "input"),
                rts = c("vrs", "crs", "nirs", "ndrs"),
                orientation = c("in", "out"), trans = NULL, \dots)

\method{print}{dea_undesirable}(x, \dots)
}
\arguments{
\item{x, y, data}{Inputs and the \emph{desirable} outputs, in any of the forms
accepted by \code{\link{dea}}. For the \code{print} method, an object of class
\code{"dea_undesirable"}.}
\item{bad}{The undesirable outputs, in any of the same forms. One column per
bad output, one row per DMU.}
\item{method}{\code{"translate"} for Seiford and Zhu (2002), \code{"input"} for
Hailu and Veeman (2001). There is no default that is right for every question;
see \strong{Which route} below.}
\item{rts, orientation}{As in \code{\link{dea}}. \code{method = "translate"}
\strong{requires} \code{"vrs"} and \code{"in"} and refuses anything else, for
the reason given below.}
\item{trans}{The translation constant for \code{method = "translate"}: one value
per bad output or a single value for all of them. The default is
\code{max(bad) + 1} per column, which is \code{deaR}'s. The score does not
depend on it.}
\item{\dots}{Passed to \code{\link{dea}} --- \code{slack}, \code{peers},
\code{scaling}, \code{multipliers}, \code{xref}/\code{yref}.}
}
\value{
An object of class \code{c("dea_undesirable", "dea")}, so every method for
\code{"dea"} applies. In addition to the usual components:
\item{bad}{The undesirable outputs as supplied.}
\item{bad_method, bad_trans}{The route used, and the translation constants
(\code{NULL} for \code{method = "input"}).}
\item{slack_bad}{One column per bad output: the amount by which it could be
reduced, in the caller's own units, under either route.}
}
\section{Which route}{
A bad output is not an output. Free disposal runs the wrong way for it, so
listing it in \code{y} and solving the ordinary program asks the technology to
let a unit produce \emph{more} of it and rewards a unit for doing so. Three
routes answer this; two need no new solver and are offered here.

\code{method = "translate"} replaces each bad \eqn{b} by \eqn{w - b}, so that
more of the transformed variable is better, and runs the ordinary program. The
transformed variable is a genuine output: a cleaner unit has more of it, and free
disposal of it correctly says a unit could always pollute more.

\code{method = "input"} puts the bad among the inputs. Free disposal then runs
the right way, it needs no constant and no invariance, and it works under every
technology and orientation. What it does not encode is \strong{null-jointness}:
the program is free to recommend less bad at no cost in good output, which the
physics of most production forbids. That objection is real, and it is why the
third route exists.

The third route, the weak-disposability technology of Fare and Grosskopf, puts
the bads on equality rows with an abatement factor and is the only one of the
three that encodes null-jointness. It needs a new program rather than a
rearrangement of this one, and is \code{\link{dea_weak}}.
}
\section{Why the translate route refuses most models}{
\eqn{w} is arbitrary, so the answer is only meaningful where the model does not
depend on it --- that is, where the model is invariant to translating that
output. Measured on a 40-DMU design, running this transformation at
\eqn{w = \max(b)+1, 20, 100} and \eqn{1000}, the largest change in \emph{any}
DMU's score:
\tabular{lrr}{
rts \tab input \tab output \cr
vrs \tab \strong{7.0e-11} \tab 2.5e-01 \cr
crs \tab 1.8e-02 \tab 2.6e-02 \cr
nirs \tab 4.0e-03 \tab 2.5e-01 \cr
ndrs \tab 1.8e-02 \tab 2.6e-02
}
It holds in exactly one cell, with eight orders of magnitude between that cell
and every other --- not a tolerance question. And not "a convex technology"
either: \code{nirs} is convex and moves by 4e-3. So the other combinations are
refused rather than served with a number, as \code{\link{dea_profit}} refuses a
cone. \code{deaR::undesirable_basic}'s help says \code{rts} must be \code{vrs}
and says nothing about the orientation, which is the half of the precondition
that bites hardest: studying outputs is exactly when a reader reaches for an
output orientation.

\strong{The second precondition is numerical, and the failure is silent.} A
large \eqn{w} makes the transformed column nearly constant, and after it is
scaled to mean one its values differ only in the last digits. Writing
\code{rel} for the column's relative spread \eqn{\mbox{range}(b)/\mbox{mean}(w-b)},
measured against the default on the same design:
\tabular{lrr}{
rel \tab largest score change \tab DMUs unsolved \cr
1.4e+00 (default) \tab --- \tab 0 \cr
7.7e-04 \tab 1.7e-08 \tab 0 \cr
7.7e-06 \tab 2.5e-06 \tab 1 \cr
7.7e-07 \tab 0.188 \tab 3 \cr
7.7e-09 \tab 0.188 \tab \strong{0}
}
The last row is why this is warned about rather than merely documented: every
program reports a clean optimum, nothing is \code{NA}, and the scores are wrong
by 0.188. A relative spread below 1e-3 warns. The default never triggers it.
}
\references{
Hailu, A. and Veeman, T. S. (2001). Non-parametric productivity analysis with
undesirable outputs: an application to the Canadian pulp and paper industry.
\emph{American Journal of Agricultural Economics} 83, 605--616.

Seiford, L. M. and Zhu, J. (2002). Modeling undesirable factors in efficiency
evaluation. \emph{European Journal of Operational Research} 142, 16--20.

Fare, R. and Grosskopf, S. (2004). Modeling undesirable factors in efficiency
evaluation: comment. \emph{European Journal of Operational Research} 157,
242--245.
}
\seealso{\code{\link{dea_weak}} for the weak-disposability technology,
\code{\link{dea}}, \code{\link{dea_ddf}} for a directional measure that takes an
arbitrary direction.}
\examples{
set.seed(4)
n <- 40
x <- matrix(runif(n * 2, 2, 10), n, 2, dimnames = list(NULL, c("x1", "x2")))
good <- matrix((x[, 1] * x[, 2])^0.4 * exp(-abs(rnorm(n, 0, 0.2))), n, 1,
               dimnames = list(NULL, "output"))
co2 <- matrix(good[, 1]^1.2 * exp(rnorm(n, 0, 0.15)), n, 1,
              dimnames = list(NULL, "CO2"))

r <- dea_undesirable(x, good, co2)
r
head(r$slack_bad)        # reducible CO2, in the caller's units

## The same data with the bad treated as an input, which needs no constant.
dea_undesirable(x, good, co2, method = "input")

## The score does not depend on the translation constant -- that is the whole
## justification for the route, and it is why the other models are refused.
all.equal(dea_undesirable(x, good, co2, trans = 1000)$eff, r$eff,
          tolerance = 1e-8)
}
