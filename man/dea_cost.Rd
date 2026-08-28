\name{dea_cost}
\alias{dea_cost}
\alias{dea_revenue}
\alias{dea_profit}
\alias{print.dea_price}
\alias{summary.dea_price}
\alias{efficiency.dea_price}
\alias{fitted.dea_price}
\alias{peers.dea_price}
\alias{nobs.dea_price}
\title{Cost, revenue and profit efficiency}
\description{
Efficiency when prices are known. Technical efficiency asks whether a DMU is on
the frontier; these ask whether it is at the right \emph{point} on it. A DMU can
be perfectly efficient technically and still be buying the wrong input mix for
the prices it faces, and only a price-based measure sees that.
}
\usage{
dea_cost(x, y, w, data = NULL,
         rts = c("vrs", "crs", "nirs", "ndrs"),
         peers = TRUE, scaling = TRUE,
         xref = NULL, yref = NULL, dataref = NULL)

dea_revenue(x, y, r, data = NULL,
            rts = c("vrs", "crs", "nirs", "ndrs"),
            peers = TRUE, scaling = TRUE,
            xref = NULL, yref = NULL, dataref = NULL)

dea_profit(x, y, w, r, data = NULL,
           direction = c("both", "in", "out", "unit", "mean"),
           rts = c("vrs", "nirs"),
           peers = TRUE, scaling = TRUE,
           xref = NULL, yref = NULL, dataref = NULL)
}
\arguments{
\item{x, y, data}{Inputs and outputs, in any of the forms accepted by
\code{\link{dea}}.}
\item{w}{Input prices. Either a length-\eqn{p} numeric vector, meaning one
price list faced by every DMU, or an \eqn{n \times p} matrix or data frame with
one row per DMU, or -- with \code{data} -- column names or a one-sided formula.
All prices must be strictly positive. A bare vector is read as a shared price
list, \emph{not} as one observation per DMU; when \eqn{n = p} makes the two
readings indistinguishable, the call is refused rather than guessed.}
\item{r}{Output prices, in the same forms, with \eqn{q} in place of \eqn{p}.}
\item{rts}{Returns to scale, as in \code{\link{dea}}. \code{dea_profit()}
accepts only \code{"vrs"} and \code{"nirs"}; see \strong{Details}.}
\item{direction}{For \code{dea_profit()}, the direction \eqn{g} the profit gap
is normalized by, in the forms \code{\link{dea_ddf}} accepts.}
\item{peers, scaling, xref, yref, dataref}{As in \code{\link{dea}}.}
}
\value{
An object of class \code{"dea_price"}, a list whose components include:
\item{eff}{Overall efficiency. For cost and revenue, a ratio in \eqn{(0, 1]}
where 1 is best. For profit, the Nerlovian gap, where \strong{0 is best and
larger is worse}.}
\item{technical, allocative}{The two factors of \code{eff}. Cost and revenue
efficiency \emph{multiply}: \code{eff == technical * allocative}. The Nerlovian
gap \emph{adds}: \code{eff == technical + allocative}.}
\item{optimal, observed}{Minimum cost and observed cost, or maximum revenue and
observed revenue. Not present for \code{dea_profit()}, which reports
\code{profit_max} and \code{profit_obs} instead.}
\item{optimal_q}{The cost-minimizing input mix, or the revenue-maximizing
output mix, in the caller's units. Also returned by \code{fitted()}.}
\item{lambda, status, dmu, ref, n, nref, p, q}{As in \code{\link{dea}}.}
}
\details{
\strong{The decomposition is the point.} Cost efficiency is
\deqn{CE_o = \frac{\min\{w_o'x : (x, y_o) \in T\}}{w_o'x_o}
           = \theta_o \times AE_o}
the first factor being the input-oriented radial score and the second the loss
from using the wrong mix at these prices. Both are in \eqn{(0, 1]}, and a DMU
is cost efficient only if it is both technically and allocatively efficient.
Revenue efficiency decomposes the same way against the output-oriented score.

\strong{Profit is not a ratio, and this implementation does not pretend it is.}
Observed profit is routinely zero and sometimes negative, so \eqn{\Pi_o/\Pi^*}
is undefined or meaningless exactly where the question matters. The measure
used is the Nerlovian one, the profit gap normalized by the value of a
direction,
\deqn{\frac{\Pi^* - \Pi_o}{w_o'g_x + r_o'g_y}}
which the Chambers-Chung-Fare duality splits exactly into the directional
distance \eqn{\beta} in the same direction plus an allocative remainder. Both
parts are non-negative and bigger is worse, matching \code{\link{dea_ddf}} and
reversing the convention of the other two functions.

\strong{Constant and non-decreasing returns are refused for profit}, and this
is a property of the problem rather than a gap in the implementation: over a
cone, one profitable reference DMU makes maximum profit unbounded, because
scaling that DMU up scales its profit up with it. Maximum profit is a
variable-returns concept.

\strong{An external reference set can make the cost program infeasible} under
variable returns, for the same reason it can for \code{\link{dea_sbm}}: the
program must match the evaluated DMU's output from a convex combination of the
reference DMUs, and an outsider producing more than any of them cannot be
matched. \code{NA} is returned with a warning naming the cause. Under constant
returns the cone always reaches, so this does not arise.

\strong{Prices are assumed exogenous to the DMU.} Nothing here tests that, and
where a DMU is large enough to move its own prices the measures answer a
question about a price-taker that the DMU is not.
}
\references{
Farrell, M. J. (1957). The measurement of productive efficiency.
\emph{Journal of the Royal Statistical Society A} 120, 253--281.

Chambers, R. G., Chung, Y. and Fare, R. (1998). Profit, directional distance
functions, and Nerlovian efficiency. \emph{Journal of Optimization Theory and
Applications} 98, 351--364.
}
\seealso{\code{\link{dea}}, \code{\link{dea_ddf}}, \code{\link{dea_cross}}}
\examples{
X <- as.matrix(charnes1981[, paste0("x", 1:5)])
Y <- as.matrix(charnes1981[, paste0("y", 1:3)])

## One price list faced by every school site.
ce <- dea_cost(X, Y, w = c(1.4, 0.9, 2.1, 1.2, 1.0), rts = "crs")
ce

## The decomposition is exact, and this is what it buys: sites that are on the
## frontier but are not buying the cheapest mix that would keep them there.
sum(ce$technical > 1 - 1e-9)      # technically efficient
sum(ce$eff == 1)                  # AND allocatively efficient

## Profit efficiency needs prices on both sides, and variable returns.
pe <- dea_profit(X, Y, w = c(1.4, 0.9, 2.1, 1.2, 1.0), r = c(3, 2.5, 1.8))
all.equal(unname(pe$eff), unname(pe$technical + pe$allocative))
}
