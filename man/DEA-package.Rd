\name{DEA-package}
\alias{DEA-package}
\docType{package}
\title{\packageTitle{DEA}}
\description{\packageDescription{DEA}}
\details{
Technical efficiency -- five estimators, one interface and one
result object:
\describe{
  \item{\code{\link{dea}}}{radial Debreu-Farrell efficiency under constant,
  variable, non-increasing or non-decreasing returns, or free disposal, in
  either orientation, with second-stage slacks and super-efficiency.}
  \item{\code{\link{dea_sbm}}}{the slacks-based measure of Tone (2001).}
  \item{\code{\link{dea_ddf}}}{the directional distance function of Chambers,
  Chung and Fare (1996), which handles zero and negative data.}
  \item{\code{\link{dea_add}}}{the additive model of Charnes et al. (1985),
  unweighted or as the Range Adjusted Measure.}
  \item{\code{\link{dea_rts}}}{scale efficiency and the returns-to-scale
  classification.}
  \item{\code{\link{dea_categorical}}, \code{\link{dea_nd}}}{environments the
  DMU did not choose: an ordered category becomes a nested reference set, and a
  fixed input stops being something theta contracts.}
  \item{\code{\link{dea_mult}}}{the multiplicative, or piecewise
  Cobb-Douglas, model -- a log-linear envelope and a geometric peer target
  rather than a linear one and an arithmetic target. No other package on CRAN
  has it.}
  \item{\code{\link{dea_weak}}}{the weak-disposability technology, where a bad
  output can be reduced only by scaling the activity back -- the one route of
  the three that encodes null-jointness.}
  \item{\code{\link{dea_undesirable}}}{efficiency when some outputs are bad,
  by the Seiford-Zhu transformation or by treating them as inputs -- with the
  models the first of those is not valid for refused rather than served.}
}
Over time, on a long-format panel:
\describe{
  \item{\code{\link{dea_panel}}}{one row per DMU per period, with the period
  order taken from the data rather than guessed.}
  \item{\code{\link{dea_malmquist}}}{the Malmquist productivity index of
  Fare, Grosskopf, Norris and Zhang (1994), split into efficiency change and
  technical change, and efficiency change into pure and scale parts.}
}
When prices are known, efficiency splits into a technical part and a part that
is about buying the wrong mix at those prices:
\describe{
  \item{\code{\link{dea_cost}}, \code{\link{dea_revenue}}}{cost and revenue
  efficiency, each factoring exactly into technical times allocative.}
  \item{\code{\link{dea_profit}}}{Nerlovian profit inefficiency, which adds
  rather than multiplies and ties back to the directional distance function.}
}
When they are not, and the efficient set is too large to rank:
\describe{
  \item{\code{\link{dea_cross}}}{cross-efficiency -- peer appraisal rather
  than self-appraisal -- with the Doyle-Green secondary goals, which bracket an
  answer that is otherwise not reproducible across solvers.}
  \item{\code{\link{multipliers}}}{the optimal weights themselves, from the
  multiplier form of any radial fit.}
  \item{\code{\link{dea_wr}}, \code{\link{dea_cone}}}{weight restrictions,
  assurance regions and the cone-ratio model, for when the freedom to choose
  any prices at all is the problem.}
}
and, because a DEA score is an estimate and not a measurement,
\describe{
  \item{\code{\link{dea_boot}}}{bias correction and confidence intervals by the
  smoothed homogeneous bootstrap of Simar and Wilson (1998).}
  \item{\code{\link{dea_reg}}}{the two-stage regression of Simar and Wilson
  (2007), for when efficiency is to be explained by covariates -- the place
  applied DEA most often goes wrong.}
  \item{\code{\link{dea_rts_test}}}{the test of constant returns for the
  sample as a whole, by the bootstrap of Simar and Wilson (2002), where
  \code{dea_rts()} only classifies.}
  \item{\code{\link{dea_rate}}}{the rate at which any of this can converge,
  which depends on the number of inputs and outputs and is never
  root-\eqn{n}.}
  \item{\code{\link{dea_sim}}}{a technology whose efficiency is known in closed
  form, for checking an estimator against a truth rather than against another
  implementation.}
}
The DESCRIPTION file:
\packageDESCRIPTION{DEA}
\packageIndices{DEA}
}
\author{\packageAuthor{DEA}}
\keyword{package}
\seealso{
  \url{https://www.davidharrybernstein.com/software}

  \url{https://github.com/davidhbernstein/DEA}
}
\examples{
sim <- dea_sim(100, p = 2, q = 1, seed = 1)

fit <- dea(sim$x, sim$y, rts = "vrs", orientation = "in")
summary(fit)

dea_sbm(sim$x, sim$y, rts = "vrs")
dea_ddf(sim$x, sim$y, direction = "both", rts = "vrs")
dea_rts(sim$x, sim$y)

## Productivity over time.
later <- dea_sim(100, p = 2, q = 1, seed = 2)
dea_malmquist(rbind(sim$x, later$x * 0.95), rbind(sim$y, later$y * 1.08),
              id = rep(1:100, 2), period = rep(1:2, each = 100))

## With prices: on the frontier, but at the wrong point on it?
dea_cost(sim$x, sim$y, w = c(1.2, 0.8), rts = "crs")

## Without prices: ranking the units that all scored 1.
dea_cross(sim$x, sim$y, secondary = "benevolent")
}
