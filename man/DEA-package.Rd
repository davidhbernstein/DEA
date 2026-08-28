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
}
and, because a DEA score is an estimate and not a measurement,
\describe{
  \item{\code{\link{dea_boot}}}{bias correction and confidence intervals by the
  smoothed homogeneous bootstrap of Simar and Wilson (1998).}
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

## With prices: on the frontier, but at the wrong point on it?
dea_cost(sim$x, sim$y, w = c(1.2, 0.8), rts = "crs")

## Without prices: ranking the units that all scored 1.
dea_cross(sim$x, sim$y, secondary = "benevolent")
}
