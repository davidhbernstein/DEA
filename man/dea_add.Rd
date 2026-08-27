\name{dea_add}
\alias{dea_add}
\title{Additive efficiency, and the weighted measures built on it}
\description{
The additive model of Charnes, Cooper, Golany, Seiford and Stutz (1985):
maximize total slack over the technology, with no orientation and no radial
contraction. Available unweighted, or normalized as the Range Adjusted Measure
or the Measure of Inefficiency Proportions.
}
\usage{
dea_add(x, y, data = NULL,
        measure = c("ram", "mip", "unweighted"),
        rts = c("vrs", "crs", "nirs", "ndrs"),
        peers = TRUE, scaling = TRUE,
        xref = NULL, yref = NULL, dataref = NULL)
}
\arguments{
\item{x, y, data}{Inputs and outputs, in any of the forms accepted by
\code{\link{dea}}. \code{measure = "mip"} additionally requires every value to
be strictly positive, because it divides each slack by the DMU's own level.}
\item{measure}{Which weights the objective carries:
\describe{
  \item{\code{"ram"}}{Range Adjusted Measure (Cooper, Park and Pastor 1999).
  Each slack is divided by the range of its variable across the reference set,
  and by \eqn{p + q}. Units invariant, and the reported score
  \eqn{\rho = 1 - \mbox{objective}} lies in \eqn{[0, 1]}, because no slack
  can exceed its variable's range. The default.}
  \item{\code{"mip"}}{Measure of Inefficiency Proportions. Each slack is
  divided by the evaluated DMU's own level of that variable, and by
  \eqn{p + q}. Units invariant, and at most 1 --- but \strong{not bounded
  below by 0}: an input slack cannot exceed the DMU's own input, so its term
  is capped, while an output slack has no such ceiling and a badly performing
  DMU can score negative. RAM is bounded on both sides, which is what it was
  designed for, and is the better default.}
  \item{\code{"unweighted"}}{The original 1985 model. Reported as the raw
  slack total, so \strong{0 is efficient and larger is worse} --- and it
  carries the units of the data. See \sQuote{Details}.}
}}
\item{rts}{Returns to scale, as in \code{\link{dea}}. There is no \code{"fdh"}
option --- the measure is defined on a convex technology.}
\item{peers, scaling, xref, yref, dataref}{As in \code{\link{dea}}, except that
\code{scaling} is ignored for \code{measure = "unweighted"}: that measure is
defined on the caller's units, so rescaling the columns would change the
answer rather than leave it alone.}
}
\value{
An object of class \code{"dea"} with \code{model = "additive"} and the chosen
\code{measure} recorded. \code{efficient} is the Pareto-Koopmans efficient set
on every scale.
}
\details{
\strong{One model, three scales.} All three options solve the same program over
the same technology and differ only in the objective weights. The literature
often presents RAM and MIP as separate models; they are one model read three
ways, and the projection and the efficient set are identical across them.

\strong{The classic version is not units invariant, and that is why the other
two exist.} With unit weights the objective adds a slack measured in
staff-hours to one measured in euros: double an input's units and its slack
halves, which can change which DMU looks worst. The efficient \emph{set} is
unaffected --- a slack is zero in any units --- but the \emph{ordering} of the
inefficient DMUs is not. \code{\link{efficiency}(fit, "score")} therefore
refuses to put the unweighted measure on a common scale, and
\code{scaling = TRUE} is ignored for it.

Solving in the caller's units also means the unweighted program is poorly
conditioned when the columns differ wildly in magnitude, and \code{dea_add()}
warns above a spread of \eqn{10^3}. On a design with one input multiplied by
1000, one DMU failed numerically outright. That is a real property of the
measure and another reason to prefer \code{"ram"}.

\strong{No orientation.} The additive model does not ask by what factor inputs
could shrink or outputs grow; it asks how far the DMU is from the frontier in
every coordinate at once. That makes it the natural tool when the question is
\emph{which} units are dominated rather than \emph{by how much}, and it is why
its efficient set coincides exactly with \code{\link{dea_sbm}}'s
\eqn{\rho = 1} set and with \code{\link{dea}}'s Pareto-Koopmans set.

\strong{Infeasibility with an external reference set.} As for
\code{\link{dea_sbm}}, the balance rows are equalities with non-negative
slacks, so a solution exists only if some convex combination of the reference
DMUs weakly dominates the evaluated point in every input and output. A point
outside the estimated technology has none; \code{dea_add()} warns and reports
\code{NA}. For the same reason \code{"ram"} can return a score below 0 when
scoring out of sample, since a slack may then exceed the reference range.
}
\references{
Charnes, A., Cooper, W. W., Golany, B., Seiford, L. and Stutz, J. (1985).
Foundations of data envelopment analysis for Pareto-Koopmans efficient
empirical production functions. \emph{Journal of Econometrics} 30, 91--107.

Cooper, W. W., Park, K. S. and Pastor, J. T. (1999). RAM: a range adjusted
measure of inefficiency for use with additive models, and relations to other
models and measures in DEA. \emph{Journal of Productivity Analysis} 11,
5--42.
}
\seealso{\code{\link{dea_sbm}}, \code{\link{dea}}, \code{\link{dea_ddf}}}
\examples{
data(charnes1981)
x <- charnes1981[, paste0("x", 1:5)]
y <- charnes1981[, paste0("y", 1:3)]

dea_add(x, y, measure = "ram", rts = "vrs")

## The efficient set is the Pareto-Koopmans set, identical to the one the
## radial model finds after its second stage and the one dea_sbm() scores 1.
a  <- dea_add(x, y, measure = "ram", rts = "vrs")
r  <- dea(x, y, rts = "vrs", orientation = "in", slack = TRUE)
identical(unname(a$efficient), unname(r$efficient))

## RAM is units invariant; the unweighted model is not.
x2 <- x; x2[, 1] <- x2[, 1] * 1000
max(abs(dea_add(x2, y, measure = "ram")$eff - dea_add(x, y, measure = "ram")$eff))
max(abs(dea_add(x2, y, measure = "unweighted")$eff -
        dea_add(x, y, measure = "unweighted")$eff))
}
