\name{charnes1981}
\alias{charnes1981}
\docType{data}
\title{Program Follow Through: the original data envelopment analysis data set}
\description{
Seventy US primary school sites from Program Follow Through, a federally
sponsored programme of remedial assistance to disadvantaged pupils. This is the
data set on which Charnes, Cooper and Rhodes introduced and applied the model
that carries their name, and it remains the standard worked example in the
field.
}
\usage{data(charnes1981)}
\format{
A data frame with 70 rows and 10 columns.
\describe{
  \item{\code{site}}{School site name. Not unique: 38 distinct names across 70
  rows, since several cities contributed more than one site --- New York and
  Philadelphia nine each.}
  \item{\code{pft}}{Logical. \code{TRUE} for the 49 sites enrolled in Program
  Follow Through (rows 1--49), \code{FALSE} for the 21 that were not (rows
  50--70).}
  \item{\code{x1}}{Input: education level of the mother.}
  \item{\code{x2}}{Input: highest occupation of a family member.}
  \item{\code{x3}}{Input: parental visits to the school.}
  \item{\code{x4}}{Input: time spent with children on school-related topics.}
  \item{\code{x5}}{Input: number of teachers at the site.}
  \item{\code{y1}}{Output: reading score.}
  \item{\code{y2}}{Output: mathematics score.}
  \item{\code{y3}}{Output: self-esteem score.}
}
}
\details{
Five inputs and three outputs on seventy units, which by the rule of thumb
\eqn{n \ge 3(p+q)} is comfortably enough --- 70 against 24 --- and is
nonetheless a problem in eight dimensions. Expect a large efficient set, and
see the discussion of the curse of dimensionality in \code{\link{dea_rate}}
before reading much into the scores of the units that reach 1.

The \code{pft} flag makes this a natural example for an explicit reference
technology: scoring the enrolled sites against the frontier spanned by the
non-enrolled ones, or the reverse, is a metafrontier comparison and is what
\code{xref}/\code{yref} are for.
}
\source{
Charnes, A., Cooper, W. W. and Rhodes, E. (1981). Evaluating program and
managerial efficiency: an application of data envelopment analysis to Program
Follow Through. \emph{Management Science} 27, 668--697.

The values here are identical to those distributed as
\code{Benchmarking::charnes1981} and \code{npsf::ccr81}; the two were checked
against each other before this copy was made and agree on every one of the
\eqn{70 \times 8} measurements. The columns have been renamed and the
programme indicator turned into a logical, but no number has been altered.
}
\note{
The site names are as distributed, including the misspellings
\dQuote{Chigago} and \dQuote{Philidelphia}. They are left alone rather than
silently corrected, so that this copy can be compared against the others
without a diff appearing where none belongs.
}
\references{
Charnes, A., Cooper, W. W. and Rhodes, E. (1978). Measuring the efficiency of
decision making units. \emph{European Journal of Operational Research} 2,
429--444.
}
\examples{
data(charnes1981)
x <- charnes1981[, c("x1", "x2", "x3", "x4", "x5")]
y <- charnes1981[, c("y1", "y2", "y3")]

## The model of the original paper: constant returns, input oriented.
ccr <- dea(x, y, rts = "crs", orientation = "in")
ccr

## Variable returns, and the scale decomposition between the two.
dea_rts(x, y, orientation = "in")

\donttest{
## Score the enrolled sites against the frontier spanned by the sites that
## were NOT enrolled -- a metafrontier comparison, not an internal ranking.
inp <- charnes1981$pft
dea(x[inp, ], y[inp, ], rts = "vrs", orientation = "in",
    xref = x[!inp, ], yref = y[!inp, ])
}
}
\keyword{datasets}
