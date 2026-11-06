\name{dea_ladder}
\alias{dea_ladder}
\alias{print.dea_ladder}
\alias{summary.dea_ladder}
\alias{plot.dea_ladder}
\title{The efficiency ladder: how much of a score one peer is worth}
\description{
Drop the peer a DMU's score leans on hardest, re-score it against what is left,
and repeat. The result says how much of the score was owed to a single
reference unit --- the practical form of DEA's sensitivity to one observation,
for which the estimator itself offers no residual and no standard error.
}
\usage{
dea_ladder(object, dmu, max_steps = NULL, threshold = 1e-8)
}
\arguments{
\item{object}{A fit from \code{\link{dea}}: the radial model, scored against
itself, with \code{peers = TRUE} and \code{super = FALSE}. See
\strong{Details} for why it takes a fit rather than data.}
\item{dmu}{The DMU to audit, by name or by position. One DMU.}
\item{max_steps}{Stop after this many rungs. The default runs the ladder out.}
\item{threshold}{A \code{lambda} above this makes a reference DMU a peer, as in
\code{\link{peers}}.}
}
\value{
An object of class \code{"dea_ladder"}. Its \code{table} has one row per rung:
\code{step} (0 is the fit as reported), \code{dropped} (the peer removed at
that rung), \code{lambda} (the weight that peer carried \emph{in the reference
set it was dropped from}, not in the original fit), \code{n_ref} (reference
DMUs remaining), \code{eff} and \code{status}. \code{print()} shows the first
ten rungs, \code{summary()} adds the single largest move and its share of the
total, and \code{plot()} draws the climb.
}
\details{
\strong{It only ever moves one way, and that is a check rather than a
coincidence.} Removing reference DMUs can only shrink the estimated
technology, so the evaluated point can only become \emph{less} enveloped:
\eqn{\theta} weakly increases toward 1 under an input orientation and
\eqn{\phi} weakly decreases toward 1 under an output one. A ladder that moved
the other way would mean the refits were not nested.

\strong{The evaluated DMU stays in the reference set.} Removing it would make
the program unbounded in the usual case and meaningless in the rest. Keeping it
is also what makes the ladder terminate: once its own peers are gone it
references itself alone and scores exactly 1. The ladder therefore ends on the
frontier by construction, and the number worth reading is how many rungs that
took and where the score moved most --- not the endpoint.

\strong{Which peer goes first is defined, not left to the solver.} Each rung
drops the remaining peer with the largest \code{lambda}, ties broken by
reference-set order. That matters because peer sets are not unique when the
program has multiple optima, which is common, so an undefined rule would give a
different ladder on a different machine. What the ladder measures is the
sensitivity of the score to the peer set \emph{this fit returned};
\code{peers()} shows it, and the caveat \code{summary()} prints about
non-uniqueness applies here in full.

\strong{Why it takes a fit rather than data.} Every setting that defines the
technology --- \code{rts}, \code{orientation}, \code{scaling} --- has to match
the fit whose score is being audited, and re-specifying them by hand invites
auditing a different model from the one fitted. The radial model only: the
non-radial fits do not carry the arguments a refit would need (a direction, a
measure), and guessing them would be the same mistake one level down.

\strong{Reading it with the isolated DMUs \code{summary()} reports.}
\code{summary()} on a \code{"dea"} fit lists the efficient DMUs that envelop no
one --- units nothing in the sample is benchmarked against, efficient by
default rather than by comparison. The ladder is the other side of that
question: for an \emph{inefficient} DMU, how much of its score rests on any one
of the units that do envelop it. A DMU whose score jumps most of the way to 1
on the first rung was being measured against a single observation.
}
\references{
Charnes, A., Cooper, W. W., Lewin, A. Y. and Seiford, L. M. (1994).
\emph{Data Envelopment Analysis: Theory, Methodology and Application}. Kluwer.
Chapter 21, section 6.

Bogetoft, P. and Otto, L. (2011). \emph{Benchmarking with DEA, SFA, and R}.
Springer, New York.
}
\seealso{\code{\link{peers}}, \code{\link{dea}},
\code{\link{dea_boot}} for sampling noise rather than peer dependence}
\examples{
s <- dea_sim(30, p = 2, q = 1, seed = 7)
fit <- dea(s$x, s$y, rts = "vrs", orientation = "in")

## The least efficient DMU: how much of that score is one peer worth?
lad <- dea_ladder(fit, dmu = which.min(fit$eff))
lad
summary(lad)

## The score can only rise as the technology shrinks.
all(diff(lad$table$eff) >= -1e-9)

## Several DMUs at once, as a table of first-rung sensitivities.
who <- order(fit$eff)[1:3]
sapply(who, function(i) {
  l <- dea_ladder(fit, dmu = i, max_steps = 1)
  diff(l$table$eff)
})
}
