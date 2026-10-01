\name{dea_categorical}
\alias{dea_categorical}
\alias{print.dea_categorical}
\alias{summary.dea_categorical}
\title{Categorical environmental variables, as nested reference sets}
\description{
Score each DMU only against units whose environment was no more favourable than
its own. The categorical model of Banker and Morey (1986b), implemented by
restricting the reference set rather than by their mixed-integer program.
}
\usage{
dea_categorical(x, y, category, data = NULL,
                rts = c("vrs", "crs", "nirs", "ndrs"),
                orientation = c("in", "out"),
                slack = TRUE, peers = TRUE, scaling = TRUE)

\method{print}{dea_categorical}(x, \dots)
\method{summary}{dea_categorical}(object, \dots)
}
\arguments{
\item{x, y, data}{Inputs and outputs, in any of the forms accepted by
\code{\link{dea}}. For the \code{print} method, an object of class
\code{"dea_categorical"}.}
\item{category}{The environment, as an \strong{ordered} factor or a numeric
rank, one entry per DMU, or a column name when \code{data} is given.
\strong{Higher means more favourable.} An unordered factor is refused.}
\item{rts, orientation, slack, peers, scaling}{As in \code{\link{dea}}.}
\item{object}{An object of class \code{"dea_categorical"}.}
\item{\dots}{Ignored.}
}
\value{
An object of class \code{c("dea_categorical", "dea")}. Beyond the usual
components:
\item{category, category_levels, n_category}{The environment as an ordered
factor, its levels from least to most favourable, and how many there are.}
\item{reference_size}{How many DMUs each unit was scored against.}
\item{pooled_eff}{What the same fit gives with no categories at all.}
\item{category_effect}{How much the DMU's efficiency improves, on the common
\eqn{(0,1]} scale, from being compared only with units no better placed. Never
negative.}
}
\details{
A school in a deprived catchment, a hospital in a rural district, a branch in a
small town: the DMU operates in an environment it did not choose and cannot
change, and scoring it against units in easier environments measures the
environment rather than the management. The environment is \emph{ordered} but is
not a quantity --- "more deprived" is a rank, not a number of units of anything
--- so it cannot simply be added as an input.

The fix is a restriction on who may be a peer, and nothing else. Order the
categories from least to most favourable; a DMU in category \eqn{k} is scored
against every DMU in categories \eqn{1..k}. That is \code{xref}/\code{yref}
applied once per category: no new program and no integer variables. Two
identities follow and are worth knowing, because they make the output readable:
the \emph{top} category's scores are exactly the pooled ones, and the
\emph{bottom} category's are exactly an ordinary fit on that subset alone.

\strong{Banker and Morey's own formulation is deliberately not built.} Theirs is
a mixed-integer program, and Charnes, Cooper, Lewin and Seiford (1994, ch. 3) are
explicit that it "incorrectly prescribes categories that are undefined or
meaningless" --- it can hand a DMU a target in a category that does not exist.
Adapting the solution procedure rather than the model is both correct and allows
several categorical variables at once, which the mixed-integer version does not.

\strong{Which end is which is asked for explicitly}, because getting it backwards
does not error, it silently answers the opposite question. \code{category} must
be an ordered factor or a numeric rank; an unordered factor is refused, since
alphabetical order is not an ordering of environments. The reference-set size per
category is printed so the direction is visible in the output rather than taken
on trust.

\strong{Nesting costs sample size, and the cost falls on category 1}, which is
scored against its own members alone. A design that is comfortable pooled can be
far below the \eqn{n \ge 3(p+q)} rule of thumb once split, with almost everything
efficient by dimension rather than by performance. The scores stay finite and
plausible, so a warning is the only notice there is, and it is raised per
category with the category named.

\strong{If the categories are not nested, this is the wrong function.} Public
versus private universities are not ordered --- neither environment's advantages
are a subset of the other's --- and inventing an order produces a number that
answers nothing. The right analysis there is a separate one per group, or a
metafrontier through \code{xref}/\code{yref}. \code{summary()} says so.
}
\references{
Banker, R. D. and Morey, R. C. (1986b). The use of categorical variables in data
envelopment analysis. \emph{Management Science} 32, 1613--1627.

Charnes, A., Cooper, W. W., Lewin, A. Y. and Seiford, L. M. (1994).
\emph{Data Envelopment Analysis: Theory, Methodology and Application}. Kluwer.
}
\seealso{\code{\link{dea_nd}} for variables that are continuous but outside the
DMU's control, \code{\link{dea}}.}
\examples{
set.seed(21)
n <- 90
x <- matrix(runif(n * 2, 2, 10), n, 2,
            dimnames = list(NULL, c("staff", "budget")))
env <- rep(1:3, length.out = n)
y <- matrix((x[, 1] * x[, 2])^0.4 * c(1, 1.15, 1.32)[env] *
            exp(-abs(rnorm(n, 0, 0.2))), n, 1, dimnames = list(NULL, "output"))
catg <- factor(c("deprived", "average", "affluent")[env],
               levels = c("deprived", "average", "affluent"), ordered = TRUE)

r <- dea_categorical(x, y, catg)
r

## The top category is scored against everyone, so it gains nothing.
tapply(r$category_effect, r$category, mean)
}
