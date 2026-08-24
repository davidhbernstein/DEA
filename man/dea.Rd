\name{dea}
\alias{dea}
\title{Radial technical efficiency}
\description{
Debreu-Farrell radial efficiency by data envelopment analysis, under any of
five returns-to-scale assumptions and either orientation, with optional
second-stage slacks and Andersen-Petersen super-efficiency.
}
\usage{
dea(x, y, data = NULL,
    rts = c("vrs", "crs", "nirs", "ndrs", "fdh"),
    orientation = c("in", "out"),
    slack = TRUE, super = FALSE, peers = TRUE, scaling = TRUE,
    xref = NULL, yref = NULL, dataref = NULL)
}
\arguments{
\item{x}{Inputs, one row per decision-making unit (DMU) and one column per
input. A matrix, a data frame, or a numeric vector for a single input. With
\code{data} supplied it may instead be a character vector of column names or a
one-sided formula such as \code{~ labour + capital}.}
\item{y}{Outputs, in the same forms as \code{x}.}
\item{data}{Optional data frame holding the columns named by \code{x} and
\code{y}.}
\item{rts}{Returns to scale, i.e. the restriction placed on the intensity
weights:
\describe{
  \item{\code{"vrs"}}{variable, \eqn{\sum\lambda = 1}; the BCC model of Banker,
  Charnes and Cooper (1984). The default.}
  \item{\code{"crs"}}{constant, \eqn{\lambda} unrestricted; the CCR model of
  Charnes, Cooper and Rhodes (1978).}
  \item{\code{"nirs"}}{non-increasing, \eqn{\sum\lambda \le 1}.}
  \item{\code{"ndrs"}}{non-decreasing, \eqn{\sum\lambda \ge 1}.}
  \item{\code{"fdh"}}{free disposal hull, \eqn{\lambda} binary; drops convexity
  altogether (Deprins, Simar and Tulkens 1984).}
}
Matching ignores case, and \code{"drs"}/\code{"irs"} are accepted as aliases
for \code{"nirs"}/\code{"ndrs"}, which is how \pkg{Benchmarking} spells them.}
\item{orientation}{\code{"in"} contracts inputs at fixed output and returns
\eqn{\theta \in (0, 1]}; \code{"out"} expands outputs at fixed input and
returns \eqn{\phi \ge 1}. Both equal 1 on the frontier. Use
\code{\link{efficiency}(fit, "score")} to put either on a common \eqn{(0, 1]}
scale.}
\item{slack}{Run the second-stage slack-maximizing program. Its cost is small
-- the constraint matrix does not depend on the DMU -- and without it the
\code{efficient} component reports radial efficiency rather than
Pareto-Koopmans efficiency. See \sQuote{Details}.}
\item{super}{Andersen-Petersen (1993) super-efficiency: score each DMU against
a technology it is excluded from, which breaks the ties among the efficient
units. See \sQuote{Details} for the infeasibility that comes with it.}
\item{peers}{Keep the \eqn{n \times n} matrix of intensity weights. Set
\code{FALSE} for large samples, where it is the dominant memory cost.}
\item{scaling}{Rescale every input and output column to mean one before
solving. Radial efficiency is invariant to the units of each variable, so this
changes no answer and conditions the linear program considerably better;
slacks are returned in the original units either way.}
\item{xref, yref}{An explicit reference technology: the DMUs that SPAN the
frontier, where \code{x} and \code{y} are the DMUs SCORED against it. Both must
be given together. Default \code{NULL} uses \code{x} and \code{y} themselves,
which is ordinary DEA.}
\item{dataref}{Data frame for \code{xref}/\code{yref} when those are given as
names or formulas. Defaults to \code{data}.}
}
\value{
An object of class \code{"dea"}: a list whose components are
\item{eff}{Named numeric vector of efficiency scores, on the orientation's own
scale.}
\item{model, rts, orientation, super}{The specification, as resolved.}
\item{lambda}{Matrix of intensity weights, evaluated DMUs by reference DMUs, or
\code{NULL} when \code{peers = FALSE}.}
\item{sum_lambda}{Row sums of \code{lambda} from the first-stage program.}
\item{slack_x, slack_y}{Input and output slacks in the original units, or
\code{NULL} when \code{slack = FALSE}.}
\item{efficient}{Logical. Pareto-Koopmans efficiency when slacks were computed
-- on the frontier radially AND with no slack left -- and radial efficiency
otherwise.}
\item{status}{Per-DMU \pkg{lpSolveAPI} return code; 0 is an optimal solve.}
\item{x, y, xref, yref}{The data, as coerced.}
\item{n, nref, p, q}{Counts of evaluated DMUs, reference DMUs, inputs and
outputs.}
\item{total_time}{Elapsed seconds.}
\item{call}{The matched call.}
}
\details{
\strong{The radial score is not the whole story.} A DMU projected onto a
vertical or horizontal face of the frontier is radially efficient and still
dominated: some input can be cut further without giving up any output. The
second stage maximizes total slack at the radial projection, and
\code{efficient} combines the two into the Pareto-Koopmans judgement. Where the
two differ, \code{print()} says so. \code{\link{dea_sbm}} builds the same
judgement into a single number.

\strong{Peer sets are not unique.} These linear programs frequently have
multiple optima, so \code{lambda} is one valid set of benchmarks rather than
the set. Efficiency scores are unique; the weights that achieve them generally
are not, and no DEA package can make them so.

\strong{Super-efficiency and infeasibility.} Removing a DMU from its own
reference set can leave the program with no solution at all: under
\code{"vrs"}, \code{"nirs"} or \code{"ndrs"} a DMU at the edge of the input
space may have no other unit able to dominate it. \code{eff} is \code{NA}
there. That is the correct answer and not a numerical failure -- a package
that returns a large finite number instead is reporting something the model
does not define. Infeasibility cannot occur for input-oriented \code{"crs"}.

\strong{The curse of dimensionality.} DEA is consistent but converges slowly,
at a rate that worsens with every variable added: see \code{\link{dea_rate}}.
With \eqn{n < 3(p+q)} most DMUs are efficient by dimension alone and
\code{dea()} warns.

\strong{Sampling variation.} Because the frontier is spanned by the observed
DMUs it lies inside the true one, so every score is biased toward 1. See
\code{\link{dea_boot}} for the bias correction and confidence intervals.
}
\references{
Andersen, P. and Petersen, N. C. (1993). A procedure for ranking efficient
units in data envelopment analysis. \emph{Management Science} 39, 1261--1264.

Banker, R. D., Charnes, A. and Cooper, W. W. (1984). Some models for estimating
technical and scale inefficiencies in data envelopment analysis.
\emph{Management Science} 30, 1078--1092.

Charnes, A., Cooper, W. W. and Rhodes, E. (1978). Measuring the efficiency of
decision making units. \emph{European Journal of Operational Research} 2,
429--444.

Deprins, D., Simar, L. and Tulkens, H. (1984). Measuring labor-efficiency in
post offices. In M. Marchand, P. Pestieau and H. Tulkens (eds),
\emph{The Performance of Public Enterprises}. North-Holland.
}
\seealso{\code{\link{dea_sbm}}, \code{\link{dea_ddf}}, \code{\link{dea_rts}},
\code{\link{dea_boot}}, \code{\link{efficiency}}, \code{\link{dea_sim}}}
\examples{
sim <- dea_sim(60, p = 2, q = 1, seed = 1)

fit <- dea(sim$x, sim$y, rts = "vrs", orientation = "in")
fit
head(efficiency(fit))
head(peers(fit))

## The same call from a data frame, naming columns.
d <- sim$data
dea(x = c("x1", "x2"), y = "y1", data = d, rts = "crs")

## Score one group against another group's frontier.
dea(sim$x[1:10, ], sim$y[1:10, , drop = FALSE], rts = "vrs",
    xref = sim$x[11:60, ], yref = sim$y[11:60, , drop = FALSE])
}
