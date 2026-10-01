\name{dea_nd}
\alias{dea_nd}
\alias{print.dea_nd}
\title{Non-discretionary inputs and outputs}
\description{
Radial efficiency when some inputs or outputs are outside the DMU's control:
they still constrain who may be a peer, but no credit is given for changing
them. Banker and Morey (1986a).
}
\usage{
dea_nd(x, y, data = NULL, nd_x = NULL, nd_y = NULL,
       rts = c("vrs", "crs", "nirs", "ndrs"),
       orientation = c("in", "out"),
       slack = TRUE, peers = TRUE, scaling = TRUE,
       xref = NULL, yref = NULL, dataref = NULL)

\method{print}{dea_nd}(x, \dots)
}
\arguments{
\item{x, y, data}{Inputs and outputs, in any of the forms accepted by
\code{\link{dea}}. For the \code{print} method, an object of class
\code{"dea_nd"}.}
\item{nd_x, nd_y}{Which inputs and outputs are non-discretionary, as a logical
mask, column positions, or column names. \code{NULL} means none, and the fit is
then exactly \code{\link{dea}}'s.}
\item{rts, orientation, slack, peers, scaling, xref, yref, dataref}{As in
\code{\link{dea}}.}
\item{\dots}{Ignored.}
}
\value{
An object of class \code{c("dea_nd", "dea")}, with \code{nd_x} and \code{nd_y}
as named logical masks and \code{n_fixed} counting them. \code{efficient} is
Pareto-Koopmans over the \emph{discretionary} coordinates only: residual slack
on a fixed variable is not something the DMU failed to do.
}
\details{
Some of a DMU's inputs are not its to choose. A school cannot vary the number of
pupils it is sent; a hospital cannot choose its catchment. The ordinary radial
program contracts every input by the same \eqn{\theta}, so it reports a unit as
inefficient for failing to do something it has no power to do.

The change is one row type:
\deqn{\sum_j \lambda_j x_{ij} \le \theta x_{io} \quad i \mbox{ discretionary}}
\deqn{\sum_j \lambda_j x_{ij} \le x_{io} \qquad\;\; i \mbox{ fixed}}
so a fixed input still constrains who may be a peer --- a unit with twice the
pupils is not a fair comparison --- while no credit is given for cutting it.

\strong{The two stages mean different things.} Stage one scales only the
oriented side, so only \code{nd_x} changes an input-oriented score and only
\code{nd_y} an output-oriented one. Stage two maximises a sum of slacks, and
there \emph{both} masks apply in \emph{both} orientations: a slack the DMU
cannot act on should earn it no credit whichever way the program is oriented.
That is Banker and Morey's own rule --- the non-discretionary slacks are dropped
from the objective and kept in the constraints --- and it is why \code{nd_y} is
not ignored under an input orientation even though it leaves \eqn{\theta} alone.
}
\section{Fixing a variable makes the measured inefficiency larger}{
Every intuition says otherwise, so it is worth stating flatly:
\eqn{\theta_{nd} \le \theta_{ordinary}}, always. Measured over 8704 DMU-fits
across four technologies and random dimensions, the largest excess of
\eqn{\theta_{nd}} over \eqn{\theta_{ordinary}} was 1.2e-12, which is solver
noise; on the example below the mean \eqn{\theta} falls from 0.901 to 0.789.
Output oriented the mirror holds, \eqn{\phi_{nd} \ge \phi_{ordinary}}.

The reason is one line about the feasible set. For \eqn{\theta < 1} the fixed row
\eqn{\sum_j \lambda_j x_{ij} \le x_{io}} is \emph{looser} than
\eqn{\le \theta x_{io}}, so the non-discretionary program's feasible set contains
the ordinary one and its minimum can only fall. Economically: \eqn{\theta} in the
ordinary model is a single proportional contraction of every input at once, so a
peer must use proportionally less of everything, pupils included --- a demanding
requirement. Here a peer need only use no more pupils, and the whole contraction
lands on the teachers.

So the two are not the same quantity measured more or less fairly. They answer
different questions --- "could all inputs have been smaller together?" against
"with these pupils, could the teachers have been fewer?" --- and a table
comparing them as though one were a corrected version of the other is comparing
two estimands.
}
\section{Non-discretionary outputs are not the mirror image}{
The symmetry looks exact and is not, and chapter 21 of Charnes, Cooper, Lewin and
Seiford (1994) flags the confusion as a recurring interpretational error.

"Non-discretionary input" has a clean meaning: a quantity the DMU is handed.
"Non-discretionary output" is usually used for something quite different --- an
output the DMU influences only \emph{indirectly}, through its inputs, such as
sales through advertising. That is not fixed; it is controlled at one remove, and
this formulation does not model it. Marking it \code{nd_y} tells the program the
DMU should get no credit for producing more of it, which for sales is plainly
wrong.

Use \code{nd_y} only for an output that is genuinely exogenous --- a statutory
quota, a count of mandated cases. If the right description is "controlled
indirectly", neither this nor any other model in this package says so, and
\code{dea_nd()} prints a message whenever \code{nd_y} is used so the distinction
is made deliberately.
}
\references{
Banker, R. D. and Morey, R. C. (1986a). Efficiency analysis for exogenously
fixed inputs and outputs. \emph{Operations Research} 34, 513--521.

Charnes, A., Cooper, W. W., Lewin, A. Y. and Seiford, L. M. (1994).
\emph{Data Envelopment Analysis: Theory, Methodology and Application}. Kluwer.
}
\seealso{\code{\link{dea_categorical}} for an environment that is ordered rather
than measured, \code{\link{dea}}.}
\examples{
set.seed(31)
n <- 60
teachers <- runif(n, 10, 40)
pupils <- runif(n, 150, 600)          # sent to the school, not chosen by it
x <- cbind(teachers = teachers, pupils = pupils)
y <- matrix(teachers^0.3 * pupils^0.6 * exp(-abs(rnorm(n, 0, 0.2))), n, 1,
            dimnames = list(NULL, "attainment"))

r <- dea_nd(x, y, nd_x = "pupils")
r

## Lower than the ordinary fit, which is the direction to expect.
c(ordinary = mean(dea(x, y)$eff), non_discretionary = mean(r$eff))
}
