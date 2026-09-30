\name{dea_wr}
\alias{dea_wr}
\alias{dea_cone}
\alias{wr_bound}
\alias{wr_ratio}
\alias{wr_linear}
\alias{print.dea_wr}
\alias{print.dea_restriction}
\alias{c.dea_restriction}
\title{Weight restrictions, assurance regions and the cone-ratio model}
\description{
An unrestricted DEA score lets every DMU pick the prices that flatter it most,
including prices that value an input or an output at zero. \code{dea_wr()}
solves the radial model with the admissible weights restricted: bounds on
individual multipliers (Dyson and Thanassoulis 1988), bounds on their ratios
--- assurance regions (Thompson et al. 1986) --- and general linear
inequalities (Wong and Beasley 1990). \code{dea_cone()} is the cone-ratio model
of Charnes et al. (1989), which needs no new linear program at all.
}
\usage{
wr_bound(on = c("v", "u"), index, lower = NA, upper = NA)
wr_ratio(on = c("v", "u"), numerator, denominator, lower = NA, upper = NA)
wr_linear(v = NULL, u = NULL, type = c("<=", ">=", "="), rhs = 0)

dea_wr(x, y, restrictions, data = NULL,
       rts = c("vrs", "crs", "nirs", "ndrs"),
       orientation = c("in", "out"),
       scaling = TRUE,
       xref = NULL, yref = NULL, dataref = NULL)

dea_cone(x, y, A = NULL, B = NULL, data = NULL,
         rts = c("crs", "vrs", "nirs", "ndrs"),
         orientation = c("in", "out"), ...)
}
\arguments{
\item{on}{Which family of weights the restriction is about: \code{"v"} for the
inputs, \code{"u"} for the outputs.}
\item{index, numerator, denominator}{Column positions, counted within that
family.}
\item{lower, upper}{Bounds. Give one or both.}
\item{v, u}{For \code{wr_linear()}, coefficient vectors on the input and the
output weights, one entry per column; either may be \code{NULL}.}
\item{type, rhs}{The sense and right-hand side of the inequality.}
\item{x, y, data, rts, orientation, scaling, xref, yref, dataref}{As in
\code{\link{dea}}.}
\item{restrictions}{A restriction set built by the three constructors above and
combined with \code{c()}, or \code{NULL} for none.}
\item{A, B}{Generator matrices for the admissible input and output weight
cones, one row per generator and one column per input or output.
\code{NULL} means the identity, which is no restriction.}
\item{\dots}{Passed to \code{\link{dea}}.}
}
\value{
\code{dea_wr()} returns an object of class \code{c("dea_wr", "dea")} carrying
\code{eff}, the weights \code{v}, \code{u} and \code{u0}, and \code{status}.
\strong{It deliberately carries no \code{lambda}, peers or slacks}; see
\strong{Details}. \code{dea_cone()} returns an ordinary \code{"dea"} object,
because it is an ordinary fit on transformed data, with the generators recorded
in \code{$cone}.
}
\details{
\strong{The score comes from the multiplier program, and with restrictions the
two forms are no longer interchangeable.} Strong duality still holds, but the
dual of the restricted multiplier program is not the ordinary envelopment
program. It acquires one extra non-negative term per restriction, entering the
constraints in exactly the place a slack does:

\deqn{X'\lambda + A_x z \le \theta x_o, \qquad Y'\lambda - A_y z \ge y_o,
      \qquad z \ge 0.}

Those \eqn{Az} terms are the \emph{residue}. They are not slacks and the
\eqn{\lambda} beside them does not project the DMU onto a point of the
technology. So \code{dea_wr()} reports the score and the weights and nothing
else, and \code{peers()} and \code{slacks()} on the result raise an error
explaining why rather than returning something plausible. That is the reason
this is a separate function and not an argument to \code{\link{dea}}.

With \code{restrictions = NULL} the two forms can still be compared, and they
agree to solver precision; the package's tests assert it, because it is what
makes the restricted results believable.

\strong{A restriction can only make a DMU look worse.} Shrinking the admissible
weight set shrinks the feasible region of a maximum, so the input-oriented
score falls; the output-oriented program is a minimum, so its score rises.
Either way the DMU loses ground. A restriction that improved a score would mean
a row written with the wrong sign.

\strong{Restrictions are written in the caller's units.} The solver divides
every column by its mean, so the weight it reports applies to \eqn{x_i/s_i} and
the caller's weight is \eqn{v_i/s_i}; the coefficients of a restriction are
divided by the same factors on the way in. Getting that wrong would not error,
it would quietly impose a different restriction, so the tests require the
restricted score to be identical with \code{scaling} on and off.

\strong{Restrictions can make the program infeasible}, and that is the one
place a DEA program here can be infeasible on ordinary data --- ask for weights
no price vector satisfies and there is no answer. Those DMUs come back
\code{NA} with a warning that points at the restrictions rather than at the
data, because an unrestricted fit of the same data cannot fail that way.

\strong{The cone-ratio model needs no new program.} When the admissible weights
are polyhedral cones in \emph{sum} form, \eqn{V = \{A'\alpha : \alpha \ge 0\}}
and \eqn{U = \{B'\gamma : \gamma \ge 0\}}, substituting \eqn{v = A'\alpha} into
\eqn{v'x} gives \eqn{\alpha'(Ax)}, so the cone-ratio model is exactly a CCR
model on transformed data \eqn{\tilde X = XA'}, \eqn{\tilde Y = YB'} --- and
\eqn{\alpha \ge 0} is what an ordinary DEA program already imposes.

\code{dea_cone()} and \code{dea_wr()} are \emph{not} wrappers for each other. A
cone in sum form becomes a set of inequalities only through its polar, which is
a separate computation and is not generally as small. Where the two can be
written both ways they agree, which the package's tests check on a case where
the restriction binds.
}
\references{
Charnes, A., Cooper, W. W., Wei, Q. L. and Huang, Z. M. (1989). Cone ratio data
envelopment analysis and multi-objective programming.
\emph{International Journal of Systems Science} 20, 1099--1118.

Dyson, R. G. and Thanassoulis, E. (1988). Reducing weight flexibility in data
envelopment analysis. \emph{Journal of the Operational Research Society} 39,
563--576.

Thompson, R. G., Singleton, F. D., Thrall, R. M. and Smith, B. A. (1986).
Comparative site evaluations for locating a high-energy physics lab in Texas.
\emph{Interfaces} 16, 35--49.

Wong, Y.-H. B. and Beasley, J. E. (1990). Restricting weight flexibility in
data envelopment analysis. \emph{Journal of the Operational Research Society}
41, 829--835.
}
\seealso{\code{\link{dea}}, \code{\link{multipliers}}, \code{\link{dea_cross}}}
\examples{
sim <- dea_sim(40, p = 2, q = 2, seed = 3)

## An assurance region: neither input may be valued more than twice the other.
rs <- c(wr_ratio("v", 2, 1, lower = 0.5, upper = 2),
        wr_ratio("u", 2, 1, lower = 0.5, upper = 2))
rs
fit <- dea_wr(sim$x, sim$y, rs, rts = "crs")
fit

## Restricting the weights can only cost a DMU ground.
free <- dea(sim$x, sim$y, rts = "crs", slack = FALSE)
all(fit$eff <= free$eff + 1e-9)

## A floor under one weight, so no input can be valued at nothing.
dea_wr(sim$x, sim$y, wr_bound("v", 1, lower = 0.05), rts = "crs")

## The cone-ratio model, as a transformation of the data.
A <- rbind(c(1, 0.4), c(0.3, 1))
dea_cone(sim$x, sim$y, A = A, rts = "crs")
}
