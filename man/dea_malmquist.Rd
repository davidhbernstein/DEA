\name{dea_malmquist}
\alias{dea_malmquist}
\alias{dea_panel}
\alias{print.dea_panel}
\alias{summary.dea_panel}
\alias{print.dea_malmquist}
\alias{summary.dea_malmquist}
\title{Panel data and the Malmquist productivity index}
\description{
\code{dea_panel()} holds a long-format panel --- one row per DMU per period ---
and is the layer every index over time is built on. \code{dea_malmquist()}
computes the Malmquist productivity index of Fare, Grosskopf, Norris and Zhang
(1994) on it, for each pair of consecutive periods, with its decomposition into
efficiency change and technical change and the further split of efficiency
change into a pure and a scale part.
}
\usage{
dea_panel(x, y, id, period, data = NULL)

dea_malmquist(x, y = NULL, id = NULL, period = NULL, data = NULL,
              rts = c("crs", "vrs", "nirs", "ndrs"),
              orientation = c("in", "out"),
              scale_decomposition = TRUE,
              scaling = TRUE)
}
\arguments{
\item{x, y}{Inputs and outputs, stacked over DMU-periods, in any of the forms
\code{\link{dea}} accepts. For \code{dea_malmquist()}, \code{x} may instead be
a \code{"dea_panel"} object, in which case \code{y}, \code{id} and
\code{period} are taken from it.}
\item{id}{Which DMU each row belongs to: a vector, or --- with \code{data} ---
a column name or a one-sided formula. Labels, not numbers: they are compared,
never arithmetic.}
\item{period}{Which period each row belongs to, in the same forms. \strong{The
order matters}; see \strong{Details}.}
\item{data}{Optional data frame the column names or formulas refer to.}
\item{rts}{Returns to scale. The default \code{"crs"} is the technology the
index is defined on; see \strong{Details} before changing it.}
\item{orientation}{\code{"in"} or \code{"out"}. Both report \eqn{M > 1} as
growth.}
\item{scale_decomposition}{Split efficiency change into a pure and a scale
part. Defined only for the \code{"crs"} index and ignored, with a warning if it
was asked for explicitly, otherwise.}
\item{scaling}{As in \code{\link{dea}}.}
}
\value{
\code{dea_panel()} returns an object of class \code{"dea_panel"} with the
stacked \code{X} and \code{Y}, the \code{id} and \code{period} of each row, the
\code{periods} in the order they will be used, and whether the panel is
\code{balanced}.

\code{dea_malmquist()} returns an object of class \code{"dea_malmquist"}. Its
\code{table} has one row per DMU per consecutive period pair, with columns
\code{id}, \code{from}, \code{to}, \code{malmquist}, \code{effch},
\code{techch} and --- under \code{rts = "crs"} --- \code{purech} and
\code{scalech}. \code{distances} carries the four distance functions each row
was built from. A value above 1 is growth in every column.
}
\details{
\strong{The four distance functions, and what their names mean.} For a DMU
observed in periods \eqn{t} (\dQuote{from}) and \eqn{t+1} (\dQuote{to}), the
index needs that DMU scored against both periods' technologies:

\tabular{ll}{
  \code{d_ff} \tab the from-period observation on the from-period technology \cr
  \code{d_tt} \tab to on to \cr
  \code{d_tf} \tab the to-period observation on the \emph{from}-period technology \cr
  \code{d_ft} \tab from on \emph{to}
}

The convention is \strong{observation first, technology second}, and it is
worth stating because the two mixed terms are easy to swap and swapping them
replaces efficiency change with its reciprocal. Each is an ordinary
\code{\link{dea}} fit with \code{xref}/\code{yref} pointing at the other
period; the separation of the evaluated and reference sets is what makes this
possible at all. Then

\deqn{EC = d_{tt}/d_{ff}, \quad
      TC = \sqrt{(d_{tf}/d_{tt})(d_{ff}/d_{ft})}, \quad M = EC \times TC.}

Every score is put on the common \eqn{(0, 1]} efficiency scale first, so the
output-oriented case converts \eqn{\phi} to \eqn{1/\phi} once at the source
rather than carrying a second set of formulas. The mixed terms routinely leave
\eqn{(0, 1]} --- an observation scored against another period's frontier need
not lie inside it --- and that is the measurement rather than an error.

\strong{A convention the field is split on.} This package follows Fare et
al.'s: the index is built from Shephard distance functions, so \eqn{M > 1} is
productivity growth whichever side it is measured from, and the two
orientations agree under \code{"crs"}. \code{productivity::malm()} reports the
same numbers in both orientations. \code{Benchmarking::malmq()} instead builds
the index from the native scale of each orientation, so \strong{its
output-oriented index is the reciprocal of this one}. Neither is wrong; they
are different conventions, and a number copied between them without checking
is.

\strong{Why the default is constant returns.} The mixed distance functions can
be \emph{infeasible} under variable returns: a convex hull cannot extrapolate
to a point outside it, while a cone can always be scaled to reach one. That is
why the index is conventionally computed on a cone. \code{"vrs"} is allowed
here; the infeasible terms come back \code{NA} and the count is reported per
period pair rather than left to be discovered.

\strong{The pure/scale split is an identity, so it is worth saying what it
does and does not establish.} With \eqn{v} the variable-returns own-period
scores,

\deqn{PEC = v_{tt}/v_{ff}, \qquad SEC = (d_{tt}/v_{tt})/(d_{ff}/v_{ff}),}

and \eqn{EC = PEC \times SEC} follows by cancellation --- so that product
cannot fail and checking it establishes nothing. What the package's tests check
is that each piece equals the same quantity computed from independent
variable- and constant-returns fits.

\strong{The order of the periods is data, not a guess.} A factor's levels are
honoured as given; anything else is sorted, which is right for numbers and for
ISO-like dates and \emph{wrong} for labels such as \code{"Q1"}, \code{"Q10"},
\code{"Q2"}. Reversing the order inverts every index, so the panel object
prints the order it settled on and a factor is how to override it.

\strong{An unbalanced panel is an ordinary case.} Each period pair uses the
DMUs observed in both; anything seen in only one has no index for that pair.
\code{summary()} on the panel reports how many DMUs each pair has and how many
observations went unmatched.

\strong{The technology is estimated per period}, from that period's DMUs alone.
Pooling the periods would define the frontier shift away.
}
\references{
Fare, R., Grosskopf, S., Norris, M. and Zhang, Z. (1994). Productivity growth,
technical progress, and efficiency change in industrialized countries.
\emph{American Economic Review} 84, 66--83.

Caves, D. W., Christensen, L. R. and Diewert, W. E. (1982). The economic theory
of index numbers and the measurement of input, output, and productivity.
\emph{Econometrica} 50, 1393--1414.
}
\seealso{\code{\link{dea}}, \code{\link{dea_rts}}}
\examples{
## A two-period panel with a frontier that moves outward.
a <- dea_sim(25, p = 2, q = 1, seed = 1)
b <- dea_sim(25, p = 2, q = 1, seed = 2)
pan <- dea_panel(rbind(a$x, b$x * 0.92), rbind(a$y, b$y * 1.12),
                 id = rep(1:25, 2), period = rep(c(2023, 2024), each = 25))
pan

m <- dea_malmquist(pan, rts = "crs")
m
head(m$table)

## The four distances the index was built from.
head(m$distances)

## Column names work too, when the panel is a data frame.
df <- data.frame(firm = rep(1:25, 2), year = rep(c(2023, 2024), each = 25),
                 k = pan$X[, 1], l = pan$X[, 2], q = pan$Y[, 1])
dea_malmquist(c("k", "l"), "q", "firm", "year", data = df, rts = "crs")
}
