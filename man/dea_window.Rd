\name{dea_window}
\alias{dea_window}
\alias{print.dea_window}
\alias{summary.dea_window}
\alias{plot.dea_window}
\title{Window analysis}
\description{
Window analysis: a DMU observed in each of \eqn{T} periods is treated as
\eqn{T} distinct DMUs inside a rolling window of \eqn{w} consecutive periods,
giving \eqn{T - w + 1} analyses of up to \eqn{n w} units each. It is the
standard answer to a panel whose individual periods are too short to
discriminate, and it measures how stable a DMU's score is across the windows it
appears in.
}
\usage{
dea_window(x, y = NULL, id = NULL, period = NULL, data = NULL,
           width = 3L,
           rts = c("vrs", "crs", "nirs", "ndrs", "fdh"),
           orientation = c("in", "out"),
           scaling = TRUE)
}
\arguments{
\item{x, y}{Inputs and outputs, stacked over DMU-periods, in any of the forms
\code{\link{dea}} accepts. \code{x} may instead be a \code{"dea_panel"}
object, in which case \code{y}, \code{id} and \code{period} are taken from it.}
\item{id}{Which DMU each row belongs to, as in \code{\link{dea_panel}}.}
\item{period}{Which period each row belongs to. \strong{The order matters} and
comes from the panel; see \strong{Details}.}
\item{data}{Optional data frame the column names or formulas refer to.}
\item{width}{The window width in periods, at least 2 and at most the number of
periods. \strong{A vector is allowed and is the intended use}; see
\strong{Details}.}
\item{rts, orientation, scaling}{As in \code{\link{dea}}, applied to every
window.}
}
\value{
An object of class \code{"dea_window"}. Its \code{table} has one row per
(width, window, DMU-period) with columns \code{width}, \code{window},
\code{win_from}, \code{win_to}, \code{id}, \code{period}, \code{eff} and
\code{status}; \code{panel} is the panel it was built on. \code{print()} gives
the mean score by window, \code{summary()} the per-DMU mean, \code{range}
and \code{ldr} (see Details), and \code{plot()} the mean score by window, one line per width.
}
\details{
\strong{The output is a spread, not a score.} A DMU-period in the interior of
the panel appears in up to \eqn{w} windows and therefore receives up to
\eqn{w} different scores, because each window's technology is estimated from a
different set of points. Reading any one of them as \dQuote{the} score throws
away what the method produces. \code{summary()} reports two spreads per DMU,
and only one of them is about the windows. \code{ldr} is the largest
difference between scores of the \emph{same} DMU-period in different windows
of the same width: the observation is held fixed and only the window moves, so
a DMU with a small \code{ldr} is being measured by a stable technology and one
with a large \code{ldr} has a score that depends on which periods its window
happens to contain. \code{range} is the spread of every score the DMU
received, over all its periods, windows and widths; it also includes the DMU's
genuine change between periods, so it can be large under a perfectly stable
technology.

\strong{Why the pooling is the point.} DEA discriminates badly when the sample
is small relative to \eqn{p + q} --- the usual heuristic is
\eqn{n \ge 3(p + q)}, which \code{\link{dea}} warns about --- and a single
period of a short panel often fails it where a window of three periods does
not. The price is that each DMU-period is treated as a separate unit, so a
window's technology is spanned by \eqn{n w} points rather than \eqn{n}.

\strong{What pooling assumes, and it is the opposite of Malmquist.} Putting
several periods into one technology assumes the frontier did not move within
the window. A DMU that looks worse in a later window may simply be facing a
frontier that advanced. Window analysis measures stability with the frontier
held fixed; \code{\link{dea_malmquist}} measures the movement of the frontier
itself. Neither substitutes for the other, and they are routinely confused.

\strong{The width is ad hoc, and the literature says so.} No theory picks
\eqn{w}. The standard advice is to run two or three widths and see whether the
answer moves, so \code{width} takes a vector and the result carries a
\code{width} column --- the advice costs an argument rather than a loop to be
reconciled by hand. A width equal to the number of periods is one window, which
is a pooled fit rather than a window analysis; it is allowed, and it is the
only case the method degenerates in.

\strong{The period order is inherited from \code{\link{dea_panel}}.} A window
is a contiguous block of that order, so a wrong order does not fail here --- it
silently groups periods that are not neighbours. The panel object prints the
order it settled on, and a factor is how to set it.

\strong{Unbalanced panels are an ordinary case.} A window contains whatever
DMU-periods were observed in its periods, so a DMU absent from one period
simply has fewer scores. The \code{n} column of \code{summary()} says how many
each DMU got.
}
\references{
Charnes, A., Clark, C. T., Cooper, W. W. and Golany, B. (1984). A developmental
study of data envelopment analysis in measuring the efficiency of maintenance
units in the U.S. air forces. \emph{Annals of Operations Research} 2, 95--112.

Charnes, A., Cooper, W. W., Lewin, A. Y. and Seiford, L. M. (1994).
\emph{Data Envelopment Analysis: Theory, Methodology and Application}. Kluwer.
Chapter 3, section 5.
}
\seealso{\code{\link{dea_panel}}, \code{\link{dea_malmquist}},
\code{\link{dea}}}
\examples{
## A four-period panel of 15 DMUs: 15 units a period is thin for p + q = 3,
## which is the situation window analysis exists for.
d <- lapply(1:4, function(t) dea_sim(15, p = 2, q = 1, seed = 100 + t))
pan <- dea_panel(do.call(rbind, lapply(d, `[[`, "x")),
                 do.call(rbind, lapply(d, `[[`, "y")),
                 id = rep(1:15, times = 4), period = rep(1:4, each = 15))

w <- dea_window(pan, width = 3)
w
head(w$table)

## Two widths, which is how the width's arbitrariness is meant to be handled:
## read as stable only what does not move between them.
w2 <- dea_window(pan, width = c(2, 3))
summary(w2)
}
