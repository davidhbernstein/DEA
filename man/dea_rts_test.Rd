\name{dea_rts_test}
\alias{dea_rts_test}
\alias{print.dea_rts_test}
\alias{summary.dea_rts_test}
\title{Nonparametric test of returns to scale}
\description{
Test whether the technology as a whole may be taken to be constant-returns (or
non-increasing-returns) against the variable-returns alternative, by the
bootstrap of Simar and Wilson (2002). Where \code{\link{dea_rts}} classifies
each DMU with no standard error attached, this attaches a p-value to a statement
about the sample.
}
\usage{
dea_rts_test(x, y, data = NULL, h0 = c("crs", "nirs", "ndrs"),
             orientation = c("in", "out"), B = 2000, alpha = 0.05,
             bw = c("silverman", "nrd0", "ucv", "sj"), seed = NULL,
             ncores = 1L, progress = interactive(), scaling = TRUE)

\method{print}{dea_rts_test}(x, \dots)
\method{summary}{dea_rts_test}(object, \dots)
}
\arguments{
\item{x, y, data}{Inputs and outputs, in any of the forms accepted by
\code{\link{dea}}. For the \code{print} method, an object of class
\code{"dea_rts_test"}.}
\item{h0}{The technology under the null: \code{"crs"}, \code{"nirs"} or
\code{"ndrs"}. The alternative is always \code{"vrs"}. Benchmarking's spellings
\code{"drs"} and \code{"irs"} are accepted, as elsewhere in this package.}
\item{orientation}{\code{"in"} or \code{"out"}, as in \code{\link{dea}}. Read
\strong{What this test cannot be asked} below before choosing \code{"in"}.}
\item{B}{Bootstrap replications. The p-value's resolution is
\eqn{1/(B+1)}, so a small \code{B} is not merely imprecise, it is a floor on
how small a p-value can be reported; \code{B} below 200 warns.}
\item{alpha}{Nominal level, used for \code{reject} and for the simultaneous
level applied to the per-DMU tests.}
\item{bw}{Bandwidth for the kernel that smooths the resampled efficiencies, as
in \code{\link{dea_boot}}, or a single positive number.}
\item{seed}{Integer seed. Set and restored, so the caller's random stream is
left as it was found.}
\item{ncores}{Replications to run in parallel, as in \code{\link{dea_boot}}.}
\item{progress}{Print a replication counter.}
\item{scaling}{As in \code{\link{dea}}. The statistic is units invariant either
way; this only conditions the programs.}
\item{object}{An object of class \code{"dea_rts_test"}.}
\item{\dots}{Ignored.}
}
\value{
An object of class \code{"dea_rts_test"}, a list whose components include
\item{statistic}{The two observed statistics, \code{ratio_of_means} and
\code{mean_of_ratios}, both in \eqn{(0, 1]}.}
\item{p_value, reject}{One-sided lower-tail p-values for each, and whether each
rejects at \code{alpha}.}
\item{table}{One row per DMU: the two scores, the per-DMU \code{scale_eff}, its
p-value, whether the DMU is scale-efficient at the simultaneous level, and
\code{reachable} (see below).}
\item{boot}{The \code{B} by 2 matrix of statistics from the null distribution.}
\item{boot_ratio}{The \code{n} by \code{B} matrix of per-DMU ratios.}
\item{alpha_individual}{The Sidak level \eqn{1 - (1-\alpha)^{1/n}}.}
}
\details{
Under the null the restricted technology is the true one, so the restricted and
the unrestricted estimator are consistent for the same distance and their ratio
tends to 1. Under the alternative the restricted technology is strictly larger
than the truth, the restricted scores are strictly worse, and the ratio tends to
something below 1. So \strong{the test is one-sided and rejects for small
values}. Two aggregations of the per-DMU ratio are reported, both from Simar and
Wilson (2002, section 4): the ratio of the two means, and the mean of the
per-DMU ratios. They are different numbers and need not move together --- the
first weights a DMU by its own level and the second does not.

\strong{The scale the statistic is written on.} Scores enter on this package's
common \eqn{(0, 1]} scale --- \eqn{\theta} input-oriented, \eqn{1/\phi}
output-oriented, which is what \code{efficiency(fit, "score")} returns. That
makes the statistic lie in \eqn{(0, 1]} and the test left-tailed in both
orientations, and makes \code{mean_of_ratios} exactly the mean of
\code{\link{dea_rts}}'s \code{scale_eff} column.

It also settles something the two packages on CRAN with this test settle
differently. Output-oriented, \code{npsf::nptestrts} aggregates on the
\eqn{\phi} scale while this package and \code{rDEA::rts.test} aggregate on the
\eqn{1/\phi} scale. On the 40-DMU design used for the recorded reference values:
\tabular{lrrr}{
                  \tab this package / rDEA \tab npsf \tab 1 / npsf \cr
ratio of means    \tab 0.858356 \tab 1.221348 \tab 0.818767 \cr
mean of ratios    \tab 0.862094 \tab 1.244737 \tab 0.803382
}
The third column is not the first. A mean does not commute with a reciprocal, so
neither aggregate transfers between the scales --- and that is true of the mean
of ratios as much as of the ratio of means, even though the \emph{per-DMU}
ratios are exact reciprocals of each other. Both are legitimate statistics for
the same hypothesis. Input-oriented the question does not arise and all three
packages agree to 1e-13.

\strong{The null distribution} comes from the smoothed homogeneous bootstrap of
\code{\link{dea_boot}}, driven by the \emph{restricted} scores rather than the
variable-returns ones: resampling and smoothing \eqn{\{s, 2-s\}} and moving each
DMU so that its distance to the estimated restricted frontier is the resampled
value gives a pseudo-sample whose own frontier satisfies the null by
construction. Both technologies are then fitted to that pseudo-sample and the
statistic recomputed. Two details of that are choices rather than consequences:

\itemize{
\item the \emph{technology} is the pseudo-sample and the points \emph{evaluated}
against it are the original observations, as in \code{\link{dea_boot}} and as in
both reference packages. One thing this buys: in each orientation the untouched
half of the data guarantees feasibility, so no replication can be infeasible.
\item one pseudo-sample per \emph{replication}, shared by all \eqn{n} evaluated
DMUs. The statistic is an average over DMUs, so averaging instead over \eqn{n}
independent pseudo-samples would shrink its spread and leave the null
distribution too narrow. \code{npsf::nptestrts} draws a fresh pseudo-sample
inside its DMU loop; \code{rDEA::rts.test} draws one per replication, as here.
}

The p-value is \eqn{(\mbox{count}+1)/(B+1)} rather than \eqn{\mbox{count}/B},
the standard Monte Carlo correction for the observed statistic being itself a
draw under the null. It is therefore never exactly zero, which is the honest
report for a quantity whose resolution is \eqn{1/(B+1)}.

\strong{Per-DMU tests come free} --- each replication already yields the whole
vector of pseudo-ratios --- and are \eqn{n} simultaneous tests, so
\code{scale_efficient} is read at the Sidak level \eqn{1-(1-\alpha)^{1/n}}
rather than at \code{alpha}. With \code{h0 = "nirs"} the per-DMU test applied to
the scale-inefficient DMUs is what distinguishes decreasing from increasing
returns.
}
\section{What this test cannot be asked}{
\strong{On a bounded input support the input-oriented version of this test is not
valid, and no bootstrap repairs it.} The variable-returns program cannot put
\eqn{\sum_i \lambda_i x_{ij}} below \eqn{\min_i x_{ij}}, because a convex
combination of the sample's \eqn{j}th input is at least its smallest value. So a
DMU whose restricted projection needs \eqn{\theta x_{oj} < \min_i x_{ij}} in any
input \eqn{j} has \eqn{\theta_{VRS} > \theta_{CRS}} for a reason that has
nothing to do with returns to scale, and no sample size removes it. This is the
same identification problem that \code{\link{dea_sim}}'s \code{x_range} argument
exists for, reaching a test instead of a convergence rate.

Measured on \code{dea_sim(returns = 1)}, where the technology \emph{is} a cone
and the null is exactly true, so a consistent statistic must approach 1:
\tabular{lrrrr}{
n                     \tab 80 \tab 200 \tab 800 \tab 3200 \cr
input, ratio of means \tab 0.887 \tab 0.903 \tab 0.894 \tab 0.898 \cr
output, ratio of means \tab 0.970 \tab 0.985 \tab 0.991 \tab 0.997
}
The input-oriented statistic does not converge, so the input-oriented test
rejects a true null with probability tending to one. Splitting the same designs
on exactly the condition above separates a consistent statistic from a constant:
at \eqn{p = 1} the reachable DMUs give 0.982, 0.991, 0.9991, 0.9995 while the
34\% that are not reachable give 0.755, 0.751, 0.771, 0.776.

So \code{reachable} is computed per DMU --- the output-oriented mirror is
\eqn{\phi y_{oj} > \max_i y_{ij}} --- reported by \code{print}, and warned about
past a quarter of the sample. The benign case is not "none flagged" but "a share
that \emph{shrinks} with n": output-oriented on the same cone that share runs
0.188, 0.020, 0.050 where input-oriented it holds at 0.362, 0.340, 0.338.
Neither \code{npsf::nptestrts} nor \code{rDEA::rts.test} reports anything of the
kind, and both will run the input-oriented test on such a design without
comment.

Separately, Simar and Wilson (2011) show that the full-sample bootstrap used
here is not consistent for tests of this kind and propose an \eqn{m} out of
\eqn{n} subsampling alternative. That is not implemented; this function is the
2002 procedure, which is what both packages on CRAN provide.
}
\references{
Simar, L. and Wilson, P. W. (1998). Sensitivity analysis of efficiency scores:
how to bootstrap in nonparametric frontier models.
\emph{Management Science} 44, 49--61.

Simar, L. and Wilson, P. W. (2002). Non-parametric tests of returns to scale.
\emph{European Journal of Operational Research} 139, 115--132.

Simar, L. and Wilson, P. W. (2011). Inference by the \eqn{m} out of \eqn{n}
bootstrap in nonparametric frontier models.
\emph{Journal of Productivity Analysis} 36, 33--53.
}
\seealso{\code{\link{dea_rts}} for the per-DMU classification and scale
efficiency, \code{\link{dea_boot}} for the resampling scheme this builds on.}
\examples{
sim <- dea_sim(80, p = 2, q = 1, returns = 0.6, seed = 2)
## Output oriented: see 'What this test cannot be asked' above.
r <- dea_rts_test(sim$x, sim$y, orientation = "out", B = 299, seed = 1,
                  progress = FALSE)
r
## The two aggregations are different tests, and on this design they disagree.
head(r$table)
}
