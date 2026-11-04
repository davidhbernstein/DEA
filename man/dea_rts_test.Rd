\name{dea_rts_test}
\alias{dea_rts_test}
\alias{print.dea_rts_test}
\alias{summary.dea_rts_test}
\title{Nonparametric test of returns to scale}
\description{
Test whether the technology as a whole may be taken to be constant-returns (or
non-increasing-returns) against the variable-returns alternative, by the
bootstrap of Simar and Wilson (2002) or the \eqn{m} out of \eqn{n}
subsampling of Simar and Wilson (2011). Where \code{\link{dea_rts}} classifies
each DMU with no standard error attached, this attaches a p-value to a statement
about the sample.
}
\usage{
dea_rts_test(x, y, data = NULL, h0 = c("crs", "nirs", "ndrs"),
             orientation = c("in", "out"),
             method = c("bootstrap", "subsample"), B = 2000, alpha = 0.05,
             bw = c("silverman", "nrd0", "ucv", "sj"),
             m = NULL, m_grid = NULL, B_select = 200L, window = 3L,
             seed = NULL, ncores = 1L, progress = interactive(),
             scaling = TRUE)

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
\item{method}{How the null distribution is obtained: \code{"bootstrap"}, the
smoothed bootstrap of Simar and Wilson (2002), or \code{"subsample"}, the
\eqn{m} out of \eqn{n} subsampling of Simar and Wilson (2011). See
\strong{Subsampling} below.}
\item{B}{Bootstrap replications, or subsamples under
\code{method = "subsample"}. The p-value's resolution is
\eqn{1/(B+1)}, so a small \code{B} is not merely imprecise, it is a floor on
how small a p-value can be reported; \code{B} below 200 warns.}
\item{alpha}{Nominal level, used for \code{reject} and for the simultaneous
level applied to the per-DMU tests.}
\item{bw}{Bandwidth for the kernel that smooths the resampled efficiencies, as
in \code{\link{dea_boot}}, or a single positive number. Bootstrap only.}
\item{m}{Subsample size, \eqn{3 \le m < n}. \code{NULL} chooses it by
minimum volatility over \code{m_grid}. Subsampling only.}
\item{m_grid, B_select, window}{The candidate sizes, the subsamples drawn at
each, and the half-width of the volatility window, as in
\code{\link{dea_subsample}}. Used only when \code{m} is \code{NULL}.}
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
\item{method}{Which null distribution was used.}
\item{kappa, m, tau, boot_tau}{Subsampling only: the rate, the subsample size,
the observed \eqn{n^\kappa(1-S_n)} and the \code{B} by 2 matrix of
\eqn{m^\kappa(1-S^*_m)}.}
\item{m_grid, critical_path, volatility}{Subsampling with \code{m} chosen:
the candidates, the \eqn{1-\alpha} quantile of \code{boot_tau} at each, and
the volatility the choice minimised.}
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
\section{How well it works, measured}{
\code{dea_sim(returns = 1)} draws a cone, so H0 is exactly true and the size of
the test can be measured rather than assumed. Rejection rate at a nominal 5\%,
output oriented, \code{B = 400}, 200 replications (100 at \eqn{n = 150}):
\tabular{lrrrr}{
                   \tab n=60, r=1 \tab n=150, r=1 \tab n=60, r=0.8 \tab n=60, r=0.6 \cr
ratio of means     \tab 0.025 \tab 0.000 \tab 0.020 \tab 0.155 \cr
mean of ratios     \tab 0.050 \tab 0.010 \tab 0.185 \tab 0.780 \cr
\code{rDEA}        \tab 0.025 \tab 0.000 \tab 0.035 \tab 0.175 \cr
\code{npsf}        \tab 0.315 \tab 0.240 \tab 0.680 \tab 0.990
}
Three things to take from it.

\strong{On this design \code{mean_of_ratios} is the more powerful statistic},
0.780 against 0.155 at \eqn{r = 0.6} with a paired \eqn{t} of \eqn{-18.2}, at a
size no worse. \strong{But in more dimensions its bootstrap version rejects a true
null far too often} --- see \strong{Size in more dimensions} below --- so read
this as a statement about \eqn{p = 2}, \eqn{q = 1}, not as a general
recommendation.

\strong{The size is conservative, not nominal, and the power is low.} At
\eqn{r = 0.8} the ratio of means rejects 2.0\% of the time, below its own size; at
\eqn{r = 0.6} --- an elasticity of scale of 0.6, which is not subtle --- 15.5\%.
\code{rDEA::rts.test} sits at 3.5\% and 17.5\% on the same datasets, so this is
the 2002 procedure and not this implementation of it. Do not read a failure to
reject as evidence for constant returns at these sample sizes.

\strong{\code{npsf::nptestrts}'s columns are not comparable.} It rejects a true
null 31.5\% of the time, because it draws a fresh pseudo-sample inside its DMU
loop, so its apparently higher power is the same miscalibration seen from the
other side.
}
\section{Subsampling}{
Simar and Wilson (2011) show that the full-sample bootstrap is not consistent
for tests of this kind, and propose \eqn{m} out of \eqn{n} subsampling instead.
\code{method = "subsample"} implements it. Draw \eqn{m < n} DMUs without
replacement, fit \emph{both} technologies to the subsample alone, and compute the
statistic \eqn{S^*_m} from it. Then compare
\deqn{\tau_n = n^\kappa (1 - S_n) \quad\mbox{with}\quad \tau^*_m = m^\kappa (1 - S^*_m),}
rejecting for large \eqn{\tau_n}. Three things about that are choices:
\itemize{
\item \strong{No recentring.} \code{\link{dea_subsample}} centres on the
full-sample estimate because it is approximating the distribution of an error
around an unknown value. A test knows the null value: under H0 the statistic
tends to 1. Under H1, \eqn{S_n} and \eqn{S^*_m} tend to the same \eqn{c < 1},
so \eqn{\tau_n} outgrows \eqn{\tau^*_m}. That growth is the power, and centring
on \eqn{S_n} would subtract it.
\item \strong{\eqn{\kappa} is the variable-returns rate}, \eqn{2/(p+q+1)},
whatever \code{h0} is. The statistic contrasts two estimators and converges no
faster than the slower one; under \code{h0 = "crs"} that is not the
restricted one.
\item \strong{\eqn{m} is chosen by minimum volatility of the critical value}, the
\eqn{1-\alpha} quantile of \eqn{\tau^*_m}, using the same rule and the same
full-window restriction as \code{\link{dea_subsample}}. The path is returned in
\code{critical_path}, so the choice can be inspected rather than trusted.
}
There are \strong{no per-DMU tests} under subsampling. A given DMU appears in only
\eqn{m/n} of the subsamples, so its own statistic would rest on a fraction of
\code{B}. Those columns come back \code{NA} rather than being computed from
whichever draws happened to include the DMU.

\strong{What it buys, measured.} Same design as above, both methods on the same
datasets in every replication, 200 replications, \code{B = 400}:
\tabular{lrrrrrr}{
                   \tab n=60, r=1 \tab r=0.8 \tab r=0.6 \tab n=200, r=1 \tab r=0.8 \tab r=0.6 \cr
bootstrap, ratio of means \tab 0.020 \tab 0.020 \tab 0.155 \tab 0.005 \tab 0.440 \tab 0.990 \cr
bootstrap, mean of ratios \tab 0.035 \tab 0.195 \tab 0.755 \tab 0.035 \tab 0.985 \tab 1.000 \cr
subsample, ratio of means \tab 0.025 \tab 0.265 \tab 0.745 \tab 0.000 \tab 0.885 \tab 1.000 \cr
subsample, mean of ratios \tab 0.020 \tab 0.275 \tab 0.760 \tab 0.000 \tab 0.925 \tab 1.000
}
The \eqn{r = 1} columns are size; the rest are power. Subsampling
\strong{repairs the ratio of means}: paired gains of 0.245, 0.590 and 0.445
(\eqn{t} = 8.0, 16.6, 12.6). After the repair the two statistics behave alike,
so the ratio of means' weakness under the bootstrap was in the null
distribution, not in the statistic. \strong{For the mean of ratios it is not an
improvement}: +0.080 at \eqn{n = 60, r = 0.8} (\eqn{t = 2.6}), nothing at
\eqn{r = 0.6}, and \eqn{-0.060} at \eqn{n = 200, r = 0.8} (\eqn{t = -3.3}). It
is also more conservative: no rejections at all under a true cone at
\eqn{n = 200}, against 3.5\% for the bootstrap's mean of ratios.

On this design, then, the bootstrap's mean of ratios is the best column. In more
dimensions it is not --- see \strong{Size in more dimensions} --- and there
subsampling is the method whose size holds. The chosen \eqn{m} sits low in its grid (median 18--29 at \eqn{n = 200}), the
same tendency \code{\link{dea_subsample}} shows.

\strong{Why it is conservative: the rate, not \eqn{m}.} Fixing \eqn{m} across the
whole grid at \eqn{n = 200} (same datasets, 200 replications) gives a size of
0.010 at \eqn{m = 18} and 0.000 from \eqn{m = 29} to 117. Power against
\eqn{r = 0.8} falls steadily as \eqn{m} grows: 0.955, 0.945, 0.895, 0.845,
0.695. So the small \eqn{m} the rule picks is the best choice on both counts. What
makes the test conservative is \eqn{\kappa}. On \code{dea_sim(returns = 1)} at
\eqn{n} from 50 to 1600, \eqn{1 - S_n} falls at \eqn{n^{-0.594}} (95\% interval
\eqn{-0.621} to \eqn{-0.566}), faster than the \eqn{n^{-0.5}} assumed here. Had
the test used 0.594, its size at \eqn{m = 29} would have been about 0.045 and its
power about 0.990. That is a recomputation from the same draws, rejecting when
\eqn{\tau_n} exceeds the 95\% quantile, which differs from the p-value rule by
up to 0.01. It is \strong{not} used, because it was measured on one design
(\eqn{p = 2}, \eqn{q = 1}, Cobb-Douglas, output oriented). Estimating the rate
from the data instead would remove the power, since under H1 the statistic
stops shrinking and the estimated rate goes to zero.

Nor is the gap a constant that could be corrected for. The same measurement
over nine designs, 100 replications each, \eqn{n} from 50 to 1600:
\tabular{lrrr}{
design \tab theory \tab measured \tab 95\% interval \cr
\eqn{p = 1, q = 1} \tab -0.667 \tab -0.759 \tab [-0.789, -0.728] \cr
\eqn{p = 2, q = 1} \tab -0.500 \tab -0.594 \tab [-0.613, -0.574] \cr
\eqn{p = 3, q = 1} \tab -0.400 \tab -0.476 \tab [-0.490, -0.460] \cr
\eqn{p = 2, q = 2} \tab -0.400 \tab -0.397 \tab [-0.412, -0.382] \cr
\eqn{p = 3, q = 2} \tab -0.333 \tab -0.303 \tab [-0.315, -0.289] \cr
CES, \eqn{\rho = -10} \tab -0.500 \tab -0.587 \tab [-0.607, -0.567] \cr
2 facets \tab -0.500 \tab -0.617 \tab [-0.635, -0.596] \cr
4 facets, \eqn{p = 3} \tab -0.400 \tab -0.455 \tab [-0.471, -0.439] \cr
half-normal \eqn{u} \tab -0.500 \tab -0.552 \tab [-0.569, -0.535]
}
With one output the statistic falls faster than theory in every design. With
two it matches theory at \eqn{p = 2} and falls slower at \eqn{p = 3}. The
inefficiency distribution matters too. So the theoretical rate is kept: it is
wrong in both directions, and no other single number is less wrong.
}
\section{Size in more dimensions}{
\strong{The bootstrap's \code{mean_of_ratios} over-rejects a true null, badly,
as the dimension grows.} Rejection rate under a true cone, output oriented,
\eqn{n = 200}, \code{B = 400}, 200 replications, nominal 5\%:
\tabular{lrrr}{
                          \tab p=2, q=1 \tab p=2, q=2 \tab p=3, q=2 \cr
bootstrap, ratio of means \tab 0.005 \tab 0.000 \tab 0.005 \cr
bootstrap, mean of ratios \tab 0.020 \tab 0.090 \tab 0.245 \cr
subsample, ratio of means \tab 0.005 \tab 0.050 \tab 0.045 \cr
subsample, mean of ratios \tab 0.005 \tab 0.040 \tab 0.035
}
Monte Carlo standard errors are 0.030 at most. This is the 2002 procedure and
not this implementation of it: \code{rDEA::rts.test} rejects 25.0\% of the time
on the same 200 datasets at \eqn{p = 3, q = 2} (\code{pvalue48}), a paired
difference of \eqn{-0.005} (se 0.015), and the two packages agree on the
observed statistic to 6e-8.

So: \strong{subsampling is the only method whose size held in every design
measured}, at or below nominal, though conservative when \eqn{q = 1}. The
bootstrap's ratio of means never over-rejects but has little power. The
bootstrap's mean of ratios is the most powerful at \eqn{p = 2}, \eqn{q = 1} and
cannot be trusted beyond it. With more than one output, or \eqn{p + q > 3},
read \code{method = "subsample"}.
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

The subsampling version shares this limitation: it changes the null
distribution, not the statistic, and the statistic is what fails to converge.
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
efficiency, \code{\link{dea_boot}} for the resampling scheme this builds on,
\code{\link{dea_subsample}} for the subsampling machinery.}
\examples{
sim <- dea_sim(80, p = 2, q = 1, returns = 0.6, seed = 2)
## Output oriented: see 'What this test cannot be asked' above.
r <- dea_rts_test(sim$x, sim$y, orientation = "out", B = 299, seed = 1,
                  progress = FALSE)
r
## The two aggregations are different tests, and on this design they disagree.
head(r$table)
\donttest{
## The m out of n version: the same statistic, a different null distribution.
dea_rts_test(sim$x, sim$y, orientation = "out", method = "subsample",
             B = 299, seed = 1, progress = FALSE)
}
}
