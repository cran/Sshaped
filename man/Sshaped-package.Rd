\name{Sshaped-package}
\alias{Sshaped-package}
\alias{Sshaped}
\docType{package}
\title{
\packageTitle{Sshaped}}
\description{
Estimation of an S-shaped function and its corresponding inflection point via a least squares approach. A sequential mixed primal-dual bases algorithm is implemented for the fast computation of the estimator. The same algorithm can also be used to solve other shape-restricted regression problems, such as convex regression. See Fraser and Massam (1989) and Feng et al. (2022).
}
\details{
Consider the nonparametric estimation of an S-shaped regression function. The least squares
estimator provides a very natural, tuning-free approach, but results in a non-convex optimisation problem, since the inflection point is unknown. Nevertheless, the estimator may  be regarded as a projection onto a finite union of convex cones, which allows us to propose a mixed primal-dual bases algorithm for its efficient, sequential computation. 

In the current version of the package, we use this algorithm to implement the least squares regression estimator under the following shape-restrictions: S-shaped functions, i.e. increasing convex to the left of the inflection point and increasing concave to the right of the inflection point; and increasing and convex functions (as a by-product of the former). The corresponding plot and predict methods are also included. In the future, we plan to also include the estimation of additive S-shaped functions, where the covariates are multivariate for the regression.
}

\references{
Fraser, D. A. S.  and Massam, H. (1989). A mixed primal-dual bases algorithm for regression under inequality constraints. Application to concave regression.  Scandinavian Journal of Statistics, Volume 16, Pages 65-74.

Feng, O. Y., Chen, Y., Han, Q., Carroll, R. J. and Samworth, R. J. (2022). Nonparametric, tuning-free estimation of S-shaped functions. Journal of the Royal Statistical Society: Series B (Statistical Methodology), Volume 84, Issue 4, Pages 1324-1352. <doi:10.1111/rssb.12481>
}

\author{
\itemize{
\item \href{mailto:oyf20@cam.ac.uk}{Oliver Y. Feng}
\item \href{mailto:y.chen101@lse.ac.uk}{Yining Chen}
\item \href{mailto:qh85@stat.rutgers.edu}{Qiyang Han}
\item \href{mailto:carroll@stat.tamu.edu}{Raymond J. Carroll}
\item \href{mailto:r.samworth@statslab.cam.ac.uk}{Richard J. Samworth}
}
}

\examples{
# Generate data
set.seed(1)
x <- seq(-1,1,0.005)
y <- sin(x*pi/2) + rnorm(length(x))

# Fit S-shape
output <- sshapedreg(x,y)

# Plot
plot(output)

# prediction at new design points
xnew=rnorm(5)
predict(output,xnew)
}