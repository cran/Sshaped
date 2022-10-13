\name{cvxreg}
\alias{cvxreg}
\title{Estimation of an increasing and convex function}
\description{
This function computes the least squares increasing and convex regression estimator by a sequential mixed primal-dual bases algorithm.
}
\usage{
cvxreg(x,y)
}
\arguments{
  \item{x}{a numeric vector that contains all the design points. NB. for the current version, we require all values to be distinct (but not necessarily in ascending order).}
  \item{y}{a numeric vector that contains the values of the response with respect to the design points.}
}

\value{
An object of class \code{sshaped}, which contains the following fields: 
  \item{x}{covariates copied from input}
  \item{y}{response copied from input}
  \item{fitted}{the fitted values of the regression function with respect to the design points.}
  \item{rss}{the value of the minimised residual sum of squares of the fit}
  \item{inflection}{the location of the inflection point, which equals \code{max(x)} here}
  \item{shape}{the shape enforced in the fit, here equals "\code{convex}"
}
}
\details{
For more details, see \code{\link{Sshaped}}. 
}

\examples{
x<-seq(0,5,0.01)
y<-x^2 + rnorm(length(x))
output<-cvxreg(x,y)
plot(output)
xnew=rexp(5); predict(output,xnew)
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

