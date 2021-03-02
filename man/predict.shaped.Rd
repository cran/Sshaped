\name{predict.sshaped}
\alias{predict.sshaped}

\title{
Predict method for a \code{sshaped} fit
}

\description{
This function obtains predictions at new design points from a fitted \code{sshaped} object.
}

\usage{
\method{predict}{sshaped}(object, xnew, ...)
}

\arguments{
  \item{object}{A fitted \code{sshaped} object produced by \code{\link{sshapedreg}} or  \code{\link{cvxreg}}.}

  \item{xnew}{An optional numeric \code{vector} specifying the location at which prediction is required. Its values can lie outside the range of the original design points. This argument can be missing, in which case predictions are made at the same values of the covariates used to compute the object.}
  
  \item{...}{Other arguments passed to \code{predict}.}
}

\value{A numeric vector containing predictions at new data points.}

\author{
\itemize{
\item \href{mailto:oyf20@cam.ac.uk}{Oliver Y. Feng}
\item \href{mailto:y.chen101@lse.ac.uk}{Yining Chen}
\item \href{mailto:qh85@stat.rutgers.edu}{Qiyang Han}
\item \href{mailto:carroll@stat.tamu.edu}{Raymond J. Carroll}
\item \href{mailto:r.samworth@statslab.cam.ac.uk}{Richard J. Samworth}
}
}

\seealso{\code{\link{sshapedreg}}, \code{\link{cvxreg}}}

\examples{
## See examples for the functions sshapedreg and cvxreg.
}
