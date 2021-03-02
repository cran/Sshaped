\name{plot.sshaped}
\alias{plot.sshaped}
\title{Plot of a \code{sshaped} object}

\description{This function takes a fitted \code{sshaped} object produced 
     by \code{\link{sshapedreg}} or  \code{\link{cvxreg}},  and plots the observations, the fitted function and the location of the estimated inflection point (where applicable).
}

\usage{
\method{plot}{sshaped}(x, ...)
}

\arguments{
     \item{x}{A fitted \code{sshaped} object produced 
       by \code{\link{sshapedreg}} or  \code{\link{cvxreg}}.}

     \item{...}{Other arguments passed to \code{plot}.}
}

\details{
     A plot that shows the observations, the fitted function and the location of the estimated inflection point (where applicable, i.e. when \code{object$shape} equals "\code{sshaped}".).
}


\value{No return value. A plot is produced.}


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
