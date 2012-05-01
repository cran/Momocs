\name{eFa}
\alias{eFa}
\title{Elliptical Fourier Analysis on Coo objects.}
\description{
Once the number of harmonics to calculate and the number of smoothing iterations to perform have been determined,
calculates elliptic Fourier transforms on the list of \eqn{(x; y)} outline coordinates included in \code{Coo}-objects. 
}
\usage{
eFa(coo, nb.h = 32, smooth.it = 0, fromrt = FALSE)
}
\arguments{
  \item{coo}{the \code{Coo}-object}
  \item{nb.h}{\code{integer} indicating how many harmonics to calculate}
  \item{smooth.it}{\code{integer} indicating how many smoothing iterations to perform}
  \item{fromrt}{\code{logical} indicating whether the position of the starting point has to be preserved or not}
}

\examples{
data(bottles.cont)
eFa(bottles.cont@coo[[1]])
}
\keyword{Misc functions}
