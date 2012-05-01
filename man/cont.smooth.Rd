\name{cont.smooth}
\alias{cont.smooth}
\title{
Smoothes a list or a matrix of \eqn{(x; y)} coordinates.
}
\description{
Applies a simple algorithm to smooth outlines, particularly to remove, if needed, digitalization artefacts.
}
\usage{
cont.smooth(M, n)
}
\arguments{
  \item{M}{\code{list} or a matrix of \eqn{(x; y)} coordinates to smooth}
  \item{n}{\code{integer} indicating how many smoothing iterations to perform}
}
\details{
The algorithm used is simplistic: the new \eqn{(x; y)_n} coordinates are calculated as:
\deqn{\frac{1}{4}\times(x; y)_{n-1} + \frac{1}{2}\times(x; y)_{n} + \frac{1}{4}\times(x; y)_{n+1}}
}
\value{
A \code{list} of smoothed coordinates.
}
\examples{
fake <- list(x=1:10, y=20:11)
cont.smooth(fake, 1)
cont.smooth(fake, 10)
}
\keyword{Misc functions}