\name{A.mshape}
\alias{A.mshape}

\title{
Calculates the mean shape of an array of landmarks.
}

\description{
\code{A.mshape} returns the mean shape of an array of landmarks when provided with landmarks as rows, \eqn{(x; y)} coordinates as columns, and individuals as the third dimension of the array.
}

\usage{
A.mshape(A)
}

\arguments{
  \item{A}{An \code{array} of \eqn{(x; y)} coordinates, typically, landmarks.}
  }

\value{
Returns a matrix of \eqn{(x; y)} coordinates.
}

\seealso{
\link{A.points}, \link{A.segments}, \link{A.plot} \link{coo.plot}, and the others coo.utilities.
}
\examples{
data(gorf.dat) # we import gorf.data from shapes
A.plot(gorf.dat, pch=20)
m <- A.mshape(gorf.dat)
points(m, pch=5)
}
\keyword{Utilities}