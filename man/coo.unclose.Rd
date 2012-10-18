\name{coo.unclose}
\alias{coo.unclose}

\title{
Uncloses a list or matrix of coordinates.
}

\description{
\code{coo.unclose} tests if \code{coo} is closed, then simply removes the last coordinate.
}

\usage{
coo.unclose(coo)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  }

\value{
Returns a matrix of \code{(x; y)}coordinates.
}

\seealso{
\link{coo.close} and \link{is.closed}
}

\examples{
data(gorf.dat)
coo <- gorf.dat[,,1]
is.closed(coo)
coo.c <- coo.close(coo)
is.closed(coo.c)
coo.cu <- coo.unclose(coo)
is.closed(coo.cu)
}

\keyword{coo Utilities}