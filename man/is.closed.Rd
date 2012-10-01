\name{is.closed}
\alias{is.closed}

\title{
Tests if a list or matrix of coordinates is closed.
}

\description{
\code{is.closed} tests if the last coordinate of \code{coo} provided as a \code{list} or a \code{matrix} is closed, \emph{i.e.} if the last coordinate is identical as the first.
}

\usage{
is.closed(coo)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  }

\value{
A \code{logical}, either \code{TRUE} or \code{FALSE}.
}

\seealso{
\link{coo.close} and \link{coo.unclose}.
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