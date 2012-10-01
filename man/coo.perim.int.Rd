\name{coo.perim.pts}
\alias{coo.perim.pts}

\title{
Calculates the euclidean distance between points of a list or matrix of coordinates.
}

\description{
\code{coo.perim.pts} returns the euclidean distances between all successive (x; y) coordinates.
}

\usage{
coo.perim.pts(coo)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  }

\value{
Returns a \code{vector} of euclidean distances.
}

\seealso{
\link{ed}, \link{coo.perim}.
}

\examples{
data(bot)
coo.perim.pts(bot@coo[[1]])
}

\keyword{coo Utilities}