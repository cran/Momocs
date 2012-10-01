\name{coo.perim}
\alias{coo.perim}

\title{
Calculates the perimeter of a list or matrix of coordinates.
}

\description{
\code{coo.perim} returns the sum of the euclidean distances between all successive (x; y) coordinates.
}

\usage{
coo.perim(coo)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  }

\value{
Returns the perimeter length.
}

\examples{
data(bot)
coo.perim(bot@coo[[1]])
}

\keyword{coo Utilities}