\name{coo.centsize}
\alias{coo.centsize}

\title{
Calculates the centroid size of a list or matrix of coordinates.
}

\description{
\code{coo.centsize} returns the centroid size of the shape, \emph{i.e.} the square root of the sum of squared distances from each point to the centroid of the shape.
}

\usage{
coo.centsize(coo)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  }

\value{
Returns the centroid size.
}

\seealso{
\link{coo.centpos}.
}

\examples{
data(bot)
coo <- bot@coo[[1]]
coo.centsize(coo)
coo.centsize(coo.scale(coo)) # by default coo.scale scales to centroid size.
}

\keyword{coo Utilities}