\name{coo.centpos}
\alias{coo.centpos}

\title{
Calculates the position of the centroid of a list or a matrix of coordinates.
}

\description{
\code{coo.centpos} returns the centroid position of the shape that is the mean of x and y coordinates.
}

\usage{
coo.centpos(coo)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  }

\value{
Returns the \code{(x; y)} coordinates of the centroid.
}

\seealso{
\link{coo.centsize}.
}

\examples{
data(bot)
coo <- bot@coo[[1]]
coo.plot(coo, cent=FALSE) # we do not plot the centroid
cent <- coo.centpos(coo)
points(cent[1], cent[2], pch=3, cex=5, col="red")
}

\keyword{coo Utilities}