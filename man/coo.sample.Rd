\name{coo.sample}
\alias{coo.sample}

\title{
Samples points along the curvilinear abscissa.
}

\description{
\code{coo.sample} samples \code{n} points in \code{coo} along the curvilinear abscissa.
}

\usage{
coo.sample(coo, n)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  \item{n}{\code{integer}. The number of points to sample.}
  }

\details{
So far points are just sampled along the \code{coo} provided) and not along the "true" curvilinear abscissa, but in most cases, differences should be very small.
}

\value{
Returns a matrix of \code{(x; y)}coordinates.
}

\seealso{
\link{coo.sample.rr}.
}

\examples{
data(bot)
coo <- bot@coo[[1]]
coo.plot(coo,  main="24 points sampled along the outline")
points(coo.sample(coo, 24), col="red", pch=20)
}
\keyword{coo Utilities}