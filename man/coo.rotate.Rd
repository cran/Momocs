\name{coo.rotate}
\alias{coo.rotate}

\title{
Rotates a list or matrix of coordinates.
}

\description{
\code{coo.rotate} rotates (counter-clockwise) \code{coo} with a \code{theta} angle (in radians) .
}

\usage{
coo.rotate(coo, theta)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  \item{theta}{\code{numeric}. The angle to rotate the shape.}
  }

\details{
For those not familiar with linear mapping the providence is there: \url{http://en.wikipedia.org/wiki/Linear_mapping}}

\value{
A \code{list} with \code{x; y} components or a a matrix of \code{(x; y)}coordinates.
}

\seealso{
\link{coo.rotate.center} if one wants to rotate shapes with another center of rotation than the origin.
}

\examples{
data(bot)
coo <- bot@coo[[1]]
coo.plot(coo,  main="The bottle revolution", xlim=c(-1.5e3, 1.5e3))
r <- seq(pi/6, 2*pi, pi/6)
cols <- col.summer(12)
for (i in seq(along=r)) {
  coo.draw(coo.rotate(coo, r[i]), border=cols[i], col=NA)  
}
}
\keyword{Utilities}