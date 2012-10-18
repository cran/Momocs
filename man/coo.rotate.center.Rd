\name{coo.rotate.center}
\alias{coo.rotate.center}

\title{
Rotates a list or matrix of coordinates with any center of symmetry.
}

\description{
\code{coo.rotate.center} rotates (counter-clockwise) \code{coo} with a \code{theta} angle (in radians) and with the center of symmetry specified by \code{center} .
}

\usage{
coo.rotate.center(coo, theta, center=c(0, 0))
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  \item{theta}{\code{numeric}. The angle to rotate the shape.}
  \item{center}{ The (x; y) coordinates of the center of rotation.}
  }

\details{
This function simply centers the shape, rotate it, and then "uncenter" it.
}

\value{
A \code{list} with \code{x; y} components or a a matrix of \code{(x; y)}coordinates.
}

\seealso{
\link{coo.rotate}.
}
\examples{
data(bot)
coo <- coo.scale(coo.center(bot@coo[[1]]))
coo.plot(coo,  main="The bottle's revolution", ylim=c(-10, 50))
r <- seq(pi/6, 2*pi, pi/6)
cols <- col.summer(12)
for (i in seq(along=r)) {
  coo.draw(coo.rotate.center(coo, r[i], center=c(0, 20)), border=cols[i], col=NA)  
  }
abline(v=0, h=0, col="grey60", lty=2)
points(0, 20, pch=3, col="red")

}
\keyword{Utilities}