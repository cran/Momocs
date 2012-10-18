\name{coo.trans}
\alias{coo.trans}

\title{
Transles a list or a matrix of coordinates.
}

\description{
\code{coo.trans} translates \code{coo} by \code{x} and {y} on the two dimensions, respectively.
}

\usage{
coo.trans(coo, x, y)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  \item{x}{\code{numeric}. The x-axis translation.}
  \item{y}{\code{numeric}. The y-axis translation.}
  }

\value{
Returns a matrix of \code{(x; y)}coordinates.
}

\examples{

data(bot)
coo <- coo.scale(coo.center(bot@coo[[1]]))

tx <- seq(0, 10*pi, length=50)
ty <- sin(tx)*5
cols <- col.summer(50)

coo.plot(xlim=c(0, 10*pi), main="The bottle's wave")
lines(tx, ty, col="grey60", lty=2)
for (i in seq(along=tx)) {
  coo.draw(coo.trans(coo, tx[i], ty[i]), col=NA, border=cols[i])
}


}
\keyword{Utilities}