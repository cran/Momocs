\name{coo.center}
\alias{coo.center}


\title{
Centers a list or matrix of coordinates.
}

\description{
\code{coo.center} centers the \code{coo}'s centroid on the origin.
}

\usage{
coo.center(coo)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  }

\value{
Returns a matrix of \code{(x; y)}coordinates.
}

\examples{
def.par <- par(no.readonly = TRUE)
layout(matrix(1:2, 1, 2))
data(bot)
coo <- bot@coo[[1]]
coo.plot(coo, main="A shape")
abline(v=0, h=0, col="grey60", lty=2)
coo.plot(coo.center(coo), border="red", col=NA, main="A centered shape")
abline(v=0, h=0, col="grey60", lty=2)
par(def.par)
}
\keyword{coo Utilities}