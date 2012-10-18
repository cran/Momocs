\name{coo.force2close}
\alias{coo.force2close}

\title{
Forces a list or matrix or coordinates to close.
}

\description{
\code{coo.force2close} force one coo to close, i.e. the difference \eqn{(x_n - x_1)/n} is added to each x-coordinates and similarly for y-coordinates.
}

\usage{
coo.force2close(coo)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  }

\value{
Returns a matrix of \code{(x; y)}coordinates.
}


\examples{
data(bot)
coo <- bot@coo[[1]]
coo.plot(coo, ylim=c(-200, 1200), xlim=c(0, 500), main="A half bottle forced to close")
coo.draw(coo[1:1500,], border="blue", col=NA)
coo.draw(coo.force2close(coo[1:1500, ]), border="red", col=NA)
}

\keyword{coo Utilities}
