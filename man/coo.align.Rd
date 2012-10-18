\name{coo.align}
\alias{coo.align}

\title{
Aligns a list or matrix of coordinates.
}

\description{
\code{coo.align} aligns \code{coo} using its best fitting ellipse.
}

\usage{
coo.align(coo)
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
coo.plot(coo.align(coo), border="red", col=NA, main="An aligned shape")
par(def.par)
}
\keyword{coo Utilities}