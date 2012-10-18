\name{coo.smooth}
\alias{coo.smooth}

\title{
Smoothes list and matrices of coordinates.
}

\description{
\code{coo.smooth} performs \code{n} smoothing iteration on \code{coo}.
}

\usage{
coo.smooth(coo, n)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  \item{n}{\code{integer}. The number of iterations to perform.}
  }


\value{
Returns a matrix of \code{(x; y)}coordinates.
}

\examples{
data(bot)
coo <- coo.sample(bot@coo[[1]], 100)
coo.plot(coo, points=TRUE, main="Take it not too smooth")
s <- seq(10, 1000, length=10)
cols <- col.summer(10)
for (i in seq(along=s)) {
  coo.draw(coo.smooth(coo, s[i]), col=NA, border=cols[i])
  }
}
\keyword{coo Utilities}