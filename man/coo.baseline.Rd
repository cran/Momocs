\name{coo.baseline}
\alias{coo.baseline}

\title{
Re-register a list or matrix of coordinates.
}

\description{
\code{coo.baseline} registers \code{coo} on a (new) baseline, \emph{e.g.} Bookstein's coordinates.
}

\usage{
coo.baseline(coo, ldk1 = 1, ldk2 = 2, t1=c(-0.5, 0), t2=c(0.5, 0))
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  \item{ldk1}{Row index of the first reference landmark.}
  \item{ldk2}{Row index of the second reference landmark.}
  \item{t1}{(x; y) coordinates of the first target point.}
  \item{t2}{(x; y) coordinates of the second target point.}
  }

\details{
\code{coo.baseline} returns by default Bookstein's coordinates.
}

\value{
Returns a matrix of \code{(x; y)}coordinates.
}

\examples{
def.par <- par(no.readonly = TRUE)
layout(matrix(1:2, 1, 2))
data(bot)
coo <- bot@coo[[1]]
coo <- coo.sample(coo, 12)
coo.plot(coo, main="Some landmarks", points=TRUE)
coo.plot(coo.baseline(coo), border="red", col=NA,
    points=TRUE, main="Reregistered using Bookstein's coordiantes")
abline(v=c(-0.5, 0.5), col="grey60", lty=2)
box()
par(def.par)
}

\keyword{coo Utilities}