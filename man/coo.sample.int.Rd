\name{coo.sample.int}
\alias{coo.sample.int}

\title{
Given an outline, interpolates points along the curvilinear abscissa.
}

\description{
\code{coo.sample.int} interpolates \code{n} points in \code{coo} along the curvilinear abscissa.
}

\usage{
coo.sample.int(coo, n)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  \item{n}{\code{integer}. The number of points to sample.}
  }

\details{
This function circumvent the problem exposed in the Details section of \link{coo.sample}, \emph{i.e.} points are here samples along the "true" curvilear abscissa and not based on \link{seq}.
}

\value{
Returns a matrix of \code{(x; y)}coordinates.
}

\seealso{
\link{coo.sample}, \link{coo.sample.rr}.
}

\examples{
data(bot)
coo <- coo.sample(bot@coo[[1]], 12)
coo.plot(coo, points=TRUE)
coo.draw(coo.sample.int(coo, 20), col="#B2222255", border=NA)
coo.draw(coo.sample(coo, 20),     col="#1874CD55", border=NA)
}
\keyword{coo Utilities}