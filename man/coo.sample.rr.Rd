\name{coo.sample.rr}
\alias{coo.sample.rr}

\title{
Samples points with "equally spaced" angles.
}

\description{
\code{coo.sample.rr} samples \code{n} points in \code{coo} so that the radii departing from the coo's centroid are equals. This function uses polar coordinates and complex algebra.
}


\usage{
coo.sample.rr(coo, n)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  \item{n}{\code{integer}. The number of points to sample.}
  }

\value{
Returns a list with components:
\itemize{
\item $\code{pixindices} vector of radii indices;
\item $\code{radii} vector of sampled radii lengths;
\item $\code{phase} vector of phases;
\item $\code{coord} coordinates of sampled points for a centered shape;
\item $\code{orig.coord} coordinates of sampled points on the original shape;
}}

\seealso{
\link{coo.sample}.
}

\examples{
data(bot)
coo <- bot@coo[[1]]
coo.plot(coo,  main="24 points with equally spaced radii")
rad <- coo.sample.rr(coo, 24)$orig.coord
cent <- coo.centpos(coo)
segments(cent[1], cent[2], rad[, 1], rad[, 2], col=col.summer(24))
points(coo.sample.rr(coo, 24)$orig.coord, pch=20)
}
\keyword{coo Utilities}