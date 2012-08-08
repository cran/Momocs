\name{Utilities: edm}
\alias{edm}

\title{
Pairwaise euclidean distance between two matrices/lists of coordinates.
}
\description{
This function calculates point-to-point euclidean distance between two matrices containing the same number of \eqn{(x; y)} coordinates. It can be used to calculate deviations between shapes that has previously been aligned.
}
\usage{
edm(m1, m2)
}
\arguments{
  \item{m1}{The first \code{list} or \code{matrix} of coordinates.}
  \item{m2}{The second \code{list} or \code{matrix} of coordinates.}
}
\details{
If one wish to align two (or more shapes) Procrustes surimposition may provide a better solution.
}
\value{
A \code{vector} of euclidean distances between pairwise coordinates in the two matrices.}

\seealso{
\link{dist}
}

\examples{
x <- matrix(1:10, nc=2)
edm(x, x)
edm(x, x+1)
}
\keyword{Utilities}