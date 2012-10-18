\name{vecs.param}
\alias{vecs.param}

\title{
Some vector utilities.
}
\description{
Returns ratio of norms and signed angle between two vectors provided as four numeric.
}
\usage{
vecs.param(r1, i1, r2, i2)
}

\arguments{
  \item{r1}{the "real" part of the first vector, i.e. difference in x-coordinates.}
  \item{i1}{the "imaginary" part of the first vector, i.e. difference in y-coordinates.}
  \item{r2}{the "real" part of the second vector, i.e. difference in x-coordinates.}
  \item{i2}{the "imaginary" part of the second vector, i.e. difference in y-coordinates.}
}

\value{
A list with two components: \code{r.norms} the ratio of (norm of vector 1)/(norm of vector 2) and \code{d.angle} the signed angle 'from' the first 'to' the second vector.
}

\examples{
vecs.param(1, 0, 0, 2)
}
\keyword{Utilities}
