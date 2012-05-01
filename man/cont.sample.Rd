\name{cont.sample}
\alias{cont.sample}
\title{
Samples points along a list of \eqn{(x; y)} coordinates.}
\description{
Given a list of \eqn{(x; y)} coordinates, samples equidistant points along this outline.
}
\usage{
cont.sample(coo, n)
}
\arguments{
  \item{coo}{the \code{Coo} object}
  \item{n}{\code{integer} indicating the number of points to sample}
}
\value{
Returns a, usually shortened, \code{list} of {eqn(x; y)} coordinates.
}
\examples{
fake <- list(x=1:100, y=101:200)
cont.sample(fake, 10)
}
\keyword{Misc functions}
