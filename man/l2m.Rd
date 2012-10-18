\name{l2m}
\alias{l2m}

\title{
Converts a list of coordinates to a matrix.}

\description{
\code{l2m} converts a \code{list} with x and y components to a 2-col \code{matrix} of coordinates.
}


\usage{
l2m(l)
}

\arguments{
  \item{l}{A \code{list} with x and y coordinates as components.}
  }

\value{
Returns a matrix of \code{(x; y)}coordinates.
}

\seealso{
\link{m2l}.
}

\examples{
l <- list(x=1:5, y=5:1)
l2m(l)

}
\keyword{coo Utilities}