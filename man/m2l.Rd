\name{m2l}
\alias{m2l}

\title{
Convert a matrix of coordinates to a list of coordinates.
}

\description{
\code{m2l} converts a matrix of \code{(x; y)}coordinates to a list with \code{x; y} components.
}

\usage{
m2l(m)
}

\arguments{
  \item{m}{A 2-columns \code{matrix} containing x and y coordinates.}
  }

\value{
Returns a \code{list} with \code{x; y} components.
}

\seealso{
\link{l2m}.
}

\examples{
	\dontrun{
data(gorf.dat)
m2l(gorf.dat[,,1])
}
}
\keyword{coo Utilities}