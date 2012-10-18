\name{coo.scale}
\alias{coo.scale}


\title{
Scales a list or matrix of coordinates.
}

\description{
\code{coo.scale} resizes \code{coo} so that it can be included in a square of \code{scale} side.
}

\usage{
coo.scale(coo, scale)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  \item{scale}{The scale to use.}
  }

\value{
Returns a matrix of \code{(x; y)}coordinates.
}

\seealso{
\link{coo.template}.
}

\keyword{coo Utilities}