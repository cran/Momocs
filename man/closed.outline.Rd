\name{closed.outline}
\alias{closed.outline}
\title{
Closes outlines.
}
\description{
Closes lists of outline coordinates.
}
\usage{
closed.outline(cont)
}
\arguments{
  \item{cont}{
  \code{list} or \code{matrix} of \eqn{(x; y)} coordinates}
}

\value{
The \code{list} of \eqn{(x; y)} coordinates provided with the first coordinates added at the end of the list.}

\examples{
fake <- list(x=1:4, y=5:8)
closed.outline(fake)
}
\keyword{Misc functions}