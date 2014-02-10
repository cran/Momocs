\name{coo.close}
\alias{coo.close}

\title{
Closes a list or matrix of coordinates.
}

\description{
\code{coo.close} closes \code{coo}, \emph{i.e.} makes the last coordinates be the first than the first.
}

\usage{
coo.close(coo)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  }


\value{
Returns a matrix of \code{(x; y)}coordinates.
}

\seealso{
\link{coo.unclose} and \link{is.closed}
}

\examples{
\dontrun{
data(gorf.dat) # we import gorf.data from shapes package
coo <- gorf.dat[,,1]
is.closed(coo)
coo.c <- coo.close(coo)
is.closed(coo.c)
coo.cu <- coo.unclose(coo)
is.closed(coo.cu)
}
}

\keyword{coo Utilities}