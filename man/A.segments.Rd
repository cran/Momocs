\name{A.segments}

\alias{A.segments}

\title{
Utilities to manipulate arrays of homologous coordinates, typically landmarks.
}

\description{
\code{A.segments} draws interlandmarks segments as an array with landmarks as rows, \eqn{(x; y)} coordinates as columns, and individuals as the third dimension of the array.
}

\usage{
A.segments(A, col=NA, border.col="grey60", border.lty=1)
}

\arguments{
  \item{A}{An \code{array} of \eqn{(x; y)} coordinates, typically, landmarks.}
  \item{col}{A vector of colors to use for drawing landmarks.}
  \item{border.col}{A color for these segments.}
  \item{border.lty}{A \code{lty} for these segments.}
}

\value{
No returned value.
}

\seealso{
\link{A.plot}, \link{A.points}, \link{A.mshape} \link{coo.plot}, and the others coo.utilities.
}

\examples{
	\dontrun{
data(gorf.dat) # we import gorf.data from shapes package
A.segments(gorf.dat, border.lty=2, border.col="grey95")
m <- A.mshape(gorf.dat)
points(m, pch=5)
	}
}
\keyword{Utilities}