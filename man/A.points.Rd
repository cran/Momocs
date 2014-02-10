\name{A.points}
\alias{A.points}

\title{
Plots the points that corresponds to landmarks provided as an array.
}

\description{
\code{A.points} displays landmarks when provided as an array with landmarks as rows, \eqn{(x; y)} coordinates as columns, and individuals as the third dimension of the array.
}


\usage{
A.points(A, col, palette=col.summer, pch=3, cex=1)
}

\arguments{
  \item{A}{An \code{array} of \eqn{(x; y)} coordinates, typically, landmarks.}
  \item{col}{A vector of colors to use for drawing landmarks.}
  \item{palette}{If \code{col} is not provided, then the color palette to use for landmarks.}
  \item{pch}{A \code{pch} for points.}
  \item{cex}{A \code{cex} for points.}
  }

\value{
No returned value.
}

\seealso{
\link{A.segments}, \link{A.mshape}, \link{A.plot} \link{coo.plot}, and the others coo.utilities.
}

\examples{
	\dontrun{
data(gorf.dat) # we import gorf.data from shapes package
A.points(gorf.dat)
	}
}
\keyword{Utilities}