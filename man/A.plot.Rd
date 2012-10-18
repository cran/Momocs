\name{A.plot}
\alias{A.plot}

\title{
Plots an array of homologous coordinates, typically landmarks.
}

\description{
\code{A.plot} is used to plot an array of landmarks when provided with landmarks as rows, \eqn{(x; y)} coordinates as columns, and individuals as the third dimension of the array.
}

\usage{
A.plot(A, col, palette = col.summer, xlim, ylim, border.col = NA, border.lty = 2, pch = 3, cex = 1)
}

\arguments{
  \item{A}{An \code{array} of \eqn{(x; y)} coordinates, typically, landmarks.}
  \item{col}{A vector of colors to use for drawing landmarks.}
  \item{palette}{If \code{col} is not provided, then the color palette to use for landmarks.}
  \item{xlim}{A vector of length 2 specifying a custom \code{xlim} of the ploting area.}
  \item{ylim}{A vector of length 2 specifying a custom \code{ylim} of the ploting area.}
  \item{border.col}{A color for these segments.}
  \item{border.lty}{A \code{lty} for these segments.}
  \item{pch}{A \code{pch} for points.}
  \item{cex}{A \code{cex} for points.}
  }

\value{
No returned value.
}


\seealso{
\link{A.points}, \link{A.segments}, \link{A.mshape} \link{coo.plot}, and the others coo.utilities.
}
\examples{
data(gorf.dat) # we import gorf.data from shapes
A.plot(gorf.dat, pch=20)
}
\keyword{Utilities}