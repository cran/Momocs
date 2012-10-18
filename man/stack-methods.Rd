\name{stack(Coo)}
\alias{stack}
\alias{stack-methods}
\alias{stack,Coo-method}
\alias{plot,Coo,ANY-method}

\title{
Plots all the outlines from a \code{Coo} on the same graph.
}

\description{
\code{stack(Coo)} plots all the outlines from a \code{Coo} on the same graph with graphical options. Note that you can also use \code{plot(Coo)} with the very same arguments.
}

\usage{
stack(x, cols, borders,
           ldk=TRUE, ldk.pch=3, ldk.col="red", ldk.cex=1, ...)
}

\arguments{
  \item{x}{The \code{Coo} object to plot.}
  \item{cols}{A \code{vector} of colors for drawing the outlines. Either a single value or of length exactly equals to the number of coordinates.}
  \item{borders}{A \code{vector} of colors for drawing the borders. Either a single value or of length exactly equals to the number of coordinates.}
  \item{ldk}{\code{logical}. Whether to display landmarks (if any).}
  \item{ldk.pch}{A \code{pch} for these landmarks.}
  \item{ldk.col}{A color for these landmarks.}
  \item{ldk.cex}{A \code{cex} for these landmarks.}
  \item{...}{Additional arguments to be passed to \link{coo.draw}.}
}

\value{
No returned value.
}

\seealso{
\link{stack}, \link{diapo}.
}


\examples{
data(mosquito)
stack(mosquito, borders="#1A1A1A22", first.point=FALSE)

data(hearts)
stack(hearts, borders="#1A1A1A22", ldk=FALSE)
stack(hearts, borders="#1A1A1A22", ldk=TRUE, ldk.col="#1A1A1A55")
stack(hearts, borders="#1A1A1A22", ldk=TRUE, ldk.col=col.summer(4), ldk.pch=20)
}
\keyword{Utilities}