\name{dev.quant}
\alias{dev.quant}
\docType{methods}
\alias{dev.quant-methods}
\alias{dev.quant,Coo-method}
\title{Calculates and plots sum of euclidean deviations between original and reconstructed outlines.}
\description{
Calculates and plots sum of euclidean deviations between one or a range of original and reconstructed shapes, normalized by the calliper length, \emph{i.e.} the longest length measured between two outlines points.
}
\usage{
dev.quant(Coo, id = 1:length(Coo@coo), nb.h = 32, smooth.it = 0, plot=TRUE)
}
\arguments{
  \item{Coo}{the \code{Coo} object}
  \item{id}{\code{integer} indicating the single or the range of outlines indices to consider}
  \item{nb.h}{\code{integer} indicating how many harmonics to calculate}
  \item{smooth.it}{\code{integer} indicating how many smoothing iterations to perform}
  \item{plot}{a \code{logical} indicating whether to plot or not the results}
}

\examples{
\dontrun{
data(bottles.cont)
dev.quant(bottles.cont, id=4)
dev.quant(bottles.cont, id=4, nb.h=12)
}
}
\keyword{Coo methods}
