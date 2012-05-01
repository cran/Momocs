\name{dev.qual}
\alias{dev.qual}
\docType{methods}
\alias{dev.qual-methods}
\alias{dev.qual,Coo-method}
\title{Calculates and plots reconstructed outlines.}
\description{
Calculates and plots inverse reconstruction of outlines based on the list of coordinates in a \code{Coo}-object and a given number of harmonics and/or smoothing iterations. This methods is the visual way to calibrate Elliptical Fourier Analysis parameters.
}
\usage{
dev.qual(Coo, id = 1:length(Coo@coo),
	 nb.h = 32, smooth.it = 0, range = seq(1, nb.h, len=4))
}
\arguments{
  \item{Coo}{the \code{Coo}-object}
  \item{id}{\code{integer} indicating the single or the range of outlines indices to consider}
  \item{nb.h}{\code{integer} indicating how many harmonics to calculate}
  \item{smooth.it}{\code{integer} indicating how many smoothing iterations to perform}
  \item{range}{\code{integer} indicating the range of harmonics orders to explore}
}
\examples{
data(bottles.cont)
\dontrun{
dev.qual(bottles.cont)
dev.qual(bottles.cont, id=24)
dev.qual(bottles.cont, id=24, nb.h=64)
dev.qual(bottles.cont, range=seq(1,16))
dev.qual(bottles.cont, smooth.it = 50)
}
}

\keyword{Coo methods}
