\name{get.Nef}
\alias{get.Nef}
\docType{methods}
\alias{get.Nef-methods}
\alias{get.Nef,Coo-method}
\title{Calculates Elliptical Fourier Analysis.}
\description{
Calculates Elliptical Fourier Analysis with specified parameters.
}
\usage{
get.Nef(Coo, nb.h=32, smooth.it=0, fromrt=FALSE)
}
\arguments{
  \item{Coo}{the \code{Coo} object}
  \item{nb.h}{\code{integer} indicating how many harmonics to calculate}
  \item{smooth.it}{\code{integer} indicating how many smoothing iterations to perform}
  \item{fromrt}{\code{logical} whether the position of the first point has to be preserved}
}
\value{
a \code{Coo}-object is returned.
}

\examples{
data(bottles.cont)
nef <- get.Nef(bottles.cont)
pca3(nef)
}

\keyword{Coo methods}