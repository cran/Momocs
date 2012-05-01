\name{pca2shp}
\alias{pca2shp}
\title{Reconstructs a shape given using PCA.}
\description{
Provided with a harmonic matrix coefficient on which to perform PCA, and given the \eqn{(PC_1; PC_2)} coordinates, it reconstructs the corresponding shape.}
\usage{
pca2shp(pc1 = 0, pc2 = 0, data, nb.h = ncol(data)/4, nb.pts = 500,
		amp = 1, col = "black", lwd = 2, plot = TRUE)}
\arguments{
  \item{pc1}{\code{numeric} indicating the position on the first PC axis}
  \item{pc2}{\code{numeric} indicating the position on the second PC axis}
  \item{data}{the harmonic coefficient \code{matrix}}
  \item{nb.h}{\code{integer} indication how many harmonics to use}
  \item{nb.pts}{\code{integer} indicating the number of points sampled from the reconstructed outlines}
  \item{amp}{\code{numeric} indicating the magnifying factor}
  \item{col}{\code{integer} or \code{character} indicating the color to ue for drawing the shape}
  \item{lwd}{\code{numeric} indicating the \code{lwd} for the shape}
  \item{plot}{code{logical} indicating whether to plot the shape}
}
\examples{
data(bottles.nef)
pca2shp(0.5, 0.5, bottles.nef@coeff, amp=2)
pca2shp(0, 0, bottles.nef@coeff) # "average" shape
}

\keyword{Misc functions}
