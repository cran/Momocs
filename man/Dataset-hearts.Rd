\name{Dataset: hearts}
\alias{hearts}
\docType{data}
\title{Hearts dataset}
\description{
The \code{hearts} dataset is a \code{Coo}-class object that contains outline coordinates of 450 hand-drawn hearts by 15 different persons.
}
\usage{
data(hearts)
}
\format{
A \code{Coo}-object that contains in the slot \code{@coo}, the lists of \eqn{(x; y)} coordinates and in the slot \code{@fac} the corresponding \code{factor}s.
}
\source{
We thank the fellows of the Ecology Department of the French Institute of Pondicherry that drawn the hearts, that then have been smoothed, scaled, centered, and reduced to 80 coordinates per outline.
}

\examples{
data(hearts)
hp <- pca(eFourier(hearts, nb.h=20, norm=FALSE, start=FALSE))
dudi.plot(hp, "aut",
		  pch=20, star=FALSE, cex.points=0.2, clabel=0.5, axesell=FALSE,
          nr.shp=15, nc.shp=15, scale.shp=0.2, pos.shp="full", col.shp=NULL)
}
\keyword{Datasets}
