\name{dataset: hearts}
\alias{hearts}
\docType{data}
\title{hearts dataset}
\description{
The \code{hearts} dataset is a \code{Coo}-class object that contains outline coordinates of 240 hand-drawn hearts by 8 different persons.
}

\usage{
data(hearts)
}
\format{
A \code{Coo}-object that contains in the slot \code{@coo}, the lists of \eqn{(x; y)} coordinates and in the slot \code{@fac} a single factor "aut" with the name of the authors. Four landmarks are also registered.
}
\source{
We thank the fellows of the Ecology Department of the French Institute of Pondicherry that drawn the hearts, that then have been smoothed, scaled, centered, and reduced to 80 coordinates per outline.
}

\examples{
data(hearts)
levels(hearts@fac[,"aut"]) # thank you guys !
panel(hearts, cols=rep(col.summer(8), each=30), names=TRUE)
}
\keyword{Datasets}
