\name{dataset: trilo}
\alias{trilo}
\docType{data}
\title{trilo dataset}
\description{
The \code{trilo} dataset is a \code{Coo}-class object that contains 64 coordinates of 50 cephalic outlines from different ontogenetic stages of  trilobite.
}

\usage{
data(trilo)
}
\format{
A \code{Coo}-object that contains in the slot \code{@coo}, the lists of \eqn{(x; y)} coordinates and ontogenetic stages in the slot \code{@fac}.
}

\source{
Hammer & Harper (2005) \emph{Paleontological Data Analysis} ed. Blackwell. See \url{http://folk.uio.no/ohammer/past/}.

Arranged from: \url{http://folk.uio.no/ohammer/past/outlines.dat}. The original data included 51 outlines and 5 ontogenetic stages, but one of them has just a single outline thas has been removed.
}


\examples{
data(trilo)
panel(trilo, names=TRUE)
}
\keyword{Datasets}
