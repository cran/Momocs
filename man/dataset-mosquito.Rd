\name{dataset: mosquito}
\alias{mosquito}
\docType{data}
\title{mosquito dataset}
\description{
The \code{mosquito} dataset is a \code{Coo}-class object that contains outline coordinates of 126 outlines of mosquito wings used in Rohl and Archie (1984).
}
\usage{
data(mosquito)
}
\format{
A \code{Coo}-object that contains in the slot \code{@coo}, the lists of \eqn{(x; y)} coordinates. No grouping factor.
}

\source{
Rohlf F, Archie J. 1984. A comparison of Fourier methods for the description of wing shape in mosquitoes (Diptera: Culicidae). \emph{Systematic Biology}: 302-317.

Arranged from: \url{http://life.bio.sunysb.edu/morph/data/RohlfArchieWingOutlines.nts}.
}

\examples{
data(mosquito)
panel(mosquito)
}
\keyword{Datasets}
