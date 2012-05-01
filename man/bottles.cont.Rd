\name{bottles dataset}
\alias{bottles.cont}
\alias{bottles.nef}
\docType{data}
\title{Two "bottles" datasets of outlines and harmonic coefficients.}
\description{
Two datasets are provided. First, \code{bottles.cont} contains 20 whisky and 20 beer bottles full outlines coordinates. Then, \code{bottles.nef} consists of a \code{Nef} object obtained with default parameters (32 harmonics, no smoothing) on \code{bottles.cont}.
}
\usage{
data(bottles.cont)
data(bottles.nef)
}
\format{
A \code{Coo}- class object that contains in the slot \code{@coo}, the lists of \eqn{(x; y)} coordinates.
A \code{Nef}-class object that contains in the slot \code{@coeff}, the matrix of harmonic coefficients.}
\source{
Images grabbed on the internet and prepared by the package's authors.
}
\examples{
data(bottles.cont)
bottles.cont

data(bottles.nef)
bottles.nef
}
\keyword{Package}
