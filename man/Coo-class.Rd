\name{Coo-class}
\Rdversion{1.1}
\docType{class}
\alias{Coo-class}
\alias{Coo}

\title{Class \code{"Coo"}.}
\description{
This class contains, so far, a single slot \code{@coo} that contains list of \eqn{(x; y)} coordinates. Calibration methods and Elliptical Fourier Analysis can be applied to \code{Coo}-objects. Additional slots will be implemented and store the dataset structure.  A \code{Coo}-object builder, \code{Coo()} is available to coerce a matrix and create a \code{Coo}-object.
}

\section{Slots}{
  \describe{
    \item{\code{coo}:}{a list of \code{(x; y)} coordinates.}
  }
}
\section{Methods}{
  \describe{
	\item{dev.qual}{Calculates and plots inverse reconstruction of outlines}
	\item{dev.quant}{Calculates and plots deviations between original and reconstructed outlines}
	\item{get.Nef}{Calculates Elliptical Fourier Analysis with specified parameters}
	\item{harm.pow}{Calculates and plots the Fourier power spectrum}
	\item{plot}{Plots a single or a range of outlines}
	 }
}

\examples{
data(bottles.cont)
\dontrun{
bottles.cont
plot(bottles.cont)
dev.qual(bottles.cont)
dev.quant(bottles.cont)
(get.Nef(bottles.cont))
bottles.cont@coo # to access the coordinates list
}
}
\keyword{Coo methods}
