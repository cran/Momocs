\name{Outline: eFourier}
\docType{methods}
\alias{eFourier}
\alias{eFourier-methods}
\alias{eFourier,Coo-method}
\title{Methods for Function \code{eFourier}}
\description{
This method performs an elliptical Fourier analysis on a \code{Coo}-class object, with the specified parameters (number of harmonics and number of smoothing iterations) and returns a \code{Nef}-class object containing harmonic coefficients (usually normalized). It accepts the same arguments as the function \link{efourier}.
}
\section{Methods}{
\describe{
\item{\code{signature(Coo = "Coo",
           nb.h      = 32,
           smooth.it = 0,
           normalize = TRUE,
           start     = FALSE)}}{
}
}}
\seealso{
The same arguments as \link{efourier} can be used with this method.
}

\keyword{Outline Analysis}

