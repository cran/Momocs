\name{eFourier}
\alias{eFourier}
\alias{eFourier-methods}
\alias{eFourier,Coo-method}

\title{
Calculates elliptical Fourier analysis on \code{Coo} objects.
}

\description{
\code{eFourier} performs an elliptical Fourier analysis on a \code{Coo}-class object, with the specified parameters (number of harmonics and number of smoothing iterations) and returns a \code{Coe}-class object containing harmonic coefficients (usually normalized). It accepts the same arguments as the function \link{efourier}.
}

\usage{
eFourier(Coo, nb.h = 32, smooth.it = 0, norm = TRUE, start = FALSE)
}

\arguments{
  \item{Coo}{A \code{Coo} object.}
  \item{nb.h}{\code{integer}. The number of harmonics to use}
  \item{smooth.it}{\code{integer}. The number of smoothing iterations to perform.}
  \item{norm}{\code{logical}. Whether to print or not diagnosis messages.}
  \item{start}{\code{logical}. Whether to define the first point of the outlines as an homologous point, and thus preserve this information. See Details and References in \link{efourier}.}
}

\seealso{
\link{efourier} and \link{efourier.norm}
}

\examples{
data(bot)
eFourier(bot)
}

\keyword{elliptical Fourier analysis}

