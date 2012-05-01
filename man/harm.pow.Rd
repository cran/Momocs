\name{harm.pow}
\alias{harm.pow}
\docType{methods}
\alias{harm.pow-methods}
\alias{harm.pow,Coo-method}
\title{Calculates and plots Fourier harmonic spectra.}
\description{
Calculates and plots Fourier power spectra calculated as: \eqn{Power_n= \frac{A^2_n+B^2_n+C^2_n+D^2_n}{2} }.
}
\arguments{
  \item{coo}{the \code{Coo} object}
  \item{nb.h}{\code{integer} indicating how many harmonics to calculate}
  \item{smooth.it}{\code{integer} indicating how many smoothing iterations to perform}
  \item{plot}{\code{logical} indicating whether to plot the results}
  \item{max.h}{\code{integer} specifying the total number of harmonics to include}
  \item{first}{\code{logical} indicating whether to include the first harmonic}
}

\keyword{Coo methods}