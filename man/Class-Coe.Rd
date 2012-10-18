\name{Class: Coe}
\Rdversion{1.1}
\docType{class}
\alias{Coe}
\alias{Coe-class}

\alias{show,Coe-method}
\alias{boxplot,Coe-method}
\alias{hist,Coe-method}

\title{Class \code{"Coe"}}
\description{Objects of the class \code{Coe} aim at storing harmonic coefficients after a Fourier analysis, either \link{eFourier}, \link{rFourier} and \link{tFourier}.
}

\section{Slots}{
  \describe{
    \item{\code{coe}:}{A \code{matrix} of harmonic coefficients.}
    \item{\code{method}:}{A \code{character} that indicates the method used.}
    \item{\code{fac}:}{\code{data.frame} defining the grouping structure, if any. Columns must be \code{factor}s.}
    \item{\code{names}:}{\code{character} \code{vector} to name every shape.}
    \item{\code{nb.h}:}{\code{numeric} giving the number of harmonics.}
  }
}
\section{Methods}{
See \code{showMethods(class="Coe")}.
}

\examples{

## Working zone...
showClass("Coe")
}
\keyword{Classes}
