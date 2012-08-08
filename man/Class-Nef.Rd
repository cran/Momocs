\name{Class: Nef}
\Rdversion{1.1}
\docType{class}
\alias{Nef}
\alias{Nef-class}

\alias{show,Nef-method}
\alias{boxplot-methods}
\alias{boxplot,Nef-method}
\alias{hist-methods}
\alias{hist,Nef-method}

\title{Class \code{"Nef"}}
\description{Objects of the class \code{Nef} aim at storing harmonic coefficients after an elliptic Fourier analysis. Usually they are normalized (nothe the "N" in "Nef"). This class will be derived in further versions for handling matrices of coefficients from other Fourier-based approaches (see \link{tfourier} and \link{rfourier}).
}

\section{Slots}{
  \describe{
    \item{\code{coeff}:}{A \code{matrix} of harmonic coefficients.}
    \item{\code{names}:}{\code{character} \code{vector} to name every shape.}
    \item{\code{fac}:}{\code{data.frame} defining the grouping structure, if any. Columns must be \code{factor}s.}
    \item{\code{nb.h}:}{\code{numeric} giving the number of harmonics.}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{Nef}-getters (will be restored in a further version)}
    \item{[<-}{\code{Nef}-setters (will be restored in a further version)}
    \item{boxplot}{displays the variation of harmonic coefficients.}
    \item{ellipse.par}{calculates ellipse parameters.}
    \item{harm.contrib}{displays the contribution of every harmonic.}
    \item{hist}{displays the distribution of of the harmonic coefficients values.}
    \item{manova.nef}{a wrapper for Multivariate ANOVA.}
    \item{pca}{a wrapper for performing a PCA with \code{ade4}.}
    \item{show}{displays a summary.}
    \item{tps}{calculates Thin Plate Splines on outlines.}
	 }
}

\examples{
showClass("Nef")
}
\keyword{Classes}
