\name{Nef-class}
\Rdversion{1.1}
\docType{class}
\alias{Nef-class}
\alias{Nef}
\title{Class \code{"Nef"}.}
\description{
A class that contains all the information to visualize and performe multivariate analysis. Contains so far a single slot \code{@coeff} containing the harmonic coefficient matrix after an Elliptical Fourier Analysis. A \code{Nef}-object builder, \code{Nef()} is available.
}
\section{Slots}{
  \describe{
    \item{\code{coeff}:}{a matrix of harmonic coefficients.}
  }
}
\section{Methods}{
  \describe{
    \item{manova.nef}{Calculates MANOVA on a harmonic coefficient matrix}
    \item{morph.sp}{Plots the morphological space}
    \item{pca.tps}{Plots a single PCA with deformation grids}
    \item{pca}{Plots a single PCA}
    \item{pca3}{Plots all the first three PCA axes}
    \item{show}{A simple object description}
    \item{tps.grid}{Thin Plate Splines deformation grids between two shapes}
    \item{tps.iso}{TPS and iso-deformation lines between two shapes}
    \item{tps.vf}{TPS and "vector field" of deformation between two shapes}
    \item{traj}{Calculates shape intermediates}
	 }
}
\keyword{Nef methods}
