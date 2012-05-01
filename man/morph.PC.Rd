\name{morph.PC}
\alias{morph.PC}
\docType{methods}
\alias{morph.PC-methods}
\alias{morph.PC,Nef-method}
\title{Plots the morphological space along PC axis.}
\description{
Given a matrix of harmonic coefficients, a \code{Nef}-object, calculates and plots morphological space \emph{i.e.} reconstructed shapes using and distributed on the orthonormal set defined by Principal Component axes.
}
\usage{
morph.PC(Nef, sd.nb=1, pca.ax=seq(1, 3))
}
\arguments{
	\item{Nef}{the \code{Nef} object}
	\item{sd.nb}{a \code{numeric} given the number of standard deviation to represent shape deviation along each PC axis}
	\item{pca.ax}{a numeric or a vector of numeric indicating on which Principal Component to display variation}}
	
\examples{
data(bottles.nef)
morph.PC(bottles.nef)
morph.PC(bottles.nef, 2, 1:3)
}
\keyword{Nef methods}


