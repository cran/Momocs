\name{manova.nef}
\alias{manova.nef}
\docType{methods}
\alias{manova.nef-methods}
\alias{manova.nef,Nef-method}
\title{Calculates MANOVA on a harmonic coefficient matrix.}
\description{
Calculates Multivariate Analysis of Variance (MANOVA) on the harmonic coefficient matrix contained in \code{Nef}-objects.
}
\usage{manova.nef(Nef, fac, harmonics.retained, drop=FALSE)}
\arguments{
	\item{Nef}{
	the \code{Nef}-object}
	\item{fac}{a factor indicating the grouping desing}
	\item{harmonics.retained}{code{integer} indicating how many harmonics to include}
	\item{drop}{code{logical} indicating whether to drop or retain the first harmonic}
}
\details{
This function is a wrapper to calculate MANOVAs \emph{i.e.} test the significance of \emph{between} vs. \emph{within} geometric differences between sets of shapes. If not specified, the number of harmonics retained is calculated so that it is lower than the number of individuals minus two.
}
\examples{
data(bottles.nef)
fac <- factor(rep(c("beer", "whisky"), each=20))
manova.nef(bottles.nef, fac=fac)
}
\keyword{Nef methods}
