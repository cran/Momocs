\name{Multivariate Analysis: manova.nef}
\alias{manova.nef}
\alias{manova.nef-methods}
\alias{manova.nef,Nef-method}

\title{
Multivariate ANOVA on Nef objects
}
\description{
A simple wrapper for Multivariate Analysis of Variance on the matrix of harmonic coefficients in a \code{Nef} object.
}
\usage{
manova.nef(Nef, fac, retain=floor((nrow(Nef@coeff)-2)/4), drop=0)
}
\arguments{
	\item{Nef}{a \code{Nef} object.}
	\item{fac}{\code{factor} defining which groups of individuals to compare.}
	\item{retain}{\code{numeric}. The number of harmonics to retain.}
	\item{drop}{\code{numeric}. The number of harmonics to drop.}
}
\details{The number of harmonic coefficient can not be higher than the number of outlines ; you can specify it with \code{retain} and \code{drop} (see \link{coeff.sel}) or let this method select the highest possible number of harmonics.}
\value{
The MANOVA summary is printed.}

\seealso{
See \link{manova} to change the defaults parameters of this method.
}
\examples{
data(bot)
botF <- eFourier(bot)
manova.nef(botF, "type")
}
\keyword{Multivariate Analysis}