\name{manova.Coe}
\alias{manova.Coe}
\alias{manova.Coe-methods}
\alias{manova.Coe,Coe-method}

\title{
Multivariate ANOVA on \code{Coe} objects
}
\description{
A simple wrapper for Multivariate Analysis of Variance on the matrix of harmonic coefficients in a \code{Coe} object.
}

\usage{
manova.Coe(Coe, fac, retain, drop=0)
}
\arguments{
	\item{Coe}{a \code{Coe} object.}
	\item{fac}{\code{factor} defining which groups of individuals to compare.}
	\item{retain}{\code{numeric}. The number of harmonics to retain.}
	\item{drop}{\code{numeric}. The number of harmonics to drop.}
}
\details{The number of harmonic coefficient can not be higher than the number of outlines ; you can specify it with \code{retain} and \code{drop} (see \link{coeff.sel}) or let this method select the highest possible number of harmonics.}

\value{
The MANOVA summary is printed and also returned invisibly.}

\seealso{
See \link{manova} to change the defaults parameters of this method.
}

\examples{
data(bot)
botF <- eFourier(bot)
m <- manova.Coe(botF, "type")
m # if you need something from this MANOVA..
}
\keyword{Multivariate Analysis}