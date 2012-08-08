\name{Outline: ellipse.par}
\alias{ellipse.par}
\alias{ellipse.par-methods}
\alias{ellipse.par,Nef-method}

\title{
Calculate ellipse parameters on a \code{Nef}-class object.
}
\description{
It uses \link{ellpar} function to calculate ellipse parameters for all individuals and all harmonics contained in a \code{Nef}-class object.}

\usage{
ellipse.par(Nef, range = 1:nrow(Nef@coeff), nb.pts = 120)
}

\arguments{
  \item{Nef}{The \code{Nef} object on which to calculate ellipses parameters.}
  \item{range}{Vector of \code{integer} that specifies the range of harmonics on which to calculate ellipse parameters. By default, all of them.}
  \item{nb.pts}{\code{integer}. How many points on the outline are used to estimate ellipse parameters (see \bold{Details} section below).}
}
\details{This method is not exact, see \link{ellpar}}

\value{
Four matrices are returned (with individual as rows and harmonic rank as columns).
	\item{a}{a \code{matrix} of all halves of the major axis.}
	\item{b}{a \code{matrix} of all halves of the minor axis.}
	\item{e}{a \code{matrix} of all eccentricities.}
	\item{phi}{a \code{matrix} of all ellipse angles.}
}
\references{
See \bold{References} in \link{ellpar}.
}

\examples{
data(bot)
botF <- eFourier(bot)
ellipse.par(botF, 1:5, 64)
}
\keyword{Outline Analysis}