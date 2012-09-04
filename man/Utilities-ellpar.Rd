\name{Utilities: ellpar}
\alias{ellpar}

\title{
Calculate ellipse parameters.
}

\description{
Given a \code{matrix} or \code{list} of coordinates that defines an ellipse, returns its geometrical parameters, namely: half-length of the major and the minor axis and excentricity.
}

\usage{
ellpar(coo)
}

\arguments{
  \item{coo}{a \code{matrix} or a \code{list} of \code{(x; y)} coordinates.}
}

\value{
	\item{a }{\code{numeric}. Half length of the major axis.}
	\item{b }{\code{numeric}. Half length of the minor axis.}
	\item{e }{\code{numeric}. Ellipse excentricity.}
}

\details{The aim of this function is to provide complementary descriptors for ellipses that are more explicit for humans than harmonic coefficients. For that reason they are also used directly instead of harmonic coefficients (See Schmittbuhl et al. 2003). The same study provided an exact approach to directly obtain ellipse parameters from harmonic coefficients but is not yet implemented in Momocs. So far, these parameters are estimated using matrices of coordinates obtained from harmonic coefficients, \emph{i.e.} by "redrawn" ellipses. Keep in mind that even if error rate is low (moreover if the number of points provided is high), the values returned must not be considered as exact.
}

\references{
Schmittbuhl M, Allenbach B, Le Minor J-M, Schaaf A, Minor J-marie L. 2003. Elliptical Descriptors: Some Simplified Morphometric Parameters for the Quantification of Complex Outlines. \emph{Mathematical Geology} \bold{35}: 853-871.
}

\seealso{
\link{ellipse.par}, the corresponding method for \code{Nef}-class objects.
}
\examples{
data(bot)
}
\keyword{Utilities}