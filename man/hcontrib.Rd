\name{hcontrib}
\alias{hcontrib}
\alias{hcontrib-methods}
\alias{hcontrib,Coe-method}
\title{
Calculates and displays the contribution of every harmonic.
}
\description{
\code{hcontrib} applies \code{amp.h} values for each of the \code{harm.range} harmonics when reconstruction shapes. This thus help to visualize the respective contribution of every harmonic, in other words their contribution in describing shapes.
}
\usage{
hcontrib(Coe, id = 1, harm.range = 1:floor(Coe@nb.h / 4),
    amp.h   = c(0, 0.5, 1, 2, 5, 10),
    palette = col.sari, title   = "Harmonic contribution")
}
\arguments{
	\item{Coe}{The \code{Coe} object}
	\item{id}{\code{integer}. The \code{id} of the shape to display. A range of \code{id}s can be passed to \code{harm.pow}}
	\item{harm.range}{A \code{vector} of \code{integer} giving the harmonic range to calculate. See \code{nb.h} for \code{harm.pow}.}
	\item{amp.h}{A \code{vector} of numeric to specify the amplification factors to use.}
	\item{palette}{A color palette such those included in Momocs or produced with \link{colorRampPalette}.}
	\item{title}{A \code{character}to use as the title.}
}

\value{
No returned value. 
}

\examples{
#FILL
}
\keyword{Calibration}