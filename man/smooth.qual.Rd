\name{smooth.qual}
\alias{smooth.qual}
\alias{smooth.qual-methods}
\alias{smooth.qual,Coo-method}

\title{
Displays effect of smoothing on shape reconstruction.
}
\description{Helps to visually estimate the number of smoothing iterations required by drawing original shape and smoothed outlines.
}
\usage{
smooth.qual(Coo, id=1,
           smooth.range = c(10, 50, 200, 500, 1000),
           palette = col.summer)
}
\arguments{
	\item{Coo}{The \code{Coo}-object}
	\item{id}{\code{integer}. The \code{id} of the shape to display.}
	\item{smooth.range}{A \code{vector} of \code{integer} theat specify the number of smoothing iteration to perform and display.}
	\item{palette}{A color palette.}
	}

\references{
See Claude and Crampton & Crampton works below for details on smoothing outlines before an outline analysis:

Haines AJ, Crampton JS. 2000. Improvements To The Method Of Fourier Shape Analysis As Applied In Morphometric Studies. \emph{Palaeontology} \bold{43}: 765-783.

Crampton JS. 1995. Elliptical Fourier shape analysis of fossil bivalves: some practical considerations. \emph{Lethaia} \bold{28}: 179-186.
}
\examples{
data(bot)
smooth.qual(bot, 15)
}
\keyword{Calibration}