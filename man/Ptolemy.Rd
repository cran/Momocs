\name{Ptolemy}
\alias{Ptolemy}
\alias{Ptolemy-methods}
\alias{Ptolemy,Coo-method}

\title{Ptolemaic ellipses and illustration of eFourier}

\description{
\code{Ptolemy} is a method for \code{Coo} objects to calculate and display Ptolemaic ellipses which illustrates intuitively the principle behing elliptical Fourier analysis.
}

\usage{
Ptolemy(Coo, id = 1, t = seq(0, 2 * pi, length = 7)[-1], nb.h = 3, 
        nb.pts = 360, palette = col.sari, legend = FALSE) 
}

\arguments{
	\item{Coo}{The \code{Coo} object on which to display Ptolemaic ellipses.}
	\item{id}{The \code{id} on which to display Ptolemaic ellipses.}
	\item{t}{A \code{vector} af angles (in radians) on which to display ellipses.}
	\item{nb.h}{\code{integer}. The number of harmonics to display.}
	\item{nb.pts}{\code{integer}. The number of points to use to display shapes.}
	\item{palette}{A color palette.}
	\item{legend}{\code{logical}. Whether to plot the legend box.}
}

\references{
This method has been inspired by the figures found in the followings papers.

Kuhl FP, Giardina CR. 1982. Elliptic Fourier features of a closed contour. \emph{Computer Graphics and Image Processing} \bold{18}: 236-258.

Crampton JS. 1995. Elliptical Fourier shape analysis of fossil bivalves: some practical considerations. \emph{Lethaia} \bold{28}: 179-186.
}

\seealso{
\link{efourier}. An intuitive explanation of elliptic Fourier analysis can be found in the \bold{Details} section of the \link{efourier} function.
}

\examples{
data(bot)
Ptolemy(bot, 1)
}

\keyword{elliptical Fourier analysis}

