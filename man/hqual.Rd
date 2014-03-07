\name{hqual}
\alias{hqual}
\alias{hqual-methods}
\alias{hqual,Coo-method}

\title{
Displays how shapes are described with a given number of harmonics.
}

\description{
\code{hqual} is used to calculate and displays reconstructed shapes using a range of harmonic number.
}

\usage{
hqual(Coo, method = c("efourier", "rfourier", "tfourier"),
                 id = sample(Coo@coo.nb, 1), smooth.it = 0,
                 harm.range = c(1, 2, 4, 8, 16, 32),
                 scale = TRUE, center = TRUE, align = TRUE,
                 plot.method = c("stack", "panel")[1], legend = TRUE,
                 legend.title = "# harmonics",
                 palette = col.summer, shp.col = "#70809033",
                 shp.border = "#708090EE")
}
\arguments{
	\item{Coo}{The \code{Coo} object}
	\item{method}{A character, either "efourier", "rfourier" or "tfourier" (partial matches are allowed) to use to calculate morphological that indicates which method to use.}
	\item{id}{\code{integer}. The \code{id} of the shape to display. A range of \code{id}s can be passed to \code{harm.pow}}
	\item{smooth.it}{\code{integer}. The number of smoothing iteration to perform.}
	\item{harm.range}{A \code{vector} of \code{integer} giving the harmonic range to calculate. See \code{nb.h} for \code{harm.pow}.}
	\item{scale}{\code{logical}. Whether to scale or not the shape.}
	\item{center}{\code{logical}. Whether to center or not the shape.}
	\item{align}{\code{logical}. Whether to align or not the shape.}
	\item{plot.method}{A plot method to use, either \link{stack} or \link{panel}.}
	\item{legend}{\code{logical}. Whether to display a legend box.}
	\item{legend.title}{\code{character}. A title for the legend box.}
	\item{palette}{A color palette such those included in Momocs or produced with \link{colorRampPalette}.}
	\item{shp.col}{A color for the shape body.}
	\item{shp.border}{A color for the shape border.}
}


\value{
No value returned.
}

\examples{
data(bot)
hqual(bot)
}
\keyword{Calibration}