\name{hquant}
\alias{hquant}
\alias{hquant-methods}
\alias{hquant,Coo-method}

\title{
Calculates deviations between reconstructed and best possible shapes using different distance methods.
}
\description{
\code{hquant} is used to estimate the number of harmonics required to describe shapes with the required accuracy. So far, two different methods can be used: \link{edm} which calculates point-to-point distance and \link{edm.nearest} which calculates for every point of the reconstructed shape the nearest point in the best possible shape. The best possible shape is calculated using the highest possible number of harmonics.
}
\usage{
hquant(Coo, method = c("efourier", "rfourier", "tfourier"),
                 id = 1, smooth.it = 0, harm.range = c(seq(4, 24, 4),
                 200), norm.centsize = TRUE, dist.method = edm.nearest,
                 plot = TRUE, dev.plot = TRUE, title = "Deviations along the outline",
                 legend = TRUE, palette = col.summer, lineat.y = c(0.5, 0.1, 0.01))
}
\arguments{
	\item{Coo}{The \code{Coo} object}
	\item{method}{A character, either "efourier", "rfourier" or "tfourier" (partial matches are allowed) to use to calculate morphological that indicates which method to use.}
	\item{id}{\code{integer}. The \code{id} of the shape to display. A range of \code{id}s can be passed to \code{harm.pow}}
	\item{smooth.it}{\code{integer}. The number of smoothing iteration to perform.}
	\item{harm.range}{A \code{vector} of \code{integer} giving the harmonic range to calculate. See \code{nb.h} for \code{harm.pow}.}
	\item{norm.centsize}{\code{logical}. Whether to normalize distances by the centroid size of every outline.}
	\item{dist.method}{A distance method to use, e.g. either \link{edm.nearest} or \link{edm}.}
	\item{plot}{\code{logical}. Whether to plot or not the shape. If \code{FALSE}, only the results are returned.}
	\item{dev.plot}{\code{logical}. Whether to plot or not the deviation plot.}
	\item{title}{A title for the plot.}
	\item{legend}{\code{logical}. Whether to display a legend box.}
	\item{palette}{A color palette such those included in Momocs or produced with \link{colorRampPalette}.}
	\item{lineat.y}{A \code{vector} of \code{numeric} to specify where to plot dashed lines on the y-axis.}
}


\value{
Returns a matrix containing deviations for each harmonic and corresponding quantiles.
}

\examples{
\dontrun{
data(bot)
hquant(bot)
}
}
\keyword{Calibration}