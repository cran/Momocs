\name{coo.plot}
\alias{coo.plot}
\title{
Plots a single outline.
}
\description{
\code{coo.plot} is a simple wrapper for plotting shapes. It basically tunes standard \link{plot} function to display single shapes as polygons, within a standardised plotting area.
}
\usage{
coo.plot(coo=NA, col="#F5F5F5",
         border="#1A1A1A", lwd=1, lty = 1, xlim=c(-1, 1), ylim=c(-1, 1),
         points=FALSE, first.point=TRUE, centroid=TRUE,
         points.col=border, pch=1, cex=0.8, main, ...)
}

\arguments{
	\item{coo}{A \code{list} or a \code{matrix} of coordinates.}
	\item{col}{A color to fill the shape polygon.}
	\item{border}{A color for the shape border.}
	\item{lwd}{The \code{lwd} for drawing shapes.}
	\item{lty}{The \code{lty} for drawing shapes.}
	\item{xlim}{If \code{coo.plot} is called and \code{coo} is missing, then a vector of length 2 specifying the \code{ylim} of the ploting area.}
	\item{ylim}{If \code{coo.plot} is called and \code{coo} is missing, then a vector of length 2 specifying the \code{ylim} of the ploting area.}
	\item{points}{\code{logical}. Whether to display points. If missing and number of points is < 100, then points are plotted.}
	\item{first.point}{\code{logical}. Whether to display the first point.}
	\item{centroid}{\code{logical}. Whether to display centroid.}
	\item{points.col}{The color for plotting points.}
	\item{pch}{The \code{pch} for points.}
	\item{cex}{The \code{cex} for points.}
	\item{main}{\code{character}. A title for the plot.}
	\item{...}{Additional parameters for drawing the first point.}
}

\value{
No returned value.
}
\seealso{
\link{coo.draw}.
}

\examples{
data(bot)
coo.plot(bot@coo[[1]])
}
\keyword{Utilities}