\name{coo.draw}
\alias{coo.draw}
\title{
Adds a single outline on the current plot.
}
\description{
\code{coo.draw} is a light version of \link{coo.plot} that simply adds a shape on the active plot.
}
\usage{
coo.draw(coo = NA, col = "#70809033",
         border = "#708090EE", lwd=1, lty = 1, points = FALSE,
         first.point=TRUE, centroid = FALSE,
         points.col = border, pch = 20, cex = 0.25,  ...)
}

\arguments{
	\item{coo}{A \code{list} or a \code{matrix} of coordinates.}
	\item{col}{A color to fill the shape polygon.}
	\item{border}{A color for the shape border.}
	\item{lwd}{The \code{lwd} for drawing shapes.}
	\item{lty}{The \code{lty} for drawing shapes.}
	\item{points}{\code{logical}. Whether to display points. If missing and number of points is < 100, then points are plotted.}
	\item{first.point}{\code{logical}. Whether to display the first point.}
	\item{centroid}{\code{logical}. Whether to display centroid.}
	\item{points.col}{The color for plotting points.}
	\item{pch}{The \code{pch} for points.}
	\item{cex}{The \code{cex} for points.}
	\item{...}{Additional parameters for drawing the first point.}
}

\value{
No returned value.
}
\seealso{
\link{coo.plot}.
}

\examples{
data(bot)
coo.plot(bot@coo[[1]])
coo.plot(bot@coo[[2]], lwd=1.2)
}
\keyword{Utilities}