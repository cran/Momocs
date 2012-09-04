\name{Utilities: coo.plot}
\alias{coo.plot}
\alias{coo.draw}
\alias{coo.list.panel}
\title{
Plots sets of shapes.
}
\description{
\code{coo.plot} is a simple wrapper for plotting shapes. It basically tunes plot standard function to display single shapes as polygons, within a standardised plotting area. \code{coo.draw} simply adds a shape on the active plot. \code{coo.list.panel} plots a list of shapes.
}
\usage{
coo.plot(coo=NA, col="#70809033",
		border="#708090EE", xlim=c(-1, 1), ylim=c(-1, 1), points=TRUE, first.point=TRUE,
		points.col=border, pch=20, cex=0.25, main, ...)

coo.draw(coo = NA, col = "#70809033",
		border = "#708090EE", points = TRUE, first.point=TRUE,
		points.col = border, pch = 20, cex = 0.25,  ...)

coo.list.panel(coo.list, dim,
		byrow = TRUE, fromtop = TRUE,
		mar = rep(0, 4), cols, borders, density = NULL, angle = 45)
}

\arguments{
	\item{coo}{A \code{list} or a \code{matrix} of coordinates.}
	\item{col}{A color to fill the shape polygon.}
	\item{border}{A color for the shape border.}
	\item{xlim}{If \code{coo.plot} is called and \code{coo} is missing, then a vector of length 2 specifying the \code{ylim} of the ploting area.}
	\item{ylim}{If \code{coo.plot} is called and \code{coo} is missing, then a vector of length 2 specifying the \code{ylim} of the ploting area.}
	\item{points}{\code{logical}. Whether to display points. If missing and number of points is < 100, then points are plotted.}
	\item{first.point}{\code{logical}. Whether to display the first point.}
	\item{points.col}{The color for plotting points.}
	\item{pch}{The \code{pch} for points.}
	\item{cex}{The \code{cex} for points.}
	\item{main}{\code{character}. A title for the plot.}
	\item{...}{Additional parameters for drawing the first point.}
	\item{coo.list}{A \code{list} of coordinates, such as those in the \code{@coo} slot of \code{Coo} objects.}
	\item{dim}{A \code{vector} of the form \code{(nb.row, nb.cols)} to specify the panel display. If missing, shapes are arranged in a square.}
	\item{byrow}{\code{logical}. Whether to succesive shape by row or by col.}
	\item{fromtop}{\code{logical}. Whether to display shapes from the top of the plotting region.}
	\item{mar}{A \code{vector} to define margins.}
	\item{cols}{A \code{vector} of colors to fill shapes.}
	\item{borders}{A \code{vector} of colors to draw shape borders.}
	\item{density}{A \code{vector} for density of shading lines. See \link{polygon}}
	\item{angle}{A \code{vector} for shading lines.}
}

\examples{
data(bot)
coo.plot(bot@coo[[1]])
coo.list.panel(bot@coo)
}
\keyword{Utilities}