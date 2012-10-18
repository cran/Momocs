\name{morpho.space}
\alias{morpho.space}

\title{
Calculate and plots morphological spaces.}

\description{
Calculates and plots morphological spaces with many options.
}

\usage{
morpho.space(dudi, xax = 1, yax = 2, xlim, ylim, nb.pts = 300,
                 pos.shp = c("li", "circle", "range")[3], nr.shp = 6,
                 nc.shp = 5, amp.shp = 1, scale.shp = 1, rotate.shp =
                 0, circle.nb.shp = 12, circle.r.shp, plot = TRUE,
                 layer = TRUE, col.shp = "#70809011", border.shp =
                 "#708090", pch.pts = 20, col.pts = "grey40",
                 first.point = FALSE)
}
\arguments{
  \item{dudi}{a dudi.pca object.}
  \item{xax}{\code{integer}. The index of the first PC axis to use.}
  \item{yax}{\code{integer}. The index of the first PC axis to use.}
  \item{xlim}{A vector of 2 \code{numeric} indicating the x-range.}
  \item{ylim}{A vector of 2 \code{numeric} indicating the y-range.}
  \item{nb.pts}{The number of points to calculate.}
  \item{pos.shp}{\code{character}, any of the \code{("li", "circle", "range")} methods. It specifies the way shapes have to be drawn. \code{"li"} draws shapes on the actual positions on the factorial map ; \code{"circle"} draws shapes on a circle with origin as the center ; \code{"range"} draws shape on a rectangle that covers PC1 and PC2 range. Alternatively, a two columns matrix of coordinates where to calculate shapes.}
  \item{nr.shp}{\code{numeric}. The number of shape rows.}
  \item{nc.shp}{\code{numeric}. The number of shape columns.}
  \item{amp.shp}{\code{numeric}. An amplifying factor for shape deformation.}
  \item{scale.shp}{\code{numeric}. The size of shapes, relatively to 1/8 of the highest range (x or y).}
  \item{rotate.shp}{\code{numeric}. If specified, the angles (in radians) to counter-clockwise rotate shapes plotted.}
  \item{circle.nb.shp}{\code{integer}. When \code{pos == "circle"}, the number of shapes on the circle.}
  \item{circle.r.shp}{\code{numeric}. When \code{pos == "circle"}, the circle radius.}
  \item{plot}{\code{logical}. Whether to plot or not the morphological space.}
  \item{layer}{\code{logical}. Whether to add calculated shapes on an existing plot. Used for instance by \link{dudi.plot}.}
  \item{col.shp}{A color string for filling the shapes.}
  \item{border.shp}{A color string for shape borders.}
  \item{pch.pts}{A \code{pch} for plotting the points from the \code{dudi}-class object. Use \code{NA} to not display these points.}
  \item{col.pts}{A color to plot these points.}
  \item{first.point}{\code{logical}. Whether to draw or not the first point of the shapes.}
}

\details{
Note that \link{tFourier} analyses since the reconstruction of shapes do not always lead to closed shapes (see Rohlf F, Archie J. 1984. A comparison of Fourier methods for the description of wing shape in mosquitoes (Diptera: Culicidae). Systematic Biology: 302-317.)
}

\seealso{
\link{dudi.plot}.
}

\examples{
data(bot)
botF <- eFourier(bot)
botD <- pca(botF)
morpho.space(botD)

}
\keyword{Multivariate Analysis}
