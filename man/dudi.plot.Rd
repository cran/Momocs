\name{dudi.plot}
\alias{dudi.plot}

\title{
A wrapper for \code{dudi.pca} graphical functions.
}
\description{
A wrapper for \code{dudi.pca} objects produced by \link{dudi.pca} in \code{ade4} or \link{pca} in \code{Momocs}.
}
\usage{
dudi.plot(dudi, fac = NULL, xax = 1, yax = 2, grid = TRUE,
	points = TRUE, pch.points = 1, col.points = "black", cex.points = 0.8,
	labels = FALSE, label = rownames(dudi$li),
	boxes = TRUE, clabel = 0.6,
	neighbors = FALSE, col.nei = "grey90", lwd.nei = 0.5,
	star = TRUE, col.star = "grey60", cstar = 1,
	ellipses = TRUE, col.ellipse = "grey30", cellipse = 1, axesell = TRUE,
	chull = FALSE, col.chull = "grey30", optchull = c(0.5, 1),
	arrows = FALSE, edge.arrow = FALSE,
	        box.arrow = TRUE, maxnb.arrow = 10, dratio.arrow = 0.2,
	shapes = TRUE, pos.shp=c("li", "circle", "range", "full")[3], 
		nr.shp = 6, nc.shp = 5, amp.shp = 1,
		scale.shp = 0.666, nb.pts.shp=300,
		first.point.shp = FALSE, rotate.shp=0,
		circle.nb.shp = 12, circle.r.shp,
		col.shp = "#70809011", border.shp = "#708090",
	rug = TRUE, rug.ticksize = 0.01, rug.col = "#708090",
	eigen = FALSE, eigen.ratio = 0.2,
	palette = col.sari, title = substitute(dudi), legend=FALSE,
	center.orig = FALSE, zoom.plot = 1)
}
\arguments{
  \item{dudi}{a dudi.pca object.}
  \item{fac}{the \code{(col)name} of the \code{@fac} slot to use as the grouping \code{factor}.}
  \item{xax}{\code{integer}. The index of the first PC axis to use.}
  \item{yax}{\code{integer}. The index of the first PC axis to use.}
  \item{grid}{\code{logical}. Whether to draw the grid.}
  \item{points}{\code{logical}. Whether to draw the points.}
  \item{pch.points}{the \code{pch} for drawing points.}
  \item{col.points}{the \code{col} for drawing points.}
  \item{cex.points}{the \code{cex} for drawing points.}
  \item{labels}{\code{logical}. Whether to draw labels.}
  \item{label}{\code{character}. The labels to draw.}
  \item{boxes}{\code{logical}. Whether to draw labels in boxes.}
  \item{clabel}{The \code{cex} for labels.}
  \item{neighbors}{\code{logical}. Whether to draw a neighboring graph.}
  \item{col.nei}{the \code{col} for drawing neighboring graph.}
  \item{lwd.nei}{the \code{lwd} for drawing neighboring graph.}
  \item{star}{\code{logical}. Whether to draw the star.}
  \item{col.star}{the \code{col} for drawing the star.}
  \item{cstar}{\code{numeric}. The size of the star.}
  \item{ellipses}{\code{logical}. Whether to draw bivariate confidence ellipses.}
  \item{col.ellipse}{The \code{col} for drawing these ellipses.}
  \item{cellipse}{\code{numeric}. The size of this ellipse.}
  \item{axesell}{\code{logical}. Whether to draw the ellipses axes.}
  \item{chull}{\code{logical}. Whether to draw a convex hull.}
  \item{col.chull}{The \code{col} for drawing convex hulls.}
  \item{optchull}{\code{numeric}. A vector of quantiles of the chulls.}
  \item{arrows}{\code{logical}. Whether to draw variables arrows.}
  \item{edge.arrow}{\code{logical}. Whether to neutralise arrows.}
  \item{box.arrow}{\code{logical}. Whether to draw boxes around arrows' labels.}
  \item{maxnb.arrow}{\code{numeric}. Only the \code{maxnb.arrow} most important variables will be draw.}
  \item{dratio.arrow}{\code{numeric}. Same idea as above but here the threshold is a relative to the size of the grid.}
  \item{shapes}{\code{logical}. Whether to plot shapes.}
  \item{pos.shp}{\code{character}, any of the \code{("li", "circle", "range")} methods. It specifies the way shapes have to be drawn. \code{"li"} draws shapes on the actual positions on the factorial map ; \code{"circle"} draws shapes on a circle with origin as the center ; \code{"full"} (not available in \link{morpho.space}) allows to cover the whole plotting region ; \code{"range"} draws shape on a rectangle that covers PC1 and PC2 range. Alternatively, a two columns matrix of coordinates where to calculate shapes.}
  \item{nr.shp}{\code{numeric}. The number of shape rows.}
  \item{nc.shp}{\code{numeric}. The number of shape columns.}
  \item{amp.shp}{\code{numeric}. An amplifying factor for shape deformation.}
  \item{scale.shp}{\code{numeric}. The size of shapes, relatively to grid size (0.5 = half of the grid size).}
  \item{nb.pts.shp}{\code{numeric}. The number of point to use to reconstruct shapes.}
  \item{first.point.shp}{\code{logical}. Whether to plot or not the first point when plotting shapes.}
  \item{rotate.shp}{\code{numeric}. If specified, the angles (in radians) to counter-clockwise rotate shapes plotted.}
  \item{circle.nb.shp}{\code{integer}. When \code{pos == "circle"}, the number of shapes on the circle.}
  \item{circle.r.shp}{\code{numeric}. When \code{pos == "circle"}, the circle radius.}
  \item{col.shp}{A color string for filling the shapes.}
  \item{border.shp}{A color string for shape borders.}
  \item{rug}{\code{logical}. Whether to add rug on axes (see \link{rug}).}
  \item{rug.ticksize}{\code{numeric}. The relative size of rug.}
  \item{rug.col}{A color string for rug ticks.}
  \item{eigen}{\code{logical}. Whether to add the eigen values barplot.}
  \item{eigen.ratio}{\code{numeric}. The relative size of the eigen values barplot.}
  \item{palette}{A color palette for group colors such as those produced by \link{colorRampPalette}.}
  \item{title}{A \code{character} string to be add on the graph.}
  \item{legend}{\code{logical}. Whether to add or not a legend.}.
  \item{center.orig}{\code{logical}. Whether to center the graphical window on the origin.}
  \item{zoom.plot}{\code{numeric}. Will help you to keep your distances, e.g. the value of magnification. Requires \code{center.orig} to be set true.}
}
\details{
This function takes the slightly modified \code{dudi.pca} object, obtained with \link{pca} on a \code{Coe} object and widely hinges on \code{ade4} plotting facilities and graphical layers. Notice that besides dedicated \emph{shapes} argument and options, it can be used for plotting factorial maps on "regular" \code{dudi.pca} objects. By default, \code{shapes} is set to \code{FALSE} for \link{tFourier} analyses since the reconstruction of shapes do not always lead to closed shapes (see Rohlf F, Archie J. 1984. A comparison of Fourier methods for the description of wing shape in mosquitoes (Diptera: Culicidae). Systematic Biology: 302-317.).
}
\references{
See the papers below that introduce \code{ade4} and also the package's homepage: \url{http://pbil.univ-lyon1.fr/ADE-4/} and particularly the file called "td83" : Lobry JR. 2010. Les fonctions graphiques 2D du paquet ade4.

Dray, S. and Dufour, A.B. (2007): The ade4 package: implementing the duality diagram for ecologists. \emph{Journal of Statistical Software}. \bold{22}(4): 1-20.

Chessel, D. and Dufour, A.B. and Thioulouse, J. (2004): The ade4 package-I- One-table methods. \emph{R News}. \bold{4}: 5-10.

Dray, S. and Dufour, A.B. and Chessel, D. (2007): The ade4 package-II: Two-table and K-table methods. \emph{R News}. \bold{7}(2): 47-52.

}
\seealso{
\link{pca}, \link{morpho.space}.
}
\examples{
data(bot)
botF <- eFourier(bot, nb.h=32)
botD <- pca(botF)
dudi.plot(botD)
dudi.plot(botD, 1, title="botD with no class but with ellipses")
dudi.plot(botD, fac=1, chull=TRUE, rug=FALSE, shape=FALSE, title="botD with convex hull")
dudi.plot(botD, fac=1, ellipses=FALSE, neighbors=TRUE, shapes=FALSE, star=FALSE,
	col.nei="black", title="botD with Gabriel's neighboring graph")
dudi.plot(botD, labels=TRUE, points=FALSE, boxes=FALSE, shapes=TRUE, pos.shp="li",
	title="botD with labels and reconstructed shapes")
dudi.plot(botD, 1, points=FALSE, labels=TRUE, boxes=FALSE, shapes=FALSE,
	title="botD with labels and ellipse")
dudi.plot(botD, 1, arrows=TRUE, dratio.arrow=0.2, shapes=FALSE,
	title="botD with harmonic correlations")
# With some fake factors
botD <- pca(botF)
dudi.plot(botD, "type", palette=col.gallus,
    rotate.shp=pi/2, title="botD with classes") # rotated shapes
dudi.plot(botD, "type", palette=col.gallus, eigen=TRUE, title="botD with eigen values")
dudi.plot(botD, "type", pos.shp="full", title="botD with shapes(1)")
dudi.plot(botD, "type", pos.shp="range", scale.shp=0.5, shapes=TRUE,
	border.shp="firebrick3", col.shp=NA, center.orig=TRUE, 
	zoom.plot=0.8, title="botD with shapes(2)")
dudi.plot(botD, "type", pos.shp="circle", center.orig=TRUE, title="botD with shapes(3)")
dudi.plot(botD, "type", pos.shp="range", scale.shp=0.5, title="botD with shapes(4)")
dudi.plot(botD, pos.shp=as.matrix(expand.grid(seq(-0.05, 0.05, 0.025),
    seq(-0.05, 0.05, 0.025)))) # an example with a matrix provided to pos.shp

}
\keyword{Multivariate Analysis}
