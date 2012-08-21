\name{Outline: tps}
\alias{tps}
\alias{tps-methods}
\alias{tps,Nef-method}
	
\title{Thin Plin Splines}

\description{
A method to compute Thin Plate Splines and display deformation grids.
}
\usage{
tps(Nef, id.fr = 1, id.to = 5, nb.pts = 64, amp = 1,
	method = c("grid", "arr", "iso")[1],
	grid.size = 20, grid.col = "grey50",
	arr.nb = 100, arr.levels = 100,
	arr.len = 0.1, arr.ang = 30,
	arr.lwd = 1, arr.col = "grey50",
	iso.levels = 10, iso.nb = 5000,
    palette = col.summer, cont = TRUE,
	cont.col = col.gallus(2), cont.lwd = c(2, 2))
}
\arguments{
\item{Nef}{the \code{Nef} class object to use.}
\item{id.fr}{\code{integer}. The \code{id} of the 'starting' shape.}
\item{id.to}{\code{integer}. The \code{id} of the 'ending' shape.}
\item{nb.pts}{\code{integer}. The number of points to use to reconstruct shapes.}
\item{amp}{\code{numeric}. An amplifying factor for differences between shapes.}
\item{method}{\code{character}; \code{"grid"}, \code{"arr"} and \code{"iso"} returns a simple deformation grid, a 'vector field' and isodeformation lines, respectively.}
\item{grid.size}{\code{numeric}. The size of the deformation grid.}
\item{grid.col}{color. The color of the grid.}
\item{arr.nb}{\code{integer}. The number of arrows to display.}
\item{arr.levels}{\code{integer}. The number of arrows levels for color purpose.}
\item{arr.len}{\code{numeric}. The length of these arrows.}
\item{arr.ang}{\code{numeric}. The angle of the arrows edges.}
\item{arr.lwd}{\code{numeric}. The arrows line width.}
\item{arr.col}{color. The color of these arrows (if palette is not supplied).}
\item{iso.levels}{\code{integer}. The number of isolines levels.}
\item{iso.nb}{\code{integer}. The number of lines to draw.}
\item{palette}{a color palette such as those produced by \link{colorRampPalette}.}
\item{cont}{\code{logical}. Whether to draw or not the shapes outlines.}
\item{cont.col}{A \code{vector} containing the two as colors for the starting and ending shapes, respectively.}
\item{cont.lwd}{A \code{vector} containing the two line widths as numeric for the outlines.}
}
\details{Amplification is calculated after on Principal Component Analysis, \emph{i.e.} the deformation are based on the mean shape and the position of the stwo shapes on the PC1 and PC2 axes.}

\examples{
data(bot)
botF <- eFourier(bot)
tps(botF, 1, 2)
tps(botF, 1, 20, method="arr", arr.len=0.05, arr.nb=300, palette=col.sari)
tps(botF, 1, 20, method="iso", palette=col.gallus, amp=1.5)
}
\keyword{Outline Analysis}