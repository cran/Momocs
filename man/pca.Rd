\name{pca}
\alias{pca}
\alias{pca3}
\alias{pca.tps}
\docType{methods}
\alias{pca-methods}
\alias{pca,Nef-method}
\alias{pca.tps-methods}
\alias{pca.tps,Nef-method}
\alias{pca3-methods}
\alias{pca3,Nef-method}
\title{Calculates and plots Principal Component Analysis.}
\description{
Calculates and plots Principal Component Analysis usng \code{\link[=prcomp]{prcomp()}}. Methods for plotting a single PCA, a triple PCA and deformation grids are detailed below.
}
\usage{
pca(Nef, fac = NA,PCa = 1, PCb = 2,
	col = "black", pch = 1, lty=1, shp.nb=NA, shp.size,
	shp.col="#00000022", shp.border="black", title = "Principal Component Analysis",
	legend = TRUE, lab = FALSE, lab.txt = rownames(Nef@coeff), lab.cex = 1, lab.box = TRUE,
	ell = TRUE, r = 1, lwd = 1, zoom.x = 0.25, zoom.y = 0.3)

pca3(Nef,  fac = NA,
	col = 1:nlevels(fac), pch = 1:nlevels(fac), lty = rep(1,nlevels(fac)),
	lab = FALSE, lab.txt = rownames(Nef@coeff), lab.cex = 1, lab.box = TRUE,
	ell = 1, r = 1, lwd = 1, zoom = 1.4, legend = FALSE)

pca.tps(Nef, fac = NA, PCa = 1, PCb = 2,
		col = "black", pch = 1, ell = TRUE, zoom = 1.4, ncells = 20,
		title = "Deformations alongs PC axes")
}
\arguments{
	\item{Nef}{the \code{Nef} object}
	\item{fac}{the grouping \code{factor}}
	\item{PCa}{\code{integer} corresponding to the \code{a^th} PCA axis to plot}
	\item{PCb}{\code{integer} corresponding to the \code{a^th} PCA axis to plot}
	\item{col}{\code{integer} or \code{character} indicating the \code{col} for each group to plot}
	\item{pch}{\code{integer} or \code{character} indicating the \code{pch} for each groups to plot}
	\item{lty}{\code{integer} indicating the \code{lty} for each confidence ellipse to plot}
	\item{shp.nb}{\code{integer} indicating how many, if any, shapes to plot}
	\item{shp.size}{\code{numeric} indicating the size of these shapes}
	\item{shp.col}{\code{integer} or a character indicating the color of these shapes}
	\item{shp.border}{\code{integer} or a character indicating the border color of these shapes}
	\item{title}{\code{character} to add a better title to the plot}
	\item{lab}{\code{logical} indicating whether to plot labels for every point}
	\item{lab.txt}{\code{character} vector containing labels names}
	\item{lab.cex}{\code{numeric} indicating the \code{cex} size of these labels}
	\item{lab.box}{\code{logical} indicating whether to draw a border for these labs}
	\item{ell}{\code{logical} indicating whether to draw confidence ellipses for every group}
	\item{r}{\code{numeric} indicating the number of standard deviations for confidence ellipses computation}
	\item{lwd}{\code{numeric} indicating the \code{lwd} for the confidence ellipses}
	\item{zoom}{\code{numeric} used to adjust the range of the plot}
	\item{zoom.x}{\code{numeric} used to adjust the x-range of the plot}	
	\item{zoom.y}{\code{numeric} used to adjust the y-range of the plot}
	\item{legend}{\code{logical} indicating whether to add a legend on the plot}
	\item{ncells}{\code{integer} indicating the number of cells for deformation grids}}
\examples{

\dontrun{
data(bottles.nef)
fac <- factor(rep(c("beer", "whisky"), each=20))

### pca
pca(bottles.nef)
pca(bottles.nef, fac=fac)
pca(bottles.nef, fac=fac, pch=c(4,5))
pca(bottles.nef, fac=fac, pch=c(4,5), ell=FALSE)
pca(bottles.nef, fac=fac,pch=c(4,5), lty=c(2,3))
pca(bottles.nef, pch=c(4,5), fac=fac, lty=c(2,3), col=c("dodgerblue","firebrick"))
pca(bottles.nef, fac=fac, lab=T)
pca(bottles.nef, fac=fac, lab=T, lab.cex=0.8)
pca(bottles.nef, fac=fac, lab=T, lab.box=FALSE)
pca(bottles.nef, fac=fac, lab=T, lab.txt=c(letters[1:20],LETTERS[1:20]))
pca(bottles.nef, fac=fac, r=0.5)
pca(bottles.nef, fac=fac, r=0.5, zoom.x=0.1, zoom.y=0.1)
pca(bottles.nef, PCa=2, PCb=3)
pca(bottles.nef, shp.nb=5)
pca(bottles.nef, shp.nb=5, shp.col="#FF660033", shp.border="#FF6600")

### pca3
pca3(bottles.nef)
pca3(bottles.nef, fac=fac)
pca3(bottles.nef, fac=fac, pch=c(4,5))
pca3(bottles.nef, fac=fac, pch=c(4,5), lty=c(2,3))
pca3(bottles.nef, fac=fac, pch=c(4,5), lty=c(2,3), col=c("dodgerblue","firebrick"))
pca3(bottles.nef, fac=fac, pch=c(4,5), lty=c(2,3), col=c("dodgerblue","firebrick"), legend=T)

### pca.tps
pca.tps(bottles.nef)
pca.tps(bottles.nef, fac=fac)
}
}
\keyword{Nef methods}
