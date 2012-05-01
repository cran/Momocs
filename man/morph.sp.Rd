\name{morph.sp}
\alias{morph.sp}
\docType{methods}
\alias{morph.sp-methods}
\alias{morph.sp,Nef-method}
\title{Plots the morphological space.}
\description{
Given a matrix of harmonic coefficients, a \code{Nef}-object, calculates and plots morphological space \emph{i.e.} reconstructed shapes using and distributed on the orthonormal set defined by Principal Component axes.
}
\usage{
morph.sp(Nef,
           PCa = 1, PCb = 2, nb.PCa = 5, nb.PCb = 6, fac = NA, 
           morph.sp.extend = 1, zoom.extend = 1.2, asp,
           pch = 20, shp.col = NA, shp.lwd = 1, shp.size, col = "grey40", 
		   ell = FALSE, r = 1, lwd = 1, title = "Morphological space")
}		   
\arguments{
	\item{Nef}{the \code{Nef} object}
	\item{PCa}{a \code{numeric} indicating the first Principal Component axis on which to reconstruct shape}
	\item{PCb}{a \code{numeric} indicating the second Pricipal Component axis on which to reconstruct shape}
	\item{nb.PCa}{a \code{numeric} indicating how many shape to reconstruct on the first PC axis considered}
	\item{nb.PCb}{a \code{numeric} indicating how many shape to reconstruct on the second PC axis considered}
	\item{fac}{a factor indicating the grouping desing}
	\item{morph.sp.extend}{\code{integer} how much to extend morphological space reconstruction beyond range on the first PC considered}
	\item{zoom.extend}{\code{integer} indicating how much to extend the graphical window}
	\item{asp}{\code{numeric} and \code{optionnal} indicating the \code{asp} of the plotting window}
	\item{pch}{\code{integer} or a character indicating the \code{pch} for each groups to plot}
	\item{shp.col}{\code{integer} or a character indicating the \code{col} of these shapes}
	\item{shp.lwd}{\code{numeric} indicating the \code{lwd} of these shapes borders}
	\item{shp.size}{\code{numeric} for fine-tuning of shapes size}
	\item{col}{\code{integer} or a character indicating the \code{col} for each confidence ellipse to plot}
	\item{ell}{\code{logical} indicating whether to draw confidence ellipses for every group}
	\item{r}{\code{numeric} indicating the number of standard deviation for confidence ellipses computation}
	\item{lwd}{\code{numeric} indicating the \code{lwd} for the confidence ellipses}
	\item{title}{\code{character} to change the title of the plot}
}
\examples{
data(bottles.nef)
fac <- factor(rep(c("beer", "whisky"), each=20))

morph.sp(bottles.nef)
morph.sp(bottles.nef, fac=fac, ell=TRUE)
morph.sp(bottles.nef, fac=fac, nb.PCa=10, nb.PCb=10, ell=TRUE)
morph.sp(bottles.nef, PCa=2, PCb=3)

}
\keyword{Nef methods}


