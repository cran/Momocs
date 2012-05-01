\name{tps.grid}
\alias{tps.grid}
\alias{tps.iso}
\alias{tps.vf}
\docType{methods}
\alias{tps.grid-methods}
\alias{tps.grid,Nef-method}
\alias{tps.iso-methods}
\alias{tps.iso,Nef-method}
\alias{tps.vf-methods}
\alias{tps.vf,Nef-method}
\title{Thin Plate Splines deformation grids between two shapes.}
\description{
Passed with a \code{Nef}-object, and two positions on the set defined by \eqn{PC_1} and \eqn{PC_2}, calculates and plots deformation grids, "vector field" or iso-deformation lines between these two shapes.
}
\usage{
tps.grid(Nef, fr, to, nb.pts = 50, amp = 1, grid.size = 50,
		 grid.col = "grey40", cont = TRUE,
		 cont.col = c("dodgerblue3", "firebrick3"), cont.lwd = rep(3,2)) 
		  
tps.iso(Nef, fr, to, nb.pts = 200, amp = 1, iso.pts = 1000,
		col.pal = topo.colors, col.lev = 500,
		cont.to = TRUE, cont.fr = TRUE,  cont.lev = 10,
		cont.col = c("dodgerblue3", "firebrick3"),  cont.lwd = rep(3,2))
			
tps.vf(Nef, fr, to, nb.pts = 100, amp = 1, arr.nb = 300, 
		arr.len = 0.05, arr.ang = 30, arr.col = "grey40",
		arr.pal = FALSE, arr.palette = topo.colors, arr.pal.lev = 20,
		arr.lwd = 1, cont.col = c("dodgerblue3", "firebrick3"),
		cont.lwd = rep(3, 2))
}
\arguments{
	\item{Nef}{the \code{Nef} object}
	\item{fr}{a vector with two \code{numeric}s indicating the \eqn{(x; y)} coordinates of the starting point. If not provided, \code{locator(1)} is called}
	\item{to}{a vector with two \code{numeric}s indicating the \eqn{(x; y)} coordinates of the ending point. If not provided, \code{locator(1)} is called}
	\item{nb.pts}{\code{integer} indicating the number of points used to calculate deformation grids}
	\item{amp}{a \code{numeric} indicating magnifying factor for deformations}
	\item{cont}{\code{logical} indicating whether to plot original and deformed shapes}
	\item{cont.to}{\code{logical} indicating whether to plot the original shape}
	\item{cont.fr}{\code{logical} indicating whether to plot the deformed shape}
	\item{cont.col}{\code{integer} or \code{character} indicating the \code{col} of the two outlines compared}
	\item{cont.lwd}{a \code{numeric} indicating a vector containing the \code{lwd} of the two outlines compared}
	\item{grid.size}{a \code{numeric} indicating deformation grid size}
	\item{grid.col}{code{integer} or a character indicatinf the deformation grid color}
	\item{iso.pts}{a \code{integer} indicating the number of iso points to use for isolines calculation}
	\item{col.pal}{a color palette such as \code{rainbow}, \code{heat.colors} or such as build by \code{colorRampPalette} for isolines drawing}
	\item{col.lev}{\code{integer} indicating how many color levels to use}
	\item{cont.lev}{\code{integer} indicating how many isolines to calculate}
	\item{arr.nb}{\code{integer} indicating how many arrows to display}
	\item{arr.len}{a \code{numeric} indicating the arrows length}
	\item{arr.ang}{a \code{numeric} indicating the arrows angle}
	\item{arr.col}{\code{integer} or a character indicating the arrows color}
	\item{arr.pal}{\code{logical} indicating whether to use or not a color palette for drawing arrows}
	\item{arr.palette}{a color palette such as \code{rainbow}, \code{heat.colors} or such as build by \code{colorRampPalette}}
	\item{arr.pal.lev}{\code{integer} specifying how many levels to use for the color palette}
	\item{arr.lwd}{a \code{numeric} indicating the arrows \code{lwd}}	
	}
\examples{
data(bottles.nef)

\dontrun{
tps.grid(bottles.nef, fr=c(-0.05, -0.05), to=c(0.15, 0.05))
tps.vf(bottles.nef, fr=c(-0.05, -0.05), to=c(0.15, 0.05))
tps.iso(bottles.nef, fr=c(-0.05, -0.05), to=c(0.15, 0.05))
}
}
\keyword{Nef methods}
