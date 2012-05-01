\name{traj}
\alias{traj}
\docType{methods}
\alias{traj-methods}
\alias{traj,Nef-method}
\title{Calculates shape intermediates.}
\description{
Given a \code{Nef} object, and two positions on the set defined by \eqn{PC_1} and \eqn{PC_2}, calculates and plots intermediate shapes along the euclidean distance between these two shapes.
}
\usage{
traj(Nef, fr = c(0, 0), to = c(1, 1), nb.int = 50, nb.pts= 500, save = FALSE, prog = TRUE, pause=TRUE)
}
\arguments{
	\item{Nef}{the \code{Nef} object}
	\item{fr}{a vector with two \code{numeric}s indicating the \eqn{(x; y)} coordinates of the starting point. If not provided, \code{locator(1) is called}}
	\item{to}{a vector with two \code{numeric}s indicating the \eqn{(x; y)} coordinates of the ending point. If not provided, \code{locator(1) is called}}
	\item{nb.int}{code{integer} giving the number of shape intermediates to calculate}
	\item{nb.pts}{code{integer} giving the number of points of the reconstructed outlines}
	\item{save}{\code{logical} indicating whether to save the images in a dedicated folder}
	\item{prog}{\code{logical} indicating whether to plot or not a progression bar}
	\item{pause}{\code{logical} indicating whether to ask the user to display successive intermediate shapes}
}
\examples{
\dontrun{
data(bottles.nef)
traj(bottles.nef)
traj(bottles.nef, fr=c(-0.05, -0.05), to=c(0.15, 0.05))
}
}
\keyword{Nef methods}
