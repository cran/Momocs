\name{Utilities: coo.oscillo}
\alias{coo.oscillo}
\alias{coo.oscillo1}
\title{
Momocs' "oscilloscope" for periodic functions.
}
\description{
Shape analysis deals with curve fitting, whether \eqn{x(t)} and \eqn{y(t)} positions along the curvilinear abscissa or radius/tangent angle variation. We may need to represent these single or double periodic functions that are ajusted by Fourier-based method. \code{coo.oscillo} and \code{coo.oscillo1} compute and provide standardized plot when given a matrix of coordinates or a vector, respectively.
}
\usage{
coo.oscillo(coo, method = c("d0", "di")[1], plot = TRUE, rug = TRUE, 
    legend = TRUE, cols = col.gallus(2), ref=FALSE, ref.nb=8, ...) 

coo.oscillo1(coo, method = c("d0", "di")[1], plot = TRUE, rug = TRUE, 
    legend = TRUE, cols = col.gallus(1), 
    xlab = "Points sampled along the outline",  ylab = "Deviation", ...) 
}
\arguments{
	\item{coo}{A list or a matrix of coordinates.}
	\item{method}{\code{character}. Whether to calculate differences with the first point (\code{"d0"}) or the previous (\code{"di"}), ie the derivate.}
	\item{plot}{Whether to plot the results.}
	\item{rug}{\code{logical}. Whether to display a pseudo rug, that indicate if the derivate is positive.}
	\item{legend}{\code{logical}. Whether to add a legend.}
	\item{cols}{A \code{vector} of two (for \code{coo.oscillo}) or a single (\code{coo.oscillo1}) color for lines.}
	\item{ref}{\code{logical}. Whether to display the original shape besides the oscillo.}
	\item{ref.nb}{\code{integer}. The number or reference points, sampled equidistantly along the curvilinear abscissa and added on the oscillo curves.}
	\item{xlab}{\code{character}. Alternate label for x-axis.}
	\item{ylab}{\code{character}. Alternate label for y-axis.}
	\item{...}{Additional parameters than can be passed to lines.}
}
\value{
Returns a \code{list} with two or one component(s), giving the difference calculated.
}

\examples{
data(bot)
coo.oscillo(bot@coo[[1]], lty=2)
coo.oscillo1(tfourier(bot@coo[[1]])$phi)
}

\keyword{Utilities}