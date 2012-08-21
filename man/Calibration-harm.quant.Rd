\name{Calibration: harm.quant}
\alias{harm.quant}
\alias{harm.quant-methods}
\alias{harm.quant,Coo-method}

\alias{harm.qual}
\alias{harm.qual-methods}
\alias{harm.qual,Coo-method}

\alias{harm.pow}
\alias{harm.pow-methods}
\alias{harm.pow,Coo-method}

\alias{harm.contrib}
\alias{harm.contrib-methods}
\alias{harm.contrib,Nef-method}
\title{
Calculate harmonic effect.
}
\description{Provide facilities for choosing the number of harmonics to retain in elliptical Fourier analysis.
}
\usage{
harm.quant(Coo, id = 1:Coo@coo.nb, smooth.it = 0, harm.range= seq(8, 32, 6),
           scale = FALSE, center = TRUE, align = FALSE,
           plot = TRUE, legend = TRUE, palette = col.summer, lineat.y=c(1, 5, 10))

harm.qual(Coo, id = 1, smooth.it = 0, harm.range= c(1, 2, 4, 8, 16, 32),
           scale = FALSE, center = TRUE, align = FALSE, method = c("stack", "panel")[1],
		   legend = TRUE,
           palette = col.summer, shp.col="#70809033", shp.border="#708090EE")		   

harm.pow(Coo, id=1:Coo@coo.nb, probs=c(0, 0.5, 1), nb.h = 24,
           drop   = 1, smooth.it = 0, plot = TRUE, 
           legend = TRUE, title="Fourier power spectrum",
           lineat.x=seq(0, nb.h, by=6), lineat.y=c(0.9, 0.99), bw=0.1)
}
\arguments{
  \item{Coo}{The \code{Coo} object}
	\item{id}{\code{integer}. The \code{id} of the shape to display. A range of \code{id}s can be passed to \code{harm.pow}}
	\item{probs}{A \code{vector} of \code{numeric}, to define quantiles to calculate ; see \link{quantile}}
	\item{smooth.it}{\code{integer}. The number of smoothing iteration to perform.}
	\item{harm.range}{A \code{vector} of \code{integer} giving the harmonic range to calculate. See \code{nb.h} for \code{harm.pow}.}
	\item{nb.h}{\code{integer}. The maximal number of harmonics to calculate.}
	\item{scale}{\code{logical}. Whether to scale or not the shape.}
	\item{center}{\code{logical}. Whether to center or not the shape.}
	\item{align}{\code{logical}. Whether to align or not the shape.}
	\item{plot}{\code{logical}. Whether to plot or not the shape. If \code{FALSE}, only the results are returned.}
	\item{method}{A plotting method for \code{harm.qual}. If \code{"stack"} outlines are plotted on the same same, if \code{"panel"}, separate reconstructions are plotted.}
	\item{legend}{\code{logical}. Whether to display a legend box.}
	\item{title}{\code{character}. The title to add.}
	\item{drop}{\code{logical}. Whether to drop the first harmonic for plotting and power calculation.}
	\item{palette}{A color palette such those included in Momocs or produced with \link{colorRampPalette}.}
	\item{shp.col}{A color for the shape body.}
	\item{shp.border}{A color for the shape border.}
	\item{lineat.x}{A \code{vector} of \code{numeric} to specify where to plot dashed lines on the x-axis.}
	\item{lineat.y}{A \code{vector} of \code{numeric} to specify where to plot dashed lines on the y-axis.}
	\item{bw}{\code{numeric}. The width of horizontal segments drawn for each harmonic.}
}
\details{
So far, these methods can only be used for elliptical Fourier analysis, not for tangent angle or radius variation pespectives. They will be included in further versions. \code{harm.quant} is based on euclidean distance between original and reconstructed outlines harm.pow returns and plot cumulated harmonic power. The power of a given harmonic \eqn{n} is calculated as follows:
\deqn{
HarmonicPower_n= \frac{A^2_n+B^2_n+C^2_n+D^2_n}{2}
}
}

\value{
\code{harm.quant} returns a matrix containing deviations for each harmonic and corresponding quantiles.
\code{harm.pow} returns a matrix containing cumulated harmonic power for each harmonic. 
}

\examples{
data(bot)
harm.quant(bot)
harm.qual(bot)
harm.pow(bot)
}
\keyword{Calibration}