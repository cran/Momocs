\name{hpow}
\alias{hpow}
\alias{hpow-methods}
\alias{hpow,Coo-method}

\title{
Calculates harmonic power.
}
\description{
\code{hpow} is used to estimate the number of harmonics required for the three Fourier methods implemented so far in Momocs: elliptical Fourier analysis (see \link{efourier}), radii variation analysis (see \link{tfourier}) and tangent angle analysis (see \link{tfourier}). It returns and can plot cumulated harmonic power whether dropping the first harmonic or not. 
}
\usage{
hpow(Coo, method = c("efourier", "rfourier", "tfourier"),
                 id = 1:Coo@coo.nb, probs = c(0, 0.5, 1), nb.h = 24,
                 drop = 1, smooth.it = 0, plot = TRUE, legend = FALSE,
                 title = "Fourier power spectrum", 
                 lineat.y = c(0.9, 0.95, 0.99, 0.999), bw = 0.1)
}
\arguments{
	\item{Coo}{The \code{Coo} object}
	\item{method}{A character, either "efourier", "rfourier" or "tfourier" (partial matches are allowed) to use to calculate morphological that indicates which method to use.}
	\item{id}{\code{integer}. The \code{id} of the shape to display. A range of \code{id}s can be passed to \code{harm.pow}}
	\item{probs}{A \code{vector} of \code{numeric}, to define quantiles to calculate ; see \link{quantile}}
	\item{nb.h}{\code{integer}. The maximal number of harmonics to calculate.}
	\item{drop}{\code{logical}. Whether to drop the first harmonic for plotting and power calculation.}
	\item{smooth.it}{\code{integer}. The number of smoothing iteration to perform.}
	\item{plot}{\code{logical}. Whether to plot or not the shape. If \code{FALSE}, only the results are returned.}
	\item{legend}{\code{logical}. Whether to display a legend box.}
	\item{title}{\code{character}. The title to add.}
	\item{lineat.y}{A \code{vector} of \code{numeric} to specify where to plot dashed lines on the y-axis.}
	\item{bw}{\code{numeric}. The width of horizontal segments drawn for each harmonic.}
}

\details{
The power of a given harmonic \eqn{n} is calculated as follows for elliptical Fourier analysis:
\deqn{
HarmonicPower_n= \frac{A^2_n+B^2_n+C^2_n+D^2_n}{2}
}
and as follows for radii variation and tangent angle:
\deqn{
HarmonicPower_n= \frac{A^2_n+B^2_n+C^2_n+D^2_n}{2}
}
}

\value{
Returns a matrix containing cumulated harmonic power for each harmonic. 
}

\examples{
data(bot)
hpow(bot)
}
\keyword{Calibration}