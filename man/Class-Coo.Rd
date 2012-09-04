\name{Class: Coo}
\Rdversion{1.1}
\docType{class}
\alias{Coo}
\alias{Coo-class}
\alias{plot-methods}
\alias{plot,Coo-method}
\alias{show,Coo-method}
\alias{sample}
\alias{sample-methods}
\alias{sample,Coo-method}
\alias{scale}
\alias{scale-methods}
\alias{scale,Coo-method}
\alias{center}
\alias{center-methods}
\alias{center,Coo-method}
\alias{align}
\alias{align-methods}
\alias{align,Coo-method}
\alias{slide}
\alias{slide-methods}
\alias{slide,Coo-method}
\alias{smooth}
\alias{smooth-methods}
\alias{smooth,Coo-method}

\title{Class \code{"Coo"}}
\description{A class for handling morphometric datasets, \emph{i.e.} lists of \eqn{(x; y)} outline coordinates with or without a grouping structure.
}
\section{Objects from the Class}{
New \code{Coo}-objects can be created with \code{new("Coo", ...)} or more directly using the \code{Coo} builder.
}
\section{Slots}{
  \describe{
    \item{\code{coo}:}{\code{list} a list of shapes arranged in two columns coordinate matrices and that can have different lengths.}
    \item{\code{names}:}{\code{character} \code{vector} to name every shape.}
    \item{\code{fac}:}{\code{data.frame} defining the grouping structure, if any. Columns must be \code{factor}s.}
    \item{\code{ldk}:}{\code{list} defining the coordinates indices that will be considered as homologous points (not yet used in current version).}
    \item{\code{scale}:}{\code{numeric} to define the scale (not yet used in current version).}
    \item{\code{coo.nb}:}{\code{numeric}. The number of shapes.}
    \item{\code{coo.len}:}{a \code{numeric} \code{vector} that contains the number of coordinates per shape.}
    \item{\code{details}:}{any additional details, such as licence and authors can be added in this \code{list}.}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{Coo}-getters (will be restored in a further version)}
    \item{[<-}{\code{Coo}-setters (will be restored in a further version)}
    \item{align}{align the outlines.}
    \item{center}{center the outlines.}
    \item{eFourier}{perform an elliptic Fourier analysis.}
    \item{harm.pow}{calculates the harmonic power.}
    \item{harm.qual}{helps to visual calibration of the number of harmonics.}
    \item{harm.quant}{uses deviations to calibrate of the number of harmonics.}
    \item{plot}{plots a set of set.}
    \item{Ptolemy}{plots Ptolemaic ellipses, the idea behind elliptic Fourier analysis.}
    \item{sample}{sample a given number of points.}
    \item{scale}{scale the outlines.}
    \item{slide}{"slides" the outlines according to the additional argument \code{ldk.id=1}.}
    \item{smooth}{smoothes the coordinates.}    
	\item{show}{displays a summary.}
    \item{smooth.qual}{helps to estimate the number of smoothing iterations.}
	 }
}
\examples{
showClass("Coo")
# an example of the Coo builder
data(bot)
list_of_coordinates <- bot@coo
Coo(list_of_coordinates)
}
\keyword{Classes}