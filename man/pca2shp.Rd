\name{pca2shp}
\alias{pca2shp}

\title{
From a factorial map to a shape.
}
\description{Converts a position in a factorial map to a shape. This function is used internally but might be interesting for a direct used on a PCA calculated on a matrix of harmonic coefficients.
}
\usage{
pca2shp(pos, rot, mean.shp, method = c("efourier", "rfourier",
                 "tfourier"), scale = 1, amp = 1, trans = TRUE, nb.pts
                 = 64, rotate.shp)
}

\arguments{
  \item{pos}{A \code{vector} of positions on a given PC axis.}
  \item{rot}{A \code{vector} of multiplying values corresponding to the PC axis passed to \code{pos}.}
  \item{mean.shp}{A \code{vector} of harmonic coefficients corresponding to the mean shape.}
  \item{method}{A character, either "efourier", "rfourier" or "tfourier" (partial matches are allowed) to use to calculate morphological space.}
  \item{scale}{\code{numeric}. To scale the shape returned.}
  \item{amp}{\code{numeric}. To amplify the amplitude of the deformation.}
  \item{trans}{\code{logical}. Whether to translate of not the coordinate of the shape reconstructed according to \code{pos}.}
  \item{nb.pts}{numeirc. The number of points to use for reconstructed shapes.}
  \item{rotate.shp}{\code{numeric}. If specified, the angles (in radians) to rotate shapes plotted.}
  }
\value{
Returns a \code{list} of coordinates.
}

\examples{
data(bot)
}
\keyword{Utilities}
