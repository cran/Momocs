\name{Abstract: Momocs}
\alias{Momocs-package}
\alias{Momocs}
\docType{package}
\title{
Outline Analysis using R.
}
\description{
Momocs is intended to ease and popularize shape analysis of outlines (especially using elliptical Fourier analysis). It mostly hinges on the functions developped in \emph{Morphometrics with R} (Claude, 2008). From outline extraction of images and elliptical Fourier calculation to multivariate analysis and the visualization of transformations within the morphological space, Momocs provides a complete and convenient toolkit to specialists within every field that are, or may be, interested in morphological comparisons of outlines.
}

\author{
\enumerate{
\item Vincent Bonhomme, French Institute of Pondicherry, India. <bonhomme.vincent@gmail.com>
\item Sandrine Picq, UMR CBAE, Montpellier, France.
\item Julien Claude UMR, ISEM, Universite de Montpellier II, France.
}}

\references{
The companion paper of Momcs is submitted.

You can have a look to Claude (2008) that greatly inspired Momocs: Claude, J. (2008) \emph{Morphometrics with R}, Use R! series, Springer 316 pp.
}

\examples{
# We load a dataset as a \code{Coo} object
data(bot)
bot
# Some calibration
harm.qual(bot, 1)
harm.quant(bot)
harm.pow(bot)
# The idea behing elliptic Fourier analysis
Ptolemy(bot, 1)
# We calculate elliptic Fourier analysis
botF <- eFourier(bot)
botF
# We explore the matrix of coefficients
boxplot(botF)
harm.contrib(botF)
# Some analysis now
manova.nef(botF, "type")
tps(botF, 1, 2, grid.size=50, amp=1.5)
# And one PCA
botD <- pca(botF)
botD
dudi.plot(botD, "type")
}
\keyword{Package}