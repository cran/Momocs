\name{PC.contrib}
\alias{PC.contrib}

\title{
Shape variation along PC axis}

\description{
\code{PC.contrib} calculates and plots shape variation along Principal Component axes.
}

\usage{
PC.contrib(dudi, PC.r = 1:dudi$nf, sd = 2, cols = rep(NA, 3), 
    borders = c("#000080", "#000000", "#EE0000"), lwd = 1, nb.pts = 300, 
    plot = TRUE, legend = TRUE) 
}

\arguments{
  \item{dudi}{a dudi.pca object.}
  \item{PC.r}{A range of \code{integer}s indicating the PC axes on which to display shape variation.}
  \item{sd}{A \code{numeric} to indicate +/- the number of standard deviations to consider.}
  \item{cols}{A color string of length 3 to use to fill the shapes.}
  \item{borders}{A color string of length 3 to use for the shape borders.}
  \item{lwd}{A \code{numeric} to specify the \code{lwd} of the shapes.}
  \item{nb.pts}{\code{integer} to specify the number of points to draw the shapes.}
  \item{plot}{\code{logical}. Whether to plot the results.}
  \item{legend}{\code{logical}. Whether to add the legend.}
}

\value{
Invisibly returns a list that contains the coordiantes of the calculated shapes.
}

\examples{
data(bot)
#botF <- eFourier(bot)
#botD <- pca(botF)
#PC.contrib(botD)
#PC.contrib(botD, sd=1) # only one sd
#PC.contrib(botD, PC.r=1:3, sd=1, cols=paste(col.sari(3), "55", sep=""),
#   borders=rep("black", 3), legend=FALSE)
# only 3 PC axis and some cosmectics.
}

\keyword{Multivariate Analysis}
