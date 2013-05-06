\name{clust}
\alias{clust}
\alias{clust-methods}
\alias{clust,Coe-method}

\title{Hierarchical clustering on a matrix of coefficients}

\description{
\code{clust} is a method for \code{Coe} objects to calculate and display hierarchical clustering of shapes based on the distance between their coefficients, after a distance matrix calculation based on \link{dist}.
}

\usage{
clust(Coe, fac, method = "euclidean", type = "fan", palette = col.summer)
}

\arguments{
	\item{Coe}{The \code{Coe} object on which to perform hierarchical clustering.}
	\item{fac}{\code{factor} defining which column of the @fac slot to use for grouping, if any.}
	\item{method}{The distance method used by \link{dist}.}
	\item{type}{The plotting method used by \link{plot.phylo}.}
	\item{palette}{A color palette to use to plot the different groups, if any.}
}

\value{
Returns a \link{dist} object.}

\seealso{
\link{dist}, \link{hclust} and \link{plot.phylo}.}

\examples{
data(bot)
botE <- eFourier(bot, nb.h=32)
clust(botE)
}

\keyword{Multivariate Analysis}

