\name{coo.list.panel}
\alias{coo.list.panel}
\title{
Plots sets of shapes.
}
\description{
\code{coo.list.panel} plots a list of shapes if passed with a list of coordinates. Outlines are templated and on the same graphical window with the help of \link{coo.template}.
}
\usage{
coo.list.panel(coo.list, dim,
               byrow = TRUE, fromtop = TRUE,
               mar = rep(0, 4), cols, borders, density = NULL, angle = 45)
}

\arguments{
	\item{coo.list}{A \code{list} of coordinates, such as those in the \code{@coo} slot of \code{Coo} objects.}
	\item{dim}{A \code{vector} of the form \code{(nb.row, nb.cols)} to specify the panel display. If missing, shapes are arranged in a square.}
	\item{byrow}{\code{logical}. Whether to succesive shape by row or by col.}
	\item{fromtop}{\code{logical}. Whether to display shapes from the top of the plotting region.}
	\item{mar}{A \code{vector} to define margins.}
	\item{cols}{A \code{vector} of colors to fill shapes.}
	\item{borders}{A \code{vector} of colors to draw shape borders.}
	\item{density}{A \code{vector} for density of shading lines. See \link{polygon}}
	\item{angle}{A \code{vector} for shading lines.}
}

\value{
Returns (invisibly) a \code{data.frame} with position of shapes that can be used for other sophisticated plotting design.
}

\seealso{
\link{coo.plot} and \link{coo.template}.
}

\examples{
data(bot)
coo.list.panel(bot@coo)
x <- coo.list.panel(bot@coo)
x # positions of shapes returned invisibly 
# axis(1) ; axis(2) # that's a single graphical window ;)
}
\keyword{Utilities}