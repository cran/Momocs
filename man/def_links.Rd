% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coo-ldk-utilities.R
\name{def_links}
\alias{def_links}
\title{Defines links between landmarks}
\usage{
def_links(x, nb.ldk)
}
\arguments{
\item{x}{Ldk, matric or array}

\item{nb.ldk}{numeric the iterative procedure is stopped when the
user click on the top of the graphical window.}
}
\value{
a Momocs object of same class
}
\description{
Works on Ldk objects, on 2-cols matrices, 3-dim arrays (\link{MSHAPES} turns it into a matrix).
}
\examples{
\dontrun{
wm <- MSHAPES(wings)
links <- def_links(wm, 3) # click to define pairs of landmarks
ldk_links(wm, links)
}
}
\seealso{
Other ldk helpers: 
\code{\link{ldk_check}()},
\code{\link{links_all}()},
\code{\link{links_delaunay}()}
}
\concept{ldk helpers}
