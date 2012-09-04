\name{Utilities: coo-family}
\alias{coo.align}
\alias{coo.center}
\alias{coo.close}
\alias{coo.rotate}
\alias{coo.sample}
\alias{coo.sample.rr}
\alias{coo.scale}
\alias{coo.slide}
\alias{coo.smooth}
\alias{coo.template}
\alias{coo.trans}
\alias{l2m}
\alias{m2l}
\alias{l2a}

\title{
Utilities to manipulate list and matrices of coordinates.
}

\description{Some utilities to manipulate \eqn{(x; y)} coordinates when arranged as matrices or lists. Most of them are used internally but they can be used directly.
\itemize{
\item \code{coo.align} aligns \code{coo} using its best fitting ellipse.
\item \code{coo.center} centers \code{coo} on the origin.
\item \code{coo.close} closes \code{coo}, \emph{i.e.} makes the last coordinates be the first than the first.
\item \code{coo.rotate} rotates \code{coo} with a \code{theta} angle.
\item \code{coo.sample} samples \code{n} points in \code{coo}.
\item \code{coo.sample.rr} does the same but based on polar coordinates.
\item \code{coo.scale} "scales" \code{coo}, \emph{i.e.} resizes \code{coo} so that it can be included in a square of \code{scale} side.
\item \code{coo.slide} "slides" \code{coo}, \emph{i.e.} makes the \code{id1} coordinate to be the first and change the others accordingly.
\item \code{coo.smooth} performs \code{n} smoothing iteration on \code{coo}. \code{coo.template} returns \code{coo} so that the shape it is centered on the origin and inscribed in a size-side square, also centered on the origin; see \link{coo.list.panel} for an illustration of this function.
\item \code{coo.trans} translates \code{coo} by \code{x} and {y} on the two dimensions, respectively.
\item \code{m2l} converts a matrix of \code{(x; y)}coordinates to a list with \code{x; y} components.
\item \code{l2m} does the exact opposite than \code{m2l}.
\item \code{l2a} converts a list of length \code{k} with identical \code{m x n} matrices to a \code{m x n x k} array.
}
}

\usage{
coo.align(coo)
coo.center(coo)
coo.close(coo)
coo.rotate(coo, theta)
coo.sample(coo, n)
coo.sample.rr(coo, n)
coo.scale(coo, scale)
coo.slide(coo, id1)
coo.smooth(coo, n)
coo.template(coo, size)
coo.trans(coo, x, y)
m2l(m)
l2m(l)	
l2a(l)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  \item{theta}{\code{numeric}. The angle to rotate the shape.}
  \item{n}{\code{integer}. The number of points to sample.}
  \item{scale}{\code{numeric}. Indicates the scaling factor of the shape.}
  \item{size}{\code{numeric}. Indicates the length of the side "inscribing" the shape.}
  \item{id1}{\code{integer}. Specifies the index of the coordinates that has to be defined as the first coordinate.}
  \item{x}{\code{numeric}. The x-axis translation.}
  \item{y}{\code{numeric}. The y-axis translation.}
  \item{m}{A 2-columns matrix containing x and y coordinates.}
  \item{l}{A \code{list} with x and y coordinates as components.}
  }

\details{
\code{coo.scale} simply returns coordinates so that the longer axis of the shape is of size "scale".
\code{coo.template} returns coordinates so that the shape is centered on the origin and is completely inscribed in a size-side square, also centered on the origin. 
}

\value{
A \code{list} with \code{x; y} components or a a matrix of \code{(x; y)}coordinates.
}

\seealso{
\link{Coo-class}, \link{coo.list.panel}.
}
\examples{
data(bot)
b <- bot@coo[[1]]
coo.plot(b)
coo.align(b)               # a bottle aligned along its longer axis
coo.center(b)              # a bottle centered on the origin
coo.close(b)               # a "closed" bottle
coo.sample(b, 60)          # 60 points sampled along the curvilinear abscissa
coo.sample.rr(b, 60)$coord # 60 points sampled along the curvilinear abscissa
coo.trans(b, x=50, y=-200) # translate bottle (x: +50 ; y: -200)
coo.template(b, 5)         # a "templated" bottle
coo.scale(b, 5)			   # a "scaled" bottle
}
\keyword{Utilities}