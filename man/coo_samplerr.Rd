% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coo-utilities.R
\name{coo_samplerr}
\alias{coo_samplerr}
\title{Samples coordinates (regular radius)}
\usage{
coo_samplerr(coo, n)
}
\arguments{
\item{coo}{\code{matrix} of \verb{(x; y)} coordinates or any \link{Coo} object.}

\item{n}{\code{integer}, the number of points to sample.}
}
\value{
a \code{matrix} of (x; y) coordinates or a Coo object.
}
\description{
Samples n coordinates with a regular angle.
}
\details{
By design, this function samples among existing points, so using
\link{coo_interpolate} prior to it may be useful to have
more homogeneous angles. See examples.
}
\examples{
stack(bot)
bot <- coo_center(bot)
stack(coo_samplerr(bot, 12))
coo_plot(bot[1])
coo_plot(rr <- coo_samplerr(bot[1], 12))
cpos <- coo_centpos(bot[1])
segments(cpos[1], cpos[2], rr[, 1], rr[, 2])

# Sometimes, interpolating may be useful:
shp <- hearts[1] \%>\% coo_center

# given a shp, draw segments from each points on it, to its centroid
draw_rads <- function(shp, ...){
 segments(shp[, 1], shp[, 2], coo_centpos(shp)[1], coo_centpos(shp)[2], ...)
}

# calculate the sd of argument difference in successive points,
# in other words a proxy for the homogeneity of angles
sd_theta_diff <- function(shp)
   shp \%>\% complex(real=.[, 1], imaginary=.[, 2]) \%>\%
   Arg \%>\% `[`(-1) \%>\% diff \%>\% sd

# no interpolation: all points are sampled from existing points but the
# angles are not equal
shp \%>\% coo_plot(points=TRUE, main="no interpolation")
shp \%>\% coo_samplerr(64) \%T>\% draw_rads(col="red") \%>\% sd_theta_diff
# with interpolation: much more homogeneous angles
shp \%>\% coo_plot(points=TRUE)
shp \%>\% coo_interpolate(360) \%>\% coo_samplerr(64) \%T>\% draw_rads(col="blue") \%>\% sd_theta_diff
}
\seealso{
Other sampling functions: 
\code{\link{coo_extract}()},
\code{\link{coo_interpolate}()},
\code{\link{coo_sample_prop}()},
\code{\link{coo_sample}()}

Other coo_ utilities: 
\code{\link{coo_aligncalliper}()},
\code{\link{coo_alignminradius}()},
\code{\link{coo_alignxax}()},
\code{\link{coo_align}()},
\code{\link{coo_baseline}()},
\code{\link{coo_bookstein}()},
\code{\link{coo_boundingbox}()},
\code{\link{coo_calliper}()},
\code{\link{coo_centdist}()},
\code{\link{coo_center}()},
\code{\link{coo_centpos}()},
\code{\link{coo_close}()},
\code{\link{coo_down}()},
\code{\link{coo_dxy}()},
\code{\link{coo_extract}()},
\code{\link{coo_flipx}()},
\code{\link{coo_force2close}()},
\code{\link{coo_interpolate}()},
\code{\link{coo_is_closed}()},
\code{\link{coo_jitter}()},
\code{\link{coo_left}()},
\code{\link{coo_likely_clockwise}()},
\code{\link{coo_nb}()},
\code{\link{coo_perim}()},
\code{\link{coo_range}()},
\code{\link{coo_rev}()},
\code{\link{coo_right}()},
\code{\link{coo_rotatecenter}()},
\code{\link{coo_rotate}()},
\code{\link{coo_sample_prop}()},
\code{\link{coo_sample}()},
\code{\link{coo_scale}()},
\code{\link{coo_shearx}()},
\code{\link{coo_slice}()},
\code{\link{coo_slidedirection}()},
\code{\link{coo_slidegap}()},
\code{\link{coo_slide}()},
\code{\link{coo_smoothcurve}()},
\code{\link{coo_smooth}()},
\code{\link{coo_template}()},
\code{\link{coo_trans}()},
\code{\link{coo_trimbottom}()},
\code{\link{coo_trimtop}()},
\code{\link{coo_trim}()},
\code{\link{coo_untiltx}()},
\code{\link{coo_up}()},
\code{\link{is_equallyspacedradii}()}
}
\concept{coo_ utilities}
\concept{sampling functions}
