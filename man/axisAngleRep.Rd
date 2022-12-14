% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SOManifold.R
\name{axisAngleRep}
\alias{axisAngleRep}
\title{Returns the angle representation of SO(3) matrices
c.f. \url{https://en.wikipedia.org/wiki/Axis-angle_representation}}
\usage{
axisAngleRep(mfd, X)
}
\arguments{
\item{mfd}{A manifold object created by \code{\link{createM}}}

\item{X}{A matrix holding a vectorized SO(3) matrix in each column}
}
\value{
A matrix with 4 rows and the same number of columns as X. The first row contains the angles in rads (theta), and the last three rows correspond to the the axes with respect to which the rotations are performed.
}
\description{
Returns the angle representation of SO(3) matrices
c.f. \url{https://en.wikipedia.org/wiki/Axis-angle_representation}
}
