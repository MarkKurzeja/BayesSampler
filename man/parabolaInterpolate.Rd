% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayes_main.R
\name{parabolaInterpolate}
\alias{parabolaInterpolate}
\title{Interpolate a function with a parabola}
\usage{
parabolaInterpolate(x, xvec, yvec)
}
\arguments{
\item{x}{The input value to the function}

\item{xvec}{A vector of length three of the x values}

\item{yvec}{A vector of length three of the y values}
}
\description{
Take in a vector of x values and a vector of y values for a parabola, three in total, and return the parabola approximation at the value of x
}
\examples{
 
parabolaInterpolate(0.5, c(1,0,4), c(5,6,2))
}
\keyword{interpolate}
\keyword{interpolation}
\keyword{parabola}

