% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ttools.R
\name{pspec}
\alias{pspec}
\title{Power spectrum}
\usage{
pspec(s)
}
\arguments{
\item{s}{vector of time series samples}
}
\value{
Power of non-negative frequencies order DC to Nyquist.
}
\description{
Direct FFT computation of time series power spectrum
}
\examples{
s <- sin( 2 * pi * ( 1:256 ) * 0.1 )
PS <- pspec( s )
f <- seq( 0, 0.5, length.out = length( PS ) )
plot( f, PS, type = "l", xlab = "Frequency", ylab = "Power" )
}
