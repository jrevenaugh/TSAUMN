% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p2fc.R
\name{p2fc}
\alias{p2fc}
\title{Convert corner periods to fc specification for cosfilt}
\usage{
p2fc(deltat = 1, pc)
}
\arguments{
\item{deltat}{sampling interval}

\item{pc}{vector of corner periods}
}
\description{
cosfilt assumes a frequency range of 0 (DC) to 0.5 (Nyquist).  p2fc converts native
unit periods to cosfilt corner frequencies.
}
\examples{
t <- seq( 0, 500, 2 )
x <- rnorm( length( t ) ) + 100 * sin( 2 * pi * t / 20 )
x.ts <- ts( x, deltat = 2 )
pc <- c( 4, 10, 30, 40 )
fc <- p2fc( 2, pc )
f.ts <- cosfilt( x.ts, fc, fw = c( 0, 1, 1, 0 ) )
plot( x.ts )
lines( f.ts, col = "red" )
}
