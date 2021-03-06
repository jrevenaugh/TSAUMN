% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phunwrap.R
\name{phunwrap}
\alias{phunwrap}
\title{Phase unwrapping}
\usage{
phunwrap(p)
}
\arguments{
\item{p}{vector of phase values}
}
\value{
unwrapped phase values
}
\description{
Unwrap a vector of phase values using Itoh's algorithm.
}
\examples{
s <- rep( 0, 124 )
s[50] <- 1
S <- fft( s )
p <- Arg( S )
up <- phunwrap( p[1:63] )
f <- seq( 0, 0.5, length.out = 63 )
plot( f, up, xlab = "Frequency", ylab = "Phase (rad)", type = "l", main = "Estimated (black) and Theory (red)" )
lines( f, -2 * pi * f * 49, col = "red" )

}
