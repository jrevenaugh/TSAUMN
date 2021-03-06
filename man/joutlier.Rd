% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/joutlier.R
\name{joutlier}
\alias{joutlier}
\title{Recursive outlier detector/remover}
\usage{
joutlier(s, order = c(2, 0, 0), k = 10, tol = 3)
}
\arguments{
\item{s}{timeseries}

\item{order}{order specification of arima model; defaults to AR(2)}

\item{k}{hampel's half-window size.  Defines window for MAD computation}

\item{tol}{MAD outlier threshold; defaults to 3 (Pearson's rule)}
}
\value{
cleaned timeseries.
}
\description{
Locate and remove outliers using an arima model to whiten time series, Hampel's
MAD criterion to identify outliers, and kalman smoothing.  Algorithm
is recursive and continues until no further outliers are found.
}
\examples{
y <- arimaSim( 200, ar = c( 0.5, -0.1 ) )
# Add a few outliers
o <- floor( runif( 5, min = 5, max = 195 ) )
y[o] <- runif( 5, min = 4, max = 8 )
yc <- joutlier( y )
plot( y )
lines( yc, col = "red" )

}
