% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/palinsol_helper.R
\name{obliquity}
\alias{obliquity}
\title{Compute June 1, 65 N obliquity}
\usage{
obliquity(time = seq(-2e+06, 0, 1000), model = "la04")
}
\arguments{
\item{time}{vector of years at which to calculate obliquity.  Negative for past times.}

\item{model}{one of la04, ber78 or ber90 to select orbital calculation, defaults to la04 (Laskar 2004)}
}
\description{
Compute June 1, 65 N obliquity with one of three orbital calculations.
}
\examples{
o <- obliquity( time = seq( -8e5, 0, 1e3 ), model = "ber78" )

}
