% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lanczos.R
\name{lanczos_lp}
\alias{lanczos_lp}
\title{Create a low-pass Lanczos filter}
\usage{
lanczos_lp(n, low)
}
\arguments{
\item{n}{length of filter; must be odd}

\item{low}{cutoff frequency of the filter (Nyquist = 0.5)}
}
\value{
non-negative time filter coefficients (filter is symmetric about 0).
}
\description{
Compute coefficients of a low-pass Lanczos filter
}
\examples{
f <- lanczos_lp( 7, 0.1 )
print( f )

}
