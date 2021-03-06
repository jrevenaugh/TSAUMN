% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lanczos.R
\name{lanczos_bp}
\alias{lanczos_bp}
\title{Create a band-pass Lanczos filter}
\usage{
lanczos_bp(n, low, high)
}
\arguments{
\item{n}{length of filter; must be odd}

\item{low}{transition frequency of the filter (Nyquist = 0.5)}

\item{high}{transition frequency of the filter (Nyquist = 0.5)}
}
\value{
non-negative time filter coefficients (filter is symmetric about 0). 
Frequencies outside of band (low, high) are attenuated.
}
\description{
Compute coefficients of a band-pass Lanczos filter
}
\examples{
f <- lanczos_bp( 7, 0.1, 0.4 )
print( f )

}
