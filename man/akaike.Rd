% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ttools.R
\name{akaike}
\alias{akaike}
\title{Akaike 2nd-order window}
\usage{
akaike(n)
}
\arguments{
\item{n}{length of the filter; number of coefficients to generate.}
}
\value{
Filter coefficients
}
\description{
Akaike's 2nd-order windowing function
}
\examples{
n <- 100
plot( akaike( n ), type = "l" )

}
