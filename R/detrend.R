#' Detrend a time series with a low-order polynomial.
#'
#' Remove low-order polynomial "trend" from a time series.  
#' @param x time series or vector of data values
#' @param n polynomial order
#' @return ts object with low-order LS polynomial removed.
#' @export
#' @examples
#' n <- 100
#' x <- arimaSim( n, ar = 0.9, int = 1 ) + 1:n / 20
#' plot( x )
#' d <- detrend( x, 1 )
#' plot( d )
#' 
detrend <- function( x, n = 1 ) {
  l <- length( x )
  r <- residuals( lm( x ~ stats::poly( 1:l, n ) ) )
  if ( is.ts( x ) ) return( ts( r, start = time( x )[1], deltat = deltat( x ) ) )
  return( r )
}