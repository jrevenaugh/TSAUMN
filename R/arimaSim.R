#' ARIMA process realization
#'
#' Produce realization of an ARIMA process using \code{\link{arima.sim}}.  
#' @param n number of output samples
#' @param ar vector of ar coefficients; defaults to NULL (ignored)
#' @param ma vector of ma coefficients; defaults to NULL (ignored)
#' @param int integration order, a non-negative integer; defaults to 0 (no integration)
#' @param ... other arguments passed to arima.sim.  sd is the most useful
#' @return ts object with n + int samples.
#' @export
#' @examples
#' a <- arimaSim( n = 100, ar = 0.9, int = 1 )
#' plot( a )
#' b <- arimaSim( n = 100, ma = c( 0.5, 0.1 ), ar = c( -0.2 ) )
#' plot( b )
#' 
arimaSim <- function( n = 64, ar = NULL, ma = NULL, int = 0, ... ) {
  # Create an order-structured model descriptor for arima.sim
  if ( int < 0 ) stop( "Integration order must be non-negative." )
  n_ar <- length( ar )
  n_ma <- length( ma )
  order <- c( n_ar, int, n_ma )
  return( arima.sim( list( order = c( n_ar, int, n_ma ), ar = ar, ma = ma ), n = n - int, ... ) )
}