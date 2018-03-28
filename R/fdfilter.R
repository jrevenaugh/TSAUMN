#' Finite difference stencil filter
#'
#' Compute a finite-difference "stencil" for numerical differentiation  
#' @param n total length of filter.  3 <= n < 19 and odd.
#' @param d order of the derivative.  Must be an integer and satisfy
#' 1 <= d < \code{n}
#' @return Derivative filter coefficients.  Use \code{\link[stats]{filter}}
#' two-sided filter.
#' @export
#' @examples
#' f <- fdfilter( 5, 1 )
#' time <- 1:100
#' period <- 20
#' phase <- 2 * pi * time / period
#' s <- sin( phase )
#' Ds <- stats::filter( s, f )
#' plot( Ds, type = "l" )
#' lines( 2 * pi / period * cos( phase ), col = "red" )
#' 
fdfilter <- function( n = 3, d = 1 ) {
  if ( ( n - 1 ) %% 2 || n < 3 ) stop( "Filter length must be odd and >= 3." )
  l <- ( n - 1 ) / 2
  s <- c( -l:l )
  if ( floor( d ) != d || d < 1 ) stop( "Derivative order must be integer >= 1.")
  M <- matrix( rep( 0, n^2 ), nrow = n )
  for ( i in 0:( n - 1 ) ) M[i+1,] <- s^i
  x <- rep( 0, n )
  x[d + 1] <- factorial( d )
  f <- solve( M, x )
  return( rev( f ) )
}