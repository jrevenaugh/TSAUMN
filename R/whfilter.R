#' Apply derived FIR Causal Wiener-Hopf Filter
#'
#' Apply the optimal (MMSE) Wiener-Hopf causal FIR filter to a known multivariate
#' input ts x.  Uses filter output from \code{whfir}.
#' @param x a ts, vector or matrix of values (columns as time series).
#' @param s as ts or vector of values
#' @param h list returned from \code{whfir}.
#' @return Wiener-Hopf prediction of time series s.
#' @import "ltsa"
#' @export
#' @examples
#' x <- rnorm( 504 )
#' h <- seq( 1, 0, length.out = 5 )
#' s <- stats::filter( x, h, sides = 1 )
#' x <- x[-(1:4)] # Remove filter roll-on
#' s <- s[-(1:4)]
#' wh <- whfir( x, s, n = c(1, 10 ) )
#' spred <- whfilt( x, s, wh )
#' plot( s, type = "l" )
#' lines( spred, col = "red" )
#'
whfilt <- function( x, s, h ) {
  if ( is.null( x ) ) stop( "Null x" )

  multi <- is.matrix( x )
  if ( !multi ) {
    x <- matrix( x, ncol = 1 )
  }
  lx <- dim( x )[1]
  m <- dim( x )[2]
  hest <- h$h
  n <- dim( hest )[1]

  sest <- stats::filter( x[,1], hest[1:n], sides = 1 )
  fs <- seq( 1, length( hest ), n )
  fe <- seq( n, length( hest ), n )
  if ( m > 1 ) for ( j in 2:m ) sest <- sest + stats::filter( x[,j], hest[fs[j]:fe[j]], sides = 1 )
  return( ts( sest, start = time( s )[1], deltat = deltat( s ) ) )
}
