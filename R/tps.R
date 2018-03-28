#' Compute theorectical AR(p) power spectrum
#'
#' Compute theorectical AR(p) power spectrum.  
#' @param p process order (defaults to 1)
#' @param a single coefficient (p = 1) or vector of coefficients, a1, a2, ..., ap.  
#' No check is made to insure stationarity.
#' @param freq either a scalar specifying how many frequencies between 0 and 0.5 to return
#' or a vector of frequencies. 
#' @return AR(p) power spectral density
#' @export
#' @examples
#' f <- seq( 0, 0.5, 0.01 )
#' s <- arps( p = 1, a = 0.5, freq = f )
#' plot( f, s, type = "l" )

# Compute theoretical power-spectrum of a AR(p) model.  
arps <- function( p = 1, a = c( 0.9 ), freq = 100 ) {
  p <- floor( p )
  if ( p <= 0 ) stop( "P must be 1 or greater." )
  if ( length(a) != p ) stop( "length( a ) must equal p." )
  if ( length( freq ) > 1 ) {
    f <- freq
    nf <- length( f )
  } else {
    freq <- floor( ( abs( freq ) ) )
    if ( freq == 0 ) {
      freq = 100
      warning("100 substituted for nfreq.")
    }
    nf <- freq
    f <- seq( 0, 0.5, length.out = nf )
  }
  z <- exp( -2 * pi * 1i * f )
  I <- sapply( z, function( x, a, p ) { ( Mod( 1.0 - a %*% x^(1:p) )^-2 ) }, a, p )
  return( I )
}

#'
#' Compute theorectical MA(q) power spectrum.  
#' @param a process order (defaults to 1)
#' @param b single coefficient (a = 1) or vector of coefficients, b1, b2, ..., bq.  
#' No check is made to insure stationarity.
#' @param freq either a scalar specifying how many frequencies between 0 and 0.5 to return
#' or a vector of frequencies. 
#' @return MA(q) power spectral density
#' @export
#' @examples
#' f <- seq( 0, 0.5, 0.01 )
#' s <- maps( p = 1, a = 0.5, freq = f )
#' plot( f, s, type = "l" )

maps <- function( q = 1, b = c( 0.9 ), freq = 100 ) {
  q <- floor( q )
  if ( q <= 0 ) stop( "q must be 1 or greater." )
  if ( length(b) != q ) stop( "length( b ) must equal q." )
  if ( length( freq ) > 1 ) {
    f <- freq
    nf <- length( f )
  } else {
    freq <- floor( ( abs( freq ) ) )
    if ( freq == 0 ) {
      freq = 100
      warning("100 substituted for nfreq.")
    }
    nf <- freq
    f <- seq( 0, 0.5, length.out = nf )
  }
  
  z <- exp( -2 * pi * 1i * f )
  I <- sapply( z, function( x, b, q ) { Mod( 1.0 + b %*% x^(1:q) )^2 }, b, q )
  return( I )
}