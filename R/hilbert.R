#' Hilbert transform
#'
#' Compute hilbert transform using FFT
#' @param s vector of timeseries samples or ts object.
#' NA values are not allowed.
#' @return Complex valued hilbert transform of x.
#' @export
#' @examples
#' s <- sin( 2 * pi * 1:256 * 0.05 )
#' H <- hilbert( s )
#' plot( Re( H ), xlab = "Time", ylab = NULL, type = "l", ylim = c( -1.1, 1.1 ) )
#' lines( Im( H ), col = "red" )

hilbert <- function( s ) {
  n <- nextn( l <- length( s ), c( 2, 3, 5 ) )
  if ( n > l ) s <- c( s, rep( 0, n - l ) )
  S <- fft( s )
  h <- rep( 0, n )
  if ( n %% 2 == 0 ) {
    h[c( 1, n / 2 + 1 )] <- 1
    h[2:(n / 2 )] <- 2
  } else {
    h[1] <- 1
    h[2:( ( n + 1 ) / 2 )] <- 2
  }
  ht <- fft( S * h, inverse = TRUE ) / n
  return( ht[1:l] )
}

#' Envelope using analytic signal
#'
#' Compute envelope (instantaneous amplitude) of a time series using hilbert transform
#' @param s vector of timeseries samples or ts object.
#' NA values are not allowed.
#' @return Real valued (non-negative) envelope of x.
#' @export
#' @examples
#' s <- sin( 2 * pi * 1:256 * 0.05 )
#' S <- envelope( s )
#' plot( S, xlab = "Time", ylab = NULL, type = "l", ylim = c( -1.1, 1.1 ) )
#' lines( s, col = "red" )
envelope <- function( s ) Mod( hilbert( s ) )
