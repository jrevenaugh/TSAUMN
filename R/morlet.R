#' Morlet wavelet
#'
#' Compute real-value morlet wavelet  
#' @param n number of output samples
#' @param center where to center the wavelet in n-sample sequence
#' @param scale dilation factor 
#' @param w0 sinusoid radial frequency (defaults to 2 * pi )
#' @return vector of wavelet coefficients
#' @export
#' @examples
#' w <- morlet( 500, 250, 10 )
#' plot( 1:500, w, type = "l" )
morlet <- function( n = 500, center = 250, scale = 10, w0 = 2 * pi ) {
  t <- 1:n
  dt <- ( t - center ) / scale
  w <- pi^(-1/4) * cos( w0 * dt ) * exp( -dt^2 / 2 )
  c <- sqrt( 1 + exp( -w0^2 ) -2 * exp( -3/4 * w0^2 ) )
  return( c * w )
}