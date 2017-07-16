#' Phase unwrapping
#'
#' Unwrap a vector of phase values using Itoh's algorithm.  
#' @param p vector of phase values
#' @return unwrapped phase values
#' @export
#' @examples
#' s <- rep( 0, 124 )
#' s[50] <- 1
#' S <- fft( s )
#' p <- Arg( S )
#' up <- phunwrap( p[1:63] )
#' f <- seq( 0, 0.5, length.out = 63 )
#' plot( f, up, xlab = "Frequency", ylab = "Phase (rad)", type = "l", main = "Estimated (black) and Theory (red)" )
#' lines( f, -2 * pi * f * 49, col = "red" )
#' 
phunwrap <- function( p ) {
  n <- length( p )
  if ( n < 2 ) {
    warning( noquote( "Series too short to unwrap." ) )
    return( p )
  }
  Dp <- c( 0, -diff( p ) )
  s <- 2 * pi * ( as.numeric( Dp >= pi ) - as.numeric( Dp < -pi ) )
  return( p + cumsum( s ) )
}