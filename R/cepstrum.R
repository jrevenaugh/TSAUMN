#' Cepstrum of univariate real timeseries
#'
#' Compute cepstrum of a univariate, real-valued timeseries  
#' @param s real-valued, univariate timeseries
#' @param taper taper applied to s prior to fft; one of "hamming" (default), "akaike" or "gausswin"
#' @param poly polynomial order of trend removed from the log of power spectral density; defaults to 5
#' @return data.frame with quefrency and cepstrum
#' @export
#' @import "signal"
#' @examples
#' t <- 1:500
#' f <- seq( 0.05, 0.25, 0.05 )
#' s <- sapply( t, function( x, y ) sum( sin( 2 * pi * x * y ) ), y = f )
#' s <- s + rnorm( length( t ), sd = 0.5 )
#' c <- cepstrum( s, taper = "akaike" )
#' ggplot( c, aes( x = quefrency, y = cepstrum ) ) + geom_line() +
#' labs( x = "Quefrency", y = "Cepstrum" )
#' # harmonics spaced at 0.05 produce peak at quefrency of 20 (1 / 20 = 0.05 )
#' 
cepstrum <- function( s, taper = c( "hamming", "akaike", "gausswin" ), poly = 5 ) {
  taper <- match.arg( taper )
  m <- length( s )
  n <- nextn( m, 2 )
  f <- switch( taper,
               hamming = signal::hamming( m ),
               akaike = akaike( m ),
               gausswin = signal::gausswin( m ) )
  s <- c( f * s, rep( 0, n - m ) )
  S <- fft( s )
  l <- n / 2 + 1
  i <- 0 + 1i
  logP <- log( Mod( S )^2 )
  if ( poly > 0 ) logP <- residuals( lm( logP ~ stats::poly( 1:n, poly ) ) )
  c <- data.frame( quefrency = 1:l - 1, 
                   cepstrum = Re( fft( logP, inverse = TRUE )[1:l] ) / n )
  return( c )
}