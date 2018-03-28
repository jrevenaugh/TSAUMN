#' Upsample timeseries
#'
#' Upsample a timeseries using FFT  
#' @param s timeseries (as a ts object or vector of samples)
#' @param nfac upsampling factor (integer greater than 1)
#' @param pct percent length of timeseries used for cosine-squared taper of ends.  
#' This reduces edge effects, but is lossy. 
#' @param deltat sampling interval; used only if s is not a ts object
#' @return resampled timeseries.
#' @export
#' @examples
#' s <- ts( sin( 2 * pi * 1:100 / 20 ), deltat = 1 )
#' us <- fftresample( s, nfac = 10 )
#' plot( s, type = "l" )
#' lines( us, col = "red" )
fftresample <- function( s, nfac = 2, pct = 0.025, deltat = 1 ) {
  if ( nfac <= 1 || ( floor( nfac ) != nfac ) ) stop( "nfac must be an integer greater than 1." )
  if ( pct < 0 ) stop( "pct must be greater than 0." )
  if ( is.ts( s ) ) {
    deltat <- deltat( s )
    start <- time( s )[1]
  } else {
    start <- 0
  }
  n <- length( s )
  m <- nextn( n, c( 2 ) )
  s <- c( s * costap( n, pct ), rep( 0, m - n ) )
  M <- m * nfac
  S <- fft( s )
  nyq <- m / 2 + 1
  H <- c( S[1:nyq], rep( 0 + 0i, M - m - 1 ), S[nyq:m])
  h <- fft( H, inverse = T ) / m
  s <- Re( h[1:(nfac * n)] )
  s <- ts( s, deltat = deltat / nfac, start = start )
  return( s )
}

#' Cosine-squared taper
#'
#' Compute cosine-squared taper  
#' @param n length of taper
#' @param pct percent of n used to roll on and off the taper on ends.
#' @return taper vector.
#' @export
#' @examples
#' t <- costap( 100, pct = 0.2 )
#' plot( t, xlab = "Time", ylab = "Taper", type = "l" )
costap <- function( n, pct = 0.05 ) {
  m <- round( n * pct )
  c2 <- 0.5 * ( 1 - cos( pi * seq( 0, m - 1 ) / m ) )
  tap <- c( c2, rep( 1, n - 2 * m ), rev( c2 ) )
  return( tap )
}

#' Downsample timeseries
#'
#' Downsample a timeseries using FFT  
#' @param s timeseries (as a ts object or vector of samples)
#' @param nfac downsampling factor (integer greater than 1)
#' @param pct percent length of timeseries used for cosine-squared taper of ends.  
#' This reduces edge effects, but is lossy. 
#' @param deltat sampling interval; used only if s is not a ts object
#' @return resampled timeseries.
#' @export
#' @examples
#' s <- ts( sin( 2 * pi * 1:1000 / 200 ), deltat = 0.1 )
#' ds <- fftdecimate( s, nfac = 10 )
#' plot( s, type = "l" )
#' lines( ds, col = "red" )
fftdecimate <- function( s, nfac = 2, pct = 0.025, deltat = 1 )
{
  if ( nfac <= 1 ) stop( "nfac must be greater than 1." )
  if ( pct < 0 ) stop( "pct must be greater than 0." )
  if ( is.ts( s ) ) {
    deltat <- deltat( s )
    start <- time( s )[1]
  } else {
    start <- 0
  }
  n <- length( s )
  m <- nextn( nextn( n, 2 ), nfac )
  s <- c( s * costap( n, pct ), rep( 0, m - n ) )
  S <- fft( s )
  nnyq <- m / ( 2 * nfac ) + 1
  H <- rep( 0 + 0i, m )
  H[1:nnyq] <- S[1:nnyq]
  H <- H + Conj( c( 0, rev( H[-1] ) ) )
  h <- fft( H, inverse = T ) / m
  ind <- seq( 1, n, nfac )
  s <- Re( h[ind] )
  s <- ts( s, deltat = deltat * nfac, start = start )
  return( s )
}



