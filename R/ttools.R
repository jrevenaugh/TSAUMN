#' Convert data.frame object to a ts
#'
#' Convert data.frame object to a ts  
#' @param df data.frame with regularly sampled time series coefficients as columns and time as an additional column.
#' @param time character vector name of time column or its numeric index.
#' @return ts object with converted data.frame object.
#' @export
#' @importFrom dplyr select
#' @examples
#' s.df <- data.frame( time = 1:100, s = sin( 2 * pi * 1:100 / 20 ) )
#' s.ts <- df2ts( s.df, "time" )
#' 
df2ts <- function( df, time ) {
  if( !is.data.frame( df ) ) stop( "Input not a data.frame" )
  t <- df[, time]
  st <- min( t )
  dt <- t[2] - t[1]
  ts( data = dplyr::select( df, -one_of( time ) ), start = st, deltat = dt )
}

#' Convert ts object to a data.frame
#'
#' Convert ts object to a data.frame  
#' @param s time series object
#' @return data.frame with time as column 1, time series as remaining columns.
#' @export
#' @examples
#' s <- ts( data = sin( 2 * pi * (1:100) * 0.1 ), start = 1, end = 100, deltat = 1 )
#' S.df <- ts2df( s )
#' head( S.df )
#' 
ts2df <- function( s ) {
  if ( !is.ts( s ) ) stop( "Input not a time-series (ts)." )
  time <- time( s )
  data.frame( time, s )
}


#' Welch window
#'
#' Welch's windowing function  
#' @param n length of the filter; number of coefficients to generate.
#' @return Filter coefficients
#' @export
#' @examples
#' n <- 100
#' plot( welch( n ), type = "l" )
#' 
welch <- function( n ) {
    n <- n - 1
    c <- 1 - ( ( 0:n - n / 2 ) / ( n / 2 ) )^2
    return( c )
}

#' Akaike 2nd-order window
#'
#' Akaike's 2nd-order windowing function  
#' @param n length of the filter; number of coefficients to generate.
#' @return Filter coefficients
#' @export
#' @examples
#' n <- 100
#' plot( akaike( n ), type = "l" )
#' 
akaike <- function( n ) {
  m <- floor( n / 2 )
  a <- c( 0.64, 0.24, -0.06 )
  t <- seq( -m, m, length.out = n )
  i <- ( 0 + 0i )
  w <- a[1] + 2 * a[2] * cos( pi * t / m ) + 2 * a[3] * cos( 2 * pi * t / m )
  return( w )
}

#' Power spectrum
#'
#' Direct FFT computation of time series power spectrum  
#' @param s vector of time series samples
#' @return  Power of non-negative frequencies order DC to Nyquist.
#' @export
#' @examples
#' s <- sin( 2 * pi * ( 1:256 ) * 0.1 )
#' PS <- pspec( s )
#' f <- seq( 0, 0.5, length.out = length( PS ) )
#' plot( f, PS, type = "l", xlab = "Frequency", ylab = "Power" )
pspec <- function( s ) return( Mod( fftnneg( s ) )^2 )

#' Zero pad
#'
#' Pad a vector with zeros  
#' @param s vector of time series samples
#' @param n expansion factor - returned vector will be at least n times original length.
#' @return  Padded input signal
#' @export
#' @examples
#' n <- 16
#' s <- rep( 1, 32 )
#' spad <- zpad( s, n )
#' S <- fft( spad )
#' # Show 
#' f <- seq( 0, to = 0.5, length.out = length( spad ) )
#' plot( f, Mod( S[1:(length( spad ) / 2 + 1 )] ), type = "l", xlab = "Frequency", ylab = "Mod(S)" )
#' 
zpad <- function( s, n ) {
    N <- n * nextn( length( s ), c( 2, 3, 5 ) )
    return( c( s, rep( 0, N - length( s ) ) ) )
}

#' Spectral leakage
#'
#' Evaluate side-lobe leakage for a given windowing function  
#' @param w windowing function coefficients
#' @return  log(modulus) of the padded fft
#' @export
#' @examples
#' s <- rep( 1, 32 )
#' S <- leakage( s )
#' f <- seq( 0, to = 0.5, length.out = length( S ) )
#' plot( f, S, type = "l", xlab = "Frequency", ylab = "log(Mod(S))" )
#' 
leakage <- function( w ) {
    w <- zpad( w, 16 )
    x <- log( pspec( w ) )
    return( x )
}


#' Non-negative frequencies of FFT
#'
#' Compute non-negative frequency coefficients of fft  
#' @param s vector of time series samples
#' @return complex valued fft of s at non-negative frequencies, ordered DC to Nyquist
#' @export
#' @examples
#' s <- sin( 2 * pi * ( 1:256 ) * 0.1 )
#' S <- fftnneg( s )
#' f <- seq( 0, to = 0.5, length.out = length( S ) )
#' plot( f, Mod( S ), type = "l", xlab = "Frequency", ylab = "Mod(S)" )
#' 
fftnneg <- function( s ) {
  n <- nextn( length( s ), c( 2, 3, 5 ) )
  m <- floor( n / 2 + 1 )
  f <- fft( s )
  return( f[1:m] )
}
