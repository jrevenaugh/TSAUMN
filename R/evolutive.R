#' Time-Frequency Analysis
#'
#' Compute evolutive spectrum using direct FFT method
#' @param s vector of timeseries samples or ts object.  If the former, deltat should be specified,
#' NA values are not allowed.
#' @param deltat sampling interval of time series.  Ignored if s is a ts object
#' @param nfft length of padding window input to fft; controls frequency resolution
#' @param nwin number of samples in fft window
#' @param nstep number of samples to shift right between spectral estimates; controls time resolution
#' @param fl lowest frequency to output
#' @param fh highest frequence to output
#' @return List containing the following items:
#'  \item{"signal"}{Original input signal}
#'  \item{"deltat"}{Sampling interval of time series}
#'  \item{"numfreqs"}{Number of frequencies retained in analysis}
#'  \item{"freqs"}{Vector of frequencies retained in analysis}
#'  \item{"wpars"}{List of analysis parameters: fl, fh and nwin}
#' @import "signal"
#' @import "RColorBrewer"
#' @export
#' @examples
#' require( RColorBrewer )
#' # create chirp signal
#' n <- 1024
#' t <- 1:n
#' a <- 4 * n
#' x <- cos( 2 * pi * t^2 / a )
#' # perform time-frequency analysis of chirp
#' evol <- evolFFT( x, deltat = 1, nwin = 200, nstep = 5 )
#' # plot result
#' plotFTSpec( evol )

evolFFT <- function ( s, deltat = 1,
                      nfft = 1024, nwin, nstep,
                      fl = 0, fh = 0.5 / deltat ) {

  n = length( s )
  if ( n < 64 ) stop( "Really?  Way too short of a time series." )
  nyquist = 0.5 / deltat
  if ( nwin > n / 2 ) {
    warning( noquote( "Adjusted sliding window length." ) )
    nwin <- floor( n / 2 )
  }
  if ( nstep >= n / 10 ) {
    warning( noquote( "Adjusted sliding window step size." ) )
    nstep <- floor( n / 10 )
  }
  if ( ( nfft %% 2 ) != 0 ) {
    warning( noquote( "Adjusted number of fft frequencies." ) )
    nfft <- nextn( nfft, 2 )
  }

  kcol <- floor( ( n - nwin ) / nstep ) + 1
  krow <- nfft
  df <- 1 / ( nfft * deltat )
  M <- matrix( rep( 0, krow * kcol ), nrow = krow, ncol = kcol )
  m = 1:(kcol)
  ibegin <- 1 + ( ( 1:kcol ) - 1 ) * nstep
  iend <- ibegin + nwin - 1
  for ( i in 1:kcol ) {
    x = s[ibegin[i]:iend[i]]
    x <- x - mean( x )
    x <- x * signal::hamming( nwin )
    M[,i] <- c( x, rep(0, krow - nwin ) )
  }
  EP = Mod( mvfft( M ) )^2
  t = ( ibegin + nwin / 2 ) * deltat
  freqs = seq( 0, nyquist, length.out = ( nfft / 2 + 1 ) )
  fout <- which( freqs >= fl & freqs <= fh )
  nout <- length( fout )
  list( signal = s, deltat = deltat, numfreqs = nout, DSPEC = EP[fout,],
             freqs = freqs[fout], tims = t, wpars = list( fl = fl, fh = fh, nwin = nwin ) )
}

#' Time-Frequency Analysis
#'
#' Compute evolutive spectrum using MTM method
#' @param s vector of timeseries samples or ts object.  If the former, deltat should be specified,
#' NA values are not allowed.
#' @param deltat sampling interval of time series.  Ignored if s is a ts object
#' @param nw frequency bandwidth parameter of tapers
#' @param k number of tapers used, usually 2 * nw
#' @param nfft length of padding window input to fft; controls frequency resolution
#' @param nwin number of samples in fft window
#' @param nstep number of samples to shift right between spectral estimates; controls time resolution
#' @param fl lowest frequency to output
#' @param fh highest frequence to output
#' @return List containing the following items:
#'  \item{"signal"}{Original input signal}
#'  \item{"deltat"}{Sampling interval of time series}
#'  \item{"numfreqs"}{Number of frequencies retained in analysis}
#'  \item{"freqs"}{Vector of frequencies retained in analysis}
#'  \item{"wpars"}{List of analysis parameters: fl, fh and nwin}
#' @import "multitaper"
#' @import "RColorBrewer"
#' @export
#' @examples
#' require( RColorBrewer )
#' # create chirp signal
#' n <- 1024
#' t <- 1:n
#' a <- 4 * n
#' x <- cos( 2 * pi * t^2 / a )
#' # perform time-frequency analysis of chirp
#' evol <- evolMTM( x, deltat = 1, nwin = 200, nstep = 5 )
#' # plot result
#' plotFTSpec( evol )
evolMTM <- function ( s, deltat = 1, nw = 4, k = 7,
                      nfft = 1024, nwin, nstep,
                      fl = 0, fh = 0.5 / deltat ) {

  n = length( s )
  if ( n < 64 ) stop( "Really?  Way too short of a time series." )
  nyquist = 0.5 / deltat
  if ( nwin > n / 2 ) {
    warning( noquote( "Adjusted sliding window length." ) )
    nwin <- floor( n / 2 )
  }
  if ( nstep >= n / 10 ) {
    warning( noquote( "Adjusted sliding window step size." ) )
    nstep <- floor( n / 10 )
  }
  if ( ( nfft %% 2 ) != 0 ) {
    warning( noquote( "Adjusted number of fft frequencies." ) )
    nfft <- nextn( nfft, 2 )
  }

  kcol <- floor( ( n - nwin ) / nstep ) + 1
  krow <- nfft + 1
  df <- 1 / ( nfft * deltat )
  m = 1:(kcol)
  ibegin <- 1 + ( ( 1:kcol ) - 1 ) * nstep
  iend <- ibegin + nwin - 1
  EP <- matrix( rep( 0, kcol * krow ), nrow = krow, ncol = kcol )
  for ( i in 1:kcol ) {
    x = s[ibegin[i]:iend[i]]
    x <- x - mean( x )
    x <- c( x, rep( 0, nfft - nwin ) )
    EP[,i] <- multitaper::spec.mtm( ts( x, deltat = deltat ), nw = nw, k = k, plot = FALSE )$spec
  }
  t = ( ibegin + nwin / 2 ) * deltat
  freqs = seq( 0, nyquist, length.out = krow )
  fout <- which( freqs >= fl & freqs <= fh )
  nout <- length( fout )
  list( signal = s, deltat = deltat, numfreqs = nout, DSPEC = EP[fout,],
        freqs = freqs[fout], tims = t, wpars = list( fl = fl, fh = fh, nwin = nwin ) )
}
