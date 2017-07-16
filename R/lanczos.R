#' Create a Lanczos filter
#'
#' Compute coefficients of a Lanczos filter  
#' @param n length of filter; must be odd
#' @param type one of "lp" (low pass), "hp" (high pass), "bp" (band pass) or "br" (band reject)
#' @param low low frequency transition of "lp", "bp" and "br" filters (Nyquist = 0.5)
#' @param high high frequency transition of "hp", "bp" and "br" filters (Nyquist = 0.5)
#' @return non-negative time filter coefficients (filter is symmetric about 0).
#' @export
#' @examples
#' # low-pass filter
#' lp <- lanczos( 7, "lp", low = 0.1 )
#' 
#' # high-pass filter
#' hp <- lanczos( 7, "hp", high = 0.4 )
#' 
#' # band-pass filter
#' bp <- lanczos( 7, "bp", low = 0.1, high = 0.4 )
#' 
#' # band-reject filter
#' br <- lanczos( 7, "br", low = 0.1, high = 0.4 )
#' 
#' plot( lp, xlab = "Index", ylab = "Filter Coefficients", type = "l", ylim = 1.5 * range( lp, hp, bp, br ) )
#' lines( hp, col = "red" )
#' lines( bp, col = "blue" )
#' lines( br, col = "green" )
lanczos <- function( n = 7, 
                     type = c( "lp", "hp", "bp", "br" ), 
                     low = 0, high = 0.5 ) {
  type <- match.arg( type )
  switch( type, 
          lp = lanczos_lp( n, low ),
          hp = lanczos_hp( n, high ),
          bp = lanczos_bp( n, low, high ),
          br = lanczos_br( n, low, high ) )
}

#' Create a low-pass Lanczos filter
#'
#' Compute coefficients of a low-pass Lanczos filter  
#' @param n length of filter; must be odd
#' @param low cutoff frequency of the filter (Nyquist = 0.5)
#' @return non-negative time filter coefficients (filter is symmetric about 0).
#' @export
#' @examples
#' f <- lanczos_lp( 7, 0.1 )
#' print( f )
#' 
lanczos_lp <- function( n, low ) {
  order <- ( ( n - 1 ) %/% 2 ) + 1
  n <- order
  w <- rep( 0, n )
  w[1] <- 2 * low
  k <- 1:( n - 1 )
  sigma <- sin( pi * k / n ) * n / ( pi * k )
  firstfactor <- sin( 2 * pi * low * k ) / ( pi * k )
  w[2:n] <- firstfactor * sigma
  return ( w )
}

#' Create a high-pass Lanczos filter
#'
#' Compute coefficients of a high-pass Lanczos filter  
#' @param n length of filter; must be odd
#' @param high cutoff frequency of the filter (Nyquist = 0.5)
#' @return non-negative time filter coefficients (filter is symmetric about 0).
#' @export
#' @examples
#' f <- lanczos_hp( 7, 0.4 )
#' print( f )
#' 
lanczos_hp <- function( n, high ) {
  w <- lanczos_lp( n, high )
  w <- -w
  w[1] <- w[1] + 1
  return( w )
}

#' Create a band-reject Lanczos filter
#'
#' Compute coefficients of a band-reject Lanczos filter  
#' @param n length of filter; must be odd
#' @param low transition frequency of the filter (Nyquist = 0.5)
#' @param high transition frequency of the filter (Nyquist = 0.5)
#' @return non-negative time filter coefficients (filter is symmetric about 0). 
#' Frequencies inside band (low, high) are attenuated.
#' @export
#' @examples
#' f <- lanczos_br( 7, 0.1, 0.4 )
#' print( f )
#' 
lanczos_br <- function( n, low, high ) {
  wl <- lanczos_lp( n, low )
  wh <- lanczos_hp( n, high )
  return( wl + wh )
}

#' Create a band-pass Lanczos filter
#'
#' Compute coefficients of a band-pass Lanczos filter  
#' @param n length of filter; must be odd
#' @param low transition frequency of the filter (Nyquist = 0.5)
#' @param high transition frequency of the filter (Nyquist = 0.5)
#' @return non-negative time filter coefficients (filter is symmetric about 0). 
#' Frequencies outside of band (low, high) are attenuated.
#' @export
#' @examples
#' f <- lanczos_bp( 7, 0.1, 0.4 )
#' print( f )
#' 
lanczos_bp <- function( n, low, high ) {
  wl <- lanczos_lp( n, high )
  wh <- lanczos_hp( n, low )
  w <- wl + wh
  w[1] <- w[1] - 1
  return( w )
}

#' Compute Lanczos filter amplitude and phase (not unwound) spectra
#'
#' Compute Lanczos filter amplitude spectrum 
#' @param w filter coefficients returned by \code{\link{lanczos}}
#' @param n zero-padded length of filter to transform
#' @return data.frame containing amplitude and phase spectrum of filter for non-negative frequencies (DC to Nyquist).
#'  \item{"freq"}{Fourier frequencies (0 to Nyquist)}
#'  \item{"amp"}{Fourier amplitude spectrum of filter}
#'  \item{"phase"}{Fourier phase spectrum of filter}
#' @export
#' @examples
#' f <- lanczos_bp( 7, 0.1, 0.4 )
#' s <- lanczos_spec( f, 512 )
#' plot( s$freq, s$amp, type = "l", xlab = "Frequency", ylab = "Amplitude" )
#' 
lanczos_spec <- function( w, nPts ) {
  n <- nextn( nPts, 2 ) + 1
  h <- rep( 0, n )
  h[1:length( w )] <- w
  h <- h + rev( h )
  nf <- ( n - 1 ) / 2 + 1
  H <- fft( h[-n] )[1:nf]
  amp <- Mod( Conj( H ) * H )
  phase <- Arg( H )
  freq <- seq( 0, 0.5, length.out = nf )
  return( data.frame( freq = freq, amp = amp, phase = phase ) )
}
  
