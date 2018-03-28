#' Cosine-squared timeseries filtering
#'
#' Filter a regularly sampled univariate or multivariate timeseries using cosine-squared,
#' symmetric filter bank.  Filtering is performed in the frequency domain.  
#' @param s univariate or multivariate (as matrix columns) real-valued timeseries.
#' @param fc vector of "corner" frequencies specifying endpoints of frequency bands.  
#' Regardless of sample rate, these are specified in terms of the Nyquist band
#' (0,0.5).  
#' @param fw amplitudes of filter at corresponding corner frequencies.  Lengths of fc and fw must match.
#' Must be 0 (reject) or 1 (pass).  
#' @details By specifying fc and fw, one can create many filters: low-pass, high-pass, band-pass
#' and band-reject as well as their combinations.  Filter response will be 0 or 1 except across
#' intervals where fw changes.  The lengths of these transitions affect the length of the effective filter impulse response.  
#' Short transitions produce long filters.  Filtering is performed in the frequency domain and is circular. 
#' Input timeseries is padded to next factor of 2 prior to filtering. If the length is near a factor of 2,
#' you may want to pad s prior to use of cosfilt, especially if filter transition intervals are short.  
#' @return filtered timeseries
#' @export
#' @examples
#' s <- 0.5 * sin( 2 * pi * 1:200 / 20 ) + rnorm( 200 )
#' fs <- cosfilt( s, fc = c( 0, 0.1, 0.2, 0.5 ), fw = c( 1, 1, 0, 0 ) )
#' plot( s, type = "l" )
#' lines( fs, col = "red" )
#' 

cosfilt <- function( s, fc = c( 0, 0.2, 0.3, 0.5), 
                     fw = c( 1, 1, 0, 0 ) ) {
  if ( any( fc < 0 ) || any( fc > 0.5 ) ) stop( "Frequencies must be in range (0, 0.5)" )
  if ( any( as.logical( fw %% 1 ) ) || any( fw < 0 ) || any( fw > 1 ) ) stop( "Weights must be in range (0, 1)" )
  if ( length( fc ) != length( fw ) ) stop( "Lenghts of fc and fw must be equal" )
  if ( is.null( dim( s ) ) ) {
    n <- length( s )
    m <- 1
    s <- ts( matrix( s, ncol = 1 ), start = time( s )[1], deltat = deltat( s ) )
  } else{
    n <- dim( s )[1]
    m <- dim( s )[2]
  }
  npad <- nextn( n, 2 )
  l <- npad / 2 + 1
  freq <- seq( 0, 0.5, length.out = l )

  # build H
  fc <- fc[o <- order(fc)]
  fw <- fw[o]
  if ( fc[1] != 0 ) {
    fc <- c( 0, fc )
    fw <- c( fw[1], fw )
  }
  p <- length( fc )
  if ( fc[p] != 0.5 ) {
    fc <- c( fc, 0.5 )
    fw <- c( fw, fw[p] )
    p <- p + 1
  }
  H <- rep( 1, length( freq ) )
  for ( i in 1:(p-1) ) {
    ind <- which( freq >= fc[i] & freq <= fc[i+1] )
    if ( fw[i] == fw[i+1] ) {
      H[ind] <- fw[i]
    } else {
      ft <- pi * ( freq[ind] - fc[i] ) / ( fc[i + 1] - fc[i] )
      T <- 0.5 * ( 1 - cos( ft ) )
      if ( fw[i] == 1 ) T <- 1 - T
      H[ind] <- T
    }
  }

  # cycle through series
  for ( i in 1:m ) {
    st <- c( s[,i], rep( 0, npad - n ) )
    S <- fft( st )[1:l] * H
    S <- c( S, Conj( rev( S[2:(l-1)] ) ) )
    st <- fft( S, inverse = TRUE ) / npad
    s[,i] <- st[1:n]
  }
  if ( m == 1 ) s <- s[,1]
  return( Re( s ) )
}