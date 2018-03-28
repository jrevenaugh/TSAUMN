#' Deconvolve a known linear filter
#'
#' Deconvolve a known linear filter from sampled univariate or multivariate timeseries using standard spectral
#' division with or without some regularization.
#' @param s real-valued univariate or multivariate (as matrix columns) timeseries.
#' @param f real-valued filter to deconvolve.
#' @param type one of water-level ("wl", default), damped "least squqres" ("dls") or no regularization ("noreg").
#' See detail for more information  
#' @param metric one of "max" (default), "median" or "percentile."  Together with \code{alpha} this determines
#' the regularization of the deconvolution.
#' @param alpha multiplier of the max or median value of the denominator for the corresponding \mode{metric},
#' If \code{metric} is "percentile", this specifies the percentile value of the denominator to use in regularization.
#' @details If the convolved filter has spectral zeros
#' (or near zeros), any noise in that band will be greatly exaggerated by deconvolution. To control this,
#' it is common to regularize the denominator during spectral division (deconvolution in the frequency domain).
#' Two methods are provided. (1) water-level ("wl") doesn't allow the value of the denominator to drop below the scaled metric.
#' Visually, the denominator plotted against frequency is filled with water up to the scaled metric, thus
#' filling in "holes."  (2) damped "least-squares" ("dls") simply adds the scaled metric to the 
#' denominator at all frequencies.  Visually, this shifts the baseline of the denominator up.  While this
#' reduces the noise-amplification of zeros and near-zeros, it biases results at all frequencies.
#' Type "noreg" applies no regularization and will produce useless results for noisy inputs.
#' 
#' The scaled metric is one of: alpha * max( D ), alpha * median( D ) or quantile( D, alpha ) where D
#' is the denominator (function of frequency).
#' 
#' Deconvolution is circular.  Long filters will cycle around.  
#' Padding \code{s} prior to deconvolution is often desirable.
#' @return Deconvolved timeseries
#' @export
#' @examples
#' # Create spike sequence
#' r <- runif( 1024, min = 0, max = 1 )
#' c <- which( r > 0.95 )
#' r <- rep( 0, 1024 )
#' r[c] <- 5 * rnorm( length( c ) )
#'
#' # Create a wavelet
#' t <- seq( 0, 1023, 1 )
#' w <- 10 * dnorm( t, sd = 10 ) * sin( 2.0 * pi * t / 20 )
#'
#' # Convolve spikes with wavelet and add white noise.
#' z <- convolve( r, w, conj = FALSE, type = "circular" ) + rnorm( 1024, sd = 0.1 )
#'
#' # Deconvolve
#' rd <- deconvolve( z, w, type = "wl", metric = "max", alpha = 0.005 )
#' plotl( r, main = "Water Level" )
#' lines( rd, col = "red" )
#' rd <- deconvolve( z, w, type = "dls", metric = "max", alpha = 0.005 )
#' plotl( r, main = "Damped LS" )
#' lines( rd, col = "red" )

deconvolve <- function( s, f, type = c( "wl", "dls", "noreg" ), 
                       alpha = 0.1, metric = c( "median", "max", "percentile" ) ) {
  # match inputs
  type <- match.arg( type )
  metric <- match.arg( metric )
  alpha <- abs( alpha )

  # establish sizes
  if ( is.null( dim( s ) ) ) {
    n <- length( s )
    m <- 1
    s <- matrix( s, ncol = 1 )
  } else{
    n <- dim( s )[1]
    m <- dim( s )[2]
  }
  if ( !is.null( dim( f ) ) ) {
    warning( "Filters past first were ignored." )
    f <- f[,1]
  }
  n <- max( n, length( f ) )
  npad <- nextn( n, 2 )
  l <- npad / 2 + 1
  
  # prepare inverse filter
  f <- c( f, rep( 0, npad - n ) )
  F <- fft( f )
  D <- Re( Conj( F ) * F )
  reg <- switch( metric,
                 "max" = alpha * max( D ),
                 "median" = alpha * median( D ),
                 "percentile" = quantile( D, alpha ) )
  
  if ( type == "wl" ) {
    ind <- which( D < reg )
    D[ind] <- reg
  } 
  if ( type == "dls" ) D <- D + reg

  # loop over input signals
  for ( i in 1:m ) {
    st <- c( s[,i], rep( 0, npad - n ) )
    S <- fft( st )
    S <- Conj( F ) * S / D
    st <- fft( S, inverse = TRUE ) / n
    s[,i] <- st[1:n]
  }
  
  return( Re( s ) )
}