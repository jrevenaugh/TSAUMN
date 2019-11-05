#' Blackman-Tukey spectral density estimation
#'
#' Perform Blackman-Tukey power-spectral density estimate of a timeseries.
#' @param s univariate real-valued timeseries.
#' @param q maximum lag retained in the auto-correlation function.
#' @param taper apply Akaike ACF taper if true (default).
#' @param n length of zero-padded ACF function.
#' @param deltat sample interval of timeseries if s passed as a vector.
#' @param demean logical.  If true, timeseries mean is removed.
#' @param detrend logical.  If true, linear trend is removed from timeseries.
#' @param plot logical.  If true, spectrum is plotted using gplot.mtm
#' @param ... Additional parameters passed to gplot.mtm
#' @details Blackman-Tukey is a parametric power-spectral density estimator.
#' It assumes an MA(q) process with known q. In practice, increasing q
#' trades off greater spectral resolution with greater estimator variance. In all cases,
#' q should be much shorter than the timeseries itself.
#' @return list with items:
#'   \item{freq}{Vector of sampled frequencies}
#'   \item{spec}{Vector of power spectral estimates}
#'   \item{series}{Name of the input time series}
#'   \item{Method}{"Blackman Tukey"}
#' @export
#' @examples
#' s <- arimaSim( 1028, ma = 0.5 )
#' spec <- spec.bt( s, q = 5 )
#'
spec.bt <- function( s = NA, q = 20, taper = TRUE, n = 2048, deltat = 1,
                     demean = TRUE, detrend = FALSE, plot = TRUE, ... ) {
  if ( is.null( s ) ) stop( "No time series provided." )
  series <- deparse( substitute( s ) )

  if ( !is.ts( s ) ) s <- ts( s, deltat = deltat )
  else deltat <- deltat(s)
  if ( q < 1 || q > length( s ) / 2 ) stop( "Process order out of bounds." )
  if ( n < q ) stop( "Increase n or lower q." )

  if ( detrend ) {
    l <- length( s )
    s <- residuals( stats::lm( s ~ stats::poly( 1:l, 1 ) ) )
  } else if ( demean ) {
    s <- s - mean( s )
  }

  # Compute the (optionally) Akaike-windowed ACF
  a <- acf( s, lag.max = q, plot = F )$acf[,1,1]
  l <- q + 1
  if ( taper ) a <- a * akaike( 2 * l )[(l+1):(2*l)]

  # Mirror the ACF with zero padding
  x1 <- rep( 0, n )
  x1[1:(q+1)] <- a
  u <- n - q + 1
  x1[n:u] <- a[-1]

  # FFT the padded, mirrored acf function.
  f <- fft( x1 )
  u <- n / 2 + 1
  y <- Mod( f[1:u] ) * var( s )

  # Compute Fourier frequencies
  freq <- seq( 0, 0.5, length.out = u ) / deltat
  spec <- list( freq = freq, spec = y, series = series, method = "Blackman Tukey" )
  if ( plot ) {
    title = sprintf( "Blackman Tukey Spectrum of %s", series )
    print( gplot.mtm( spec, ... ) + labs( title = title ) )
  }
  return( spec )
}
