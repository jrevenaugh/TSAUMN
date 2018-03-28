#' Compute Multi-taper Averaged Lomb-Scargle periodogram of (un)evenly spaced data
#'
#' Computes the Lomb-Scargle periodogram for a time series with irregular (or regular) sampling intervals using
#' a series of "Slepian"-like tapers (the tapers are not orthogonal on the irregular sampling grid).  The ensemble
#' is averaged to produce a consistent estimator.
#' @param x The data to be analysed. x can be either a two-column numerical dataframe or matrix, with sampling times in columnn 1 and measurements in column 2, a single numerical vector containing measurements, or a single vector ts object (which will be converted to a numerical vector).
#' @param t If x is a single vector, t can be provided as a numerical vector of equal length containing sampling times. If x is a vector and times is NULL, the data are assumed to be equally sampled and times is set to 1:length(x).
#' @param nw Thompson's frequency bandwidth parameter (>= 1)
#' @param k Number of tapers, usually 2nw or 2nw - 1 (defaults to 2 nw)
#' @param demean remove mean from timeseries prior to spectral estimation
#' @param detrend remove linear trend from timeseries prior to spectral estimation
#' @param plot Logical.  If plot = TRUE, the spectrum is plotted.
#' @param ... Additional arguments passed to gplot.mtm.
#' @return object of class spec with the following list items:
#'   \item{"freq"}{A vector with spectrum frequencies}
#'   \item{"spec"}{A vector with spectral power estimates corresponding to "freq"}
#'   \item{"series"}{Name of input time series}
#'   \item{"method"}{Method name: "MTLS"}
#' @export
#' @import "lomb"
#' @import "multitaper"
#' @examples
#' x <- rnorm( 256 )
#' s <- spec.ls( x, 1:256, plot = TRUE )

spec.mtls <- function( x, t = NULL, nw = 4, k = NULL, demean = TRUE, detrend = FALSE, plot = TRUE, ... ) {
  if ( nw < 1 ) {
    nw <- 1
    warning( "nw coerced to 1" )
  }
  if ( is.null( k ) ) k = 2 * nw
  if ( k < 1 ) {
    warning( "k coerced to 1" )
    k <- 1
  }
  
  if ( !is.null( t ) ) {
    if ( !is.vector( t ) ) stop( "No multivariate methods available" )
    if ( length( x ) != length( t) ) stop( "Length of x and t vectors must be equal" )
    series <- deparse( substitute( x ) )
  }
  
  if ( is.null( t ) && is.null( ncol( x ) ) ) {
    series <- deparse( substitute( x ) )
    t <- 1:length( x )
  }
  if ( is.matrix( x ) || is.data.frame( x )) {
    if ( ncol( x ) > 2) stop( "No multivariate methods available" )
    if ( ncol( x ) == 2) {
      series <- colnames( x )[2]
      t <- x[, 1]
      x <- x[, 2]
    }
  }
  t <- t[!is.na( x )]
  x <- x[!is.na( x )]
  
  if ( detrend ) {
    l <- length( x )
    x <- residuals( stats::lm( x ~ stats::poly( t, 1 ) ) )
  } else if ( demean ) {
    x <- x - mean( x )
  }
  
  tapers <- multitaper::dpss( length( t ), k, nw )
  tapx <- tapers$v[,1] * x
  ss <- lomb::lsp( tapx, t, plot = FALSE )
  raw <- data.frame( freq = ss$scanned, spec = ss$power )
  spec <- list( freq = ss$scanned, spec = ss$power )
  freq <- ss$scanned
  power <- ss$power
  if ( k > 1 ) {
    for ( i in 2:k ) {
      tapx <- tapers$v[,i] * x
      ss <- lomb::lsp( tapx, t, plot = FALSE )
      power <- power + ss$power
    }
  }
  power <- power / k
  spec <- list( freq = freq, spec = power, series = series, method = "MTLS" )
  if ( plot ) {
    title = sprintf( "MTLS Spectrum of %s", series )
    print( gplot.mtm( spec, ... ) + labs( title = title ) )
  }
  return( spec )
}