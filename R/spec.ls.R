#' Compute Lomb-Scargle periodogram of (un)evenly spaced data
#'
#' Computes the Lomb-Scargle periodogram for a time series with irregular (or regular) sampling intervals.   
#' @param x The data to be analysed. x can be either a two-column numerical dataframe or matrix, with sampling times in columnn 1 and measurements in column 2, a single numerical vector containing measurements, or a single vector ts object (which will be converted to a numerical vector).
#' @param t If x is a single vector, t can be provided as a numerical vector of equal length containing sampling times. If x is a vector and times is NULL, the data are assumed to be equally sampled and times is set to 1:length(x).
#' @param over The oversampling factor. Must be an integer >= 1. Larger values of over lead to finer scanning of frequencies. 
#' @param demean remove mean from timeseries prior to spectral estimation
#' @param detrend remove linear trend from timeseries prior to spectral estimation
#' @param plot Logical.  If plot = TRUE, the spectrum is plotted.
#' @param ... Additional arguments passed to gplot.mtm.
#' @return object of class spec with the following list items:
#'   \item{"freq"}{A vector with spectrum frequencies}
#'   \item{"spec"}{A vector with spectral power estimates corresponding to "freq"}
#'   \item{"series"}{Name of input time series}
#'   \item{"method"}{Method name: "Lomb-Scargle"}
#' @export
#' @import "lomb"
#' @examples
#' x <- rnorm( 256 )
#' s <- spec.ls( x, 1:256, plot = TRUE )

spec.ls <- function( x, t = NULL, over = 1, demean = TRUE, detrend = FALSE, plot = TRUE, ... ) {
  if ( over != floor( over ) ) {
    over <- floor( over )
    warning( "over coerced to integer" )
  }
  if ( over < 1 ) {
    over <- 1
    warning( "over must be integer >=1. Set to 1" )
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
  
  s <- lomb::lsp( x, t, ofac = over, plot = FALSE )
  spec <- list( freq = s$scanned, spec = s$power, series = series, method = "Lomb-Scargle" )
  if ( plot ) {
    title = sprintf( "Lomb-Scargle Periodogram of %s", series )
    print( gplot.mtm( spec, ... ) + labs( title = title ) )
  }
  return( spec )
}