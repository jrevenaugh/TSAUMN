#' Plot timeseries
#'
#' Attractive \link{ggplot2} plot of an univariate timeseries  
#' @param s a timeseries
#' @return a ggplot2 object
#' @export
#' @import "ggplot2"
#' @examples
#' Timeseries <- ts( rnorm( 100 ), deltat = 1 )
#' gtsplot( Timeseries )
#' 
gplot.ts <- function( s ) {
  series <- deparse( substitute( s ) )
  if ( !is.ts( s ) ) stop( "s must be a timeseries" )
  S <- data.frame( as.numeric( time( s ) ), as.numeric( s ) )
  colnames( S ) <- c( "Time", series )
  g <- ggplot2::ggplot( S, ggplot2::aes_string( x = "Time", y = series  ) ) + 
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous( )
  return( g )
}