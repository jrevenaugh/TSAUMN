#' Plot autocorrelation function
#'
#' Attractive \link{ggplot2} plot of \code{\link{acf}} output  
#' @param x a timeseries
#' @param lag.max maximum time lag to display
#' @param type one of "correlation" (default) or "covariance."  
#' The latter is not implemented. 
#' @param demean remove timeseries mean prior to computation of ACF
#' @param ci confidence interval to display
#' @return a ggplot2 object
#' @export
#' @import "ggplot2"
#' @examples
#' x <- rnorm( 1000 )
#' plot.acf( x, lag.max = 20 )

gplot.acf <- function( x,
                      lag.max = NULL, 
                      type = c( "correlation", "covariance" ),
                      demean = TRUE,
                      ci = 0.95 ) {
  
  type <- match.arg( type )
  r <- acf( x, lag.max = lag.max, type = type, demean = demean, plot = FALSE )
  R <- data.frame( Lag = r$lag[,1,1], ACF = r$acf[,1,1] )
  n <- length( R$Lag )
  R$angle <- rep( pi / 2, n )
  R$y <- rep( 0, n )
  
  s <- qnorm( ( 1 + ci ) / 2 ) / sqrt( r$n.used )
  R$psigma <- rep( s, n )
  R$nsigma <- rep( -s, n )
  
  g <- ggplot2::ggplot( R, aes( x = Lag, y = y ) ) + 
         ggplot2::geom_spoke( aes( angle = angle, radius = ACF ) ) +
         ggplot2::geom_point( aes( y = ACF ) ) +
         ggplot2::geom_line( aes( y = psigma ), linetype = 2, color = "dodgerblue3" ) +
         ggplot2::geom_line( aes( y = nsigma ), linetype = 2, color = "dodgerblue3" ) +
         ggplot2::labs( y = "ACF" )
  
  return( g )
}

#' Plot cross-correlation function
#'
#' Attractive \link{ggplot2} plot of \code{\link{ccf}} output  
#' @param x a timeseries
#' @param y a timeseries
#' @param lag.max maximum time lag to display
#' @param type one of "correlation" (default) or "covariance."  
#' The latter is not implemented. 
#' @param demean remove timeseries mean prior to computation of ACF
#' @param ci confidence interval to display
#' @return a ggplot2 object
#' @export
#' @import "ggplot2"
#' @examples
#' x <- rnorm( 1000 )
#' plot.acf( x, lag.max = 20 )

gplot.ccf <- function( x, y,
                       lag.max = NULL, 
                       type = c( "correlation", "covariance" ),
                       demean = TRUE,
                       ci = 0.95 ) {
  
  type <- match.arg( type )
  r <- ccf( x, y, lag.max = lag.max, type = type, demean = demean, plot = FALSE )
  R <- data.frame( Lag = r$lag[,1,1], ACF = r$acf[,1,1] )
  n <- length( R$Lag )
  R$angle <- rep( pi / 2, n )
  R$y <- rep( 0, n )
  
  s <- qnorm( ( 1 + ci ) / 2 ) / sqrt( r$n.used )
  R$psigma <- rep( s, n )
  R$nsigma <- rep( -s, n )
  
  g <- ggplot2::ggplot( R, aes( x = Lag, y = y ) ) + 
    ggplot2::geom_spoke( aes( angle = angle, radius = ACF ) ) +
    ggplot2::geom_point( aes( y = ACF ) ) +
    ggplot2::geom_line( aes( y = psigma ), linetype = 2, color = "dodgerblue3" ) +
    ggplot2::geom_line( aes( y = nsigma ), linetype = 2, color = "dodgerblue3" ) +
    ggplot2::labs( y = "CCF" )
  
  return( g )
}