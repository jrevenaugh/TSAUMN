#' Plot short convolution filter
#'
#' Attractive \link{ggplot2} plot of short convolution filter coefficients
#' @param h a filter
#' @param type one of "centered" (default), "causal" or "symmetric".
#' If "centered" filter is shifted left to center on origin.  If even length, it
#' will extend further to the right of 0.  If symmetric, h is assumed to provide only
#' the time >= 0 coefficients.  Negative time will be obtained by mirroring.
#' @return a ggplot2 object
#' @export
#' @import "ggplot2"
#' @examples
#' h <- fdfilter( 5 )
#' gplot.filter( h )


gplot.filter <- function( h, type = c( "centered", "causal", "symmetric" ) ) {
  type <- match.arg( type )
  if ( !is.null( dim( h ) ) ) h <- h[,1]
  if ( type == "symmetric" ) h <- c( rev( h ), h[-1] )
  n <- length( h )
  if ( type == "causal" ) {
    t0 <- 0
    tn <- n - 1
  } else {
    t0 <- -( n - 1 ) / 2
    tn <- t0 + n - 1
  }
  t <- t0:tn
  H <- data.frame( Time = t, H = h, y = rep( 0, n ),
                   angle = pi * ( 1 - sign( h ) / 2 ) )

  g <- ggplot2::ggplot( H, aes( x = Time, y = H ) ) +
    ggplot2::geom_spoke( aes( y = y, angle = angle, radius = abs( H ) ) ) +
    ggplot2::geom_point( aes( y = H ), color = "orangered3", size = 2 ) +
    annotate( "line", x = range( t ), y = c( 0, 0 ), linetype = 3 ) +
    ggplot2::labs( y = "Filter" )

  return( g )
}
