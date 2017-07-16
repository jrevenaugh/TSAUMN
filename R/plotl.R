#' Shorthand line plotting
#'
#' Equivalent to plot( ..., type = "l" )  
#' @param ... Arguments passed to \code{\link{plot}}, including \code{\link{graphical parameters}} (see \code{\link{par}} ).
#' @export
#' @examples
#' x <- 1:50
#' y <- sin( 2 * pi * x / 10 ) + rnorm( 50 )
#' plotl( x, y, col = "blue" )

plotl <- function( ... ) { plot( ..., type = "l" ) }