#' Rauch–Tung–Striebel Kalman smoothing
#'
#' RTS kalman filter (smooth) a univariate timeseries with F and H = 1.  Simple function
#' to demonstrate kalman smoothing.
#' @param z univariate timeseries (ts) or single vector of values
#' @param q estimated process noise variance (V in StructTS$model)
#' @param r estimated measurement noise variance (h in StructTS$model)
#' @return smoothed timeseries.
#' @export
#' @examples
#' z <- ts( sin( 2 * pi * 1:200 / 20 ) + rnorm( 200 ), deltat = 1 )
#' s <- RTSsmooth( z, q = 0.1, r = 1 )
#' plot( z, type = "l" )
#' lines( s, col = "red" )
RTSsmooth <- function( z, q = 0.5, r = 0.5 ) { 
  # Figure out what type of object was passed
  if ( is.ts( z ) ) {
    deltat <- deltat( z )
    st <- time( z )[1]
  } else {
    deltat <- 1
    st <- 0
  }
  # Initiate state variables.
  n <- length( z )
  x <- rep( 0, n )
  p <- rep( 0, n )
  xh <- rep( 0, n )
  ph <- rep( 0, n )
  xn <- rep( 0, n )
  pn <- rep( 0, n )
  x[1] <- z[1]
  xh[1] <- x[1]
  p[1] <- q
  ph[1] <- p[1]
  
  # Kalman filter the series forward in time.
  for ( k in 2:n ) {
    xh[k] <- x[k - 1]
    ph[k] <- p[k - 1] + q
    y <- z[k] - xh[k]
    s <- ph[k] + r
    g <- ph[k] / s
    x[k] <- xh[k] + g * y
    p[k] <- ( 1 - g ) * ph[k]
  }
  pn[n] <- p[n]
  xn[n] <- x[n]
  
  # Now RTS filter the state space variables backward in time.
  for ( k in (n-1):1 ) {
    c = p[k] / ph[k+1]
    xn[k] <- x[k] + c * ( xn[k+1] - xh[k+1] )
    pn[k] <- p[k] + c^2 * ( pn[k+1] - ph[k+1])
  }
  
  return( ts( xn, deltat = deltat, start = st ) )
}
