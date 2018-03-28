#' Newton-Cotes timeseries integration
#'
#' Integrate a timeseries using stable Newton-Cotes filters.  
#' @param s timeseries
#' @param order Newton-Cotes filter order ( 2 <= order <= 7 ).  Beware that
#' high-frequence content can become unstable for orders > 2 (order = 2 is the trapezoid rule).
#' @param deltat sampling interval; used only if s is not a "ts" object.  Defaults to 1.
#' @param nfac oversampling factor.  Must be integer >= 1. Timeseries s is
#' resampled at sampling interval deltat/nfac prior to integration.  Does not
#' affect timing or duration of returned integrated series.  nfac > 1 provides greater
#' accuracy for timeseries with energetic high frequency components.
#' @return For s(t0 + i*deltat) for 0 <= i < length( s ), returns the definite integral
#' from t0 to t0 + i*deltat.  This is equivalent to the indefinite integral of s
#' with an arbitrary baseline.
#' @export
#' @examples
#' t <- 1:200
#' s <- sin( 2 * pi * t / 20 )
#' ints.true <- -cos( 2 * pi * t / 20 ) / ( 2 * pi / 20 )
#' ints <- ncint( s, order = 2, deltat = 1 )
#' offset <- ints.true[1] - ints[1]
#' plot( t, ints.true, type = "l" )
#' lines( t, ints + offset, col = "red" )
#' 
ncint <- function( s, order = 2, deltat = 1, nfac = 1 ) {
  if ( !is.ts( s ) ) s <- ts( s, deltat = deltat )
  if ( order < 2 || order > 7 ) 
    stop( "Order out of range 2:7." )
  if ( floor( nfac ) != nfac || nfac < 1 ) 
    stop( "Oversample factor must be integral and >= 1" )
  order <- floor( order )
  dt <- deltat( s ) / nfac
  ns <- length( s )
  t <- time( s )
  tr <- seq( t[1], t[ns] + 1 - dt, dt )
  if ( nfac != 1 ) s <- signal::resample( s, deltat( s ), dt )
  n <- length( s )
  fMat <- newcotes( order )
  intg <- rep( 0, n )
  if ( order > 2 ) {
    for ( i in 2:(order-1) ) {
      intg[i] <- t( fMat[i,1:i] ) %*% s[1:i]
    }
  }
  intg[order] <- t( fMat[order,] ) %*% s[1:order]
  for ( i in (order + 1):n ) {
    intg[i] <- t( fMat[order,] ) %*% s[(i-order+1):i] + intg[i - order + 1]
  }
  intg <- intg * dt
  indx <- seq( 1, by = nfac, length.out = ns )
  return( ts( intg[indx], start = t[1], deltat = deltat( s ) ) )
}