#' Recursive outlier detector/remover
#'
#' Locate and remove outliers using an arima model to whiten time series, Hampel's
#' MAD criterion to identify outliers, and kalman smoothing.  Algorithm
#' is recursive and continues until no further outliers are found. 
#' @param s timeseries
#' @param order order specification of arima model; defaults to AR(2)
#' @param k hampel's half-window size.  Defines window for MAD computation
#' @param tol MAD outlier threshold; defaults to 3 (Pearson's rule)
#' @return cleaned timeseries.
#' @export
#' @import "pracma"
#' @import "imputeTS"
#' @examples
#' s.df <- data.frame( time = 1:100, s = sin( 2 * pi * 1:100 / 20 ) )
#' s.ts <- df2ts( s.df, "time" )
#' 

joutlier <- function( s, order = c( 2, 0, 0 ), k = 10, tol = 3 ) {
  m1 <- arima( s, order )
  l1 <- pracma::hampel( x = residuals( m1 ), k = k, t0 = tol )
  if ( length( l1$ind ) == 0 ) {
    return( s )
  }
  s[l1$ind] <- NA
  s <- imputeTS::na.kalman( s, smooth = TRUE )
  joutlier( s, order, k )
}