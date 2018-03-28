#' Newton-Cotes integration filters
#'
#' Compute  matrix of stable Newton-Cotes integration filters up to order 7  
#' @return filter matrix
#' 
newcotes <- function( n = 2 ) {
  if ( floor( n ) != n || n < 2 || n > 7 ) stop( "Order n must be an integer in [2,7]" )
  m <- 7
  f <- matrix( rep( 0, m^2 ), nrow = m )
  f[2,] <- 1 / 2 * c( 1, 1, 0, 0, 0, 0, 0 )
  f[3,] <- 1 / 3 * c( 1, 4, 1, 0, 0, 0, 0 )
  f[4,] <- 3 / 8 * c( 1, 3, 3, 1, 0, 0, 0 )
  f[5,] <- 4 / 104 * c( 11, 26, 31, 26, 11, 0, 0 )
  f[6,] <- 5 / 336 * c( 31, 61, 76, 76, 61, 31, 0 )
  f[7,] <- 1 / 14 * c( 7, 12, 15, 16, 15, 12, 7 )
  return( f[1:n,1:n] )
}