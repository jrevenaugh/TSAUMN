#' Fast Hadamard Transform
#'
#' "Natural" order Hadamard transform.  
#' @param x vector of univariate sequence values; zero padded to next multiple of 2
#' @return vector of Hadamard transform coefficients
#' @export
#' @examples
#' a <- arimaSim( n = 124, ar = 0.9, int = 1 )
#' H <- fht( a )
#' plot( H, type = "l", xlab = "Index", ylab = "Coefficient" )
#' 
fht <- function ( x ) {
  # If x isn't a power of two in length, pad it out.
  s <- as.numeric( x )
  l <- length( s )
  n <- nextn( l, 2 )
  l2n <- log2( n )
  if ( l < n ) s[( l + 1 ):n] <- 0
  ld <- n / 2 - 1
  l1 <- 1 + 0:ld * 2
  l2 <- l1 + 1
  nrm <- 1 / sqrt( 2 ) # some fht ignore this.  
  for ( i in 1:l2n ) {
    p <- s[l1] + s[l2]
    m <- s[l1] - s[l2]
    s <- c( p, m ) * nrm
  }
  return( s )
}

#' Reorder vector into sequency order
#'
#' Reorder vector into sequency order for Walsh transform using Fast Hadamard Transform \code{link{fht}}.  
#' @param n number of samples in sequence
#' @return permutation index
#' @keywords internal
#' @export
#' 
sequency <- function( n ) {
  n <- log2( length( x ) )
  if ( floor( n ) != n ) {
    print( noquote( "Length not a factor of 2." ) )
  }
  if ( n == 0 ) return( 0 )
  else if ( n == 1 ) return( c( 0, 1) )
  l <- 2^n
  s <- 0:( l - 1 )
  index <- rep( 0, 1 )
  for ( i in 1:l ) {
    b <- decimal2binary( s[i], n )
    g <- gray2binary( rev( b ) )
    d <- binary2decimal( g )
    index[i] <- d
  }
  return( index )
}

#' Fast Walsh Transform
#'
#' Fast Walsh transform (aka "Sequency" order Hadamard)  
#' @param x vector of univariate sequence values; zero padded to next multiple of 2
#' @return vector of Walsh transform coefficients
#' @export
#' @examples
#' a <- arimaSim( n = 124, ar = 0.9, int = 1 )
#' H <- fwt( a )
#' plot( H, type = "l", xlab = "Index", ylab = "Coefficient" )
#' 
fwt <- function( x ) {
  l <- length( x )
  n <- nextn( l, 2 )
  if ( l < n ) x[( l + 1 ):n] <- 0
  n <- as.integer( log2( n ) )
  p <- sequency( n ) + 1
  y <- fht( x[p] )
  return( y )
}

#' decimal2binary
#'
#' Convert non-negative integer to binary
#' @param x non-negative integer to convert
#' @param length length of binary sequence to output
#' @return binary representation as vector of digits (0,1)
#' @keywords internal
#' @export
#' 
decimal2binary <- function( x, length ) 
{
  x <- as.integer( x )
  b <- if ( missing( length ) ) 
    NULL
  else rep( 0, length )
  i <- 0
  while ( x >= 1 ) {
    i <- i + 1
    b[i] <- x%%2
    x <- x%/%2
  }
  return( rev( b ) )
}

#' gray2binary
#'
#' Convert gray code to binary
#' @param x gray code as vector of digits (0,1)
#' @return binary representation as vector of digits (0,1)
#' @keywords internal
#' @export
#' 
gray2binary <- function( x ) 
{
  x <- as.logical( x )
  n <- length( x )
  b <- vector( mode = "logical", length = n )
  b[1] <- value <- x[1]
  if ( n > 1 ) 
    for ( i in 2:n ) {
      if ( x[i] ) 
        value <- !value
      b[i] <- value
    }
  b <- as.numeric( b )
  return( b )
}

#' binary2decimal
#'
#' Convert binary digits vector to non-negative integer value
#' @param x binary representation as vector of digits (0,1)
#' @return integer value
#' @keywords internal
#' @export
#' 
binary2decimal <- function( x ) 
{
  sum( x * 2^( rev( seq( along = x ) ) - 1 ) )
}