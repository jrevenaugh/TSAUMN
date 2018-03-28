#' FIR Causal Wiener-Hopf Filter Estimation
#'
#' Obtain the optimal (MMSE) Wiener-Hopf causal FIR filter for an known multivariate
#' input ts x and known target ts s.  Uses the Trench algorithm to compute the log-likelihood
#' of the residual ts and AIC to choose optimal filter length (within a range provided by the user).
#' Full AIC calculation not practical for signals longer than 10,000 samples.
#' @param x a ts, vector or matrix of values (columns as time series).
#' @param s as ts or vector of values
#' @param n filter length if integer or range of possible lengths if length two vector.
#' @param AIC use full AIC calculation (default) or minimum squared error to choose filter.
#' @return List with elements h, the matrix of optimal filter values (stored as columns); n, the (range of) filter length(s) tested;
#' and AIC the AIC of tested filters (corresponds to n).
#' @import "ltsa"
#' @export
#' @examples
#' x <- rnorm( 500 )
#' h <- seq( 1, 0, length.out = 5 )
#' s <- stats::filter( x, h, sides = 1 )
#' x <- x[-(1:4)] # Remove filter roll-on
#' s <- s[-(1:4)]
#' wh <- whfir( x, s, n = c(1, 10 ) )
#' plot( h, ylab = "Filter Coefficients" )
#' lines( wh$h[,1], col = "red" )
#'
whfir <- function( x, s, n = NULL, AIC = TRUE ) {
  if ( is.null( x ) ) stop( "Null x" )
  if ( is.null( s ) ) stop( "Null s" )

  multi <- is.matrix( x )
  if ( !multi ) {
     x <- matrix( x, ncol = 1 )
  }
  lx <- dim( x )[1]
  m <- dim( x )[2]

  nMin <- min( lx, length( s ) )
  x <- x[1:nMin,]
  x <- matrix( x, nrow = nMin )
  s <- s[1:nMin]
  if ( anyNA( x ) || anyNA( s ) ) stop( "NAs found." )

  if ( is.null( n ) ) {
    n <- c( 1, floor( nMin / 5 ) )
  }
  if ( length( n ) == 1 ) n <- c( n, n )
  n <- sort( n )
  if ( n[2] > nMin / 5 ) {
    n[2] <- floor( nMin / 5 )
    message( "Maximum filter length (n) coerced to ", format( n[2] ) )
    n[1] <- min( n[1], n[2] )
  }
  XX <- array( data = NA, dim = c( n[2] + 1, m, m ) )
  SX <- array( data = NA, dim = c( n[2] + 1, m ) )
  l <- seq( n[2] + 1, 2 * n[2] + 1 )

  for ( i in 1:m ) {
    for ( j in 1:m ) {
      CF <- ccf( x[,i], x[,j], type = "covariance", lag.max = n[2], plot = FALSE )
      XX[,i,j] <- CF$acf[l,1,1]
    }
    CF <- ccf( s, x[,i], type = "covariance", lag.max = n[2], plot = FALSE )
    SX[,i] <- CF$acf[l,1,1]
  }

  n_ind <- n[1]:n[2]
  Loglike <- rep( 0, length( n_ind ) )
  nbest <- n[1]
  Loglikemin <- .Machine$double.xmax
  lskp <- c( 1:n[2] - 1 )

  for ( i in n[1]:n[2] ) {
    X <- NULL
    S <- c( SX[1:i,] )
    for ( j in 1:m ) {
      Cxx <- toeplitz( XX[1:i,j,1] )
      if ( m > 1 ) for ( k in 2:m ) Cxx <- cbind( Cxx, toeplitz( XX[1:i, j, k] ) )
      X <- rbind( X, Cxx )
    }

    hest <- solve( X, S )
    sest <- stats::filter( x[,1], hest[1:i], sides = 1 )
    fs <- seq( 1, length( hest ), i )
    fe <- seq( i, length( hest ), i )
    if ( m > 1 ) for ( j in 2:m ) sest <- sest + stats::filter( x[,j], hest[fs[j]:fe[j]], sides = 1 )

    l <- which( !is.na( sest ) )
    ds <- s[-lskp] - sest[-lskp]
    j <- i - n[1] + 1

    if ( AIC ) {
      ads <- acf( ds, lag.max = length( ds ) - 1, type = "covariance", plot = FALSE )$acf[,1,1]
      Loglike[j] <- 2 * i * m - 2 * ltsa::TrenchLoglikelihood( ads, ds )
    } else {
      Loglike[j] <- sum( ds^2 )
    }
    if ( Loglike[j] < Loglikemin ) {
      h <- hest
      Loglikemin <- Loglike[j]
      nbest <- i
    }
  }

  return( list( h = matrix( h, nrow = nbest ), n = n_ind, AIC = Loglike ) )
}
