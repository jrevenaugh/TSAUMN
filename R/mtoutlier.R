#' Martin-Thomson outlier removal
#'
#' Martin-Thomson filter-cleaning of a (presumed) AR timeseries  
#' @param y timeseries; no NA allowed
#' @param ar c(min, max) AR order in modeling; defaults to c(1, 10)
#' @param tol outlier threshold value; defaults to 3 (Pearson's Rule)
#' @param set0 replace outlier with 0 (TRUE, default) or scaled tol (FALSE)
#' @param ... other arguments passed to \code{\link{arima}}
#' @return cleaned timeseries
#' @export
#' @examples
#' y <- arimaSim( 200, ar = c( 0.5, -0.1 ) )
#' # Add a few outliers
#' o <- floor( runif( 5, min = 5, max = 195 ) )
#' y[o] <- runif( 5, min = 4, max = 8 )
#' yc <- mtoutlier( y )
#' plot( y )
#' lines( yc, col = "red" )
#' 
mtoutlier <- function( y, ar = c( 1, 10 ), tol = 3, set0 = TRUE, ...) {
  
  # Find best-fit AR(p) process in given range
  prange <- ar[1]:ar[2]
  np <- length( prange )
  mod <- vector( "list", np )
  aic <- rep( 0, np )
  for ( i in prange ) {
    mod[[i]] <- arima( y, order = c( i, 0, 0 ), ... )
    aic[i] <- AIC( mod[[i]] )
  }
  p <- which.min( aic )
  
  # First MT pass
  yc <- mtfc( y, p, mod[[p]], tol, set0 )
  
  # Refit arima to cleaned series
  for ( i in prange ) {
    mod[[i]] <- arima( yc, order = c( i, 0, 0 ), ... )
    aic[i] <- AIC( mod[[i]] )
  }
  p <- which.min( aic )
  
  # Second MT pass
  yc <- mtfc( y, p, mod[[p]], tol, set0 )

  return( yc )
}


Psi <- function( tau, tol = 3, set0 = TRUE ) {
  if ( abs( tau ) > tol ) {
    if ( set0 ) tau <- 0
    else tau <- sign( tau ) * tol
  }
  return( tau )
}

w <- function( tau, tol = 3 ) {
  if ( abs( tau ) < tol ) w <- 1
  else w <- sign( tau ) * tol / tau
  return( w )
}

mtfc <- function( y, p, mod, tol = 3, set0 = TRUE ) {
  Phi <- diag( p - 1 )
  Phi <- cbind( Phi, rep( 0, p - 1 ) )
  Phi <- rbind( mod$coef[1:p], Phi )
  ymean <- mod$coef[p + 1]
  Q <- matrix( rep( 0, p^2), nrow = p )
  Q[1,1] <- mod$sigma2
  
  M <- Q
  Xt <- y[p:1]
  y <- y - ymean
  yc <- y
  P <- 0 * M
  
  n <- length( y )
  for ( i in (p + 1):n ) {
    yhat <- ( Phi %*% Xt )[1]
    st <- sqrt( M[1,1] )
    tau <- ( y[i] - yhat ) / st
    Xt <- Phi %*% Xt + M[,1] / st * Psi( tau, tol, set0 )
    P <- M - w( tau, tol ) * M[,1] %o% M[,1] / st^2
    M <- Phi %*% P %*% t( Phi ) + Q
    yc[i] <- Xt[1]
  }
  return( yc + ymean )
}

