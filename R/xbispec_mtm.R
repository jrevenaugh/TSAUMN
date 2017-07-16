#' Cross-bispectrum
#'
#' Compute cross-bispectrum with the multitaper method
#' @param x regularly sampled univariate time series.  Can be passed as a ts 
#' object, a data.frame or simple vector of values.  In the case of the latter, deltat must be specified.
#' If passed as a data.frame, the first column must provide timing, the second the data itself.
#' @param y regularly sampled univariate time series.  Can be passed as a ts 
#' object, a data.frame or simple vector of values.  In the case of the latter, deltat must be specified.
#' If passed as a data.frame, the first column must provide timing, the second the data itself.  Must share sample times with x.
#' @param deltat time interval between samples, used only if x is passed as a vector.  Defaults to 1.
#' @param nw multitaper frequency bandwidth parameter, defaults to 4.
#' @param k number of tapers, typically 2 * nw - 1, defaults to 7.
#' @param nPoly order of the polynomial used for DC/trend removal, defaults to 3.
#' @return list containing frequency, bispectrum, bispectral power, and bicoherence.
#'  \item{"freq1"}{f1 frequencies (0 to Nyquist)}
#'  \item{"freq2"}{f12 frequencies (-Nyquist to Nyquist / 2)}
#'  \item{"bispec"}{Complex bispectrum [f1, f2]}
#'  \item{"coher"}{bicoherence [f1, f2]}
#'  \item{"bpspec"}{bispectral power [f1,f2]}
#' @import "multitaper"
#' @export
#' @examples
#' t <- 1:256
#' f <- c( 1/20, 1/5 )
#' s <- sin( 2 * pi * t %o% f )
#' x <- apply( s, 1, sum ) + rnorm( length( t ), sd = 0.25 )
#' y <- sin( 2 * pi * t / 4 ) + rnorm( length( t ), sd = 0.25 )
#' b <- xbispec.mtm( x = x, y = y, deltat = 1 )
#' image( b$freq1, b$freq2, t( b$coher ), useRaster = TRUE )
xbispec.mtm <- function( x, y, deltat = 1, nw = 4, k = 7, nPoly = 3 ) {
  
  # Determine what was passed for x and establish times and dt.
  if ( is.data.frame( x ) ) {
    t <- x[,1]
    s <- x[,2]
    dt <- t[2] - t[1]
  }
  else if ( is.ts( x ) ) {
    s <- as.numeric( x )
    dt <- deltat( x )
  } else if ( is.numeric( x ) ) {
    dt <- deltat
    s <- x
  }
  nPts <- length( s )
  s <- residuals( lm( s ~ poly( seq_along( s ), nPoly ) ) )
  
  # Determine what was passed for y and insure it matches x.
  if ( is.data.frame( y ) ) {
    sy <- y[,2]
  }
  else if ( is.ts( y ) ) {
    sy <- as.numeric( y )
  } else if ( is.numeric( y ) ) {
    sy <- y
  }
  if ( length( sy ) != nPts ) stop( "Time series lengths are not equal." )
  nPts <- length( s )
  if ( length( sy ) != nPts ) stop( "Time series lengths are not equal." )
  s <- residuals( lm( s ~ poly( seq_along( s ), nPoly ) ) )
  sy <- residuals( lm( sy ~ poly( seq_along( sy ), nPoly ) ) )
  
  if ( k > ( 2 * nw + 1 ) ) print( noquote( "Consider fewer tapers or larger time/bandwidth parameter." ) )
  if ( nPts * k >= 2000 ) print( noquote( "Be patient--this will take a while." ) )
  
  # Generate the DPSS tapers.
  tap <- multitaper::dpss( n = nPts, k = k, nw = nw )
  taper <- tap$v
  eign <- tap$eigen
  
  # Loop over tapers, computing the FFT of each tapering of the time series
  # Time series are padded to the next multiple of 2.
  nPad <- nextn( nPts, 2 )
  spec <- matrix( rep( complex( 0, 0 ), k * nPad ), ncol = k )
  ypec <- matrix( rep( complex( 0, 0 ), k * nPad ), ncol = k )
  
  stap <- rep( 0.0, nPad )
  for ( i in 1:k ) {
    stap[1:nPts] <- s * taper[,i]
    spec[,i] <- fft( stap )
    stap[1:nPts] <- sy * taper[,i]
    ypec[,i] <- fft( stap )
  }

  # Create the various returned matrices
  nFl <- nPad / 2 + 1
  nFs <- nPad / 4 + 1
  l2 <- c( -( nFl:2), 1:nFs )
  m2 <- seq_along( l2 )
  lFs <- length( l2 )
  freqs <- seq( 0, 0.5 / dt, length.out = nFl )
  bcoher <- matrix( rep( -Inf, lFs * nFl ), nrow = lFs )
  bpspec <- bcoher
  bispec <- matrix( rep( complex( 0, 0 ), lFs * nFl ), nrow = lFs )
  
  Plist <- vector( "list", k )
  Pmat <- matrix( rep( 0, k^2 ), nrow = k )
  Psum <- 0
  for ( k1 in 1:k ) {
    Plist[[k1]] <- Pmat
    for ( k2 in 1:k ) {
      for ( k3 in 1:k ) {
        Plist[[k1]][k2,k3] <- sum( taper[,k1] * taper[,k2] * taper[,k3] )
        Psum <- Psum + Plist[[k1]][k2,k3]^2
      }
    }
  }  

  # Loop over frequencies (2)
  fp <- matrix( rep( 0, k^2 ), nrow = k )
  M21 <- rep( 0, nFl )
  for ( i1 in 1:nFl ) {
    M21[i1] <- sum( Mod( spec[i1, ] )^2 ) / ( k * nPts )
    for ( i in 1:lFs ) {
      i2s <- l2[i]
      md2 <- m2[i]
      i2 <- abs( i2s )
      i12 <- i1 + i2s - 1
      if ( i12 > nFl || i12 < 1 || i2 > i1 ) next()
      M22 <- M21[i2]
      M2s <- sum( Mod( spec[i12,] )^2 ) / ( k * nPts )
      M3 <- 0
      fp <- spec[i1,] %o% spec[i2,]
      for ( k3 in 1:k ) {
         M3 <- M3 + Conj( ypec[i12,k3] ) * sum( Plist[[k3]] * fp )
      }
      M3 <- M3 / ( nPts^2 * Psum )
      bispec[md2,i1] <- M3
      bcoher[md2,i1] <- Mod( M3 / sqrt( M21[i1] * M22 * M2s ) )^2
      bpspec[md2,i1] <- Mod( M3 )^2
    }
  }
  
  return( list( bspec = bispec, 
                coher = bcoher, 
                bpspec = bpspec, 
                freq1 = freqs,
                freq2 = c( -rev( freqs ), freqs[2:(( nFl - 1 ) / 2 + 1 )]) ) )
  
}