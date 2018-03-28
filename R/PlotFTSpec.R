#' Plot fourier/time spectral analysis
#'
#' Plot output of \code{\link{evolFFT}} and \code{\link{evolMTM}} in \code{\link{TSAUMN}}
#' @param evol return value of \code{\link{evolFFT}} or \code{\link{evolMTM}}
#' @param Log plot log10 of power spectrum, defaults to FALSE
#' @param fscale one of "Frequency" (default), "Period" or "Octave" The latter two display period on a log10 scale and frequency on a log2 scale respectively. plot log2 of frequency, defaults to FALSE
#' @param pal color palette, defaults to a dark spectrum
#' @param xlab label for time axis, defaults to "Time"
#' @param ylab label for frequency/period axis, defaults to "Frequency" or "Period" as appropriate
#' @param dynRange orders of magnitude beneath peak to include in color scale.  Defaults to full range.
#' @export
#' @import "RColorBrewer"
#' @import "scales"
#' @examples
#' require( RColorBrewer )
#' # create chirp signal
#' n <- 1024
#' t <- 1:n
#' a <- 4 * n
#' x <- cos( 2 * pi * t^2 / a )
#' # perform time-frequency analysis of chirp
#' evol <- evolFFT( x, deltat = 1, nwin = 200, nstep = 5 )
#' # plot result
#' plotFTSpec( evol )


plotFTSpec <- function( evol, Log = F, fscale = c( "Frequency", "Octave", "Period" ),
                        pal = colorRampPalette( RColorBrewer::brewer.pal( 11, "Spectral" ) )(256),
                        xlab = "Time", ylab = "Frequency",
                        dynRange = 0 )
{
    fscale <- match.arg( fscale )
    aspect <- 1 / 1.2
    par( mar = c( 3, 2.5, 0, 5 ), cex = 1.3 )

    s <- evol$signal
    TF <- TRUE
    n <- length( s )
    dur <- ( n + 1 ) * evol$deltat
    dtinc <- ( evol$tims[2] - evol$tims[1] ) / 2
    evol$tims <- evol$tims - dtinc
    evol$tims <- c( evol$tims, evol$tims[length(evol$tims)] + 2 * dtinc )
    x <- evol$tims / dur
    if ( fscale == "Period" ) {
        temp <- 1 / evol$wpars$fh
        if ( evol$wpars$fl == 0 ) {
          evol$numfreqs <- evol$numfreqs - 1
          evol$freqs <- evol$freqs[-1]
          evol$DSPEC <- evol$DSPEC[-1,]
          evol$wpars$fl <- evol$freqs[1]
        }
        evol$wpars$fh <- 1 / evol$wpars$fl
        evol$wpars$fl <- temp
        evol$freqs <- 1 / evol$freqs
        if ( ylab == "Frequency" ) ylab = "Period"
        TF <- FALSE
    } else if ( fscale == "Octave" ) {
      evol$numfreqs <- evol$numfreqs - 1
      evol$freqs <- evol$freqs[-1]
      evol$DSPEC <- evol$DSPEC[-1,]
      evol$wpars$fl <- evol$freqs[1]
      TF <- FALSE
    }
    fmax <- evol$wpars$fh
    fmin <- evol$wpars$fl
    frange <- evol$wpars$fh - evol$wpars$fl
    lf <- which( evol$freqs >= fmin & evol$freqs <= fmax )
    if ( fscale == "Period" ) {
      y <- log10( evol$freqs[lf] )
      prange <- max( y ) - min( y )
      y <- ( y - min( y ) ) / prange
      y <- y * 0.8 * aspect
      y <- rev( y )
      y <- c( y[1] - ( y[2] - y[1] ) / 2, y )
    } else if ( fscale == "Octave" ) {
      freqs <- evol$freqs[lf]
      dfinc <- ( freqs[2] - freqs[1] ) / 2
      freqs <- c( freqs - dfinc, freqs[length(freqs)] + dfinc )
      freqs <- log2( freqs )
      ymin <- min( freqs )
      ymax <- max( freqs )
      yrange <- ymax - ymin
      y <- ( freqs - ymin ) / yrange * 0.8 * aspect
    } else {
      freqs <- evol$freqs[lf]
      dfinc <- ( freqs[2] - freqs[1] ) / 2
      freqs <- c( freqs - dfinc, freqs[length(freqs)] + dfinc )
      y <- ( freqs - freqs[1] ) / frange * 0.8 * aspect
    }
    s <- s - mean( s )
    s <- s / max( abs( s ) ) * 0.07 * aspect + 0.925 * aspect
    nf <- length( lf )
    dspec <- evol$DSPEC
    if ( fscale == "Period" ) dspec <- dspec[rev( lf ), ]
    else dspec <- dspec[lf, ]

    if ( dynRange > 0 ) {
        Min <- max( dspec ) / dynRange
        ll <- which( dspec < Min, arr.ind = TRUE )
        dspec[ll] <- Min
    }
    if ( Log ) dspec = log10( dspec )

    plot( c(0, 1), c(0, 1), axes = F, type = "n", xlab = NA, ylab = NA )
    rect( 0, 0.0, 1, 0.8 * aspect, col = "lightgray", border = NA )
    image( x, y, t( dspec ), col = pal, axes = F, add = T, xlab = NA, ylab = NA, useRaster = TF )
    lines( seq( 0, 1, length.out = n ), s, col = "black" )
    rect( 0, 0.85 * aspect, 1, 1 * aspect, col = NA, border = "black" )
    rect( 0, 0.0, 1, 0.8 * aspect, col = NA, border = "black" )
    if ( evol$wpars$nwin > 0 ) lines( c( 0, evol$wpars$nwin / n ), c( 0, 0 ), lwd = 5, xpd = T )

# Add axes

    nT <- 5
    xTic <- pretty( evol$tims, n = nT )
    trange <- range( evol$tims )
    lX <- which( xTic >= trange[1] & xTic <= trange[2] )
    xTic <- xTic[lX]
    nX <- length( xTic )

    if ( fscale == "Period" ) {
      nY <- floor( log10( fmax ) - log10( fmin ) )
      nT <- nY
      yTic <- scales::trans_breaks( "log10", function(x) 10^x, n = nT )(c( fmin, fmax ))
      nY <- length( which( yTic >= fmin & yTic <= fmax ) )
      while( nY > 5 ) {
        yTic <- scales::trans_breaks( "log10", function(x) 10^x, n = nT )(c( fmin, fmax ))
        nY <- length( which( yTic >= fmin & yTic <= fmax ) )
        nT <- nT - 1
      }
    } else if ( fscale == "Octave" ) {
      nY <- floor( log2( fmax ) - log2( fmin ) )
      nT <- nY
      yTic <- scales::trans_breaks( "log2", function(x) 2^x, n = nT )(c( fmin, fmax ))
      nY <- length( which( yTic >= fmin & yTic <= fmax ) )
      while( nY > 5 ) {
        yTic <- scales::trans_breaks( "log2", function(x) 2^x, n = nT )(c( fmin, fmax ))
        nY <- length( which( yTic >= fmin & yTic <= fmax ) )
        nT <- nT - 1
      }
    } else {
      yTic <- pretty( freqs, n = nT )
    }
    yTic <- signif( yTic, digits = 1 )
    lY <- which( yTic >= fmin & yTic <= fmax )
    yTic <- yTic[lY]
    nY <- length( yTic )

    xval <- rep( "text", nX )
    xat <- rep( 0, nX )
    for ( i in 1:nX ) {
      xval[i] = sprintf( "%g", xTic[i] )
      xat[i] <- ( xTic[i] - evol$tims[1] ) / diff( trange )
    }
    axis( 1, at = xat, labels = xval, pos = 0 )
    mtext( side = 1, at = 0.5, text = xlab, line = 2, cex = 1.3 )

    yval <- rep( "text", nY )
    yat <- rep( 0, nY )
    if ( fscale == "Period" ) {
      fmin <- log10( fmin )
      fmax <- log10( fmax )
      prange <- fmax - fmin
    } else if ( fscale == "Octave" ) {
      fmin <- log2( fmin )
      fmax <- log2( fmax )
      prange <- fmax - fmin
    }
    for ( i in 1:nY ) {
        yval[i] = sprintf( "%g", yTic[i] )
        if ( fscale == "Period" ) {
          yat[i] <- ( log10( yTic[i] ) - fmin ) / prange * 0.8 * aspect
        } else if ( fscale == "Octave" ) {
          yat[i] <- ( log2( yTic[i] ) - fmin ) / prange * 0.8 * aspect
        } else {
          yat[i] <- 0.8 * yTic[i] / fmax * aspect
        }
    }

    mtext( side = 2, at = 0.4 * aspect, text = ylab, line = 1.5, cex = 1.3 )
    axis( 2, at = yat, labels = yval, pos = 0 )

# Add legend

    xl <- 1.02
    xu <- 1.09
    yl <- 0.00
    yu <- 0.8 * aspect
    i <- seq( along = pal )
    dy <- ( yu - yl ) / length( i )
    y1 <- yl + ( i - 1 ) * dy
    y2 <- y1 + dy
    rect( xl, y1, xu, y2, col = pal, xpd = T, border = NA )
    rect( xl, yl, xu, yu, col = NA, border = "black", xpd = T )
    la <- range( dspec )
    rlab = paste(sep = " ", format.default(la[1], digits = 2, scientific = !Log ) )
    if ( nchar( rlab ) > 4 ) rlab = paste(sep = " ", format.default(la[1], digits = 2, scientific = T ) )
    text( xu, yl, labels = rlab, xpd = T, pos = 4, cex = 0.75 )
    rlab = paste(sep = " ", format.default(la[2], digits = 2, scientific = !Log  ) )
    if ( nchar( rlab ) > 4 ) rlab = paste(sep = " ", format.default(la[2], digits = 2, scientific = T ) )
    text( xu, yu, labels = rlab, xpd = T, pos = 4, cex = 0.75 )
    rlab = "Power"
    if ( Log ) rlab = "Log Power"
    text( xu + 0.05, 0.4 * aspect, rlab, xpd = T, srt = -90 )

    invisible(  )
}
