#' Plot fourier/time spectral analysis
#'
#' Plot output of \code{\link{evolfft}} and \code{\link{evolMTM}} in \code{\link{RSEIS}}  
#' @param evol return value of \code{\link{evolfft}} or \code{\link{evolMTM}}
#' @param Log plot log10 of power spectrum, defaults to FALSE
#' @param Period plot spectrum against period rather than frequency, defaults to FALSE
#' @param pal color palette, defaults to a dark spectrum
#' @param xlab label for time axis, defaults to "Time"
#' @param ylab label for frequency/period axis, defaults to "Frequency" or "Period" as appropriate
#' @param dynRange orders of magnitude beneath peak to include in color scale.  Defaults to full range.
#' @export
#' @import "RColorBrewer"
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


plotFTSpec <- function( evol, Log = F, Period = F, 
                        pal = colorRampPalette( RColorBrewer::brewer.pal( 11, "Spectral" ) )(256), 
                        xlab = "Time", ylab = "Frequency", 
                        dynRange = 0 )
{
    aspect <- 1 / 1.2
    par( mar = c( 3, 2.5, 0, 5 ), cex = 1.3 )

    s <- evol$signal
    TF <- TRUE
    n <- length( s )
    dur <- n * evol$deltat
    x <- evol$tims / dur
    if ( Period ) {
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
    }
    fmax <- evol$wpars$fh
    fmin <- evol$wpars$fl
    frange <- fmax - fmin
    lf <- which( evol$freqs >= fmin & evol$freqs <= fmax )
    y <- ( evol$freqs[lf] - fmin ) / frange * 0.8 * aspect
    if ( Period ) y <- rev( y )
    s <- s - mean( s )
    s <- s / max( abs( s ) ) * 0.07 * aspect + 0.925 * aspect
    nf <- length( lf )
    dspec <- evol$DSPEC
    if ( Period ) dspec <- dspec[rev( lf ), ]
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
    xInc <- ticInt( dur, nTic = nT )
    yInc <- ticInt( frange, nTic = nT )
    
    nT <- dur / xInc + 1
    if ( ( nT - 1 ) * xInc > 1.01 * dur ) nT <- nT - 1
    xval <- rep( "text", nT )
    xat <- rep( 0, nT )
    for ( i in 1:nT ) { 
      xval[i] = sprintf( "%g", 0 + ( i - 1 ) * xInc )
      xat[i] <- ( i - 1 ) / ( nT - 1 )
    }
    axis( 1, at = xat, labels = xval, pos = 0 )
    mtext( side = 1, at = 0.5, text = xlab, line = 2, cex = 1.3 )

    nT <- frange / yInc + 1
    if ( ( nT - 1 ) * yInc > 1.01 * fmax ) nT <- nT - 1
    yval <- rep( "text", nT )
    yat <- rep( 0, nT )
    for ( i in 1:nT ) { 
        yval[i] = sprintf( "%g", fmin + ( i - 1 ) * yInc )
        yat[i] <- 0.8 * ( i - 1 ) / ( nT - 1 ) * aspect
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

ticInt <- function( range, nTic = 5 )
{
    uTS <- range / ( nTic - 1 )
    x <- ceiling( log10( uTS ) - 1 )
    pow10x <- 10^x
    rTR <- floor( uTS / pow10x ) * pow10x
    return( rTR )
}