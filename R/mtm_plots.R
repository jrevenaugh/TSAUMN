#' Plot results from \code{\link{spec.mtm}}
#'
#' Attractive \link{ggplot2} plot of \code{\link{spec.mtm}} output  
#' @param mspec \code{\link{spec.mtm}} output
#' @param jack include jackknife confidence intervals if TRUE.  These are computed by \code{\link{spec.mtm}}
#' only if "jackknife = TRUE".  The confidence level defaults to 95%.  This can be changed in \code{\link{spec.mtm}} 
#' with the jkCIProb argument.
#' @param period display spectrum against period, rather than frequency, if TRUE 
#' @param trans apply a transform to the spectral power; defaults to no transform
#' @return a ggplot2 object
#' @export
#' @import "ggplot2"
#' @import "scales"
#' @examples
#' s <- ts( sin( 2 * pi * 1:512 / 32 ) + rnorm( 512, sd = 0.25 ), deltat = 1 )
#' mspec <- spec.mtm( s, jackknife = TRUE, plot = FALSE )
#' mtm.plot( mspec, jack = TRUE, period = TRUE, trans = "log10" )
mtm.plot <- function( mspec, jack = FALSE, period = F, trans = c( "identity", "log10", "sqrt" ) )
{
  pTrans <- match.arg( trans )
  f <- mspec$freq
  s <- mspec$spec
  if ( period ) {
      f <- 1.0 / mspec$freq[-1]
      s <- mspec$spec[-1]
  }
  s.max <- max( s )
  s.min <- min( s )
  spec.df <- data.frame( Frequency = f, Power = s )
  if ( is.null( mspec$mtm$jk ) ) {
    if ( jack ) warning( noquote( "Must set jackknife = TRUE in spec.mtm" ) )
    jack <- FALSE
  }
  if ( jack ) {
    if ( period ) {
       uci <- mspec$mtm$jk$upperCI[-1]
       lci <- mspec$mtm$jk$lowerCI[-1]
    }
    else {
      uci <- mspec$mtm$jk$upperCI
      lci <- mspec$mtm$jk$lowerCI
    }
    s.max <- max( mspec$mtm$jk$maxVal, uci )
    s.min <- min( mspec$mtm$jk$minVal, lci )
    jk.df <- data.frame( UCI = uci, LCI = lci )
    spec.df <- cbind( spec.df, jk.df )
  }
  g <- ggplot2::ggplot( spec.df, aes( x = Frequency, y = Power ) ) 
  
  if ( jack ) g <- g + ggplot2::geom_ribbon( aes( ymin = LCI, ymax = UCI ), fill = "gray80" )
  g <- g + ggplot2::geom_line()
  
  if ( pTrans == "identity" ) {
       ys <- NULL 
       ll <- labs( y = "Power" )
  }
  else if ( pTrans == "log10" ) {
       ys <- ggplot2::scale_y_continuous(trans = pTrans, limits = c( s.min, s.max ), 
                         breaks = scales::trans_breaks(pTrans, function(x) 10^x),
                         labels = scales::trans_format(pTrans, scales::math_format(10^.x) ) )
       ll <- ggplot2::labs( y = "Log Power" )
}
  else {
       ys <- ggplot2::scale_y_continuous(trans = pTrans, limits = c( s.min, s.max ),
                             breaks = scales::trans_breaks(pTrans, function(x) x^2 ),
                             labels = scales::trans_format(pTrans, scales::math_format(.x) ) )
       ll <- ggplot2::labs( y = "Root Power" )
       
  }                    
  g <- g + ys + ll      
  if ( period ) {
    xs <- ggplot2::scale_x_continuous(trans = "log10",
                             breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10", scales::math_format(10^.x) ) )
    g <- g + xs + labs( x = "Period" ) + ggplot2::annotation_logticks( sides = "tb" )
  }
  return( g )
}

#' Plot F-test results from \code{\link{spec.mtm}}
#'
#' Attractive \link{ggplot2} plot of \code{\link{spec.mtm}} F-test output  
#' @param mspec \code{\link{spec.mtm}} output ("Ftest = TRUE" )
#' @param uci upper confidence interval to highlight in plot
#' @param period display F-test against period, rather than frequency, if TRUE 
#' @param trans apply a transform to the spectral power; defaults to no transform
#' @return a ggplot2 object
#' @export
#' @import "ggplot2"
#' @import "scales"
#' @examples
#' s <- ts( sin( 2 * pi * 1:512 / 32 ) + rnorm( 512 ), deltat = 1 )
#' mspec <- spec.mtm( s, Ftest = TRUE, plot = FALSE )
#' mtm.ftest( mspec )
#' 
mtm.ftest <- function( mspec, uci = 0.95, period = F, trans = c( "identity", "log10", "sqrt" ) )
{
  if( is.null( mspec$mtm$Ftest ) ) {
    stop( "Must include 'Ftest = TRUE' in spec.mtm call.")
  }
  pTrans <- match.arg( trans )
  freq <- mspec$freq
  if ( period ) freq <- 1.0 / freq[-1]
  n <- length( freq )
  f.max <- 1.05 * max( freq )
  if ( period ) f.max <- max( freq ) ^ 1.02
  f <- mspec$mtm$Ftest
  if ( period ) f <- f[-1]
  k <- mspec$mtm$k
  fc <- qf( uci, 2, 2 * k - 2 )
  fq <- c( 0.5, 0.9, 0.95, 0.99 )
  fl <- qf( fq, 2, 2 * k - 2 )
  spec.df <- data.frame( freq, f, rep( fc, n ) ) 
  nms <- c( "Frequency", "FStat", "Crit" )
  for ( i in 1:length( fl ) ) {
    spec.df <- cbind( spec.df, rep( fl[i], n ) )
    nms <- c( nms, sprintf( "FL%d", i ) )
  }
  colnames( spec.df ) = nms
  s.max <- max( spec.df[,-1] )
  s.min <- 0.01 * mean( spec.df$FStat )
  g <- ggplot2::ggplot( spec.df, aes( x = Frequency ) ) + 
    ggplot2::labs( y = "F Statistic" )
  g <- g + ggplot2::geom_line( aes( y = FStat ) )
  g <- g + ggplot2::geom_line( aes( y = Crit ), color = "red" )
  g <- g + ggplot2::geom_ribbon( aes( ymin = s.min, ymax = Crit ), alpha = 0.1 )
  for ( i in 1:length( fl ) ) {
    g <- g + ggplot2::geom_line( aes_string( y = nms[i + 3] ), lty = 3 )
    l <- format( fq[i] )
    g <- g + ggplot2::annotate( "text", x = f.max, y = fl[i], label = l )
  }

  if ( pTrans == "sqrt" ) {
    ys <- ggplot2::scale_y_continuous(trans = "sqrt" )
    g <- g + ys
  } else if ( pTrans == "log10" ) {
    ys <- ggplot2::scale_y_continuous(trans = "log10", limits = c( s.min, s.max ),
                                      breaks = scales::trans_breaks("log10", function(x) 10^x),
                                      labels = scales::trans_format("log10", scales::math_format(10^.x) ) )
    g <- g + ys
  }
  if ( period ) {
    xs <- ggplot2::scale_x_continuous(trans = "log10", 
                             breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10", scales::math_format(10^.x) ) )
    g <- g + xs + ggplot2::labs( x = "Period" ) + ggplot2::annotation_logticks( sides = "tb" )
  }  
  
  return( g )
}