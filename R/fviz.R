#' Plot amplitude and phase of linear filter
#'
#' Plot amplitude and phase of a linear filter  
#' @param h filter coefficients; see details below
#' @param type one of "centered", "causal", "symmetric" or "arma."
#' The first three types assume \code{h} is a vector of filter coefficients.
#' If type is "centered" these are centered around zero.  If the filter length is even,
#' the filter extends further to positive time.  Type "causal" assumes coefficients start
#' with time 0.  Type "symmetric" assumes \code{h} provides the time >= 0 coefficients.  
#' These are mirrored for negative times.
#' Type "arma" is for filters created in the \code{\link{signal}} package.
#' @param phase include filter phase spectrum facet if true
#' @param hide if true, function returns a list with with filter z transform (see below)
#' and ggplot object
#' @param twoway if true and filter type is ARMA, response is computed for two-way filter
#' @return List with frequency "f" and z representation of filter "H"
#' @import tidyr
#' @export
#' @examples
#' f <- fdfilter( 5, 1 )
#' H <- fviz( f )
#' 
fviz <- function( h, 
                  type = c( "centered", "causal", "symmetric", "arma" ), 
                  phase = TRUE,
                  hide = FALSE,
                  twoway = FALSE ) {
  nf <- 256
  type <- match.arg( type )
  f <- seq( 0, 0.5, length.out = nf )
  if ( type != "arma" ) {
    if ( type == "symmetric" ) h <- c( rev( h ), h[-1] )
    n <- length( h )
    if ( type == "causal" ) {
      t0 <- 0
      tn <- n
    } else {
      t0 <- -( n - 1 ) / 2
      tn <- t0 + n - 1
    }
    H <- rep( 0 + 0i, nf )
    z <- exp( 2 * pi * 1i * f )
    for ( k in 1:n ) {
      H <- H + z^( k - 1 + t0 ) * h[k]
    }
  } else {
    h <- as.Arma( h )
    n <- length( h$b )
    d <- length( h$a )
    l <- max( n, d )
    Z <- rep( 1, nf )
    z <- exp( pi * 1i * f * 2 )
    for ( i in 2:l) Z <- cbind( Z, z^(i-1) )
    Hn <- Z[,1:n] %*% as.matrix( h$b )
    Hd <- Z[,1:d] %*% as.matrix( h$a )
    H <- Hn / Hd
    if ( twoway ) H <- Conj( H ) * H
  }
  Amplitude <- Mod( H )
  Phase <- phunwrap( Arg( H ) )
  H.df <- data.frame( Frequency = f, Amplitude, Phase )
  H.gr <- tidyr::gather( H.df, key = Key, value = Value, -Frequency )
  if ( !phase ) H.gr <- dplyr::filter( H.gr, Key == "Amplitude" )
  g <- ggplot( H.gr, aes( x = Frequency, y = Value ) ) +
    geom_line( aes( group = Key ) ) +
    facet_grid( Key ~ ., scales = "free_y" )
  if ( !hide ) {
    print( g )
    return( list( f = f, H = H ) )
  } else {
    return( list( f = f, H = H, plot = g ) )
  }
}