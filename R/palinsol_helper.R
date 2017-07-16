#' Compute June 1, 65 N insolation
#'
#' Compute June 1, 65 N insolation with \code{\link{Insol}} and one of three orbital calculations.
#' @param time vector of years at which to calculate insolation.  Negative for past times.
#' @param model one of la04, ber78 or ber90 to select orbital calculation, defaults to la04 (Laskar 2004)
#' @export
#' @import "palinsol"
#' @examples
#' s <- insolation( time = seq( -8e5, 0, 1e3 ), model = "ber78" )
#' 

insolation <- function( time = seq( -2e6, 0, 1e3 ), model = "la04" ) {
  model <- match.arg( model, c( "la04", "ber78", "ber90" ) )
  insol <- switch( model,
                   la04 =  sapply( time, function(tt) palinsol::Insol( orbit = palinsol::la04( tt ) ) ),
                   ber78 = sapply( time, function(tt) palinsol::Insol( orbit = palinsol::ber78( tt ) ) ),
                   ber90 = sapply( time, function(tt) palinsol::Insol( orbit = palinsol::ber90( tt ) ) ) )
  
  return( data_frame( time = time, insolation = insol ) )
}

#' Compute June 1, 65 N eccentricity
#'
#' Compute June 1, 65 N eccentricity with one of three orbital calculations.
#' @param time vector of years at which to calculate eccentricity.  Negative for past times.
#' @param model one of la04, ber78 or ber90 to select orbital calculation, defaults to la04 (Laskar 2004)
#' @export
#' @import "palinsol"
#' @examples
#' ecc <- eccentricity( time = seq( -8e5, 0, 1e3 ), model = "ber78" )
#' 
eccentricity <- function( time = seq( -2e6, 0, 1e3 ), model = "la04" ) {
  model <- match.arg( model, c( "la04", "ber78", "ber90" ) )
  eccen <- switch( model,
                   la04 =  sapply( time, function(tt) palinsol::la04( tt )['ecc'] ),
                   ber78 = sapply( time, function(tt) palinsol::ber78( tt )['ecc'] ),
                   ber90 = sapply( time, function(tt) palinsol::ber90( tt )['ecc'] ) )
  return( data_frame( time = time, eccentricity = eccen ) )
}
  
#' Compute June 1, 65 N obliquity
#'
#' Compute June 1, 65 N obliquity with one of three orbital calculations.
#' @param time vector of years at which to calculate obliquity.  Negative for past times.
#' @param model one of la04, ber78 or ber90 to select orbital calculation, defaults to la04 (Laskar 2004)
#' @export
#' @import "palinsol"
#' @examples
#' o <- obliquity( time = seq( -8e5, 0, 1e3 ), model = "ber78" )
#' 
obliquity <- function( time = seq( -2e6, 0, 1e3 ), model = "la04" ) {
  model <- match.arg( model, c( "la04", "ber78", "ber90" ) )
  obliq <- switch( model,
                   la04 =  sapply( time, function(tt) palinsol::la04( tt )['eps'] ),
                   ber78 = sapply( time, function(tt) palinsol::ber78( tt )['eps'] ),
                   ber90 = sapply( time, function(tt) palinsol::ber90( tt )['eps'] ) )
  return( data_frame( time = time, obliquity = obliq ) )
}

#' esin
#'
#' Compute esin term used in precession index computation.  There should be no need to call this separate from \code{\link{precession}}
#' @param orb orbital calculation output from one of \code{\link{la04}}, \code{\link{ber78}} or \code{\link{ber90}}
#' @keywords internal
#' @import "palinsol"
#' @export
esin <- function( orb ) orb['eps'] * sin( orb['varpi'] )

#' Compute June 1, 65 N precession index
#'
#' Compute June 1, 65 N precession index with one of three orbital calculations.
#' @param time vector of years at which to calculate precession.  Negative for past times.
#' @param model one of la04, ber78 or ber90 to select orbital calculation, defaults to la04 (Laskar 2004)
#' @export
#' @import "palinsol"
#' @examples
#' p <- precession( time = seq( -8e5, 0, 1e3 ), model = "ber78" )
#' 
precession <- function( time = seq( -2e6, 0, 1e3 ), model = "la04" ) {
  model <- match.arg( model, c( "la04", "ber78", "ber90" ) )
  p <- switch( model,
                   la04 =  sapply( time, function(tt) esin( palinsol::la04( tt ) ) ),
                   ber78 = sapply( time, function(tt) esin( palinsol::ber78( tt ) ) ),
                   ber90 = sapply( time, function(tt) esin( palinsol::ber90( tt ) ) ) )
  return( data_frame( time = time, precession = p ) )
}
