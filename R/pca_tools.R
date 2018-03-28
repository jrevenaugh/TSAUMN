#' Cumulative variance reduction
#'
#' Compute and (optionally) plot cumulative variance reduction of PCs returned from \code{\link{prcomp}}  
#' @param p "prcomp" list from \code{\link{prcomp}}
#' @param plot produce a ggplot2 plot
#' @param n starting with PC #1, the number of PCs that should be included in results, defaults to all
#' @return data.frame with PC index, the associated eigenvalue, variance (squared eigenvalue) and cumulative variance.
#' @export
#' @import "ggplot2"
#' @examples
#' p <- prcomp(USArrests, scale = TRUE)
#' pcaVar( p )
#' 
pcaVar <- function( p, plot = TRUE, n = length( p$sdev ) ) {
  
  # Compute variance and cumulative variance reduction
  eig <- p$sdev^2
  variance <- eig * 100 / sum( eig )
  cumvar <- cumsum( variance )
  var <- data.frame( PC = 1:n, Eigen = eig, Var = variance, CumVar = cumvar )
  
  # Plot if requested
  if ( plot ) {
    g <- ggplot2::ggplot( data = var[1:n,], aes( PC, Var ) ) + 
          ggplot2::geom_bar( stat = "identity", fill = "steelblue", color = "black" ) + 
          ggplot2::labs( x = "Eigenvector", y = "Percentage of Variance Explained", 
                         title = "Variance Distribution" )
    print( g )
  }
  
  return( var )
}

#' Principal component plot
#'
#' Plot (optionally) scaled principal components produced by \code{\link{prcomp}}  
#' @param p "prcomp" list from \code{\link{prcomp}}
#' @param which a vector with index vector of which PCs to show
#' @param scale if TRUE, scale PCs by the reciprocal of standard deviation.  If TRUE
#' PCs share a common scale.  If FALSE (default) scale depends on variance reduction. 
#' @export
#' @import "ggplot2"
#' @import "tidyr"
#' @examples
#' # Not the best example--abscissa axes are very crowded.
#' p <- prcomp(USArrests, scale = TRUE)
#' pcaPPC( p, scale = TRUE )
#' 
pcaPPC <- function( p, which = 1:length( p$sdev ), scale = FALSE  ) {
  # Build data.frame with (optionally) scaled PCs and index
  t <- rownames( p$x )
  if ( scale ) p$x[,which] <- scale( p$x[,which], scale = p$sdev[which], center = FALSE )
  x.df <- data.frame( Time = t, p$x[,which] )
  
  # Gather x.df into glyph-ready form
  x.gr <- tidyr::gather( x.df, key = variable, value = value, -Time )
  
  # Plot
  g <- ggplot2::ggplot( data = x.gr, aes( x = Time ) ) + 
    ggplot2::geom_line( aes( y = value, group = variable ) ) + 
    ggplot2::facet_wrap( "variable" ) +
    ggplot2::labs( x = "Sample Index ", y = "Value", title = "Principal Components" )
  
  return( g )
}