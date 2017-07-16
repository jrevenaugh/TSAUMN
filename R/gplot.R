gplot <- function( x, y = NULL, ... ) {
  if ( is.null( y ) ) {
    y <- x
    x <- seq_along( y )
  }
  df <- data_frame( x = x, y = y )
  ggplot( df, aes( x = x, y = y ) ) + geom_line(...)
}