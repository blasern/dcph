#' Plot persistence
#' 
#' Plot persistence diagram or barcode
#' 
#' @param pers A persistence data.frame
#' @param mode Either "diag" for persistence diagram or "bars" for 
#' barcodes
#' @param relative Should relative diameters be used 
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line geom_segment theme theme_bw element_blank
#' @export
plot_persistence <- function(pers, mode = c("diag", "bars"), relative = TRUE){
  mode = match.arg(mode)
  pers$Dimension <- as.factor(pers$Dimension)
  if (relative){
    maxdiam <- attr(pers, "maxdiam")
    pers$Birth <- pers$Birth/maxdiam
    pers$Death <- pers$Death/maxdiam
    attr(pers, "maxdiam") <- 1
  }
  if (mode == "diag"){
    p <- plot_diag(pers)
  }
  if (mode == "bars"){
    p <- plot_bars(pers)
  }
  return(p)
}

plot_diag <- function(pers){
  maxdiam <- attr(pers, "maxdiam")
  p <- ggplot2::ggplot(pers, ggplot2::aes_string(x = "Birth", y = "Death")) + 
    ggplot2::geom_point(ggplot2::aes_string(color = "Dimension")) + 
    ggplot2::geom_line(data = data.frame(Birth = c(0, maxdiam), Death = c(0, maxdiam))) +
    ggplot2::theme_bw(base_size = 14) 
  return(p)
}

plot_bars <- function(pers){
  maxdiam <- attr(pers, "maxdiam")
  pers <- pers[order(pers$Dimension, pers$Birth, pers$Death), ]
  pers$y <- 1:nrow(pers)
  p <- ggplot2::ggplot(pers) +
    ggplot2::geom_segment(ggplot2::aes_string(x = "Birth", xend = "Death", color = "Dimension", y = "y", yend = "y")) +
    ggplot2::theme_bw(base_size = 14) + 
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous("", labels = NULL) +
    ggplot2::scale_x_continuous("Diameter", limits = c(0, maxdiam)) 
  p
  return(p)
}