#' Anchor functions for divisive cover
#' 
#' @param points indices of points in patch 
#' @param data data
#' @param distance_fct distance function
#' @name anchor_fct
#' @export
anchor_extremal <- function(points, data, distance_fct){
  points[arrayInd(which.max(distance_fct(data, points, points)), rep(length(points), 2))]
}
#' @rdname anchor_fct
#' @export
anchor_heuristic_extremal <- function(points, data, distance_fct){
  p1 <- which.max(distance_fct(data, points[1], points))
  p2 <- which.max(distance_fct(data, p1, points))
  points[c(p1, p2)]
}
