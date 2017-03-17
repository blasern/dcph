#' Distance functions for divisive cover
#' 
#' @param metric metric in `cdist`
#' @param data data for divisive cover
#' @param X,Y indices
#' @name distance_fct
#' @export
distance_cdist <- function(metric = "euclidean"){
  function(data, X, Y){
    cdist(data[X, , drop = FALSE], data[Y, , drop = FALSE], metric = metric)
  }
}
#' @rdname distance_fct
#' @export
distance_euclidean <- function(data, X, Y){
  cdist(data[X, , drop = FALSE], data[Y, , drop = FALSE], "euclidean")
}
#' @rdname distance_fct
#' @export
distance_matrix <- function(data, X, Y){
  data[X, Y, drop = FALSE]
}