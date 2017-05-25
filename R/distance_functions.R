#' Distance functions for divisive cover
#' 
#' Distance functions that can be used in \code{\link{divisive_cover}}. These functions
#' have arguments \code{data}, \code{X} and \code{Y}, where \code{data} is the data 
#' matrix and \code{X} and \code{Y} are indices. 
#' 
#' @details 
#' The function \code{distance_euclidean} returns the euclidean distance of 
#' \code{data[X, ]} and \code{data[Y, ]}. 
#' 
#' The function \code{distance_matrix} assumes that the data matrix is a distance 
#' matrix and returns \code{data[X, Y]}. 
#' 
#' The function \code{distance_cdist} takes a metric as an argument and returns a
#' distance function with arguments \code{data}, \code{X} and \code{Y}. See
#' \code{\link[rdist]{cdist}}. 
#' 
#' @param metric metric in `cdist`
#' @param data data for divisive cover
#' @param X,Y indices
#' @name distance_fct
#' @export
#' @importFrom rdist cdist
distance_cdist <- function(metric = "euclidean"){
  function(data, X, Y){
    rdist::cdist(data[X, , drop = FALSE], data[Y, , drop = FALSE], metric = metric)
  }
}
#' @rdname distance_fct
#' @export
distance_euclidean <- function(data, X, Y){
  rdist::cdist(data[X, , drop = FALSE], data[Y, , drop = FALSE], "euclidean")
}
#' @rdname distance_fct
#' @export
distance_matrix <- function(data, X, Y){
  data[X, Y, drop = FALSE]
}