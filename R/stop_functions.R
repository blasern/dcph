#' Stop Functions
#' 
#' Possible stop criteria are based on the relative filter value 
#' or the maximum number of nodes. 
#' 
#' @param relative_filter maximal relative filter value of final cover
#' @param max_nodes maximum number of nodes the algorithm should generate
#' 
#' @name stop_fct
#' @export
stop_relative_filter <- function(relative_filter){
  stopifnot(relative_filter >= 0)
  function(cover, next_division){
    if (length(cover@subsets[[next_division]]@indices) == 1) {
      return(TRUE)
    }
    else {
      return(cover@subsets[[next_division]]@filter_value <= cover@subsets[[1]]@filter_value * relative_filter)
    }
  }
}

#' @rdname stop_fct
#' @export
stop_max_nodes <- function(max_nodes){
  stopifnot(max_nodes >= 1)
  function(cover, next_division){
    if (length(cover@subsets[[next_division]]@indices) == 1) {
      return(TRUE)
    }
    else {
      return(length(cover@subsets) >= max_nodes)
    }
  }
}