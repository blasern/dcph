#' Stop Functions
#' 
#' Possible stop criteria are based on the relative diameter 
#' or the maximum number of nodes. 
#' 
#' @param relative_diameter maximal diameter of the final cover
#' @param max_nodes maximum number of nodes the algorithm should generate
#' 
#' @name stop_fct
#' @export
stop_relative_diameter <- function(relative_diameter){
  stopifnot(relative_diameter >= 0)
  function(cover, next_division){
    if (length(cover@subsets[[next_division]]@indices) == 1) {
      return(TRUE)
    }
    else {
      return(cover@subsets[[next_division]]@diameter < cover@diameter * relative_diameter)
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
