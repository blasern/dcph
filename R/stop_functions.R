#' Stop Functions
#' 
#' Possible stop criteria are based on the relative filter value 
#' or the maximum number of nodes. 
#' 
#' @details 
#' These functions return stop functions that can be used in \code{\link{divisive_cover}}. The stop functions
#' return \code{TRUE} if the division should be stoped and \code{FALSE} otherwise.
#' 
#' The function \code{stop_relative_filter} returns a function that stops the division when the subset that should 
#' next be divided has a relative filter value of less than \code{relative_filter}. 
#' 
#' The function \code{stop_max_nodes} returns a function that stops the division when \code{max_nodes} nodes have 
#' been created. 
#' 
#' The function \code{stop_relative_filter_max_nodes} returns a function which stops the division if one of the 
#' above criteria is satisfied.  
#' 
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

#' @rdname stop_fct
#' @export
stop_relative_filter_max_nodes <- function(relative_filter, max_nodes){
  function(cover, next_division){
    if (length(cover@subsets[[next_division]]@indices) == 1) {
      return(TRUE)
    }
    else {
      return(length(cover@subsets) >= max_nodes || cover@subsets[[next_division]]@filter_value <= cover@subsets[[1]]@filter_value * relative_filter)
    }
  }
}