#' Division functions
#' 
#' This function does the actual division in divisive cover. 
#' 
#' @details 
#' The function \code{relative_factor_division} returns a function that divides
#' a patch into subpatches. All points with distance from one anchor point 
#' larger than \code{relative_factor} times the distance from the other anchor 
#' point come in one subpatch and vice versa. 
#' 
#' The function \code{relative_gap_division} uses the \code{relative_factor_division}
#' with \code{relative_factor = (1-2 * relative_gap)/(1+2 * relative_gap)}
#' 
#' @param relative_gap delta-parameter for delta-filtered cover (0 < delta <= 1/2)
#' @param relative_factor corresponds to (1-2 * delta)/(1+2 * delta)
#' @name division_fct
#' @export
relative_factor_division <- function(relative_factor){
  function(data, patch, distance_fct){
  
    # base points
    a <- patch@anchor_points[1]
    b <- patch@anchor_points[2]
    
    # distances to base points
    dist_a <- distance_fct(data, a, patch@indices)
    dist_b <- distance_fct(data, b, patch@indices)
    
    # indices 
    A <- patch@indices[dist_b >= dist_a * relative_factor]
    B <- patch@indices[dist_a >= dist_b * relative_factor]
    
    return(list(patch(indices = A), patch(indices = B)))
  }
}

#' @rdname division_fct
#' @export
relative_gap_division <- function(relative_gap){
  relative_factor_division((1-2 * relative_gap)/(1+2 * relative_gap))
}