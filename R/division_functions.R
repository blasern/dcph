#' Division functions
#' 
#' @param relative_gap delta-parameter for delta-filtered cover (0 < delta <= 1/2)
#' @param relative_factor corresponds to (1-2 * delta)/(1+2 * delta)
relative_factor_division <- function(relative_factor){
  function(data, patch, distance_fct){
  
    # base points
    a <- patch@anchor_points[1]
    b <- patch@anchor_points[2]
    
    # distances to base points
    dist_a <- distance_fct(data, a, patch@indices)
    dist_b <- distance_fct(data, b, patch@indices)
    
    # indices 
    A <- patch@indices[dist_b / dist_a >= relative_factor]
    B <- patch@indices[dist_a / dist_b >= relative_factor]
    
    return(list(patch(indices = A), patch(indices = B)))
  }
}

relative_gap_division <- function(relative_gap){
  relative_factor_division((1-2 * relative_gap)/(1+2 * relative_gap))
}