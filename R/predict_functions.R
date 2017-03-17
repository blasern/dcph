#' Prediction functions
#' 
#' @param relative_gap delta-parameter for delta-filtered cover (0 < delta <= 1/2)
#' @param relative_factor corresponds to (1-2 * delta)/(1+2 * delta)
#' @name prediction_fct
#' @export
relative_factor_prediction <- function(relative_factor){
  function(data, newdata, patch, anchor, distance_fct){
    
    # base points
    a <- anchor
    b <- patch@anchor_points[patch@anchor_points != anchor]
    
    # distances to base points
    dist_a <- distance_fct(rbind(data[a, ], newdata), 1, 1 + patch@predicted)
    dist_b <- distance_fct(rbind(data[b, ], newdata), 1, 1 + patch@predicted)
    
    # indices 
    A <- patch@predicted[dist_b / dist_a >= relative_factor]
    
    return(A)
  }
}

#' @rdname prediction_fct
#' @export
relative_gap_prediction <- function(relative_gap){
  relative_factor_prediction((1-2 * relative_gap)/(1+2 * relative_gap))
}