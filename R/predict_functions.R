#' Prediction functions
#' 
#' The prediction function should usually take the same arguments as
#' the division function. 
#' 
#' @param relative_factor corresponds to (1-2 * delta)/(1+2 * delta)
#' @param relative_gap delta-parameter for delta-filtered cover (0 < delta <= 1/2)
#' @param euclidean logical, if TRUE then (1-delta)/(1+delta) is used instead of 
#' (1-2 * delta)/(1+2 * delta)
#' 
#' @export
#' @rdname prediction_fct
relative_factor_prediction <- function(relative_factor){
  function(data, newdata, patch, anchor, distance_fct){
    
    # base points
    a <- anchor
    b <- patch@anchor_points[patch@anchor_points != anchor]
    
    stopifnot(length(a) == 1, length(b) == 1)
    
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
relative_gap_prediction <- function(relative_gap, euclidean = FALSE){
  if (euclidean){
    return(relative_factor_prediction((1-relative_gap)/(1+relative_gap)))
  }
  else {
    return(relative_factor_prediction((1-2 * relative_gap)/(1+2 * relative_gap)))
  }
}