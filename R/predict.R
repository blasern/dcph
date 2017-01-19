#' Predict
#' 
#' Predict in which subsets new data come to lie when given a previous
#' divisive cover 
#' @param object a divisive cover
#' @param newdists a matrix of distances between the new points and the 
#' points used to train the divisive cover
#' @param ... further arguments passed fromother methods
predict.cover <- function(object, newdists, ...){
  if (nrow(newdists) != nrow(object@distance_matrix)){
    stop("newdists must contain distances to all points of the original data points")
  }
  relative_distance <- object@parameters[["relative_distance"]]
  
  newcover <- object
  newcover@type = "predict"
  newcover@subsets[[1]]@predicted <- 1:ncol(newdists)
  for (ix in 2:length(object@subsets)){
    newcover <- divide_pred(cover = newcover, 
                            index = ix, 
                            distance_matrix = newdists)
  }
  newcover
}

divide_pred <- function(cover, index, distance_matrix){
  parent_patch <- cover@subsets[[cover@subsets[[index]]@parent]]
  
  # base points
  basepoints <- parent_patch@basepoints
  bp <- basepoints %in% cover@subsets[[index]]@indices
  a <- basepoints[bp]
  b <- basepoints[!bp]
  dist_a <- distance_matrix[a, parent_patch@predicted]
  dist_b <- distance_matrix[b, parent_patch@predicted]
  
  # indices 
  A <- parent_patch@predicted[dist_b / dist_a >= 1 - cover@parameters[["relative_distance"]]]
  
  # update 
  cover@subsets[[index]]@predicted <- A
  
  # return cover
  return(cover)
}