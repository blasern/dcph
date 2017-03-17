#' Predict
#' 
#' Predict in which subsets new data come to lie when given a previous
#' divisive cover 
#' @param object a divisive cover
#' @param newdists a matrix of distances between the new points and the 
#' points used to train the divisive cover
#' @param newdata new data in case of fast divisive cover
#' @param ... further arguments passed from other methods
#' @examples 
#' # generate data
#' rcircle <- function(N, r, sd){
#'   radius <- rnorm(N, r, sd)
#'   angle <- runif(N, 0, 2 * pi)
#'   data.frame(x = radius * cos(angle), 
#'              y = radius * sin(angle))
#' }
#' data_matrix <- rcircle(200, 1, .1)
#' 
#' # run divisive cover
#' dc <- divisive_cover(distance_matrix = dist(data_matrix), 
#'                      delta = 0.1, 
#'                      stop_fct = stop_relative_diameter(relative_diameter = 0.2))
#'                      
#' # predict
#' pred <- predict(dc, newdists = as.matrix(dist(data_matrix))[, 1:10])
#' 
#' # plot
#' sc <- subcover(pred, 0.6, "snapshot")
#' col <- predict_coloring(sc)
#' \dontrun{
#' plot(sc, coloring = col)
#' }
#' @export
predict.cover <- function(object, newdists = NULL, newdata = NULL, ...){
  if (!(object@type %in% c("divisive", "fast_divisive"))){
    stop("can only predict from divisive cover")
  }
  if (!is.null(newdists) && nrow(newdists) != nrow(object@distance_matrix)){
    stop("newdists must contain distances to all points of the original data points")
  }
  if (!is.null(newdata) && ncol(newdata) != ncol(object@data)){
    stop("newdata must be of same dimension as original data")
  }
  
  newcover <- object
  newcover@type = "predict"
  
  if (object@type == "divisive") {
    newcover@subsets[[1]]@predicted <- 1:ncol(newdists)
    for (ix in 2:length(object@subsets)){
      newcover <- divide_pred(cover = newcover, 
                              index = ix, 
                              distance_matrix = newdists)
    }
  }
  if (object@type == "fast_divisive") {
    newcover@subsets[[1]]@predicted <- 1:nrow(newdata)
    distance_function <- object@parameters[["distance_function"]]
    for (ix in 2:length(object@subsets)){
      newcover <- fast_divide_pred(cover = newcover, 
                                   index = ix, 
                                   data = newdata, 
                                   distance_function = distance_function)
    }
  }

  return(newcover)
}

divide_pred <- function(cover, index, distance_matrix){
  parent_patch <- cover@subsets[[cover@subsets[[index]]@parent]]
  
  # find factor
  delta <- cover@parameters[["delta"]]
  relative_distance <- (1-2 * delta)/(1+2 * delta)
  
  # base points
  basepoints <- parent_patch@basepoints
  bp <- basepoints %in% cover@subsets[[index]]@indices
  a <- basepoints[bp]
  b <- basepoints[!bp]
  dist_a <- distance_matrix[a, parent_patch@predicted]
  dist_b <- distance_matrix[b, parent_patch@predicted]
  
  # indices 
  A <- parent_patch@predicted[dist_b / dist_a >= relative_distance]
  
  # update 
  cover@subsets[[index]]@predicted <- A
  
  # return cover
  return(cover)
}

fast_divide_pred <- function(cover, index, data, distance_function){
  parent_patch <- cover@subsets[[cover@subsets[[index]]@parent]]
  
  # find factor
  delta <- cover@parameters[["delta"]]
  relative_distance <- (1-2 * delta)/(1+2 * delta)
  
  # base points
  basepoints <- parent_patch@basepoints
  bp <- basepoints %in% cover@subsets[[index]]@indices
  a <- basepoints[bp]
  b <- basepoints[!bp]
  dist_a <- as.matrix(distance_function(cover@data[a, , drop = FALSE], data[parent_patch@predicted, , drop = FALSE]))
  dist_b <- as.matrix(distance_function(cover@data[b, , drop = FALSE], data[parent_patch@predicted, , drop = FALSE]))
  
  # indices 
  A <- parent_patch@predicted[dist_b / dist_a >= relative_distance]

  # update 
  cover@subsets[[index]]@predicted <- A
  
  # return cover
  return(cover)
}