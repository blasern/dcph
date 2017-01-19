#' Predict
#' 
#' Predict in which subsets new data come to lie when given a previous
#' divisive cover 
#' @param object a divisive cover
#' @param newdists a matrix of distances between the new points and the 
#' points used to train the divisive cover
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
#'                      relative_diameter = 0.2, relative_distance = 0.3)
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