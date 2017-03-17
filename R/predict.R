#' Predict
#' 
#' Predict in which subsets new data come to lie when given a previous
#' divisive cover 
#' @param object a divisive cover
#' @param newdata new data 
#' @param predict_fct predict function (see below)
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
#' dc <- divisive_cover(data = data_matrix, 
#'                      division_fct = relative_gap_division(0.1), 
#'                      stop_fct = stop_relative_filter(relative_filter = 0.2))
#'                      
#' # predict
#' pred <- predict(dc, newdata = rcircle(10, 1, .1), 
#'                 predict_fct = relative_gap_prediction(0.1))
#' 
#' # plot
#' sc <- subcover(pred, 0.6, "snapshot")
#' col <- predict_coloring(sc)
#' \dontrun{
#' plot(sc, coloring = col)
#' }
#' @export
predict.cover <- function(object, newdata, predict_fct, ...){
  if (object@type != "divisive"){
    stop("can only predict from divisive cover")
  }
  # initialize new cover
  newcover <- object
  newcover@type = "predict"
  newcover@subsets[[1]]@predicted <- 1:nrow(newdata)
  
  # divide
  for (ix in 2:length(newcover@subsets)){
    newcover <- divide_pred(cover = newcover, 
                            index = ix, 
                            newdata = newdata, 
                            predict_fct = predict_fct)
  }

  return(newcover)
}

divide_pred <- function(cover, index, newdata, predict_fct){
  parent_patch <- cover@subsets[[cover@subsets[[index]]@parent]]
  
  # anchor points
  anchor_points <- parent_patch@anchor_points
  bp <- anchor_points %in% cover@subsets[[index]]@indices
  
  # predict
  cover@subsets[[index]]@predicted <- 
    predict_fct(data = cover@data, newdata = newdata, patch = parent_patch, 
                anchor = anchor_points[bp], 
                distance_fct = cover@parameters[["distance_fct"]])
  
  # return cover
  return(cover)
}