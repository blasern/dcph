#' Predict
#' 
#' Predict in which subsets new data come to lie when given a previous
#' divisive cover 
#' @param object a divisive cover
#' @param newdata new data 
#' @param predict_fct predict function (see below)
#' @param ... further arguments passed from other methods
#' @export
#' @name predict_cover
predict.cover <- function(object, newdata, predict_fct, ...){
  if (!(object@type %in% c("divisive", "skeleton"))){
    stop("can only predict from divisive or skeleton cover")
  }
  # initialize new cover
  newcover <- object
  newcover@type = "predict"
  newcover@subsets[[1]]@predicted <- 1:nrow(newdata)
  
  # divide
  for (ix in 1:length(newcover@subsets)){
    newcover <- divide_pred(cover = newcover, 
                            index = ix, 
                            newdata = newdata, 
                            predict_fct = predict_fct)
  }

  return(newcover)
}

divide_pred <- function(cover, index, newdata, predict_fct){
  # current patch
  parent_patch <- cover@subsets[[index]]
  
  # children 
  children <- cover@subsets[[index]]@children
  
  # anchor points
  anchor_points <- parent_patch@anchor_points
  
  if (length(children == 2)){
    # predict
    cover@subsets[[children[1]]]@predicted <- 
      predict_fct(data = cover@data, newdata = newdata, patch = parent_patch, 
                  anchor = anchor_points[1], 
                  distance_fct = cover@parameters[["distance_fct"]])
    
    cover@subsets[[children[2]]]@predicted <- 
      predict_fct(data = cover@data, newdata = newdata, patch = parent_patch, 
                  anchor = anchor_points[2], 
                  distance_fct = cover@parameters[["distance_fct"]])
  }

  # return cover
  return(cover)
}

#' Extract group from predict cover
#' 
#' @param pc prediction cover
#' @param group specify the group of the training set (not necessary if used in divisive_cover)
#' @export
#' @importFrom dplyr group_by_ select_ do right_join ungroup
#' @examples
#' iris <- datasets::iris
#' 
#' ttsplit <- sample(1:2, nrow(iris), replace = TRUE, prob = c(2, 1))
#' train <- iris[ttsplit == 1, ]
#' test <- iris[ttsplit == 2, ]
#' 
#' dc <- divisive_cover(data = train[, 1:4], 
#'                      group = train[, 5], 
#'                      distance_fct = distance_cdist("euclidean"), 
#'                      stop_fct = stop_relative_filter(relative_filter = 0), 
#'                      anchor_fct = anchor_extremal, 
#'                      filter_fct = entropy_filter, 
#'                      division_fct = relative_gap_division(relative_gap = 0.01)
#' )
#' 
#' pc <- predict(object = dc, newdata = test[, 1:4], 
#'               predict_fct = relative_gap_prediction(relative_gap = 0.01))
#' # plot
#' sc <- subcover(pc, method = "snapshot", relative_filter = 0.6)
#' 
#' predicted <- group_from_predict_cover(pc)
#' observed <- test[, 5]
#' 
#' table(observed, predicted)
#' @rdname predict_cover
group_from_predict_cover <- function(pc, group = NULL, cover_mat = NULL){
  # use external nodes
  sc <- subcover(pc, method = "external")
  # check group
  if (is.null(group)){
    group <- pc@parameters$group
  }
  # cover matrix
  cpm <- cover_predicted_matrix(sc)
  if (is.null(cover_mat)){
    cover_mat <- cover_matrix(sc)
  }
  # cover distance between predicted points and original points
  distance <- rdist::cdist(cpm, cover_mat, metric = "minkowski", p = 1)
  # get unique groups
  unique_groups <- unique(group)
  # calculate group probabilities
  probabilities <- t(apply(distance, 1, function(x, group, all_groups){
    indices <- which(x == min(x))
    possible_groups <- group[indices]
    table(factor(possible_groups, all_groups))/length(possible_groups)
  }, group = group, all_groups = unique_groups))
  # get group
  predicted <- unique_groups[apply(probabilities, 1, which.max)]
  attr(predicted, "probabilities") <- probabilities
  return(predicted)
}

cover_predicted_matrix <- function (cover) 
{
  npoints <- length(unique(unlist(lapply(cover@subsets, slot, "predicted"))))
  sapply(cover@subsets, function(x, npoints) {
    vec <- rep(0, npoints)
    vec[x@predicted] <- 1
    vec
  }, npoints = npoints)
}
