#' Divisive classification ensemble
#' 
#' A divisive classification ensemble is an algorithm similar to random forests with
#' divisive covers instead of decission trees. 
#' 
#' @param train Training data
#' @param group Class of training data
#' @param ntree Number of trees to grow
#' @param depth Maximal depth of trees
#' @param delta_range Range of delta parameters to consider
#' @param anchor_fct Which anchor function should be used (see \code{\link{anchor_fct}})
#' @param dc_ensemble Output of \code{divisive_classification_ensemble}
#' @param test Test data
#' 
#' @examples
#' # data
#' sampled_index <- sample(150)
#' train <- iris[sampled_index[1:100], 1:4]
#' group <- iris[sampled_index[1:100], 5]
#' test <- iris[sampled_index[101:150], 1:4]
#' observed <- iris[sampled_index[101:150], 5]
#' 
#' # fit divisive classification ensemble
#' fit <- divisive_classification_ensemble(train, group, ntree = 10, depth = 10, delta_range = c(0, 0))
#' # get predictions
#' predictions <- predict_divisive_classification_ensemble(fit, test)
#' 
#' # calculate misclassification error
#' table(observed, predictions)
#' tab <- table(observed, predictions)
#' misclassification_error <- 1 - sum(diag(tab)) / sum(tab)
#' misclassification_error
#' @export
divisive_classification_ensemble <- function(train, group, ntree = 100, depth = 1000, delta_range = c(0, 0.2), 
                                             anchor_fct = anchor_random_classify){
  replicate(ntree, random_division_classification(train, group, depth, delta_range, anchor_fct), simplify = FALSE)
}

random_division_classification <- function(train, group, depth, delta_range, anchor_fct){
  # sample data, features and relative gap
  rows <- sample(nrow(train), replace = TRUE)
  features <- sample(1:ncol(train), sqrt(ncol(train)))
  relative_gap <- runif(1, delta_range[1], delta_range[2])
  
  dc <- divisive_cover(data = train[rows, features, drop = FALSE], 
                       group = group[rows], 
                       distance_fct = distance_cdist("euclidean"), 
                       stop_fct = stop_relative_filter_max_nodes(relative_filter = 0, max_nodes = depth), 
                       anchor_fct = anchor_fct, 
                       filter_fct = classification_filter, 
                       division_fct = relative_gap_division(relative_gap = relative_gap))
  
  return(list(features = features, rows = rows, relative_gap = relative_gap, dc = dc))
}

#' @rdname divisive_classification_ensemble
#' @export
predict_divisive_classification_ensemble <- function(dc_ensemble, test){
  individual_predictions <- lapply(dc_ensemble, predict_dc_probabilities, test = test)
  # predictions
  predictions <- sapply(individual_predictions, identity)
  # voting
  majority_voting <- apply(predictions, 1, function(x){
    uniqv <- unique(x)
    tab <- tabulate(match(x, uniqv))
    uniqv[which.max(tab)]
  })
  
  # return
  majority_voting
}

predict_dc_probabilities <- function(dc_list, test){
  features <- dc_list$features
  relative_gap <- dc_list$relative_gap
  dc <- dc_list$dc
  pc <- predict(object = dc, newdata = test[, features], 
                predict_fct = relative_gap_prediction(relative_gap = dc_list$relative_gap))
  dc_pred <- group_from_predict_cover(pc)
}
