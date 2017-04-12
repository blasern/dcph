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
#' @param method Method for randomly changing the data. The \code{"weights"} method (default) 
#' randomly assigns weights to the columns and the \code{"features"} method randomly selects 
#' columns (cf. \code{\link[randomForest]{randomForest}}). 
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
#' fit <- divisive_classification_ensemble(train = train, 
#'                                         group = group, 
#'                                         ntree = 2, 
#'                                         depth = 10, 
#'                                         delta_range = c(0, 0.1))
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
                                             anchor_fct = anchor_random_classify, 
                                             method = c("weights", "features")){
  method <- match.arg(method)
  if (method == "features"){
    res <- replicate(ntree, random_division_classification(train, group, depth, delta_range, anchor_fct), simplify = FALSE)
  } 
  if (method == "weights"){
    res <- replicate(ntree, weighted_random_division_classification(train, group, depth, delta_range, anchor_fct), simplify = FALSE)
  } 
  return(res)
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
                       division_fct = relative_gap_division(relative_gap = relative_gap, euclidean = TRUE))
  
  return(list(features = features, rows = rows, relative_gap = relative_gap, dc = dc))
}
weighted_random_division_classification <- function(train, group, depth, delta_range, anchor_fct){
  # sample data, features and relative gap
  rows <- sample(nrow(train), replace = TRUE)
  weights <- runif(ncol(train))
  relative_gap <- runif(1, delta_range[1], delta_range[2])
  
  weighted_train <- t(t(train[rows, ]) * weights)
  
  dc <- divisive_cover(data = weighted_train, 
                       group = group[rows], 
                       distance_fct = distance_cdist("euclidean"), 
                       stop_fct = stop_relative_filter_max_nodes(relative_filter = 0, max_nodes = depth), 
                       anchor_fct = anchor_fct, 
                       filter_fct = classification_filter, 
                       division_fct = relative_gap_division(relative_gap = relative_gap, euclidean = TRUE))
  
  skelet <- cover_skeleton(dc)
  sc <- subcover(dc)
  pred_mat <- cover_prediction_matrix(sc)
  
  return(list(weights = weights, rows = rows, relative_gap = relative_gap, skelet = skelet, pred_mat = pred_mat))
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
  weights <- dc_list$weights
  used_test <- test
  if (!is.null(weights)) {
    used_test <- t(t(used_test) * weights)
  }
  if (!is.null(features)) {
    used_test <- used_test[, features, drop = FALSE]
  }
  relative_gap <- dc_list$relative_gap
  skelet <- dc_list$skelet
  pc <- predict(object = skelet, newdata = used_test, 
                predict_fct = relative_gap_prediction(relative_gap = dc_list$relative_gap, euclidean = TRUE))
  group_from_predict_cover_pred(pc, dc_list$pred_mat)
}
