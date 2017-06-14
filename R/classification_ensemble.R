#' Divisive classification ensemble
#' 
#' A divisive classification ensemble is an algorithm similar to random 
#' forests with divisive covers instead of decission trees. Currently uses
#' \code{\link{relative_gap_division}} with random gap as a division 
#' function. 
#' 
#' @param train Training data
#' @param group Class of training data
#' @param ntree Number of trees to grow
#' @param depth Maximal depth of trees
#' @param delta_range Range of delta parameters to consider
#' @param anchor_fct Which anchor function should be used (see \code{\link{anchor_fct}})
#' @param distance_fct Which distance function should be used (see \code{\link{distance_fct}})
#' @param stop_fct Which stop function should be used (see \code{\link{stop_fct}})
#' @param filter_fct Which filter function should be used (see \code{\link{filter_fct}})
#' @param dc_ensemble Output of \code{divisive_classification_ensemble}
#' @param test Test data
#' @param voting Use average probability voting or majority voting
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
#' # confusion matrix
#' table(observed, predictions)
#' @export
#' @importFrom stats predict rnorm runif
divisive_classification_ensemble <- function(train, group, ntree = 100, depth = 1000, delta_range = c(0, 0.2), 
                                             anchor_fct = anchor_extremal, 
                                             distance_fct = distance_cdist("euclidean"), 
                                             stop_fct = stop_relative_filter_max_nodes(relative_filter = 0, max_nodes = depth), 
                                             filter_fct = entropy_filter){
  res <- replicate(ntree, 
                   random_division_classification(train = train, group = group, depth = depth, 
                                                  delta_range = delta_range, anchor_fct = anchor_fct, 
                                                  distance_fct = distance_fct, stop_fct = stop_fct, 
                                                  filter_fct = filter_fct), 
                   simplify = FALSE)
  
  # class and attributes
  class(res) <- "divisive_ensemble"
  attr(res, "call") <- match.call()
  return(res)
}

random_division_classification <- function(train, group, depth, delta_range, anchor_fct, 
                                           distance_fct, stop_fct, filter_fct){
  # sample data, features and relative gap
  rows <- sample(nrow(train), replace = TRUE)
  weights <- stats::runif(ncol(train))
  relative_gap <- stats::runif(1, delta_range[1], delta_range[2])
  
  rotation <- random_rotation_matrix(ncol(train))

  rotated_train <-  as.matrix(train[rows, ]) %*% rotation
  weighted_train <- t(t(rotated_train) * weights)
  
  dc <- divisive_cover(data = weighted_train, 
                       group = group[rows], 
                       distance_fct = distance_fct, 
                       stop_fct = stop_fct, 
                       anchor_fct = anchor_fct, 
                       filter_fct = filter_fct, 
                       division_fct = relative_gap_division(relative_gap))
                       
  skelet <- cover_skeleton(dc)
  sc <- subcover(dc)
  cover_mat <- cover_matrix(sc)
  
  return(list(rotation = rotation, weights = weights, rows = rows, 
              relative_gap = relative_gap, skelet = skelet, 
              cover_mat = cover_mat, group = group[rows]))
}

#' @rdname divisive_classification_ensemble
#' @export
predict_divisive_classification_ensemble <- function(dc_ensemble, test, 
                                                     voting = c("average_probability", "majority")){
  individual_predictions <- lapply(dc_ensemble, predict_dc_probabilities, test = test)
  # voting
  voting <- match.arg(voting)
  # class prediction
  predicted_class <- switch(voting, 
                            majority = majority_voting(individual_predictions), 
                            average_probability = average_probability_voting(individual_predictions))
  # return
  predicted_class
}

average_probability_voting <- function(predictions){
  unique_groups <- levels(factor(unlist(predictions)))
  predictions <- lapply(predictions, function(x){
    attr(x, "probabilities")[, unique_groups]
  })
  average_probabilities <- Reduce("+", predictions) / length(predictions)
  predicted <- unique_groups[apply(average_probabilities, 1, which.max)]
  return(predicted)
}

majority_voting <- function(predictions){
  predictions <- sapply(predictions, identity)
  apply(predictions, 1, function(x){
    uniqv <- unique(x)
    tab <- tabulate(match(x, uniqv))
    uniqv[which.max(tab)]
  })
}

predict_dc_probabilities <- function(dc_list, test){
  # rotate and weigh  
  rotated_test <- as.matrix(test) %*% dc_list$rotation
  used_test <- t(t(rotated_test) * dc_list$weights)
  
  # predict
  pc <- predict(object = dc_list$skelet, newdata = used_test, 
                predict_fct = relative_gap_prediction(relative_gap = dc_list$relative_gap))
  group_from_predict_cover(pc, group = dc_list$group, cover_mat = dc_list$cover_mat)
}

#' Print divisive ensemble
#' 
#' @param x divisive ensemble
#' @param ... ignored
#' @export
print.divisive_ensemble <- function(x, ...){
  # call
  cat("\nCall:\n")
  print(attr(x, "call"))
  cat("\n")
  # number of covers
  cat("Number of covers: ", length(x), "\n", sep = "")
  # mean cover size
  sizes <- sapply(x, function(y) length(y$skelet@subsets))
  cat("Mean cover size: ", round(mean(sizes), 2), "; Range: ", min(sizes), "-", max(sizes), "\n", sep = "")
  # mean subset size
  subset_sizes <- unlist(lapply(x, function(y) rowSums(y$pred_mat)))
  cat("Mean subset size: ", round(mean(subset_sizes), 2), "; Range: ", min(subset_sizes), "-", max(subset_sizes), "\n", sep = "")
}

# random rotation matrix (source: Blaser & Fryzlewicz. Random Rotation Ensembles. 2016.)
random_On <- function(d){
  QR <- qr(matrix(stats::rnorm(d^2), ncol = d))
  M <- qr.Q(QR) %*% diag(sign(diag(qr.R(QR))))
  return(M)
}
random_rotation_matrix <- function(d){
  M <- random_On(d)
  if (det(M) < 0) M[, 1] <- -M[, 1]
  return(M)
}
