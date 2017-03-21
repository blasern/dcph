#' Filter functions
#' 
#' Filter after diameter or cardinality
#' 
#' @param patches Patch to apply filter
#' @param data Data used in distance function
#' @param distance_fct The used distance function
#' @param group the group for classification
#' @param ... ignored arguments
#' 
#' @name filter_fct
#' @export
diameter_filter <- function(patches, data, distance_fct, ...){
  sapply(patches, function(patch) diameter(patch, data, distance_fct))
}

diameter <- function(patch, data, distance_fct){
  distance_fct(data, patch@anchor_points[1], patch@anchor_points[2])
}

#' @rdname filter_fct
#' @export
cardinality_filter <- function(patches, ...){
  sapply(patches, function(x) length(x@indices))
}

#' @rdname filter_fct
#' @export
classification_filter <- function(patches, group, ...){
  sapply(patches, function(x) missclassification_rate(group[x@indices]))
}

missclassification_rate <- function(group){
  sum(head(sort(table(group)/length(group)), -1))
}