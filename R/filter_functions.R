#' Filter functions
#' 
#' Filter after diameter or cardinality
#' 
#' @param patches Patch to apply filter
#' @param data Data used in distance function
#' @param distance_fct The used distance function
#' @name filter_fct
#' @export
diameter_filter <- function(patches, data, distance_fct){
  sapply(patches, function(patch) diameter(patch, data, distance_fct))
}

diameter <- function(patch, data, distance_fct){
  distance_fct(data, patch@anchor_points[1], patch@anchor_points[2])
}

#' @rdname
#' @export
cardinality_filter <- function(patches, ...){
  sapply(patches, function(x) length(x@indices))
}
