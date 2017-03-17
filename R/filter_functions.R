diameter_filter <- function(patches, data, distance_fct){
  sapply(patches, function(patch) diameter(patch, data, distance_fct))
}

diameter <- function(patch, data, distance_fct){
  distance_fct(data, patch@anchor_points[1], patch@anchor_points[2])
}

cardinality_filter <- function(patches, ...){
  sapply(patches, function(x) length(x@indices))
}
