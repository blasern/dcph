cover_skeleton <- function(dc){
  # find anchors
  anchors <- sapply(dc@subsets, slot, "anchor_points")
  unique_anchors <- unique(as.integer(anchors))
  anchor_data_index <- matrix(match(anchors, unique_anchors), ncol = 2, byrow = TRUE)
  # create anchor_data
  anchor_data <- dc@data[unique_anchors, ]
  # create divisive cover skeleton
  dc@data <- anchor_data
  dc@type <- "skeleton"
  # what to do with the parameters?
  dc@parameters <- list(distance_fct = dc@parameters$distance_fct)
  # fix subsets
  anchor_data_index
  dc@subsets <- mapply(function(subset, anchor_point){
    subset@anchor_points <- anchor_point
    subset@indices <- unique(anchor_point)
    subset@predicted <- integer(0)    
    subset
  }, subset = dc@subsets, anchor_point = as.list(as.data.frame(t(anchor_data_index))))
  # return skeleton
  return(dc)
}