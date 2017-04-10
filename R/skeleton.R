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
  dc@parameters <- NULL
  # fix subsets
  dc@subsets <- apply(anchor_data_index, 1, 
                      function(x) patch(anchor_points = x, indices = unique(x)))
  # return skeleton
  return(dc)
}