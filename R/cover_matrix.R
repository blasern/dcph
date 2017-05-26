#' Cover matrix
#' 
#' Calculate the n * m matrix C, where n is the number of data points and m is the 
#' number of nodes, where \deqn{C_{ij}} is 1 if the the i-th data point is in the j-th 
#' subset and 0 otherwise.
#' 
#' @param cover A cover
#' @export
cover_matrix <- function(cover){
  patch_matrix(cover@subsets)
}

patch_matrix <- function(patches){
  indices_matrix(lapply(patches, slot, "indices"))
}

indices_matrix <- function(indices){
  # get number of points
  npoints <- length(unique(unlist(indices)))
  # calculate cover matrix
  sapply(indices, function(x, npoints){
    vec <- rep(0, npoints)
    vec[x] <- 1
    vec
  }, npoints = npoints)
}

#' @rdname cover_matrix
#' @export
predict_cover_matrix <- function(cover){
  indices_matrix(lapply(cover@subsets, slot, "predicted"))
}