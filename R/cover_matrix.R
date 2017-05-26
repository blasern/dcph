#' Cover matrix
#' 
#' Calculate the n * m matrix C, where n is the number of data points and m is the 
#' number of nodes, where \deqn{C_{ij}} is 1 if the the i-th data point is in the j-th 
#' subset and 0 otherwise.
#' 
#' @param cover A cover
#' @param n Data size
#' @export
cover_matrix <- function(cover, n = NULL){
  patch_matrix(cover@subsets, n = n)
}

patch_matrix <- function(patches, n = NULL){
  indices_matrix(lapply(patches, slot, "indices"), n = n)
}

indices_matrix <- function(indices, n = NULL){
  if (is.null(n)){
    # get number of points
    n <- length(unique(unlist(indices)))
  }
  # calculate cover matrix
  sapply(indices, function(x, npoints){
    vec <- rep(0, npoints)
    vec[x] <- 1
    vec
  }, npoints = n)
}

#' @rdname cover_matrix
#' @export
predict_cover_matrix <- function(cover, n = NULL){
  indices_matrix(lapply(cover@subsets, slot, "predicted"), n = n)
}