#' Cover matrix
#' 
#' Calculate the n * m matrix C, where n is the number of data points and m is the 
#' number of nodes, where \deqn{C_{ij}} is 1 if the the i-th data point is in the j-th 
#' subset and 0 otherwise.
#' 
#' @param cover A cover
#' @export
cover_matrix <- function(cover){
  # get number of points
  npoints <- length(unique(unlist(lapply(cover@subsets, slot, "indices"))))
  # calculate cover matrix
  sapply(cover@subsets, function(x, npoints){
    vec <- rep(0, npoints)
    vec[x@indices] <- 1
    vec
  }, npoints = npoints)
}

#' @rdname cover_matrix
#' @export
predict_cover_matrix <- function(cover){
  # get number of points
  npoints <- length(unique(unlist(lapply(cover@predict, slot, "predicted"))))
  # calculate cover matrix
  sapply(cover@subsets, function(x, npoints){
    vec <- rep(0, npoints)
    vec[x@predicted] <- 1
    vec
  }, npoints = npoints)
}