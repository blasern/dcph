#' patch distances
#' 
#' Calculate the distances between patches
#' 
#' @param data Original data
#' @param X,Y Patches to calculate distance from
#' @param metric Used metric
#' 
#' @examples
#' # generate sample data
#' rcircle <- function(N, r, sd){
#'   radius <- rnorm(N, r, sd)
#'   angle <- runif(N, 0, 2 * pi)
#'   data.frame(x = radius * cos(angle), 
#'              y = radius * sin(angle))
#' }
#' data_matrix <- rcircle(200, 1, .1)
#' 
#' # calculate divisive cover
#' dc <- divisive_cover(data = data_matrix,
#'                      division_fct = relative_gap_division(0.05), 
#'                      stop_fct = stop_relative_filter(0.5))
#' 
#' # external patches
#' patches <- slot(subcover(dc), "subsets")
#' 
#' # distance
#' patch_cdist(data_matrix, patches, patches)
#' patch_hausdorff_dist(data_matrix, patches, patches)
#' patch_centroid_dist(data_matrix, patches, patches)
#' 
#' @importFrom rdist cdist
#' @export
#' @name patch_dist
patch_cdist <- function(data, X, Y, metric = "jaccard"){
  rdist::cdist(t(patch_matrix(X, n = nrow(data))), 
               t(patch_matrix(Y, n = nrow(data))), 
               metric = metric)
}

#' @rdname patch_dist
#' @export
patch_hausdorff_dist <- function(data, X, Y, metric = "euclidean"){
  matrix(mapply(pointwise_patch_hausdorff, 
                rep(X, times = length(Y)), 
                rep(Y, each = length(X)), 
                MoreArgs = list(data = data, metric = metric)), 
         nrow = length(X), 
         ncol = length(Y))
}
  
pointwise_patch_hausdorff <- function(data, P, Q, metric = "euclidean"){
  u <- data[P@indices, , drop = FALSE]
  v <- data[Q@indices, , drop = FALSE]
  hausdorff_dist(u, v, metric = metric)
}

hausdorff_dist <- function (P, Q, metric = "euclidean") 
{
  # inspired by pracma::hausdorff_dist
  D <- cdist(P, Q, metric = metric)
  dhd_PQ <- max(apply(D, 1, min))
  dhd_QP <- max(apply(D, 2, min))
  return(max(dhd_PQ, dhd_QP))
}
  
#' @rdname patch_dist
#' @export
patch_centroid_dist <- function(data, X, Y, metric = "euclidean"){
  cdist(centroid(data, X),
        centroid(data, Y),
        metric=metric)
}

centroid <- function(data, X){
  do.call(rbind, 
          lapply(X, function(x){
            colMeans(data[slot(x, "indices"), ])
          }))
}