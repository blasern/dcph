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
#' # distance of external nodes
#' patches <- slot(subcover(dc), "subsets")
#' patch_cdist(data_matrix, patches, patches)
#' 
#' @importFrom rdist cdist
#' @export
patch_cdist <- function(data, X, Y, metric = "jaccard"){
  rdist::cdist(t(patch_matrix(X, n = nrow(data))), 
               t(patch_matrix(Y, n = nrow(data))), 
               metric = metric)
}

