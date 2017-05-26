#' patch distances
#' 
#' Calculate the distances between patches
#' 
#' @param patches Patches to calculate distance from
#' @param X,Y Indices
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
#' patch_cdist(slot(subcover(dc), "subsets"))
#' 
#' @importFrom rdist cdist
#' @export
patch_cdist <- function(patches, X = 1:length(patches), Y = 1:length(patches), metric = "jaccard"){
  rdist::cdist(t(patch_matrix(patches[X])), 
               t(patch_matrix(patches[Y])), 
               metric = metric)
}