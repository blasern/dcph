#' Radius cover function
#' 
#' Find the divisive cover that corresponds to a certain radius
#' 
#' @param cover A divisive cover
#' @param relative_radius The radius to extract
#' @examples
#' rcircle <- function(N, r, sd){
#'    radius <- rnorm(N, r, sd)
#'    angle <- runif(N, 0, 2 * pi)
#'    data.frame(x = radius * cos(angle), 
#'               y = radius * sin(angle))
#' }
#' data_matrix <- rcircle(200, 1, .1)
#' 
#' dc <- divisive_cover(data_matrix, distance_object = dist, 
#'                      relative_radius = 0.5, relative_distance = 0.2)
#' rc <- radius_cover(dc, relative_radius = 0.7)
#' 
#' \dontrun{
#' plot(data_matrix)
#' plot(rc)
#' }
#' @export
radius_cover <- function(cover, relative_radius){
  # absolute radius
  radius <- relative_radius * cover@subsets[[1]]@radius
#   birth <- sapply(cover@subsets, slot, "birth")
#   death <- sapply(cover@subsets, slot, "death")
  radii <- sapply(cover@subsets, slot, "radius")
  last <- max(which(radii > radius))
  # survivors
  surv <- cover@subsets[[last]]@survivors
  # new cover
  cv <- cover
  cv@subsets <- cv@subsets[surv]
  cv
}