#' Diameter cover function
#' 
#' Find the divisive cover that corresponds to a certain diameter
#' 
#' @param cover A divisive cover
#' @param relative_diameter The diameter to extract
#' @examples
#' rcircle <- function(N, r, sd){
#'    radius <- rnorm(N, r, sd)
#'    angle <- runif(N, 0, 2 * pi)
#'    data.frame(x = radius * cos(angle), 
#'               y = radius * sin(angle))
#' }
#' data_matrix <- rcircle(200, 1, .1)
#' 
#' dc <- divisive_cover(distance_matrix = dist(data_matrix), 
#'                      relative_diameter = 0.5, relative_distance = 0.2)
#' rc <- diameter_cover(dc, relative_diameter = 0.7)
#' 
#' \dontrun{
#' plot(data_matrix)
#' plot(rc)
#' }
#' @export
diameter_cover <- function(cover, relative_diameter){
  # absolute diameter
  diameter <- relative_diameter * cover@subsets[[1]]@diameter
  birth <- sapply(cover@subsets, slot, "birth")
  death <- sapply(cover@subsets, slot, "death")
  # survivors
  surv <- diameter > death & diameter < birth
  # new cover
  cv <- cover
  cv@subsets <- cv@subsets[surv]
  cv
}