#' Persistent homology
#' 
#' Calculate the persistent homology from the divisive cover complex
#' 
#' @param cover The divisive cover
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
#' pers <- persistent_homology(dc)
#' \dontrun{
#' plot_persistence(pers, mode = "diag")
#' plot_persistence(pers, mode = "bars")
#' }
#' @export
persistent_homology <- function(cover){
  # remove duplicates
  cover <- simplify_cover(cover, method = "duplicates")
  # calculate persistent homology
  pers <- persistence_from_cover(cover)
  # add column names 
  pers <- as.data.frame(pers)
  colnames(pers) <- c("Dimension", "Birth", "Death")
  # minimum diameter
  max_diameter <- cover@subsets[[1]]@diameter
  min_diameter <- cover@parameters$relative_diameter * max_diameter
  # add survivor
  pers <- rbind(pers, data.frame(Dimension = 0, 
                                 Birth = min_diameter, 
                                 Death = max_diameter))
  # make values compatible with minimum diameter
  pers$Birth[pers$Birth < min_diameter] <- min_diameter
  pers$Death[pers$Death < min_diameter] <- min_diameter
  # remove diagonal
  pers <- pers[pers$Birth != pers$Death, ]
  rownames(pers) <- NULL
  attr(pers, "maxdiam") <- max_diameter
  attr(pers, "mindiam") <- min_diameter
  return(pers)
}