#' Compare covers using persistent homology
#' 
#' Compare two covers using persistent homology
#' 
#' @param cover1,cover2 covers to compare 
#' @examples 
#' # generate data
#' rcircle <- function(N, r, sd){
#'   radius <- rnorm(N, r, sd)
#'   angle <- runif(N, 0, 2 * pi)
#'   data.frame(x = radius * cos(angle), 
#'              y = radius * sin(angle))
#' }
#' data_matrix <- rcircle(200, 1, .1)
#' # calculate covers
#' dc1 <- divisive_cover(distance_matrix = dist(data_matrix), 
#'                       relative_distance = 0.2,
#'                       relative_diameter = 0.1)
#' dc2 <- divisive_cover(distance_matrix = dist(data_matrix), 
#'                       relative_distance = 0.1,
#'                       relative_diameter = 0.3)
#' cover1 <- subcover(dc1, 0.1, "snapshot")
#' cover2 <- subcover(dc2, 0.5, "snapshot")
#' 
#' compare_covers_ph(cover1, cover2)
#' 
#' # plot covers
#' \dontrun{
#' plot(cover1)
#' plot(cover2)
#' plot(c(cover1, cover2))
#' }
#' @export
compare_covers_ph <- function(cover1, cover2){
  # recalculate cover diameters
  adjust_diameter <- function(cover, diameter){
    adjust_subset_diam <- function(subset, id, diameter){
      subset@diameter <- diameter
      subset@birth <- diameter
      subset@death <- diameter
      subset
    }
    cover@subsets <- mapply(adjust_subset_diam, 
                            subset = cover@subsets, 
                            diameter = diameter)
    cover
  }

  contract <- cover(
    distance_matrix = cover1@distance_matrix,
    subsets = list(patch(indices = unique(unlist(lapply(cover1@subsets, slot, "indices"))), 
                         death = 2, 
                         birth = 2, 
                         diameter = 2)))

  pers21 <- persistent_homology(c(contract, adjust_diameter(cover1, 1), adjust_diameter(cover2, 0)))
  pers12 <- persistent_homology(c(contract, adjust_diameter(cover2, 1), adjust_diameter(cover1, 0)))
  list(`1->2` = pers12, `2->1` = pers21)
}