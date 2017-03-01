#' Distance between covers
#' 
#' Compare two covers using persistent homology or entropy
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
#'                       delta = 0.1,
#'                       relative_diameter = 0.1)
#' dc2 <- divisive_cover(distance_matrix = dist(data_matrix), 
#'                       delta = 0.05,
#'                       relative_diameter = 0.3)
#' cover1 <- subcover(dc1, 0.1, "snapshot")
#' cover2 <- subcover(dc2, 0.5, "snapshot")
#' 
#' cover_distance_ph(cover1, cover2)
#' cover_distance_vi(cover1, cover2)
#' 
#' # plot covers
#' \dontrun{
#' plot(cover1)
#' plot(cover2)
#' plot(c(cover1, cover2))
#' }
#' @name cover_distance
#' @export
cover_distance_ph <- function(cover1, cover2){
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

#' @rdname cover_distance
#' @importFrom infotheo condentropy
#' @export
cover_distance_vi <- function(cover1, cover2){
  # define cover matrices
  x_mat <- cover_matrix(cover1)
  y_mat <- cover_matrix(cover2)
  
  cond1 <- infotheo::condentropy(x_mat, y_mat)
  cond2 <- infotheo::condentropy(y_mat, x_mat)
  
  vi <- cond1 + cond2
  attr(vi, "H(1|2)") <- cond1
  attr(vi, "H(2|1)") <- cond2
  return(vi)
}
