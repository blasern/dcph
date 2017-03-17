#' Persistent homology
#' 
#' Calculate the persistent homology from the divisive cover complex
#' 
#' @details The calculation of persistent homology is done with PHAT. 
#' @references Bauer U., Kerber M., Reininghaus J., Wagner H. (2014) PHAT - Persistent Homology Algorithms Toolbox.
#' @references https://bitbucket.org/phat-code/phat
#' @param cover The divisive cover
#' @param max_dim The maximal dimension to calculate
#' @param representation Boundary matrix representation from PHAT
#' @param reduction Reduction method from PHAT
#' @examples 
#' rcircle <- function(N, r, sd){
#'    radius <- rnorm(N, r, sd)
#'    angle <- runif(N, 0, 2 * pi)
#'    data.frame(x = radius * cos(angle), 
#'               y = radius * sin(angle))
#' }
#' data_matrix <- rcircle(200, 1, .1)
#' # get divisive cover
#' dc <- divisive_cover(data = data_matrix, 
#'                      division_fct = relative_gap_division(0.1), 
#'                      stop_fct = stop_relative_filter(relative_filter = 0.5))
#' # calculate persistent homology
#' pers <- persistent_homology(dc)
#' \dontrun{
#' # plot persistence diagram and barcode
#' plot_persistence(pers, mode = "diag")
#' plot_persistence(pers, mode = "bars")
#' }
#' @export
persistent_homology <- function(cover, max_dim = 1, 
                                representation = c("bit_tree_pivot_column", "vector_vector", "sparse_pivot_column", "heap_pivot_column", "full_pivot_column", "vector_heap", "vector_set", "vector_list"),
                                reduction = c("twist", "standard", "chunk", "row", "spectral_sequence")){
  # get representation/reduction right
  representation <- match.arg(representation)
  reduction <- match.arg(reduction)
  # remove duplicates
  cover <- simplify_cover(cover, method = "duplicates")
  # calculate persistent homology
  pers <- persistence_from_cover(cover, max_dim = max_dim, representation = representation, reduction = reduction)
  # add column names 
  pers <- as.data.frame(pers)
  colnames(pers) <- c("Dimension", "Birth", "Death")
  # minimum diameter
  max_diameter <- cover@subsets[[1]]@filter_value
  # min_diameter <- cover@parameters$relative_diameter * max_diameter
  min_diameter <- min(sapply(cover@subsets[cover@internal_nodes], slot, "filter_value"))
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
  attr(pers, "max_diameter") <- max_diameter
  attr(pers, "min_diameter") <- min_diameter
  attr(pers, "max_dimension") <- max_dim
  return(pers)
}