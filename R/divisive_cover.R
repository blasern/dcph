#' divisive cover 
#'
#' This function divides the data into a cover using the distance matrix. 
#'
#' @param cover an optional divisive cover used to update
#' @param data a data matrix
#' @param distance_fct the distance function used (see details below) 
#' @param stop_fct function that determines when to stop dividing (see details below) 
#' @param anchor_fct function that determined how anchor points are found (see details below)
#' @param filter_fct function that filters the cover (see details below)
#' @param division_fct function that divides the patches (see details below)
#' @param group group for classification
#' @details 
#' TODO... write detailed documentation
#' The function \code{stop_fct} is a function of ... that returns \code{TRUE} if the division should
#' be stoped and \code{FALSE} otherwise. Examples are \code{\link{stop_relative_filter}} and 
#' \code{\link{stop_max_nodes}}.
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
#' # calculate and update divisive cover
#' dc1 <- divisive_cover(data = data_matrix,
#'                       division_fct = relative_gap_division(0.05), 
#'                       stop_fct = stop_relative_filter(0.5))
#' dc2 <- divisive_cover(cover = dc1, 
#'                       data = data_matrix,
#'                       division_fct = relative_gap_division(0.05),
#'                       stop_fct = stop_max_nodes(20))
#' 
#' # get one snapshot for plotting
#' ddc1 <- subcover(dc1, relative_filter = 0, method = "snapshot")
#' ddc2 <- subcover(dc2, relative_filter = 0, method = "snapshot")
#' 
#' \dontrun{
#' plot(data_matrix)
#' plot(ddc1)
#' plot(ddc2)
#' }
#' @export
divisive_cover <- function(cover = NULL,
                           data, 
                           distance_fct = distance_cdist("euclidean"), 
                           stop_fct = stop_relative_filter(relative_filter = .5), 
                           anchor_fct = anchor_extremal, 
                           filter_fct = diameter_filter,
                           division_fct = relative_gap_division(relative_gap = 0.1), 
                           group = NULL
                           ){
  # generate initial cover
  if (is.null(cover)){
    initial_patch <- patch(id = 1L, 
                           indices = 1:nrow(data))
    initial_patch@anchor_points <- anchor_fct(points = 1:nrow(data), data = data, distance_fct = distance_fct, group = group)
    initial_patch@filter_value <- filter_fct(patches = list(initial_patch), data = data, distance_fct = distance_fct, group = group)
    initial_patch@parent_filter <- initial_patch@filter_value
    cover <- cover(data = data, 
                   subsets = list(initial_patch), 
                   external_nodes = 1L,
                   data_filter_value = initial_patch@filter_value,
                   parameters = list(distance_fct = distance_fct, 
                                     stop_fct = stop_fct, 
                                     anchor_fct = anchor_fct, 
                                     filter_fct = filter_fct,
                                     division_fct = division_fct, 
                                     stop_fct = stop_fct, 
                                     group = group), 
                   type = "divisive")
    filter_values <- initial_patch@filter_value
    next_division <- 1L
  }
  else {
    filter_values <- filter_fct(patches = cover@subsets, data = data, distance_fct = distance_fct, group = group)
    filter_values[cover@internal_nodes] <- -Inf
    next_division <- which.max(filter_values)
  }
  
  # divide into pieces
  while (!stop_fct(cover, next_division) && distance_fct(data, 
                                                         cover@subsets[[next_division]]@anchor_points[1], 
                                                         cover@subsets[[next_division]]@anchor_points[2]) != 0){
    # divide original patch into new ones
    new_patches <- division_fct(data = data, 
                                patch = cover@subsets[[next_division]], 
                                distance_fct = distance_fct)
    if (any(is.na(c(new_patches[[1]]@indices, new_patches[[2]]@indices)))) stop("Sometinng went wrong... check input")
    # update information in new patches
    new_patches <- lapply(new_patches, function(patch) {
      patch@anchor_points <- anchor_fct(points = patch@indices, data = data, distance_fct = distance_fct, group = group)
      patch@filter_value <- filter_fct(patches = list(patch), data = data, distance_fct = distance_fct, group = group)
      patch@parent <- next_division
      patch@parent_filter <- cover@subsets[[next_division]]@filter_value
      patch
    })
    new_patches[[1]]@id <- length(cover@subsets) + 1L
    new_patches[[2]]@id <- length(cover@subsets) + 2L
    # update cover
    cover@internal_nodes <- c(cover@internal_nodes, next_division)
    cover@external_nodes <- c(setdiff(cover@external_nodes, next_division), length(cover@subsets) + 1:2)
    cover@subsets <- c(cover@subsets, new_patches)
    # update filter values
    new_filters <- sapply(tail(cover@subsets, 2), slot, "filter_value")
    filter_values[next_division] <- -Inf
    filter_values <- c(filter_values, new_filters)
    next_division <- which.max(filter_values)
  }
  cover
}