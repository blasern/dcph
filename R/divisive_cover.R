#' divisive cover 
#'
#' This function divides the data into a cover using the distance matrix. 
#'
#' @param distance_matrix an n x n matrix of pairwise dissimilarities
#' @param delta delta-parameter for delta-filtered cover (0 < delta <= 1/2)
#' @param stop_fct function that determines when to stop dividing (see details below) 
#' @param cover a divisive cover used to update
#' @details 
#' The function \code{stop_fct} is a function of ... that returns \code{TRUE} if the division should
#' be stoped and \code{FALSE} otherwise. Examples are \code{\link{stop_relative_diameter}} and 
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
#' dc1 <- divisive_cover(distance_matrix = dist(data_matrix),
#'                       delta = 0.05, 
#'                       stop_fct = stop_relative_diameter(0.5))
#' dc2 <- divisive_cover(cover = dc1, 
#'                       distance_matrix = dist(data_matrix),
#'                       delta = 0.05, 
#'                       stop_fct = stop_max_nodes(20))
#' 
#' # get one snapshot for plotting
#' ddc1 <- subcover(dc1, relative_diameter = 0, method = "snapshot")
#' ddc2 <- subcover(dc2, relative_diameter = 0, method = "snapshot")
#' 
#' \dontrun{
#' plot(data_matrix)
#' plot(ddc1)
#' plot(ddc2)
#' }
#' @export
divisive_cover <- function(cover = NULL,
                           distance_matrix, 
                           delta, 
                           stop_fct = stop_relative_diameter(.5)){
  # check input
  stopifnot(0 <= delta && delta < 0.5)

  # find factor
  relative_distance <- (1-2 * delta)/(1+2 * delta)
  
  # data size
  distance_matrix <- as.matrix(distance_matrix)
  N <- nrow(distance_matrix)
  
  # generate initial cover
  if (is.null(cover)){
    basepoints <- as.integer(arrayInd(which.max(distance_matrix), dim(distance_matrix)))
    diameter <- distance_matrix[basepoints[1], basepoints[2]]
    cover <- cover(distance_matrix = distance_matrix, 
                   subsets = list(patch(1:N, basepoints = basepoints, id = 1L, diameter = diameter, birth = diameter)), 
                   parameters = list(delta = delta, stop_fct = stop_fct), 
                   type = "divisive")
    # minimal radius
    diams <- diameter
    next_division <- 1
  }
  else {
    diams <- sapply(cover@subsets, slot, "diameter")
    diams[!sapply(cover@subsets, slot, "id") %in% sapply(subcover(cover, 0, "snapshot")@subsets, slot, "id")] <- -Inf
    next_division <- which.max(diams)
  }
  
  # divide into pieces
  while (!stop_fct(cover, next_division)){
    cover <- divide(cover = cover, 
                    index = next_division, 
                    distance_matrix = distance_matrix,                     
                    relative_distance = relative_distance)
    new_diams <- sapply(tail(cover@subsets, 2), slot, "diameter")
    diams[next_division] <- -Inf
    diams <- c(diams, new_diams)
    next_division <- which.max(diams)
  }
  cover
}

divide <- function(cover, index, distance_matrix, relative_distance){
  divide_patch <- cover@subsets[[index]]
  
  # base points
  a <- divide_patch@basepoints[1]
  b <- divide_patch@basepoints[2]
  dist_a <- distance_matrix[a, divide_patch@indices]
  dist_b <- distance_matrix[b, divide_patch@indices]
  
  # indices 
  A <- divide_patch@indices[dist_b / dist_a >= relative_distance]
  B <- divide_patch@indices[dist_a / dist_b >= relative_distance]
  
  # basepoints
  basepoints_a <- A[as.integer(arrayInd(which.max(distance_matrix[A, A]), rep(length(A), 2)))]
  basepoints_b <- B[as.integer(arrayInd(which.max(distance_matrix[B, B]), rep(length(B), 2)))]
  
  # update divide_patch
  divide_patch@children <- length(cover@subsets) + 1:2
  divide_patch@death <- divide_patch@diameter
  cover@subsets[[index]] <- divide_patch
  # add new patches
  cover@subsets[length(cover@subsets) + 1:2] <- list(patch(A, 
                                                           basepoints = basepoints_a, 
                                                           diameter = distance_matrix[basepoints_a[1], basepoints_a[2]], 
                                                           birth = divide_patch@death,
                                                           parent = divide_patch@id, 
                                                           ancestor = c(divide_patch@id, divide_patch@ancestor), 
                                                           id = length(cover@subsets) + 1L), 
                                                     patch(B, 
                                                           basepoints = basepoints_b, 
                                                           diameter = distance_matrix[basepoints_b[1], basepoints_b[2]], 
                                                           birth = divide_patch@death, 
                                                           parent = divide_patch@id, 
                                                           ancestor = c(divide_patch@id, divide_patch@ancestor), 
                                                           id = length(cover@subsets) + 2L))
  # return cover
  return(cover)
}