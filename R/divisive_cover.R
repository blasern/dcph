#' divisive cover 
#'
#' This function divides the data into a cover using the distance matrix. 
#'
#' @param distance_matrix an n x n matrix of pairwise dissimilarities
#' @param relative_distance numeric to specify the amount of overlap in each division
#' @param relative_diameter numeric to specify the maximal diameter of the final cover
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
#' ddc <- subcover(dc, relative_diameter = 0.7, method = "snapshot")
#' 
#' \dontrun{
#' plot(data_matrix)
#' plot(ddc)
#' }
#' @export
divisive_cover <- function(distance_matrix, 
                           relative_distance = 0.2, 
                           relative_diameter = 0.7){
  # data size
  distance_matrix <- as.matrix(distance_matrix)
  N <- nrow(distance_matrix)
  
  # generate initial cover
  basepoints <- as.integer(arrayInd(which.max(distance_matrix), dim(distance_matrix)))
  diameter <- distance_matrix[basepoints[1], basepoints[2]]
  cover <- cover(distance_matrix = distance_matrix, 
                 subsets = list(patch(1:N, basepoints = basepoints, id = 1L, diameter = diameter, birth = diameter)), 
                 parameters = list(relative_distance = relative_distance, relative_diameter = relative_diameter), 
                 type = "divisive")
  
  # minimal radius
  min_diam <- diameter * relative_diameter
  diams <- diameter
  # divide into pieces
  while (any(diams > min_diam)){
    index <- which.max(diams)
    cover <- divide(cover = cover, 
                    index = index, 
                    distance_matrix = distance_matrix,                     
                    relative_distance = relative_distance)
    new_diams <- sapply(tail(cover@subsets, 2), slot, "diameter")
    diams[index] <- -Inf
    diams <- c(diams, new_diams)
    data.frame(diams = diams, sapply(lapply(cover@subsets, slot, "indices"), length))
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
  A <- divide_patch@indices[dist_b / dist_a >= 1 - relative_distance]
  B <- divide_patch@indices[dist_a / dist_b >= 1 - relative_distance]
  
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
                                                           id = length(cover@subsets) + 1L), 
                                                     patch(B, 
                                                           basepoints = basepoints_b, 
                                                           diameter = distance_matrix[basepoints_b[1], basepoints_b[2]], 
                                                           birth = divide_patch@death, 
                                                           parent = divide_patch@id, 
                                                           id = length(cover@subsets) + 2L))
  # return cover
  return(cover)
}