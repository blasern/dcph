#' fast divisive cover 
#'
#' This function divides the data into a cover using the a distance function to compute
#' the distance between two matrices. Instead of using the full distance matrix to 
#' calculate the maximal distance in a cover, this maximal distance is approximated. It  
#' uses an initial point and takes the point with maximal distance from the initial point
#' as a first basepoint and the point with maximal distance from the first basepoint 
#' as the second basepoint. 
#'
#' @param data data that should be divided
#' @param distance_function function that compute distance matrix between two matrices
#' @param delta delta-parameter for delta-filtered cover (0 < delta <= 1/2)
#' @param relative_diameter maximal diameter of the final cover
#' @param max_nodes maximum number of nodes the algorithm should generate
#' @examples
#' rcircle <- function(N, r, sd){
#'   radius <- rnorm(N, r, sd)
#'   angle <- runif(N, 0, 2 * pi)
#'   data.frame(x = radius * cos(angle), 
#'              y = radius * sin(angle))
#' }
#' data_matrix <- rcircle(200, 1, .1)
#' # calculate fast divisive cover
#' dc <- fast_divisive_cover(data = data_matrix,
#'                           distance_function = pdist::pdist,
#'                           delta = 0.1, 
#'                           relative_diameter = 0.5)
#' \dontrun{
#' plot(subcover(dc, .5, "snapshot"))
#' }
#' @export
fast_divisive_cover <- function(data,
                                distance_function,
                                delta, 
                                relative_diameter = 0, 
                                max_nodes = Inf){
  # check input
  if (relative_diameter == 0 && max_nodes == Inf){
    stop("Either relative diameter or max_nodes has to be specified.")
  }
  stopifnot(0 <= delta && delta < 0.5)
  
  # find factor
  relative_distance <- (1-2 * delta)/(1+2 * delta)
  
  # data size
  N <- nrow(data)
  data_matrix <- as.matrix(data)
  
  # generate initial cover
  basepoints <- get_basepoints_fast(data_matrix, distance_function)
  diameter <- as.numeric(as.matrix(distance_function(data_matrix[basepoints[1], ], data_matrix[basepoints[2], ])))
  cover <- cover(data = as.matrix(data),
                 subsets = list(patch(1:N, basepoints = basepoints, id = 1L, diameter = diameter, birth = diameter)), 
                 parameters = list(delta = delta, relative_diameter = relative_diameter, distance_function = distance_function), 
                 type = "fast_divisive")
  
  # minimal radius
  min_diam <- diameter * relative_diameter
  diams <- diameter
  # divide into pieces
  while (any(diams > min_diam) && sum(diams > -Inf) < max_nodes){
    index <- which.max(diams)
    cover <- divide_fast(cover = cover, 
                         index = index, 
                         data_matrix = data_matrix,
                         distance_function = distance_function,                     
                         relative_distance = relative_distance)
    new_diams <- sapply(tail(cover@subsets, 2), slot, "diameter")
    diams[index] <- -Inf
    diams <- c(diams, new_diams)
  }
  cover
}

get_basepoints_fast <- function(mat, distance_function){
  a <- which.max(as.matrix(distance_function(mat[1, ], mat)))
  b <- which.max(as.matrix(distance_function(mat[a, ], mat)))
  return(c(a, b))
}

divide_fast <- function(cover, index, data_matrix, distance_function, relative_distance){
  divide_patch <- cover@subsets[[index]]
  
  # base points
  a <- divide_patch@basepoints[1]
  b <- divide_patch@basepoints[2]
  dist_a <- as.matrix(distance_function(data_matrix[a, , drop = FALSE], data_matrix[divide_patch@indices, , drop = FALSE]))
  dist_b <- as.matrix(distance_function(data_matrix[b, , drop = FALSE], data_matrix[divide_patch@indices, , drop = FALSE]))
  
  # indices 
  A <- divide_patch@indices[dist_b / dist_a >= relative_distance]
  B <- divide_patch@indices[dist_a / dist_b >= relative_distance]
  
  # basepoints
  basepoints_a <- A[get_basepoints_fast(data_matrix[A, , drop = FALSE], distance_function)]
  basepoints_b <- B[get_basepoints_fast(data_matrix[B, , drop = FALSE], distance_function)]
  
  stopifnot(basepoints_a %in% A, basepoints_b %in% B)
  
  # update divide_patch
  divide_patch@children <- length(cover@subsets) + 1:2
  divide_patch@death <- divide_patch@diameter
  cover@subsets[[index]] <- divide_patch
  # add new patches
  cover@subsets[length(cover@subsets) + 1:2] <- list(patch(A, 
                                                           basepoints = basepoints_a, 
                                                           diameter = as.numeric(as.matrix(distance_function(data_matrix[basepoints_a[1], ], data_matrix[basepoints_a[2], ]))), 
                                                           birth = divide_patch@death,
                                                           parent = divide_patch@id, 
                                                           ancestor = c(divide_patch@id, divide_patch@ancestor), 
                                                           id = length(cover@subsets) + 1L), 
                                                     patch(B, 
                                                           basepoints = basepoints_b, 
                                                           diameter = as.numeric(as.matrix(distance_function(data_matrix[basepoints_b[1], ], data_matrix[basepoints_b[2], ]))), 
                                                           birth = divide_patch@death, 
                                                           parent = divide_patch@id, 
                                                           ancestor = c(divide_patch@id, divide_patch@ancestor), 
                                                           id = length(cover@subsets) + 2L))
  # return cover
  return(cover)
}