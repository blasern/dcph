#' divisive cover 
#'
#' This function divides the data into a cover using the distance matrix. 
#'
#' @param data_matrix An data.frame or matrix (required if distance_object is a function)
#' @param distance_object an n x n matrix of pairwise dissimilarities or a distance function 
#' such as \code{dist} or \code{flexclust::dist2}
#' @param distance_mode If \code{distance_object} is a matrix, \code{distance_mode} should be 
#' \code{"matrix"}, if \code{distance_object} is a function that takes 1 argument as a 
#' input and returns pairwise distances, then  \code{distance_mode} should be \code{"dist"}, 
#' and if \code{distance_object} is a function that takes 2 matrices as input and returns 
#' the distances between these two matrices, then \code{distance_mode} should be \code{"dist2"}.
#' @param gap numeric to specify the amount of overlap in each division
#' @param  numeric to specify the maximal diameter of the final cover
divisive_cover <- function(data_matrix = NULL, 
                           distance_object = NULL, 
                           distance_mode = guess_distance_mode(distance_object),
                           relative_distance = 0.2, 
                           relative_radius = 0.7, 
                           base_point = NULL){
  # data size
  N <- ifelse(distance_mode == "matrix", 
              nrow(as.matrix(distance_object)), 
              nrow(data_matrix))
  # base point
  if (is.null(base_point)){
    base_point <- sample.int(N, 1) 
  }
  base_point <- as.integer(base_point)
  
  # generate initial cover
  radius <- get_radius(patch(1:N, basepoint = base_point), 
                       data_matrix = data_matrix, distance_object = distance_object, distance_mode = distance_mode)
  cover <- cover(list(patch(1:N, basepoint = base_point, id = 1L, radius = radius)))
  
  # minimal radius
  min_radius <- radius * relative_radius
  radii <- radius
  
  # divide into pieces
  while (any(radii > min_radius)){
    index <- which.max(radii)
    cover <- divide(cover = cover, 
                    index = index, 
                    data_matrix = data_matrix,
                    distance_object = distance_object, 
                    distance_mode = distance_mode,                        
                    relative_distance = relative_distance)
    new_radii <- sapply(tail(cover@subsets, 2), slot, "radius")
    radii[index] <- -Inf
    radii <- c(radii, new_radii)
  }
  cover
}

distance_function <-  function(x, y, data_matrix, distance_object, distance_mode){
  if (distance_mode == "matrix"){
    d <- as.matrix(distance_object)[x, y]
  }
  if (distance_mode == "dist"){
    d <- as.matrix(distance_object(data_matrix))[x, y]
  }
  if (distance_mode == "dist2"){
    d <- distance_object(data_matrix[x, ], data_matrix[y, ])
  }
  return(d)
} 

get_radius <- function(patch, data_matrix, distance_object, distance_mode){
  max(distance_function(patch@basepoint, patch@indices, data_matrix, distance_object, distance_mode))
}

divide <- function(cover, index, data_matrix, distance_object, distance_mode, relative_distance){
  divide_patch <- cover@subsets[[index]]
  
  # base points
  a <- divide_patch@basepoint
  dist_a <- distance_function(a,  divide_patch@indices, 
                              data_matrix = data_matrix, distance_object = distance_object, distance_mode = distance_mode)
  b <- divide_patch@indices[which.max(dist_a)]
  dist_b <- distance_function(b,  divide_patch@indices, 
                              data_matrix = data_matrix, distance_object = distance_object, distance_mode = distance_mode)
  
  # indices 
  A <- divide_patch@indices[dist_b / dist_a >= 1 - relative_distance]
  B <- divide_patch@indices[dist_a / dist_b >= 1 - relative_distance]
  
  # update divide_patch
  divide_patch@children <- length(cover@subsets) + 1:2
  cover@subsets[[index]] <- divide_patch
  # add new patches
  cover@subsets[length(cover@subsets) + 1:2] <- list(patch(A, 
                                                           basepoint = a, 
                                                           radius = max(dist_a[match(divide_patch@indices, A, nomatch = 0)]), 
                                                           parent = divide_patch@id, 
                                                           id = length(cover@subsets) + 1L), 
                                                     patch(B, 
                                                           basepoint = b, 
                                                           radius = max(dist_b[match(divide_patch@indices, B, nomatch = 0)]),
                                                           parent = divide_patch@id, 
                                                           id = length(cover@subsets) + 2L))
  # return cover
  return(cover)
}

guess_distance_mode <- function(distance_object){
  if (is.matrix(distance_object)) dmode <- "matrix"
  else if (is.function(distance_object)) {
    dm <- suppressWarnings(try(as.matrix(distance_object(diag(c(0, 0)), diag(c(0, 0)))), silent = TRUE))
    dmode <- ifelse(is.matrix(dm), "dist2", "dist")
    }
  else {
    stop("Distance object must be either a function or a matrix")
  }
  if (dmode == "dist"){
    dm <- suppressWarnings(try(as.matrix(distance_object(diag(c(0, 0)))), silent = TRUE))
    if (!is.matrix(dm)){
      stop("Could not determine the mode of distance_object")
    }
  }
  return(dmode)
}