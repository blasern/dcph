#' divisive cover 
#'
#' This function divides the data into a cover using the distance matrix. 
#'
#' @param distance_matrix an n x n matrix of pairwise dissimilarities
#' @param delta delta-parameter for delta-filtered cover (0 < delta <= 1/2)
#' @param relative_diameter maximal diameter of the final cover
#' @param max_nodes maximum number of nodes the algorithm should generate
#' @param cover a divisive cover used to update
#' @examples
#' rcircle <- function(N, r, sd){
#'    radius <- rnorm(N, r, sd)
#'    angle <- runif(N, 0, 2 * pi)
#'    data.frame(x = radius * cos(angle), 
#'               y = radius * sin(angle))
#' }
#' data_matrix <- rcircle(200, 1, .1)
#' 
#' # calculate and update divisive cover
#' dc1 <- divisive_cover(distance_matrix = dist(data_matrix),
#'                       delta = 0.1, 
#'                       relative_diameter = 0.7)
#' dc2 <- update_divisive_cover(dc1, relative_diameter = 0.5)
#' 
#' # get one snapshot for plotting
#' ddc <- subcover(dc2, relative_diameter = 0.7, method = "snapshot")
#' 
#' \dontrun{
#' plot(data_matrix)
#' plot(ddc)
#' }
#' @export
divisive_cover <- function(distance_matrix, 
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
  distance_matrix <- as.matrix(distance_matrix)
  N <- nrow(distance_matrix)
  
  # generate initial cover
  basepoints <- as.integer(arrayInd(which.max(distance_matrix), dim(distance_matrix)))
  diameter <- distance_matrix[basepoints[1], basepoints[2]]
  cover <- cover(distance_matrix = distance_matrix, 
                 subsets = list(patch(1:N, basepoints = basepoints, id = 1L, diameter = diameter, birth = diameter)), 
                 parameters = list(delta = delta, relative_diameter = relative_diameter), 
                 type = "divisive")
  
  # minimal radius
  min_diam <- diameter * relative_diameter
  diams <- diameter
  # divide into pieces
  while (any(diams > min_diam) && sum(diams > -Inf) < max_nodes){
    index <- which.max(diams)
    cover <- divide(cover = cover, 
                    index = index, 
                    distance_matrix = distance_matrix,                     
                    relative_distance = relative_distance)
    new_diams <- sapply(tail(cover@subsets, 2), slot, "diameter")
    diams[index] <- -Inf
    diams <- c(diams, new_diams)
  }
  cover
}

#' @rdname divisive_cover
#' @export
update_divisive_cover <- function(cover, 
                         delta = cover@parameters[["delta"]], 
                         relative_diameter = 0.1, 
                         max_nodes = Inf){
  if (cover@type != "divisive"){
    stop("can only update divisive covers")
  }
  # find factor
  relative_distance <- (1-2 * delta)/(1+2 * delta)
  # extract diameters
  diameter <- cover@subsets[[1]]@diameter
  min_diam <- diameter * relative_diameter
  survivors <- sapply(cover@subsets, slot, "death") == 0
  diams <- sapply(cover@subsets, slot, "diameter")
  diams[!survivors] <- -Inf
  # update division
  # divide into pieces
  while (any(diams > min_diam) && sum(diams > -Inf) < max_nodes){
    index <- which.max(diams)
    cover <- divide(cover = cover, 
                    index = index, 
                    distance_matrix = cover@distance_matrix,                     
                    relative_distance = relative_distance)
    new_diams <- sapply(tail(cover@subsets, 2), slot, "diameter")
    diams[index] <- -Inf
    diams <- c(diams, new_diams)
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