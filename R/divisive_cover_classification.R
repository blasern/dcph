#' divisive cover 
#'
#' This function divides the data into a cover using the distance matrix. 
#'
#' @param distance_matrix an n x n matrix of pairwise dissimilarities
#' @param response a length n response vector
#' @param delta delta-parameter for delta-filtered cover (0 < delta <= 1/2)
#' @param max_nodes maximum number of nodes the algorithm should generate
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
divisive_cover_classification <- function(distance_matrix,
                                          response, 
                                          delta, 
                                          max_nodes = Inf){
  stopifnot(0 <= delta && delta < 0.5)
  # use response as factor
  response <- as.factor(response)

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
                 parameters = list(delta = delta), 
                 type = "divisive")
  
  # minimal radius
  diams <- diameter
  # divide into pieces
  while (any(diams > 0) && sum(diams > -Inf) < max_nodes){
    index <- which.max(diams)
    cover <- divide_classify(cover = cover, 
                             index = index, 
                             distance_matrix = distance_matrix,                     
                             relative_distance = relative_distance, 
                             response = response)
    new_diams <- sapply(tail(cover@subsets, 2), slot, "diameter")
    diams[index] <- -Inf
    diams <- c(diams, new_diams)
    same_class <- sapply(tail(cover@subsets, 2), function(x){
      all(duplicated(response[x@indices])[-1])
    })
    diams[seq(length(diams)-1, by = 1, length = 2)[same_class]] <- -Inf
  }
  cover
}

divide_classify <- function(cover, index, distance_matrix, relative_distance, response = response){
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
  
  # responses
  responseA <- factor(levels(response)[which.max(table(response[A]))], levels = levels(response))
  responseB <- factor(levels(response)[which.max(table(response[B]))], levels = levels(response))
  
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
                                                           id = length(cover@subsets) + 1L, 
                                                           response = responseA), 
                                                     patch(B, 
                                                           basepoints = basepoints_b, 
                                                           diameter = distance_matrix[basepoints_b[1], basepoints_b[2]], 
                                                           birth = divide_patch@death, 
                                                           parent = divide_patch@id, 
                                                           ancestor = c(divide_patch@id, divide_patch@ancestor), 
                                                           id = length(cover@subsets) + 2L, 
                                                           response = responseB))
  # return cover
  return(cover)
}

#' Classify
#' 
#' Classify from divisive cover
#' 
#' @param pred a prediction cover
divisive_cover_classify <- function(pred){
  pred <- subcover(pred, 0, "snapshot")
  predicted_class <- factor(rep(NA, length = length(unique(unlist(lapply(pred@subsets, slot, "predicted"))))), 
                            levels = levels(pred@subsets[[1]]@response))
  for (i in seq_along(pred@subsets)){
    x <- pred@subsets[[i]]
    predicted_class[x@predicted] <- x@response
  }
  predicted_class 
}