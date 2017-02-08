#' Agglomerative cover
#' 
#' Get agglomerative cover
#' 
#' @param cover original cover
#' 
#' @examples 
#' # generate data
#' rcircle <- function(N, r, sd){
#'   radius <- rnorm(N, r, sd)
#'   angle <- runif(N, 0, 2 * pi)
#'   data.frame(x = radius * cos(angle), 
#'              y = radius * sin(angle))
#' }
#' data_matrix <- rbind(rcircle(200, 1, .1), 
#'                      data.frame(x = rnorm(100, 0, .1), 
#'                                 y = rnorm(100, 0, .1)))
#'
#' # run divisive cover
#' dc <- divisive_cover(distance_matrix = dist(data_matrix), 
#'                      delta = 0.1, 
#'                      relative_diameter = 0.2)
#' # get subcover
#' sc <- subcover(dc, 0.3, "snapshot")
#' 
#' # agglomerative cover from subcover
#' ac <- agglomerative_cover(sc)
#' @export
agglomerative_cover <- function(cover){
  k <- length(cover@subsets) + 1L
  cover@subsets <- lapply(cover@subsets, function(x){
    x@birth <- NA_real_
    x@death <- NA_real_
    x
  })
  subcover <- cover
  cov_dist <- cover_dist(subcover)
  while (length(subcover@subsets) > 1){
    new_cov_dist <- cover_dist(subcover)
    cov_dist[1:2] <- new_cov_dist[1:2]
    cov_dist[3] <- ifelse(new_cov_dist[3] == 0, cov_dist[3], new_cov_dist[3])
    print(cov_dist)
    cover@subsets[[k]] <- 
      patch(indices = union(cover@subsets[[which(sapply(cover@subsets, slot, "id") == cov_dist[1])]]@indices, 
                            cover@subsets[[which(sapply(cover@subsets, slot, "id") == cov_dist[2])]]@indices), 
            id = max(sapply(subcover@subsets, slot, "id")) + 1L, 
            birth = cov_dist[3], 
            death = NA_real_)
    cover@subsets[[which(sapply(cover@subsets, slot, "id") == cov_dist[1])]]@death <- cov_dist[3]
    cover@subsets[[which(sapply(cover@subsets, slot, "id") == cov_dist[2])]]@death <- cov_dist[3]
    subcover <- cover
    subcover@subsets <- cover@subsets[is.na(sapply(cover@subsets, slot, "death"))]
    k <- k + 1
  }
  cover@type <- "agglomerative"
  cover
}

cover_dist <- function(cover){
  dist_mat <- outer(lapply(cover@subsets, slot, "indices"), 
                   lapply(cover@subsets, slot, "indices"), 
                   FUN = function(x, y, distance_matrix) 
                     mapply(function(x, y) min(distance_matrix[setdiff(x, y), setdiff(y, x)]), x, y), 
                   distance_matrix = cover@distance_matrix)
  dist_mat[dist_mat == Inf] <- 0
  diag(dist_mat) <- Inf
  indices <- as.numeric(arrayInd(which.min(dist_mat), .dim = dim(dist_mat)))
  ids <- sapply(cover@subsets[indices], slot, "id")
  d <- dist_mat[indices[1], indices[2]]
  c(ids, d)
}
