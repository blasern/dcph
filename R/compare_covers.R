#' Distance between covers
#' 
#' Compare two covers using persistent homology or entropy
#' 
#' @param ... covers to compare 
#' @param metric which metric should be calculated
#' @param dim dimension(s) for bottleneck distance
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
#' dc1 <- divisive_cover(data = data_matrix,
#'                       stop_fct = stop_relative_filter(relative_filter = 0.1))
#' dc2 <- divisive_cover(data = data_matrix, 
#'                       stop_fct = stop_relative_filter(relative_filter = 0.3))
#' cover1 <- subcover(dc1, relative_filter = 0.1, method = "snapshot")
#' cover2 <- subcover(dc2, relative_filter = 0.5, method = "snapshot")
#' 
#' cover_distance(cover1, cover2, metric = "vi")
#' cover_distance(cover1, cover2, metric = "bottleneck", dim = 0)
#' cover_distance(cover1, cover2, metric = "bottleneck", dim = 1)
#' cover_distance(cover1, cover2, metric = "interleaving")
#' 
#' # plot covers
#' \dontrun{
#' plot(cover1)
#' plot(cover2)
#' plot(c(cover1, cover2))
#' }
#' @importFrom infotheo condentropy
#' @importFrom rdist cdist
#' @importFrom utils combn
#' @importFrom combinat permn
#' @export
cover_distance <- function(..., metric = c("vi", "bottleneck", "interleaving"), dim = 1){
  metric <- match.arg(metric)
  covers <- unlist(list(...))
  switch(metric, 
         "vi" = cover_distance_vi(covers), 
         "bottleneck" = cover_distance_bottleneck(covers, dim = dim), 
         "interleaving" = cover_distance_interleaving(covers))
}

cover_distance_bottleneck <- function(covers, dim){
  # compute persistence diagrams
  pers <- lapply(covers, snapshot_persistence, max_dim = max(dim))
  # compute bottleneck distance
  dist_obj <- apply(combn(1:length(covers), 2), 2, 
                    function(index){
                      bottleneck_distance(pers[[index[1]]], pers[[index[2]]], dim = dim)
                    })
  # matrix output
  class(dist_obj) <- "dist"
  attr(dist_obj, "Size") <- length(covers)
  as.matrix(dist_obj)
}

snapshot_persistence <- function(cover, max_dim){
  # remove duplicates
  cover <- simplify_cover(cover, method = "duplicates")
  
  # change death 
  cover@subsets <- c(list(patch(indices = 1:length(unique(unlist(lapply(cover@subsets, slot, "indices")))), 
                                id = 0L, 
                                filter_value = cover@data_filter_value)),
                     cover@subsets)
                     
  # calculate persistent homology  
  pers <- persistence_from_cover(cover, max_dim = max_dim, 
                                 representation =  "vector_vector", reduction = "twist")

  # add column names 
  pers <- as.data.frame(pers)
  colnames(pers) <- c("Dimension", "Birth", "Death")
  
  # remove diagonal
  pers <- pers[pers$Birth != pers$Death, ]
  rownames(pers) <- NULL
  
  pers
}

# has to be optimized
bottleneck_distance <- function(pers1, pers2, dim){
  dim_dist <- sapply(dim, function(loc_dim){
    # get dimension right
    diag1 <- pers1[pers1[, "Dimension"] == loc_dim, c("Birth", "Death"), drop = FALSE]
    diag2 <- pers2[pers2[, "Dimension"] == loc_dim, c("Birth", "Death"), drop = FALSE]
    # make diag2 smaller than diag1
    if (nrow(diag1) < nrow(diag2)){
      temp <- diag1
      diag1 <- diag2
      diag2 <- temp
    }
    
    # point distances
    point_distances <- rdist::cdist(diag1, diag2, "maximum")
    
    # diagonal distances
    diag1_distances <- (diag1[, "Death"] - diag1[, "Birth"]) / 2
    diag2_distances <- (diag2[, "Death"] - diag2[, "Birth"]) / 2
    
    # adjust point differences by diagonal 
    point_diag_distances <- matrix(apply(expand.grid(diag1_distances, diag2_distances), 1, max), 
                                   nrow = nrow(diag1), 
                                   ncol = nrow(diag2))
    adjusted_point_distances <- pmin(point_distances, point_diag_distances)

    # possible combinations
    if (nrow(diag2) > 0){
      df1 <- do.call(rbind.data.frame, combinat::permn(1:nrow(diag2))) 
      df2 <- data.frame(t(combn(1:nrow(diag1), nrow(diag2))))
      all_combinations <- do.call(rbind, apply(df1, 1, function(x) cbind(t(x), df2)))
      colnames(all_combinations) <- c(paste0("first_", 1:ncol(df1)), paste0("second_", 1:ncol(df2)))
      # combination distance
      comb_distances <- apply(all_combinations, 1, function(x){
        point_d <- max(sapply(seq(length(x)/2), 
                              function(index, x){
                                adjusted_point_distances[x[paste0("second_", index)], x[paste0("first_", index)]]
                              }, x = x))
        
        non_point_d <- suppressWarnings(
          max(c(diag1_distances[-x[grepl("second", names(x))]], 
                diag2_distances[-x[grepl("first", names(x))]])))
        c(point_d, non_point_d)
        
        max(point_d, non_point_d)
      })
    } 
    else {
      comb_distances <- max(diag1_distances)
    }
    
    # minimum combination
    min(comb_distances)
  })
  
  # return sum of distances by dimension
  sum(dim_dist)
}

cover_distance_interleaving <- function(covers){
  # compute interleaving distance
  dist_obj <- apply(combn(1:length(covers), 2), 2, 
                    function(index){
                      interleaving_epsilon(covers[[index[1]]], covers[[index[2]]]) + 
                      interleaving_epsilon(covers[[index[2]]], covers[[index[1]]])
                    })
  # matrix output
  class(dist_obj) <- "dist"
  attr(dist_obj, "Size") <- length(covers)
  as.matrix(dist_obj)
}

interleaving_epsilon <- function(cover1, cover2){
  # get filter values
  diams1 <- sapply(cover1@subsets, "slot", "filter_value")
  diams2 <- sapply(cover2@subsets, "slot", "filter_value")
  
  # calculate minimal filter value of superset
  superdiams <- diams2[superset(cover1, cover2)]
  
  # maximal filter value
  superdiams[is.na(superdiams)] <- cover2@data_filter_value
  
  # calculate epsilon 
  epsilon <- max(superdiams - diams1)
  return(epsilon)
}

# @param x,y covers
superset <- function(x, y){
  # lengths of subsets
  lx <- length(x@subsets)
  ly <- length(y@subsets)
  # initialize result
  minimal_superset <- rep(NA_integer_, lx)
  adjmat <- as.adjacency(c(x, y))[1:lx,(lx + 1):(lx+ly)]
  dimmat <- matrix(rep(sapply(lapply(x@subsets, slot, "indices"), length), ly), nrow = lx) 
  #
  indices <- which(adjmat == dimmat, arr.ind = TRUE)
  if (nrow(indices) > 0){
    # get index with minimal filter value
    min_indices <- sapply(split(indices[, 2], indices[, 1]), function(index){
      index[which.min(sapply(y@subsets[index], slot, "filter_value"))]
    })
    replace_indices <- sapply(split(indices[, 1], indices[, 1]), "[", 1)
    minimal_superset[replace_indices] <- min_indices
  }  
  return(minimal_superset)
}

cover_distance_vi <- function(covers){
  # define cover matrices
  cover_mats <- lapply(covers, cover_matrix)
  
  # compute interleaving distance
  dist_obj <- apply(combn(1:length(covers), 2), 2, 
                    function(index){
                      infotheo::condentropy(cover_mats[[index[1]]], cover_mats[[index[2]]]) + 
                        infotheo::condentropy(cover_mats[[index[2]]], cover_mats[[index[1]]])
                    })
  # matrix output
  class(dist_obj) <- "dist"
  attr(dist_obj, "Size") <- length(covers)
  as.matrix(dist_obj)
}
