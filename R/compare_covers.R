#' Distance between covers
#' 
#' Compare two covers using persistent homology or entropy
#' 
#' @param cover1,cover2 covers to compare 
#' @param metric which metric should be calculated
#' @param dim dimension(s) for nerve distance
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
#' dc1 <- divisive_cover(distance_matrix = dist(data_matrix), 
#'                       delta = 0.1,
#'                       relative_diameter = 0.1)
#' dc2 <- divisive_cover(distance_matrix = dist(data_matrix), 
#'                       delta = 0.05,
#'                       relative_diameter = 0.3)
#' cover1 <- subcover(dc1, 0.1, "snapshot")
#' cover2 <- subcover(dc2, 0.5, "snapshot")
#' 
#' cover_distance(cover1, cover2, metric = "vi")
#' 
#' # plot covers
#' \dontrun{
#' plot(cover1)
#' plot(cover2)
#' plot(c(cover1, cover2))
#' }
#' @importFrom infotheo condentropy
#' @importFrom rdist cdist
#' @export
cover_distance <- function(cover1, cover2, metric = c("vi", "nerve", "intertwining"), dim = 1){
  metric <- match.arg(metric)
  switch(metric, 
         "vi" = cover_distance_vi(cover1, cover2), 
         "nerve" = cover_distance_nerve(cover1, cover2, dim = dim), 
         "intertwining" = cover_distance_intertwining(cover1, cover2))
}

cover_distance_nerve <- function(cover1, cover2, dim){
  # compute persistence diagrams
  pers1 <- snapshot_persistence(cover1, max_dim = max(dim))
  pers2 <- snapshot_persistence(cover2, max_dim = max(dim))
  # compute bottleneck distance
  bottleneck_distance(pers1, pers2, dim = dim)
}

snapshot_persistence <- function(cover, max_dim){
  # remove duplicates
  cover <- simplify_cover(cover, method = "duplicates")
  
  # change death 
  cover@subsets <- c(list(patch(1:length(unique(unlist(lapply(cover@subsets, slot, "indices")))), 
                                id = 0L, 
                                diameter = max(cover@distance_matrix), 
                                death = max(cover@distance_matrix),
                                birth = max(cover@distance_matrix))),
    lapply(cover@subsets, function(x) {
    x@death <- x@diameter
    x
  }))
  
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
    diag1 <- rbind(data.frame(Birth = 0, Death = 0), diag1)
    diag2 <- rbind(data.frame(Birth = 0, Death = 0), diag2)
    # make diag2 smaller than diag1
    if (nrow(diag1) < nrow(diag2)){
      temp <- diag1
      diag1 <- diag2
      diag2 <- temp
    }
    
    # point distances
    point_distances <- rdist::cdist(diag1, diag2, "maximum")
    
    # diagonal distances
    diag1_distances <- (diag1[, "Birth"] + diag1[, "Death"]) / 2
    diag2_distances <- (diag2[, "Birth"] + diag2[, "Death"]) / 2
    
    # adjust point differences by diagonal 
    point_diag_distances <- matrix(rowSums(expand.grid(diag1_distances, diag2_distances)), 
                                   nrow = nrow(diag1), 
                                   ncol = nrow(diag2))
    adjusted_point_distances <- pmin(point_distances, point_diag_distances)
    
    # possible combinations
    df1 <- data.frame(t(c(1:nrow(diag2))))
    colnames(df1) <- paste0("first_", 1:ncol(df1))
    df2 <- data.frame(t(combn(1:nrow(diag1), nrow(diag2))))
    colnames(df2) <- paste0("second_", 1:ncol(df2))
    all_combinations <- cbind(df1, df2)
    
    # combination distance
    comb_distances <- apply(all_combinations, 1, function(x){
      point_d <- sum(sapply(seq(length(x)/2), 
                            function(index, x){
                              adjusted_point_distances[x[paste0("second_", index)], x[paste0("first_", index)]]
                            }, x = x))
      
      non_point_d <- sum(diag1_distances[-x[grepl("second", names(x))]]) + 
        sum(diag2_distances[-x[grepl("first", names(x))]])
      
      point_d + non_point_d
    })
    
    # minimum combination
    min(comb_distances)
  })
  
  # return sum of distances by dimension
  sum(dim_dist)
}

cover_distance_intertwining <- function(cover1, cover2){
  stop("Intertwining distance not defined")
  eps1 <- intertwining_epsilon(cover1, cover2)
  eps2 <- intertwining_epsilon(cover2, cover1)
  return(eps1 + eps2)
}

intertwining_epsilon <- function(cover1, cover2){
  # maximal diameter
  maxdiam <- cover1@diameter
  stopifnot(isTRUE(all.equal(maxdiam, cover2@diameter)))
  # get diameters
  diams1 <- sapply(cover1@subsets, "slot", "diameter")
  diams2 <- sapply(cover2@subsets, "slot", "diameter")
  
  # calculate epsilon 
  epsilon <- 1
  # subsets
  return(epsilon)
}

cover_distance_vi <- function(cover1, cover2){
  # define cover matrices
  x_mat <- cover_matrix(cover1)
  y_mat <- cover_matrix(cover2)
  
  cond1 <- infotheo::condentropy(x_mat, y_mat)
  cond2 <- infotheo::condentropy(y_mat, x_mat)
  
  vi <- cond1 + cond2
#   attr(vi, "H(1|2)") <- cond1
#   attr(vi, "H(2|1)") <- cond2
  return(vi)
}

cover_distance_ph <- function(cover1, cover2){
  # recalculate cover diameters
  adjust_diameter <- function(cover, diameter){
    adjust_subset_diam <- function(subset, id, diameter){
      subset@diameter <- diameter
      subset@birth <- diameter
      subset@death <- diameter
      subset
    }
    cover@subsets <- mapply(adjust_subset_diam, 
                            subset = cover@subsets, 
                            diameter = diameter)
    cover
  }
  
  contract <- cover(
    distance_matrix = cover1@distance_matrix,
    subsets = list(patch(indices = unique(unlist(lapply(cover1@subsets, slot, "indices"))), 
                         death = 2, 
                         birth = 2, 
                         diameter = 2)))
  
  pers21 <- persistent_homology(c(contract, adjust_diameter(cover1, 1), adjust_diameter(cover2, 0)))
  pers12 <- persistent_homology(c(contract, adjust_diameter(cover2, 1), adjust_diameter(cover1, 0)))
  list(`1->2` = pers12, `2->1` = pers21)
}