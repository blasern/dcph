#' Anchor functions for divisive cover
#' 
#' @param points indices of points in patch 
#' @param data data
#' @param distance_fct distance function
#' @param group the group for classification
#' @param ... ignored arguments
#' 
#' @name anchor_fct
#' @export
anchor_extremal <- function(points, data, distance_fct, ...){
  points[arrayInd(which.max(distance_fct(data, points, points)), rep(length(points), 2))]
}

#' @rdname anchor_fct
#' @export
anchor_heuristic_extremal <- function(points, data, distance_fct, ...){
  p1 <- points[which.max(distance_fct(data, points[1], points))]
  p2 <- points[which.max(distance_fct(data, p1, points))]
  c(p1, p2)
}

#' @rdname anchor_fct
#' @export
anchor_classify <- function(points, data, distance_fct, group, ...){
  # if only one group, use extremal 
  if (length(unique(group[points])) == 1){
    return(anchor_extremal(points = points, data = data, distance_fct = distance_fct) )
  }
  # find biggest groups
  tab <- table(group[points])
  biggest_groups <- tail(names(tab)[order(tab)], 2)
  # calculate distances between biggest groups
  points1 <- points[group[points] == biggest_groups[1]]
  points2 <- points[group[points] == biggest_groups[2]]
  distances <- distance_fct(data, points1, points2)
  # find largest distance between biggest groups
  ai <- arrayInd(which.max(distances), c(length(points1), length(points2)))
  c(points1[ai[1]], points2[ai[2]])
}

#' @rdname anchor_fct
#' @export
anchor_heuristic_classify <- function(points, data, distance_fct, group, ...){
  # if only one group, use extremal 
  if (length(unique(group[points])) == 1){
    return(anchor_heuristic_extremal(points = points, data = data, distance_fct = distance_fct) )
  }
  # find biggest groups
  tab <- table(group[points])
  biggest_groups <- tail(names(tab)[order(tab)], 2)
  # calculate distances between biggest groups
  points1 <- points[group[points] == biggest_groups[1]]
  points2 <- points[group[points] == biggest_groups[2]]
  # approximate largest distance between biggest groups
  p1 <- points1[which.max(distance_fct(data, points2[1], points1))]
  p2 <- points2[which.max(distance_fct(data, p1, points2))]
  c(p1, p2)
}