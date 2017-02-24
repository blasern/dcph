#' Extract data from cover
#' 
#' Functions to extract relevant data from covers
#' 
#' @param cover A cover
#' @param ids Cover ids, for example extracted from plot
#' @param weights Should cover weights be calculated
#' @param data Data from which to extract subset
#' @name extract_cover
NULL

#' @rdname extract_cover
#' @export
index_from_cover_id <- function(cover, ids, weights = FALSE){
  cover_indices <- which(sapply(cover@subsets, slot, "id") %in% ids)
  data_indices <- unique(unlist(lapply(cover@subsets[cover_indices], slot, "indices")))
  if (weights){
    tab <- table(unlist(sapply(cover@subsets, slot, "indices")))[as.character(data_indices)]
    tab_int <- table(unlist(sapply(cover@subsets[cover_indices], slot, "indices")))[as.character(data_indices)]
    attr(data_indices, "weights") <- 1 / as.numeric(tab - tab_int + 1)
  }
  return(data_indices)
}

#' @rdname extract_cover
#' @export
data_from_cover_id <- function(data, cover, ids, weights = FALSE){
  data_indices <- index_from_cover_id(cover, ids, weights = weights)
  res <- data[data_indices, ]
  if (weights){
    res[, "weights"] <- attr(data_indices, "weights")
  }
  return(res)
}