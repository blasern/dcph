#' Extract data from cover
#' 
#' Functions to extract relevant data from covers
#' 
#' @param cover A cover
#' @param ids Cover ids, for example extracted from plot
#' @param data Data from which to extract subset
#' @name extract_cover
NULL

#' @rdname extract_cover
index_from_cover_id <- function(data, cover, ids){
  cover_indices <- which(sapply(cover@subsets, slot, "id") %in% ids)
  data_indices <- unique(unlist(lapply(cover@subsets[cover_indices], slot, "indices")))
  data_indices
}

#' @rdname extract_cover
data_from_cover_id <- function(data, cover, ids){
  data_indices <- index_from_cover_id()
  data[data_indices, ]
}