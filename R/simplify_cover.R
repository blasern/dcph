#' Simplify cover
#' 
#' Simplify cover by removing duplicates or subsets
#' 
#' @param x cover
#' @param method simplification method
#' @export
simplify_cover <- function(x, method = c("none", "duplicates", "subsets")){
  # apply simplification method
  method <- match.arg(method)
  switch(method, 
         'none' = x, 
         'duplicates' = remove_duplicates(x), 
         'subsets' = remove_subsets(x))
}

remove_duplicates <- function(x){
  simple <- which(!duplicated(lapply(sapply(x@subsets, slot, "indices"), sort)))
  x@subsets <- x@subsets[simple]
  attr(x, "simple") <- simple
  x
}

remove_subsets <- function(x){
  # which are the subsets
  is_subset <- rep(FALSE, length(x@subsets))
  adjmat <- as.adjacency(x)
  dimmat <- matrix(rep(sapply(sapply(x@subsets, slot, "indices"), length), length(x@subsets)), nrow = length(x@subsets)) 
  #
  is_subset[which(adjmat == dimmat, arr.ind = TRUE)[, 1]] <- TRUE
  x@subsets <- x@subsets[!is_subset]
  attr(x, "simple") <- which(!is_subset)
  x
}