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
  simple <- which(!duplicated(lapply(lapply(x@subsets, slot, "indices"), sort)))
  x@subsets <- x@subsets[simple]
  # update internal and external nodes
  internal_nodes <- match(x@internal_nodes, simple)
  x@internal_nodes <- internal_nodes[!is.na(internal_nodes)]
  external_nodes <- match(x@external_nodes, simple)
  x@external_nodes <- external_nodes[!is.na(external_nodes)]
  # return 
  attr(x, "simple") <- simple
  x
}

remove_subsets <- function(x){
  # remove duplicates first
  x <- remove_duplicates(x)
  # which are the subsets
  is_subset <- rep(FALSE, length(x@subsets))
  adjmat <- as.adjacency(x)
  dimmat <- matrix(rep(sapply(lapply(x@subsets, slot, "indices"), length), length(x@subsets)), nrow = length(x@subsets)) 
  #
  is_subset[which(adjmat == dimmat, arr.ind = TRUE)[, 1]] <- TRUE
  x@subsets <- x@subsets[!is_subset]
  # update internal and external nodes
  simple <- which(!is_subset)
  internal_nodes <- match(x@internal_nodes, simple)
  x@internal_nodes <- internal_nodes[!is.na(internal_nodes)]
  external_nodes <- match(x@external_nodes, simple)
  x@external_nodes <- external_nodes[!is.na(external_nodes)]
  # return 
  attr(x, "simple") <- attr(x, "simple")[!is_subset]
  x
}