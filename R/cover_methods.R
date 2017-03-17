#' @importFrom methods slot new
#' @export
c.cover <- function(..., recursive = FALSE){
  covers <- list(...)
  data <- lapply(covers, slot, "data")
  
  # check that the dimensions are the same 
  if (!all(sapply(data, nrow) == nrow(data[[1]]))){
    stop("Only covers of the same space can be concatinated.")
  }
  res <- cover(data = data[[1]], 
               subsets = do.call(c, lapply(covers, slot, "subsets")),
               type = "concat")
  res
}

#' @export
print.cover <- function(x, ...){
  cat(paste0("Cover of length ", length(x@subsets), " and covering dimension ", cover_dim(x), "\n")) 
  n <- min(10, length(x@subsets))
  print(head(x, n = n))
  if (length(x@subsets) > n) cat(paste0("... with ", length(x@subsets) - 10, " more subsets."))
}

#' @export
#' @importFrom methods show
#' @rdname cover-class
setMethod(show, signature = "cover", function(object){ print.cover(object) })

#' @export
print.patch <- function(x, ...){
  print(patch_df(x), ...)
}

patch_df <- function(x){
  data.frame(id = x@id, size = length(x@indices), filter_value = x@filter_value, parent_filter = x@parent_filter)
}

#' @export
#' @importFrom methods show
#' @rdname cover-class
setMethod(show, signature = "patch", function(object){ print.patch(object) })

#' @importFrom utils head
#' @export
head.cover <- function(x, n = 6L, ...){
  x <- x@subsets[1:min(n, length(x@subsets))]
  df <- do.call(rbind, 
                lapply(x, patch_df))
  df
}

#' @importFrom utils tail
#' @export
tail.cover <- function(x, n = 6L, ...){
  n <- min(n, length(x@subsets))
  x <- x@subsets[(length(x@subsets)-n+1):length(x@subsets)]
  df <- do.call(rbind, 
                lapply(x, patch_df))
  df
}

# covering dimension
cover_dim <- function(x){
  max(table(unlist(lapply(x@subsets, slot, "indices"))))
}

# adjacency matrix
#' @importFrom methods slot new
as.adjacency <- function(cover){
  adj_mat <- outer(lapply(cover@subsets, slot, "indices"), 
                   lapply(cover@subsets, slot, "indices"), 
                   FUN = function(x, y) mapply(function(x, y) length(intersect(x, y)), x, y))
  diag(adj_mat) <- 0
  adj_mat
}
