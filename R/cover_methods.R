#' Use as cover
#' 
#' Convert a matrix or a clustering to a cover
#' 
#' @param x covering matrix or list
#' @export
as.cover <- function(x){
  UseMethod("as.cover")
}

setGeneric("as.cover")

as.cover.cover <- function(x) {
  x
}

as.cover.data.frame <- function(x){
  stopifnot(all(unlist(x) %in% c(0, 1)))
  lst <- lapply(x, function(y) which(y == 1))
  as.cover(lst)
}

as.cover.matrix <- function(x){
  as.cover(as.data.frame(x))
}

as.cover.list <- function(x){
  lst <- lapply(x, patch)
  cover(subsets = lst)
}

#' @export
c.cover <- function(..., recursive = FALSE){
  covers <- list(...)
  indices_lst <- lapply(covers, slot, "indices")
  ix <- seq.int(length(indices_lst))
  
  # check that the indices are the same 
  if (!all(outer(indices_lst, 
                 indices_lst, 
                 FUN = function(x, y) mapply(function(x, y) setequal(x, y), x, y)))){
    stop("Only covers of the same space can be concatinated.")
  }
  res <- new("cover", 
             indices = indices_lst[[1]], 
             subsets = do.call(c, lapply(covers, slot, "subsets")))
  res
}

#' @export
print.cover <- function(x, ...){
  cat(paste0("Cover of length ", length(x@subsets), " and covering dimension ", cover_dim(x), "\n")) 
  n <- min(10, length(x@subsets))
  head(x, n = n)
  if (length(x@subsets) > n) cat(paste0("... with ", length(x@subsets) - 10, " more subsets."))
}

#' @export
#' @importFrom methods show
#' @rdname cover-class
setMethod(show, signature = "cover", function(object){ print.cover(object) })

#' @importFrom utils head
#' @export
head.cover <- function(x, n = 6L, ...){
  x <- x@subsets
  n <- min(n, length(x))
  cat_subset(x[1:n])
  invisible(x[1:n])
}

#' @importFrom utils tail
#' @export
tail.cover <- function(x, n = 6L, ...){
  x <- x@subsets
  n <- min(n, length(x))
  cat_subset(x[(length(x)-n+1):length(x)], 
             index = (length(x)-n+1):length(x))
  invisible(x[(length(x)-n+1):length(x)])
}

cat_subset <- function(lst, index = 1:length(lst)) {
  lines <- paste0("Subset ", index, " consisting of rows: ", sapply(lapply(lst, slot, "indices"), paste, collapse = " "))
  lines_trunc <- paste0(sapply(lines, str_trunc, max_width = getOption("width")))
  cat(lines_trunc, sep = "\n")
}

# from tibble
str_trunc <- function(x, max_width) 
{
  width <- nchar(x)
  for (i in seq_along(x)) {
    if (width[i] <= max_width[i]) 
      next
    x[i] <- paste0(substr(x[i], 1, max_width[i] - 3), "...")
  }
  x
}

# covering dimension
cover_dim <- function(x){
  max(table(unlist(lapply(x@subsets, slot, "indices"))))
}

# adjacency matrix
as.adjacency <- function(cover){
  adj_mat <- outer(lapply(cover@subsets, slot, "indices"), 
                   lapply(cover@subsets, slot, "indices"), 
                   FUN = function(x, y) mapply(function(x, y) length(intersect(x, y)), x, y))
  diag(adj_mat) <- 0
  adj_mat
}
