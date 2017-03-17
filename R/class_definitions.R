#' @useDynLib dcph, .registration = TRUE 
#' @importFrom methods setClass 
#' @importFrom Rcpp sourceCpp
NULL
.onUnload <- function (libpath) {
  library.dynam.unload("dcph", libpath)
}

# function to check that object is a cover
check_cover <- function(object){
  # error message
  errors <- character()
  
  # check that all subsets are patches
  if (any(!sapply(object@subsets, is.patch))) {
    msg <- paste0("subsets ", 
                  which(!sapply(object@subsets, is.patch)), 
                  " are not patches")
    errors <- c(errors, msg)
  }
  
  # check that the subsets cover the data
  N_data <- nrow(object@data)
  
  if (!setequal(1:N_data, 
                unique(unlist(lapply(object@subsets, slot, "indices"))))) {
    msg <- paste0("subsets do not cover the data")
    errors <- c(errors, msg)
  }
  
  if (length(errors) == 0) TRUE else errors
}

#' An S4 class to represent a cover.
#'
#' @slot data A data.frame
#' @slot subsets A list of patches
#' @param object A cover
setClass(Class = "cover", 
                  representation = representation(data = "matrix", 
                                                  subsets = "list",
                                                  parameters = "list", 
                                                  type = "character"), 
                  validity = check_cover)

cover <- function(data, subsets, 
                  parameters = list(), 
                  type = c("divisive", "fast_divisive", "snapshot", "predict", "concat")){
  new("cover", 
      data = as.matrix(data), 
      subsets = subsets, 
      parameters = parameters, 
      type = match.arg(type))
}
is.cover <- function (x){
  inherits(x, "cover")
}

#' An S4 class to represent a patch.
#' @slot id Patch id
#' @slot indices Indices of data points in patch
#' @slot predicted Indices of predicted points
#' @slot anchor_points Anchor points
#' @slot filter_value Filter value
#' @slot parent_filter Filter value of parent
#' @rdname cover
setClass(Class = "patch", 
                  representation = representation(id = "integer", 
                                                  indices = "integer",
                                                  predicted = "integer",
                                                  anchor_points = "integer", 
                                                  filter_value = "numeric", 
                                                  parent_filter = "numeric"), 
                  prototype = prototype(id = NA_integer_,
                                        indices = integer(0),
                                        predicted = integer(0), 
                                        anchor_points = integer(0), 
                                        filter_value = NA_real_, 
                                        parent_filter = NA_real_))

# patch
patch <- function(id = NA_integer_,
                  indices = integer(0),
                  predicted = integer(0), 
                  anchor_points = integer(0), 
                  filter_value = NA_real_, 
                  parent_filter = NA_real_){
  new("patch", id = id, indices = indices, predicted = predicted, anchor_points = anchor_points, 
      filter_value = filter_value, parent_filter = parent_filter)
}

is.patch <- function (x){
  inherits(x, "patch")
}