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
#' @slot internal_nodes internal nodes
#' @slot external_nodes external nodes
#' @slot data_filter_values the filter value of the entire data set
#' @slot parameter list of parameters
#' @slot type type of cover
#' @param object A cover
#' @name cover-class
setClass(Class = "cover", 
                  representation = representation(data = "matrix", 
                                                  subsets = "list",
                                                  internal_nodes = "integer", 
                                                  external_nodes = "integer",
                                                  data_filter_value = "numeric", 
                                                  parameters = "list", 
                                                  type = "character"), 
                  validity = check_cover)

cover <- function(data, subsets, 
                  internal_nodes = integer(0), 
                  external_nodes = integer(0), 
                  data_filter_value = numeric(0),
                  parameters = list(), 
                  type = c("divisive", "snapshot", "predict", "concat")){
  new("cover", 
      data = as.matrix(data), 
      subsets = subsets, 
      internal_nodes = internal_nodes, 
      external_nodes = external_nodes,
      data_filter_value = data_filter_value,
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
#' @slot children Indices of children
#' @slot filter_value Filter value
#' @slot parent Index of parent 
#' @slot parent_filter Filter value of parent
#' @rdname cover-class
setClass(Class = "patch", 
                  representation = representation(id = "integer", 
                                                  indices = "integer",
                                                  predicted = "integer",
                                                  anchor_points = "integer", 
                                                  children = "integer", 
                                                  filter_value = "numeric", 
                                                  parent = "integer",
                                                  parent_filter = "numeric"), 
                  prototype = prototype(id = NA_integer_,
                                        indices = integer(0),
                                        predicted = integer(0), 
                                        anchor_points = integer(0), 
                                        children = integer(0),
                                        filter_value = NA_real_, 
                                        parent = NA_integer_,
                                        parent_filter = NA_real_))

# patch
patch <- function(id = NA_integer_,
                  indices = integer(0),
                  predicted = integer(0), 
                  anchor_points = integer(0), 
                  children = integer(0),
                  filter_value = NA_real_, 
                  parent = NA_integer_,
                  parent_filter = NA_real_){
  new("patch", id = id, indices = indices, predicted = predicted, anchor_points = anchor_points, 
      children = children, filter_value = filter_value, parent = parent, parent_filter = parent_filter)
}

is.patch <- function (x){
  inherits(x, "patch")
}