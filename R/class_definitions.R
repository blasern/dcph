#' @useDynLib dcph
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
  N_data <- max(nrow(object@data), nrow(object@distance_matrix), na.rm = TRUE)
  
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
                  representation = representation(distance_matrix = "matrix", 
                                                  data = "matrix", 
                                                  diameter = "numeric",
                                                  subsets = "list",
                                                  parameters = "list", type = "character"), 
                  validity = check_cover)

cover <- function(distance_matrix = matrix(numeric(0)), data = matrix(numeric(0)), 
                  subsets, parameters = list(), diameter = max(distance_matrix),
                  type = c("divisive", "fast_divisive", "snapshot", "predict", "concat")){
  new("cover", 
      distance_matrix = as.matrix(distance_matrix), 
      data = as.matrix(data), 
      diameter = diameter,
      subsets = subsets, 
      parameters = parameters, 
      type = match.arg(type))
}
is.cover <- function (x){
  inherits(x, "cover")
}

# function to check that object is a patch
check_patch <- function(object){
  # error message
  errors <- character()
  
  # check that death occurs after birth
  if (!is.na(object@birth) & !is.na(object@death) & object@death > object@birth) {
    msg <- "death cannot be larger than birth"
    errors <- c(errors, msg)
  }
  
  # # check that basepoints are in indices
  # if (!all(object@basepoints %in% object@indices)) {
  #   msg <- "the basepoints need to be in the patch"
  #   errors <- c(errors, msg)
  # }
  
  if (length(errors) == 0) TRUE else errors
}

#' An S4 class to represent a patch.
#' @slot indices Indices of data points in patch
#' @slot predicted Indices of predicted points
#' @slot children Indices of children in cover
#' @slot parent Indices of parent in cover
#' @slot birth Birth time
#' @slot death Death time
#' @rdname cover
setClass(Class = "patch", 
                  representation = representation(id = "integer", 
                                                  indices = "integer",
                                                  predicted = "integer",
                                                  basepoints = "integer", 
                                                  children = "integer", 
                                                  parent = "integer",
                                                  ancestor = "integer",
                                                  birth = "numeric", 
                                                  death = "numeric", 
                                                  diameter = "numeric"), 
                  prototype = prototype(indices = integer(0),
                                        predicted = integer(0), 
                                        basepoints = integer(0), 
                                        id = NA_integer_,
                                        children = integer(0), 
                                        parent = integer(0), 
                                        ancestor = integer(0),
                                        birth = NA_real_, 
                                        death = 0, 
                                        diameter = NA_real_),
                  validity = check_patch)

# patch
patch <- function(indices, predicted = integer(0), basepoints = integer(0), id = NA_integer_, children = integer(0), parent = integer(0), ancestor = integer(0), birth = NA_real_, death = 0, diameter = NA_real_){
  new("patch", indices = indices, predicted = predicted, basepoints = basepoints, id = id, children = children, parent = parent, ancestor = ancestor, birth = birth, death = death, diameter = diameter)
}

is.patch <- function (x){
  inherits(x, "patch")
}
  