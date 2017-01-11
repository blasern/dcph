#' @importFrom methods setClass 
NULL

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
  if (!setequal(seq.int(nrow(object@data)), 
                unique(unlist(lapply(object@subsets, slot, "index"))))) {
    msg <- paste0("subsets do not cover the data")
    errors <- c(errors, msg)
  }
  
  if (length(errors) == 0) TRUE else errors
}

#' An S4 class to represent a cover.
#'
#' @slot data A data.frame
#' @slot subsets A list of patches
setClass(Class = "cover", 
                  representation = representation(data = "data.frame", subsets = "list"), 
                  prototype = prototype(data = data.frame(), subsets = list()), 
                  validity = check_cover)

is.cover <- function (x){
  inherits(x, "cover")
}

# function to check that object is a patch
check_patch <- function(object){
  # error message
  errors <- character()
  
  # check that death occurs after birth
  if (is.na(object@birth)) {
    msg <- "no birth defined" 
    errors <- c(errors, msg)
  }
  if (!is.na(object@death) & object@death < object@birth) {
    msg <- "death cannot be smaller than birth"
    errors <- c(errors, msg)
  }
  
  # check that patches have only one parent 
  if (length(object@parent) != 1){
    msg <- "patches can only have 1 parent"
    errors <- c(errors, msg)
  }
  
  if (length(errors) == 0) TRUE else errors
}

#' An S4 class to represent a patch.
#' @slot index Indices of data points in patch
#' @slot children Indices of children in cover
#' @slot parent Indices of parent in cover
#' @slot birth Birth time
#' @slot death Death time
#' @rdname cover
setClass(Class = "patch", 
                  representation = representation(index = "integer",
                                                  basepoint = "integer", 
                                                  children = "integer", 
                                                  parent = "integer", 
                                                  birth = "numeric", 
                                                  death = "numeric"), 
                  prototype = prototype(index = integer(0), 
                                        basepoint = integer(0), 
                                        children = NA_integer_, 
                                        parent = NA_integer_, 
                                        birth = NA_real_, 
                                        death = NA_real_),
                  validity = check_patch)

is.patch <- function (x){
  inherits(x, "patch")
}
  