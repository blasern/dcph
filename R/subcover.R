#' Subcover
#' 
#' Find the (divisive) cover that corresponds to a certain filter
#' 
#' @param cover A divisive cover
#' @param method Either \code{"external"} to get only external nodes or 
#' \code{"divisive"} to get a divisive cover or \code{"snapshot"}
#' to get only the the nodes that were external nodes at filter
#' value \code{relative_filter}, or \code{"range"} for all subsets 
#' within a range of filter values. 
#' @param relative_filter The filter to extract
#' @examples
#' rcircle <- function(N, r, sd){
#'    radius <- rnorm(N, r, sd)
#'    angle <- runif(N, 0, 2 * pi)
#'    data.frame(x = radius * cos(angle), 
#'               y = radius * sin(angle))
#' }
#' data_matrix <- rcircle(200, 1, .1)
#' 
#' dc <- divisive_cover(data = data_matrix, 
#'                      division_fct = relative_gap_division(0.1), 
#'                      stop_fct = stop_relative_filter(relative_filter = 0.5))
#' dc_0.7 <- subcover(dc, method = "divisive", relative_filter = 0.7)
#' sc_0.7 <- subcover(dc, method = "snapshot", relative_filter = 0.7)
#' sc_ext <- subcover(dc, method = "external")
#' 
#' \dontrun{
#' plot(data_matrix)
#' plot(sc_0.7)
#' plot(sc_ext)
#' }
#' @export
subcover <- function(cover, method = c("external", "divisive", "snapshot", "range"), relative_filter){
  method <- match.arg(method)
  if (method == "external"){
    return(external_subcover(cover))
  }
  else {
    return(relative_filter_subcover(cover = cover, relative_filter = relative_filter, method = method))
  }
}

relative_filter_subcover <- function(cover, relative_filter, method = c("divisive", "snapshot", "range")){
  method <- match.arg(method)
  # absolute filter
  filter <- relative_filter * cover@subsets[[1]]@filter_value
  birth <- sapply(cover@subsets, slot, "parent_filter")
  death <- sapply(cover@subsets, slot, "filter_value")
  death[cover@external_nodes] <- 0
  # survivors
  surv <- switch(method, 
                 "divisive" = filter <= birth, 
                 "snapshot" = filter >= death & filter <= birth, 
                 "range" = filter[1] <= birth & filter[2] >= death) 
  # new cover
  cv <- cover
  cv@parameters$relative_filter <- relative_filter
  cv@subsets <- cv@subsets[surv]
  cv@type <- method
  cv
}

external_subcover <- function(cover){
  cover@subsets <- cover@subsets[cover@external_nodes]
  cover@type <- "external"
  cover
}