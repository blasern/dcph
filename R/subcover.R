#' Subcover
#' 
#' Find the (divisive) cover that corresponds to a certain filter
#' 
#' @param cover A divisive cover
#' @param relative_filter The filter to extract
#' @param method Either "divisive" to get a divisive cover or "snapshot" 
#' to get only the lowest level
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
#' dc_0.7 <- subcover(dc, relative_filter = 0.7, method = "divisive")
#' sc_0.7 <- subcover(dc, relative_filter = 0.7, method = "snapshot")
#' 
#' \dontrun{
#' plot(data_matrix)
#' plot(sc_0.7)
#' }
#' @export
subcover <- function(cover, relative_filter, method = c("divisive", "snapshot", "range")){
  method <- match.arg(method)
  # absolute filter
  filter <- relative_filter * cover@subsets[[1]]@filter_value
  birth <- sapply(cover@subsets, slot, "parent_filter")
  death <- sapply(cover@subsets, slot, "filter_value")
  death[cover@external_nodes] <- 0
  # survivors
  surv <- switch(method, 
         "divisive" = filter < birth, 
         "snapshot" = filter >= death & filter < birth, 
         "range" = filter[1] < birth & filter[2] >= death) 
  # new cover
  cv <- cover
  cv@parameters$relative_filter <- relative_filter
  cv@subsets <- cv@subsets[surv]
  cv@type <- method
  cv
}