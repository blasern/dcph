#' Subcover
#' 
#' Find the (divisive) cover that corresponds to a certain diameter
#' 
#' @param cover A divisive cover
#' @param relative_diameter The diameter to extract
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
#' dc <- divisive_cover(distance_matrix = dist(data_matrix), 
#'                      delta = 0.1, 
#'                      stop_fct = stop_relative_diameter(relative_diameter = 0.5))
#' dc_0.7 <- subcover(dc, relative_diameter = 0.7, method = "divisive")
#' sc_0.7 <- subcover(dc, relative_diameter = 0.7, method = "snapshot")
#' 
#' \dontrun{
#' plot(data_matrix)
#' plot(sc_0.7)
#' }
#' @export
subcover <- function(cover, relative_diameter, method = c("divisive", "snapshot", "range")){
  method <- match.arg(method)
  # absolute diameter
  diameter <- relative_diameter * cover@subsets[[1]]@diameter
  birth <- sapply(cover@subsets, slot, "birth")
  death <- sapply(cover@subsets, slot, "death")
  # survivors
  surv <- switch(method, 
         "divisive" = diameter < birth, 
         "snapshot" = diameter >= death & diameter < birth, 
         "range" = diameter[1] < birth & diameter[2] >= death) 
  # new cover
  cv <- cover
  cv@parameters$relative_diameter <- relative_diameter
  cv@subsets <- cv@subsets[surv]
  cv@type <- method
  cv
}