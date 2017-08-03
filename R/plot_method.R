#' Plot cover output
#' 
#' \code{plot.cover} generates R plots of covers
#' 
#' @param x A cover
#' @param coloring color values. See \code{\link{cover_coloring}}
#' @param simplify simplification method. See \code{\link{simplify_cover}}
#' @param label how should the nodes be labeled. Set NA for no labels. 
#' @param legend should color legend be displayed
#' @param seed random seed used to make plots look the same 
#' @param device use igraph::plot.igraph or igraph::tkplot
#' @param ... additional plotting parameters. See \code{\link[igraph]{igraph.plotting}} for the complete list.
#' @author Nello Blaser
#' @keywords plotting
#' @seealso \code{\link{cover-class}}, \code{\link{subcover}}
#' @examples
#' # generate sample data
#' rcircle <- function(N, r, sd){
#'   radius <- rnorm(N, r, sd)
#'   angle <- runif(N, 0, 2 * pi)
#'   data.frame(x = radius * cos(angle), 
#'              y = radius * sin(angle))
#' }
#' data_matrix <- rcircle(200, 1, .1)
#' 
#' # calculate divisive cover
#' dc <- divisive_cover(data = data_matrix,
#'                      division_fct = relative_gap_division(0.05), 
#'                      stop_fct = stop_relative_filter(0.5))
#' 
#' # external nodes
#' sc <- subcover(dc)
#' 
#' # distance
#' patches <- slot(sc, "subsets")
#' d_centroid <- patch_centroid_dist(data_matrix, patches, patches)
#' 
#' # colors 
#' coloring <- cover_coloring(sc, color_values = data_matrix[, 1], FUN = median)
#' 
#' # plot
#' \dontrun{
#'   plot(sc, coloring = coloring)
#'   # custom layout
#'   plot(sc, coloring = coloring, layout = igraph::layout_with_kk)
#'   plot(sc, coloring = coloring, layout = igraph::norm_coords(cmdscale(d_centroid)))
#' }
#' @importFrom igraph graph_from_adjacency_matrix plot.igraph E
#' @method plot cover
#' @export
plot.cover <- function(x, coloring = NULL, simplify = c("none", "duplicates", "subsets"), 
                       label = sapply(x@subsets, slot, "id"), legend = TRUE, device = c("plot", "tkplot"), seed = 1, ...){
  force(label)
  simplify <-  match.arg(simplify)
  x <- simplify_cover(x, method = simplify)
  if (simplify != "none"){
    simple <- attr(x, "simple")
    label <- label[simple]
    color_legend <- attr(coloring, "legend")
    coloring <- coloring[simple]
    attr(coloring, "legend") <- color_legend 
  }
  plot_core(cover = x, adjacency = as.adjacency(x), coloring = coloring, 
            label = label, legend = legend, device = device, seed = seed, ...)
}

plot_core <- function(cover, adjacency, label = NA, coloring = NULL, legend = TRUE, device = c("plot", "tkplot"), seed = 1, ...){
  # set random seed for repetability
  old_seed <- .Random.seed
  on.exit(.Random.seed <<- old_seed)
  set.seed(seed)
  
  # set device
  device <- match.arg(device)
  
  # to igraph object
  g1 <- igraph::graph_from_adjacency_matrix(adjacency, 
                                            mode="undirected", 
                                            weighted = TRUE)
  
  # in case no coloring specified 
  if (is.null(coloring)) {
    coloring <- rep(1, length(cover))
  }
  
  # plot
  if (device == "plot"){
    igraph::plot.igraph(g1, 
                        vertex.label = label, 
                        vertex.color = coloring,
                        edge.width = scales::rescale(igraph::E(g1)$weight, to = c(1, 5)), 
                        vertex.size = scales::rescale(sapply(lapply(cover@subsets, slot, "indices"), length), to = c(10, 20)), 
                        ...
    )
    if (legend && !is.null(attr(coloring, "legend"))){
      legend('topright', 
             legend = round(attr(coloring, "legend")[, "value"], 2), 
             fill = attr(coloring, "legend")[, "color"])
    }
  }
  else {
    igraph::tkplot(g1, 
                   vertex.label = label, 
                   vertex.color = coloring,
                   edge.width = scales::rescale(igraph::E(g1)$weight, to = c(1, 5)), 
                   vertex.size = scales::rescale(sapply(lapply(cover@subsets, slot, "indices"), length), to = c(10, 20)), 
                   ...
    )
  }
}

#' Generate coloring for cover
#' 
#' Generates a coloring for cover to be used with \code{\link{plot.cover}}
#' 
#' @param x cover
#' @param color_values value for each row
#' @param FUN summary function 
#' @param color_value_range range of color values
#' @param ... other orguments passed to FUN
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRamp rgb
#' @importFrom stats quantile var
#' @importFrom scales rescale
#' @export
cover_coloring <- function(x, color_values, FUN = mean, color_value_range = range(color_values), ...){
  node_values <- sapply(lapply(x@subsets, slot, "indices"), 
                        function(point) FUN(color_values[point], ...))
  resized_node_values <- scales::rescale(node_values, to = c(0, 1), from = color_value_range)
  
  legend_values <- seq(min(node_values, na.rm = TRUE), max(node_values, na.rm = TRUE), length = 5)
  
  cols <- get_color(resized_node_values)
  legend_cols <- get_color(scales::rescale(legend_values, to = c(0, 1), from = color_value_range))
  
  attr(cols, "legend") <- data.frame(value = legend_values, 
                                     color = legend_cols, 
                                     stringsAsFactors = FALSE)
  return(cols)
}

get_color <- function(values){
  rgb_col <- grDevices::colorRamp(rev(RColorBrewer::brewer.pal(9, "RdYlGn")))(values) / 255
  rgb_col[is.na(rgb_col)] <- 1
  apply(rgb_col, 1, function(col) grDevices::rgb(col[1], col[2], col[3]))
}

#' @rdname cover_coloring
#' @export
predict_coloring <- function(x){
  node_legend_values <- sapply(lapply(x@subsets, slot, "predicted"), length)
  max_value <- length(unique(unlist(lapply(x@subsets, slot, "predicted"))))
  node_values <- scales::rescale(node_legend_values, to = c(0, 1), from = c(0, max_value))
  legend_values <- seq(0, 1, length = 5)
  
  cols <- get_color(node_values)
  legend_cols <- get_color(legend_values)
  
  attr(cols, "legend") <- data.frame(value = legend_values, 
                                     color = legend_cols, 
                                     stringsAsFactors = FALSE)
  
  return(cols)
}
