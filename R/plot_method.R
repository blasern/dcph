#' Plot cover output
#' 
#' \code{plot.cover} generates R plots of covers
#' 
#' @param x A cover
#' @param coloring color values. See \code{\link{cover_coloring}}
#' @param simplify simplification method. See \code{\link{simplify_cover}}
#' @param legend should color legend be displayed
#' @param seed random seed used to make plots look the same 
#' @param device use igraph::plot.igraph or igraph::tkplot
#' @param ... additional plotting parameters. See \code{\link[igraph]{igraph.plotting}} for the complete list.
#' @author Nello Blaser
#' @keywords plotting
#' @seealso \code{\link{as.cover}}
#' @importFrom igraph graph_from_adjacency_matrix plot.igraph E
#' @export 
plot.cover <- function(x, coloring = NULL, simplify = c("none", "duplicates", "subsets"), legend = TRUE, device = c("plot", "tkplot"), seed = 1, ...){
  x <- simplify_cover(x, method = match.arg(simplify))
  plot_core(cover = x, adjacency = as.adjacency(x), coloring = coloring, legend = legend, device = device, seed = seed, ...)
}

plot_core <- function(cover, adjacency, coloring = NULL, legend = TRUE, device = c("plot", "tkplot"), seed = 1, ...){
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
                        vertex.label = NA, 
                        vertex.color = coloring,
                        edge.width = resize(igraph::E(g1)$weight, min = 1, max = 5), 
                        vertex.size = resize(sapply(lapply(cover@subsets, slot, "indices"), length), min = 10, max = 20), 
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
                   vertex.label=NA, 
                   vertex.color = coloring,
                   edge.width = resize(igraph::E(g1)$weight, min = 1, max = 5), 
                   vertex.size = resize(sapply(cover, length), min = 10, max = 20), 
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
#' @importFrom RColorBrewer brewer.pal
#' @export
cover_coloring <- function(x, color_values){
  resized_values <- resize(color_values, 0, 1)
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  node_color_value <- function(values){
    if (is.numeric(values)) return(mean(values))
    else return(as.character(getmode(values)))
  }
  
  node_legend_values <- sapply(x$points_in_vertex, 
                               function(point) node_color_value(color_values[point]))
  node_values <- sapply(x$points_in_vertex, 
                        function(point) node_color_value(resized_values[point]))
  legend_values <- quantile(node_legend_values, c(0, .25, .5, .75, 1))
  
  cols <- apply(colorRamp(RColorBrewer::brewer.pal(9, "RdYlGn"))(node_values) / 255, 
                1, function(col) rgb(col[1], col[2], col[3]))
  legend_cols <- apply(colorRamp(RColorBrewer::brewer.pal(9, "RdYlGn"))(resize(legend_values, 0, 1)) / 255, 
                       1, function(col) rgb(col[1], col[2], col[3]))
  
  attr(cols, "legend") <- data.frame(value = legend_values, 
                                     color = legend_cols, 
                                     stringsAsFactors = FALSE)
  
  return(cols)
}

# function to resize nodes, edges, colors
resize <- function(x, min, max){
  if (length(x) == 0) return(NULL)
  if (var(x) == 0) return((max + min)/2)
  lambda <- (max - min)/ diff(range(x, na.rm = TRUE))
  min + (x - min(x)) * lambda
}
