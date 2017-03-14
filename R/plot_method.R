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
                   vertex.label = label, 
                   vertex.color = coloring,
                   edge.width = resize(igraph::E(g1)$weight, min = 1, max = 5), 
                   vertex.size = resize(sapply(lapply(cover@subsets, slot, "indices"), length), min = 10, max = 20), 
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
#' @importFrom grDevices colorRamp rgb
#' @importFrom stats quantile var
#' @export
cover_coloring <- function(x, color_values){
  resized_values <- resize(color_values, 0, 1)
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  node_color_value <- function(values){
    if (is.numeric(values)) return(mean(values, na.rm = TRUE))
    else return(as.character(getmode(values)))
  }
  
  node_legend_values <- sapply(lapply(x@subsets, slot, "indices"), 
                               function(point) node_color_value(color_values[point]))
  node_values <- sapply(lapply(x@subsets, slot, "indices"), 
                        function(point) node_color_value(resized_values[point]))
  legend_values <- quantile(node_legend_values, c(0, .25, .5, .75, 1), na.rm = TRUE)
  
  cols <- get_color(node_values)
  legend_cols <- get_color(resize(legend_values, 0, 1))
  
  attr(cols, "legend") <- data.frame(value = legend_values, 
                                     color = legend_cols, 
                                     stringsAsFactors = FALSE)
  
  return(cols)
}

get_color <- function(values){
  rgb_col <- grDevices::colorRamp(RColorBrewer::brewer.pal(9, "RdYlGn"))(values) / 255
  rgb_col[is.na(rgb_col)] <- 0
  apply(rgb_col, 1, function(col) grDevices::rgb(col[1], col[2], col[3]))
}

#' @rdname cover_coloring
#' @export
predict_coloring <- function(x){
  node_legend_values <- sapply(lapply(x@subsets, slot, "predicted"), length)
  node_values <- resize(node_legend_values, 0, 1)
  legend_values <- quantile(node_legend_values, c(0, .25, .5, .75, 1))
  
  cols <- get_color(node_values)
  legend_cols <- get_color(resize(legend_values, 0, 1))
  
  attr(cols, "legend") <- data.frame(value = legend_values, 
                                     color = legend_cols, 
                                     stringsAsFactors = FALSE)
  
  return(cols)
}


# function to resize nodes, edges, colors
resize <- function(x, min, max){
  if (length(x) == 0) return(NULL)
  if (length(x) == 1) return((max + min)/2)
  if (all(is.na(x))) return(x)
  if (stats::var(x, na.rm = TRUE) == 0) return((max + min)/2)
  lambda <- (max - min)/ diff(range(x, na.rm = TRUE))
  min + (x - min(x, na.rm = TRUE)) * lambda
}
