#' Persistent homology
#' 
#' Calculate the persistent homology from the divisive cover complex
#' 
#' @param cover The divisive cover
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
#'                      relative_diameter = 0.5, relative_distance = 0.2)
#' pers <- persistent_homology(dc)
#' \dontrun{
#' maxpers <- max(pers[, c("Birth", "Death")])
#' plot(pers$Birth, pers$Death, col = pers$Dimension, 
#'      xlim = c(0, maxpers), ylim = c(0, maxpers))
#' legend("bottomright", pch = 1, col = sort(unique(pers$Dimension)),
#'        legend = sort(unique(pers$Dimension)))
#' lines(c(0, maxpers), c(0, maxpers))
#' }
#' @export
persistent_homology <- function(cover){
  # calculate the filtration
  filtration <- get_filtration(cover) 
  # calculate persistent homology
  pers <- persistence_from_filtration(filtration)
  # add column names 
  pers <- as.data.frame(pers)
  colnames(pers) <- c("Birth", "Death")
  pers$Dimension <- filtration[pers$Birth, "dim"]
  pers <- pers[, c("Dimension", "Birth", "Death")]
  # use filter value instead of index
  pers[, "Birth"] <- filtration[pers[, "Birth"], "filter_value"]
  pers[, "Death"] <- filtration[pers[, "Death"], "filter_value"]
  # remove what dies imediately
  pers <- pers[pers[, "Birth"] != pers[, "Death"], ]
  return(pers)
}

get_filtration <- function(cover){
  # 3-fold overlap
  x <- lapply(cover@subsets, slot, "indices")
  arr <- sapply(x, function(a){
    sapply(x, function(a, b){
      sapply(x, function(a, b, d){
        length(Reduce(intersect, list(a, b, d)))
      }, a = a, b = b)
    }, a = a)
  })
  dim(arr) <- rep(length(x), 3)
  
  # diameters
  diams <- sapply(cover@subsets, slot, "diameter")
  parent <- sapply(cover@subsets, slot, "parent")
  diams <- diams[parent]
  diams[1] <- diams[2]
  
  # filtration
  filt <- as.data.frame(which(arr > 0, arr.ind = TRUE))
  filt <- t(apply(filt, 1, sort))
  filt <- filt[!duplicated(filt), ]
  colnames(filt) <- c("vertex_1", "vertex_2", "vertex_3")
  filt <- as.data.frame(filt)
  # dimension and filter value
  filt[, "dim"] <- 2 - apply(filt, 1, function(x){sum(duplicated(x))})
  filt[, "filter_value"] <- diams[pmin(filt[, "vertex_1"], filt[, "vertex_2"], filt[, "vertex_3"])]
  # order
  filt <- filt[order(filt$filter_value, filt$dim), ]
  filt <- filt[!(filt$dim == 1 & filt$vertex_1 == filt$vertex_2), ]
  # faces
  filt[, c("face_1", "face_2", "face_3")] <- NA_integer_
  for (i in 1:nrow(filt)){
    if (filt[i, "dim"] >= 1){
      filt[i, "vertex_1"] <- which((filt[i, "vertex_1"] == filt[, "vertex_1"]) & filt[, "dim"] == 0)
      filt[i, "vertex_2"] <- which((filt[i, "vertex_2"] == filt[, "vertex_1"]) & filt[, "dim"] == 0)
      filt[i, "vertex_3"] <- which((filt[i, "vertex_3"] == filt[, "vertex_1"]) & filt[, "dim"] == 0)
    }
    if (filt[i, "dim"] == 2){
      filt[i, "face_1"] <- which((filt[i, "vertex_1"] == filt[, "vertex_1"]) & (filt[i, "vertex_2"] == filt[, "vertex_2"]))[1]
      filt[i, "face_2"] <- which((filt[i, "vertex_1"] == filt[, "vertex_1"]) & (filt[i, "vertex_3"] == filt[, "vertex_2"]))[1]
      filt[i, "face_3"] <- which((filt[i, "vertex_2"] == filt[, "vertex_1"]) & (filt[i, "vertex_3"] == filt[, "vertex_2"]))[1]
    }
  }
  # sort vertices
  filt[, c("vertex_1", "vertex_2", "vertex_3")] <- t(apply(filt[, c("vertex_1", "vertex_2", "vertex_3")], 1, sort))
  filt[, c("face_1", "face_2", "face_3")] <- t(apply(filt[, c("face_1", "face_2", "face_3")], 1, sort, na.last = TRUE))
  # convert to matrix
  filt <- as.matrix(filt)
  filt[, c("dim", "vertex_1", "vertex_2", "vertex_3", "face_1", "face_2", "face_3", "filter_value")]
}