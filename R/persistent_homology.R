persistent_homology <- function(cover){
  filtration <- get_filtration(cover) 
  
  # calculate persistent homology
  pers <- persistence_from_filtration(filtration)
  pers <- as.data.frame(pers)
  colnames(pers) <- c("Birth", "Death")
  pers$Dimension <- filtration[pers$Birth, "dim"]
  pers <- pers[, c("Dimension", "Birth", "Death")]
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
  
  # radii
  radii <- sapply(cover@subsets, slot, "radius")
  radii <- 
    cummin(radii)
  radii
  
  # filtration
  filt <- as.data.frame(which(arr > 0, arr.ind = TRUE))
  filt <- t(apply(filt, 1, sort))
  filt <- filt[!duplicated(filt), ]
  colnames(filt) <- c("vertex_1", "vertex_2", "vertex_3")
  filt <- as.data.frame(filt)
  # dimension and filter value
  filt[, "dim"] <- 2 - apply(filt, 1, function(x){sum(duplicated(x))})
  filt[, "filter_value"] <- radii[pmax(filt[, "vertex_1"], filt[, "vertex_2"], filt[, "vertex_3"])]
  # order
  filt <- filt[order(pmax(filt$vertex_1, filt$vertex_2, filt$vertex_3), filt$dim), ]
  filt <- filt[!(filt$dim == 1 & filt$vertex_1 == filt$vertex_2), ]
  # faces
  filt[, c("face_1", "face_2", "face_3")] <- NA_integer_
  for (i in 1:nrow(filt)){
    if (filt[i, "dim"] == 2){
      filt[i, "face_1"] <- which((filt[i, "vertex_1"] == filt[, "vertex_1"]) & (filt[i, "vertex_2"] == filt[, "vertex_2"]))[1]
      filt[i, "face_2"] <- which((filt[i, "vertex_1"] == filt[, "vertex_1"]) & (filt[i, "vertex_3"] == filt[, "vertex_2"]))[1]
      filt[i, "face_3"] <- which((filt[i, "vertex_2"] == filt[, "vertex_1"]) & (filt[i, "vertex_3"] == filt[, "vertex_2"]))[1]
    }
  }
  # convert to matrix
  filt <- as.matrix(filt)
  filt[, c("dim", "vertex_1", "vertex_2", "vertex_3", "face_1", "face_2", "face_3", "filter_value")]
}