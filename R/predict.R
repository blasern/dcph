#' Predict
#' 
#' Predict in which subsets new data come to lie when given a previous
#' divisive cover 
#' @param object a divisive cover
#' @param newdata new data 
#' @param predict_fct predict function (see below)
#' @param ... further arguments passed from other methods
#' @export
#' @name predict_cover
predict.cover <- function(object, newdata, predict_fct, ...){
  if (object@type != "divisive"){
    stop("can only predict from divisive cover")
  }
  # initialize new cover
  newcover <- object
  newcover@type = "predict"
  newcover@subsets[[1]]@predicted <- 1:nrow(newdata)
  
  # divide
  for (ix in 2:length(newcover@subsets)){
    newcover <- divide_pred(cover = newcover, 
                            index = ix, 
                            newdata = newdata, 
                            predict_fct = predict_fct)
  }

  return(newcover)
}

divide_pred <- function(cover, index, newdata, predict_fct){
  parent_patch <- cover@subsets[[cover@subsets[[index]]@parent]]
  
  # anchor points
  anchor_points <- parent_patch@anchor_points
  bp <- anchor_points %in% cover@subsets[[index]]@indices
  
  # predict
  cover@subsets[[index]]@predicted <- 
    predict_fct(data = cover@data, newdata = newdata, patch = parent_patch, 
                anchor = anchor_points[bp], 
                distance_fct = cover@parameters[["distance_fct"]])
  
  # return cover
  return(cover)
}

#' Extract group from predict cover
#' 
#' @param pc prediction cover
#' @param group specify the group of the training set (not necessary if used in divisive_cover)
#' @export
#' @importFrom dplyr group_by_ select_ do right_join ungroup
#' @examples
#' iris <- datasets::iris
#' 
#' ttsplit <- sample(1:2, nrow(iris), replace = TRUE, prob = c(2, 1))
#' train <- iris[ttsplit == 1, ]
#' test <- iris[ttsplit == 2, ]
#' 
#' dc <- divisive_cover(data = train[, 1:4], 
#'                      group = train[, 5], 
#'                      distance_fct = distance_cdist("euclidean"), 
#'                      stop_fct = stop_relative_filter(relative_filter = 0), 
#'                      anchor_fct = anchor_classify, 
#'                      filter_fct = classification_filter, 
#'                      division_fct = relative_gap_division(relative_gap = 0.01)
#' )
#' 
#' pc <- predict(object = dc, newdata = test[, 1:4], 
#'               predict_fct = relative_gap_prediction(relative_gap = 0.01))
#' # plot
#' sc <- subcover(pc, 0.6, "snapshot")
#' col <- predict_coloring(sc)
#' \dontrun{
#' plot(sc, coloring = col)
#' }
#' 
#' predicted <- group_from_predict_cover(pc)
#' observed <- test[, 5]
#' 
#' table(observed, predicted)
#' @rdname predict_cover
group_from_predict_cover <- function(pc, group = NULL){
  `%>%` <- dplyr::`%>%`
  . <- NULL
  # use lowest level
  sc <- subcover(pc, 0, "snapshot")
  # check group
  if (is.null(group)){
    group <- pc@parameters$group
  }
  # cover matrix
  cpm <- cover_predict_matrix(sc)
  # find groups
  cpm_df <- as.data.frame(cpm)
  cpm_df %>% 
    dplyr::group_by_(.dots = colnames(cpm_df)) %>% 
    dplyr::do({
      subsets <- unlist(.[1, ])
      possible_groups <- group[unlist(lapply(sc@subsets[which(subsets == 1)], slot, "indices"))]
      uniqv <- unique(possible_groups)
      data.frame(group = uniqv[which.max(tabulate(match(possible_groups, uniqv)))])
    }) %>% 
    dplyr::right_join(cpm_df) %>%
    dplyr::ungroup() %>% 
    dplyr::select_("group") %>%
    unlist %>%
    unname
}

cover_predict_matrix <- function (cover) 
{
  npoints <- length(unique(unlist(lapply(cover@subsets, slot, 
                                         "predicted"))))
  sapply(cover@subsets, function(x, npoints) {
    vec <- rep(0, npoints)
    vec[x@predicted] <- 1
    vec
  }, npoints = npoints)
}
