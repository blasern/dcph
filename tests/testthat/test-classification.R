context("classification")

test_that("classification works", {
  iris <- datasets::iris
  dc <- divisive_cover(data = iris[, 1:4], 
                       group = iris[, 5], 
                       distance_fct = distance_cdist("euclidean"), 
                       stop_fct = stop_relative_filter(relative_filter = 0), 
                       anchor_fct = anchor_classify, 
                       filter_fct = classification_filter, 
                       division_fct = relative_gap_division(relative_gap = 0.01)
                       )
  
  sc <- subcover(dc, method = "external")
  
  expect_true(all(sapply(sc@subsets, function(x){
    length(unique(iris[x@indices, 5]))
  }) == 1))
})

test_that("classification heuristic works", {
  iris <- datasets::iris
  dc <- divisive_cover(data = iris[, 1:4], 
                       group = iris[, 5], 
                       distance_fct = distance_cdist("euclidean"), 
                       stop_fct = stop_relative_filter(relative_filter = 0), 
                       anchor_fct = anchor_heuristic_classify, 
                       filter_fct = classification_filter, 
                       division_fct = relative_gap_division(relative_gap = 0.01)
  )
  
  sc <- subcover(dc, method = "external")
  
  expect_true(all(sapply(sc@subsets, function(x){
    length(unique(iris[x@indices, 5]))
  }) == 1))
})
