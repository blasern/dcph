context("simplify")

test_that("simplification method works as expected", {
  x <- cover(data = matrix(1:3, nrow = 3), 
             subsets = lapply(list(1L, 1L, 2L, 3L, c(1L, 2L), c(2L, 1L), c(1L, 2L, 3L)), patch, id = NA_integer_))
  no_duplicates <- cover(data = matrix(1:3, nrow = 3), 
                         subsets = lapply(list(1L, 2L, 3L, c(1L, 2L), c(1L, 2L, 3L)), patch, id = NA_integer_))
  attr(no_duplicates, "simple") <- c(1L, 3L, 4L, 5L, 7L)
  no_subsets <- cover(data = matrix(1:3, nrow = 3), 
                      subsets = lapply(list(c(1L, 2L, 3L)), patch, id = NA_integer_))
  attr(no_subsets, "simple") <- 7L
  
  expect_equal(simplify_cover(x, method = "none"), simplify_cover(x))
  expect_equal(simplify_cover(x, method = "duplicates"), no_duplicates)
  expect_equal(simplify_cover(x, method = "subsets"), no_subsets)  
})
