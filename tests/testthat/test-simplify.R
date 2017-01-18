context("simplify")

test_that("simplification method works as expected", {
  x <- cover(dist(1:3), 
             lapply(list(1L, 1L, 2L, 3L, c(1L, 2L), c(2L, 1L), c(1L, 2L, 3L)), patch))
  no_duplicates <- cover(dist(1:3), 
                         lapply(list(1L, 2L, 3L, c(1L, 2L), c(1L, 2L, 3L)), patch))
  no_subsets <- cover(dist(1:3), 
                      lapply(list(c(1L, 2L, 3L)), patch))
  
  expect_equal(simplify_cover(x, method = "none"), simplify_cover(x))
  expect_equal(simplify_cover(x, method = "duplicates"), no_duplicates)
  expect_equal(simplify_cover(x, method = "subsets"), no_subsets)  
})
