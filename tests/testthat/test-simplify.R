context("simplify")

test_that("simplification method works as expected", {
  x <- as.cover(list(1L, 1L, 2L, 3L, c(1L, 2L), c(2L, 1L), c(1L, 2L, 3L)))
  no_duplicates <- as.cover(list(1L, 2L, 3L, c(1L, 2L), c(1L, 2L, 3L)))
  no_subsets <- as.cover(list(c(1L, 2L, 3L)))
  
  expect_equal(simplify_cover(x, method = "none"), simplify_cover(x))
  expect_equal(simplify_cover(x, method = "duplicates"), no_duplicates)
  expect_equal(simplify_cover(x, method = "subsets"), no_subsets)  
})
