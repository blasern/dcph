context("classes")

test_that("can create simple patches", {
  expect_true(is.patch(patch(1L)))
  expect_true(is.patch(patch(c(1L, 2L))))
})

test_that("can create simple covers", {
  expect_true(is.cover(as.cover(list(1L, 2L, 3L))))
  expect_true(is.cover(as.cover(matrix(c(0, 1, 1, 0, 1, 1), nrow = 2))))
})

test_that("cover dimension is correct", {
  expect_equal(cover_dim(as.cover(list(1L, 2L, 3L))), 1)
  expect_equal(cover_dim(as.cover(list(1L, 1L, 3L))), 2)
  expect_equal(cover_dim(as.cover(list(1L, 1L, 1L))), 3)
})

test_that("covers are concatinated correctly", {
  expect_equal(c(as.cover(list(1L, 2L, 3L)), as.cover(list(c(1L, 2L, 3L)))), 
               as.cover(list(1L, 2L, 3L, c(1L, 2L, 3L))))
})
