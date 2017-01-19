context("classes")

test_that("can create simple patches", {
  expect_true(is.patch(patch(1L)))
  expect_true(is.patch(patch(c(1L, 2L))))
})

test_that("can create simple covers", {
  expect_true(is.cover(cover(diag(1:3), list(patch(indices = 1:3)))))
})

test_that("cover dimension is correct", {
  expect_equal(cover_dim(cover(diag(1:3), list(patch(indices = 1:3)))), 1)
  expect_equal(cover_dim(cover(diag(1:3), list(patch(indices = 1:3), patch(indices = 1:3)))), 2)
})

test_that("covers are concatinated correctly", {
  c1 <- cover(diag(1:3), list(patch(indices = 1:3)))
  c2 <- cover(diag(1:3), list(patch(indices = 1:3), patch(indices = 1:3)))
  expect_equal(c(c1, c1)@subsets, 
               c2@subsets)
})
