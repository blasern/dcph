context("classes")


test_that("empty patch can be created", {
  expect_true(is.patch(new("patch", 
                           birth = 0 
  )))
})

test_that("empty cover can be created", {
  expect_true(is.cover(new("cover", 
                           data = data.frame(x = 1), 
                           subsets = list(new("patch", 
                                              index = 1L,
                                              birth = 0 
                           ))
                           )))
})
