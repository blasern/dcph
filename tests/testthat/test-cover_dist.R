context("cover_dist")

test_that("cover distances work", {
  data_matrix <- structure(c(-1.020, -0.559, 0.702, 0.543, 1.032, -0.908, -0.447, -0.158, 
                             -0.984, -0.013, -0.565, 0.88, -0.937, 0.726, -0.538, -0.651, 
                             -0.797, -0.870, 0.443, 0.811, -0.246, -0.927, 1.035, -0.993, 
                             0.100, -0.669, -0.803, 1.143, 0.356, 0.883, -0.940, -0.410, 
                             0.407, -0.529, -0.772, 0.737, -0.034, -0.694, -0.914, 0.377), 
                           .Dim = c(20L,2L), 
                           .Dimnames = list(NULL, c("x", "y")))
  dc <- divisive_cover(data = data_matrix, stop_fct = stop_relative_filter(relative_filter = 0.1))
  cover1 <- subcover(dc, 0.3, "snapshot")
  cover2 <- subcover(dc, 0.8, "snapshot")
  cover3 <- subcover(dc, 0.5, "snapshot")
  
  expect_equal(cover_distance(cover1, cover2, metric = "vi"), as.matrix(dist(c(0, 0.333044760052531))))
  expect_equal(cover_distance(cover1, cover2, metric = "bottleneck", dim = 0), as.matrix(dist(c(0, 5.25019588636283))))
  expect_equal(cover_distance(cover1, cover2, metric = "bottleneck", dim = 1), as.matrix(dist(c(0, 2.07656484777033))))
  expect_equal(cover_distance(cover1, cover2, metric = "interleaving"), as.matrix(dist(c(0, 2.76307092273138))))
  
  # can compare more than two covers  
  expect_equal(dim(cover_distance(cover1, cover2, cover3, metric = "vi")), c(3, 3))
  expect_equal(cover_distance(cover1, cover2, cover3, metric = "vi"), cover_distance(list(cover1, cover2, cover3), metric = "vi"))
  expect_equal(dim(cover_distance(cover1, cover2, cover3, metric = "bottleneck", dim = 0)), c(3, 3))
  expect_equal(dim(cover_distance(cover1, cover2, cover3, metric = "bottleneck", dim = 1)), c(3, 3))
  expect_equal(dim(cover_distance(cover1, cover2, cover3, metric = "interleaving")), c(3, 3))
})
