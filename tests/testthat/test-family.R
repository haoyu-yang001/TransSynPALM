test_that("families are the connected components of the relatedness graph", {
  grm <- Matrix::sparseMatrix(i = c(1:6, 1, 2, 4, 5), j = c(1:6, 2, 1, 5, 4),
                              x = c(rep(1, 6), .5, .5, .25, .25), dims = c(6, 6))
  f <- family_labels(grm)
  expect_equal(length(unique(f)), 4)
  expect_equal(f[1], f[2])
  expect_equal(f[4], f[5])
  expect_false(f[3] %in% f[c(1, 4, 6)])
})

test_that("a target-only relatedness matrix can join families", {
  grm <- Matrix::Diagonal(4)
  grm_tar <- Matrix::sparseMatrix(i = c(1, 2, 1, 2), j = c(1, 2, 2, 1), x = c(1, 1, .5, .5))
  f <- family_labels(grm, grm_tar = grm_tar, tar_pos = c(2, 4))
  expect_equal(f[2], f[4])
  expect_equal(length(unique(f)), 3)
})
