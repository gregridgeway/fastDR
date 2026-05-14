test_that("ks is zero for identical weighted distributions", {
  expect_equal(
    ks(x = c(1, 2, 1, 2), z = c(0, 0, 1, 1), w = rep(1, 4)),
    0
  )
})

test_that("ks detects separated distributions", {
  expect_equal(
    ks(x = c(1, 2, 3, 4), z = c(0, 0, 1, 1), w = rep(1, 4)),
    1
  )
})

test_that("ks is invariant to within-group weight scaling", {
  x <- c(1, 2, 3, 4, 2, 3, 4, 5)
  z <- rep(0:1, each = 4)
  w <- c(1, 2, 1, 2, 3, 1, 2, 1)

  expect_equal(
    ks(x, z, w),
    ks(x, z, ifelse(z == 1, 10 * w, 4 * w))
  )
})

test_that("ks returns zero when all x values are tied", {
  expect_equal(ks(rep(1, 6), rep(0:1, each = 3), rep(1, 6)), 0)
})
