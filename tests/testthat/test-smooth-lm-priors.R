test_that("smooth.lm prior rows are omitted when penalty is zero", {
  mx <- matrix(
    c(1, 1, 0, 0, 0.5),
    nrow = 1,
    dimnames = list(NULL, c("Intercept", "w", "y", "treat", "x1"))
  )

  expect_identical(.add.smooth.lm.prior.rows(mx, "gaussian", 0), mx)
})

test_that("gaussian smooth.lm prior rows penalize covariates only", {
  mx <- matrix(
    c(1, 1, 0, 0, 0.5, -0.5),
    nrow = 1,
    dimnames = list(NULL, c("Intercept", "w", "y", "treat", "x1", "x2"))
  )

  aug <- .add.smooth.lm.prior.rows(mx, "gaussian", 2)

  expect_equal(nrow(aug), 3)
  expect_equal(aug[2:3, "Intercept"], c(0, 0))
  expect_equal(aug[2:3, "w"], c(2, 2))
  expect_equal(aug[2:3, "y"], c(0, 0))
  expect_equal(aug[2:3, "treat"], c(0, 0))
  expect_equal(aug[2:3, c("x1", "x2")], diag(2), ignore_attr = TRUE)
})

test_that("binomial smooth.lm prior rows represent log-F(m,m)", {
  mx <- matrix(
    c(1, 1, 0, 0, 0.5, -0.5),
    nrow = 1,
    dimnames = list(NULL, c("Intercept", "w", "y", "treat", "x1", "x2"))
  )

  aug <- .add.smooth.lm.prior.rows(mx, "quasibinomial", 2)

  expect_equal(nrow(aug), 5)
  expect_equal(aug[2:5, "w"], rep(1, 4))
  expect_equal(aug[2:5, "y"], c(0, 0, 1, 1))
  expect_equal(aug[2:5, "treat"], rep(0, 4))
  expect_equal(aug[2:3, c("x1", "x2")], diag(2), ignore_attr = TRUE)
  expect_equal(aug[4:5, c("x1", "x2")], diag(2), ignore_attr = TRUE)
})

test_that("poisson smooth.lm prior rows represent log-gamma(m,m)", {
  mx <- matrix(
    c(1, 1, 0, 0, 0.5, -0.5),
    nrow = 1,
    dimnames = list(NULL, c("Intercept", "w", "y", "treat", "x1", "x2"))
  )

  aug <- .add.smooth.lm.prior.rows(mx, "quasipoisson", 2)

  expect_equal(nrow(aug), 3)
  expect_equal(aug[2:3, "w"], c(2, 2))
  expect_equal(aug[2:3, "y"], c(1, 1))
  expect_equal(aug[2:3, "treat"], c(0, 0))
  expect_equal(aug[2:3, c("x1", "x2")], diag(2), ignore_attr = TRUE)
})
