test_that("sampled DR covariance decomposition matches exact covariance with all cases", {
  skip_if_not_installed("survey")

  set.seed(41)
  n <- 40
  dat <- data.frame(
    y = rnorm(n),
    treat = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    x2 = factor(sample(letters[1:3], n, replace = TRUE)),
    w = runif(n, 0.5, 2)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = dat)
  fit <- survey::svyglm(y ~ treat + x1 + x2, design = des, rescale = FALSE)

  pred_data <- dat
  pred_rows <- rbind(pred_data, pred_data)
  pred_rows[1:n, "treat"] <- 0
  pred_rows[(n + 1):(2 * n), "treat"] <- 1

  pred_exact <- predict(fit, newdata = pred_rows, type = "response", vcov = TRUE)
  u <- cbind(
    rep(1:0, each = n),
    rep(0:1, each = n),
    rep(c(-1, 1), each = n)
  ) / n
  exact <- sqrt(diag(t(u) %*% vcov(pred_exact) %*% u))

  pred_diag <- predict(fit, newdata = pred_rows, type = "response")
  Vdiag <- as.numeric(vcov(pred_diag))
  VdiagE0 <- sum(Vdiag[1:n])
  VdiagE1 <- sum(Vdiag[(n + 1):(2 * n)])
  VdiagY0Y1 <- .fastDR_paired_prediction_cov_sum(fit, pred_data, "treat")
  VdiagTE <- VdiagE0 + VdiagE1 - 2 * VdiagY0Y1

  V <- vcov(pred_exact)
  V00 <- V[1:n, 1:n]
  V11 <- V[(n + 1):(2 * n), (n + 1):(2 * n)]
  V01 <- V[1:n, (n + 1):(2 * n)]
  VoffdiagE0 <- sum(V00) - sum(diag(V00))
  VoffdiagE1 <- sum(V11) - sum(diag(V11))
  VoffdiagY0Y1 <- sum(V01) - sum(diag(V01))
  approx <- sqrt(c(
    VdiagE0 + VoffdiagE0,
    VdiagE1 + VoffdiagE1,
    VdiagTE + VoffdiagE0 + VoffdiagE1 - 2 * VoffdiagY0Y1
  )) / n

  expect_equal(approx, exact, tolerance = 1e-10)
})

test_that("sampled DR covariance approximation is close to exact covariance", {
  skip_if_not_installed("survey")
  skip_if_not(
    identical(Sys.getenv("FASTDR_HEAVY_TESTS"), "true"),
    "Set FASTDR_HEAVY_TESTS=true to run the n=5000 exact covariance check."
  )

  set.seed(20260514)
  n <- 5000
  n0 <- 2500
  x1 <- rnorm(n)
  x2 <- runif(n, -1, 1)
  x3 <- factor(sample(letters[1:4], n, replace = TRUE))
  treat <- rbinom(n, 1, plogis(-0.2 + 0.3 * x1 - 0.2 * x2))
  y <- 1 + 0.7 * treat + 0.5 * x1 - 0.25 * x2 + rnorm(n)
  pred_data <- data.frame(
    y = y,
    treat = treat,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    w = runif(n, 0.5, 2)
  )
  fit <- survey::svyglm(
    y ~ treat + x1 + x2 + x3,
    design = survey::svydesign(ids = ~1, weights = ~w, data = pred_data),
    rescale = FALSE
  )

  pred_rows <- rbind(pred_data, pred_data)
  pred_rows[1:n, "treat"] <- 0
  pred_rows[(n + 1):(2 * n), "treat"] <- 1
  pred_exact <- predict(fit, newdata = pred_rows, type = "response", vcov = TRUE)
  u <- cbind(
    rep(1:0, each = n),
    rep(0:1, each = n),
    rep(c(-1, 1), each = n)
  ) / n
  exact <- sqrt(diag(t(u) %*% vcov(pred_exact) %*% u))

  pred_diag <- predict(fit, newdata = pred_rows, type = "response")
  Vdiag <- as.numeric(vcov(pred_diag))
  VdiagE0 <- sum(Vdiag[1:n])
  VdiagE1 <- sum(Vdiag[(n + 1):(2 * n)])
  VdiagY0Y1 <- .fastDR_paired_prediction_cov_sum(fit, pred_data, "treat")
  VdiagTE <- VdiagE0 + VdiagE1 - 2 * VdiagY0Y1

  set.seed(99)
  sampled_data <- pred_data[sample.int(n, n0), ]
  sampled_rows <- rbind(sampled_data, sampled_data)
  sampled_rows[1:n0, "treat"] <- 0
  sampled_rows[(n0 + 1):(2 * n0), "treat"] <- 1
  pred_sample <- predict(
    fit,
    newdata = sampled_rows,
    type = "response",
    vcov = TRUE
  )
  V <- vcov(pred_sample)
  V00 <- V[1:n0, 1:n0]
  V11 <- V[(n0 + 1):(2 * n0), (n0 + 1):(2 * n0)]
  V01 <- V[1:n0, (n0 + 1):(2 * n0)]
  VoffdiagE0 <- sum(V00) - sum(diag(V00))
  VoffdiagE1 <- sum(V11) - sum(diag(V11))
  VoffdiagY0Y1 <- sum(V01) - sum(diag(V01))
  VoffdiagE0 <- n * (n - 1) * VoffdiagE0 / (n0 * (n0 - 1))
  VoffdiagE1 <- n * (n - 1) * VoffdiagE1 / (n0 * (n0 - 1))
  VoffdiagY0Y1 <- n * (n - 1) * VoffdiagY0Y1 / (n0 * (n0 - 1))
  approx <- sqrt(c(
    VdiagE0 + VoffdiagE0,
    VdiagE1 + VoffdiagE1,
    VdiagTE + VoffdiagE0 + VoffdiagE1 - 2 * VoffdiagY0Y1
  )) / n

  expect_lt(max(abs((approx - exact) / exact)), 0.01)
})
