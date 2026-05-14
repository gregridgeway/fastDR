test_that("fastDR can return a propensity-score-only fit", {
  skip_if_not_installed("gbm3")
  skip_if_not_installed("survey")

  dat <- make_fastdr_data()
  fit <- fastDR(
    make_fastdr_forms(),
    dat,
    n.trees = 14,
    interaction.depth = 1,
    shrinkage = 0.05,
    ps.only = TRUE,
    par_details = gbm3::gbmParallel(1, 1024)
  )

  expect_s3_class(fit, "fastDR")
  expect_null(fit$effects)
  expect_named(fit$w, as.character(dat$id))
  expect_named(fit$p, as.character(dat$id))
  expect_true(all(fit$p >= 0 & fit$p <= 1))
  expect_equal(colnames(fit$balance.tab), c("control", "treatment", "KS"))
  expect_true(is.finite(fit$ESS))
})

test_that("fastDR fits a complete gaussian ATT analysis", {
  skip_if_not_installed("gbm3")
  skip_if_not_installed("survey")

  dat <- make_fastdr_data()
  fit <- fastDR(
    make_fastdr_forms(),
    dat,
    n.trees = 14,
    interaction.depth = 1,
    shrinkage = 0.05,
    keepGLM = FALSE,
    par_details = gbm3::gbmParallel(1, 1024)
  )

  expect_s3_class(fit, "fastDR")
  expect_named(fit$effects, "y")
  expect_equal(rownames(fit$effects$y), c("un", "ps", "dr"))
  expect_equal(colnames(fit$effects$y), c("E.y1", "E.y0", "se.y1", "se.y0", "TE", "se.TE", "p"))
  expect_equal(fit$estimand, "ATT")
  expect_equal(fit$y.dist, "gaussian")
  expect_false("glm.dr" %in% names(fit))
  expect_true(all(is.finite(as.matrix(fit$effects$y))))
  expect_true(all(fit$effects$y$p >= 0 & fit$effects$y$p <= 1))
  expect_lt(abs(fit$effects$y["dr", "TE"] - 0.45), 0.15)
})

test_that("fastDR supports ATE and weighted inputs", {
  skip_if_not_installed("gbm3")
  skip_if_not_installed("survey")

  dat <- make_fastdr_data()
  fit <- suppressWarnings(
    fastDR(
      make_fastdr_forms(weights = TRUE),
      dat,
      estimand = "ATE",
      n.trees = 14,
      interaction.depth = 1,
      shrinkage = 0.05,
      keepGLM = FALSE,
      par_details = gbm3::gbmParallel(1, 1024)
    )
  )

  expect_equal(fit$estimand, "ATE")
  expect_s3_class(fit, "fastDR")
  expect_true(all(fit$w > 0))
  expect_true(is.finite(fit$effects$y["dr", "TE"]))
})
