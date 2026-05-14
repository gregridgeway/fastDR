test_that("fastDR requires all formula components", {
  dat <- make_fastdr_data()
  forms <- make_fastdr_forms()

  expect_error(fastDR(forms[names(forms) != "y.form"], dat), "missing y.form")
  expect_error(fastDR(forms[names(forms) != "t.form"], dat), "missing t.form")
  expect_error(fastDR(forms[names(forms) != "x.form"], dat), "missing x.form")
  expect_error(fastDR(forms[names(forms) != "key.form"], dat), "missing key.form")
})

test_that("fastDR rejects malformed formulas", {
  dat <- make_fastdr_data()
  forms <- make_fastdr_forms()

  forms$y.form <- ~.
  expect_error(fastDR(forms, dat), "not allowed in y.form", fixed = TRUE)

  forms <- make_fastdr_forms()
  forms$t.form <- ~treat + x1
  expect_error(fastDR(forms, dat), "only one treatment indicator")

  forms <- make_fastdr_forms()
  forms$key.form <- ~id + x1
  expect_error(fastDR(forms, dat), "key.form should have only one term")

  forms <- make_fastdr_forms(weights = TRUE)
  forms$weights.form <- ~sw + x1
  expect_error(fastDR(forms, dat), "weights.form should have only one term")
})

test_that("fastDR rejects invalid distributions and estimands", {
  dat <- make_fastdr_data()
  forms <- make_fastdr_forms()

  expect_error(fastDR(forms, dat, y.dist = gaussian()), "should be given as a character string")
  expect_error(fastDR(forms, dat, y.dist = c("gaussian", "quasipoisson")), "Length of y.dist")
  expect_error(fastDR(forms, dat, estimand = "ATC"), "estimand must either")
})

test_that("fastDR protects reserved data column names", {
  dat <- make_fastdr_data()
  forms <- make_fastdr_forms()

  dat$w <- 1
  expect_error(fastDR(forms, dat), "variable names 'w'")

  dat <- make_fastdr_data()
  dat$samp.w <- 1
  expect_error(fastDR(forms, dat), "variable names 'samp.w'")
})

test_that("fastDR rejects overlapping outcome, treatment, and covariates", {
  dat <- make_fastdr_data()

  forms <- make_fastdr_forms()
  forms$y.form <- ~y + treat
  expect_error(fastDR(forms, dat), "Treatment indicator in t.form")

  forms <- make_fastdr_forms()
  forms$x.form <- ~x1 + treat
  expect_error(fastDR(forms, dat), "should not be included in x.form")

  forms <- make_fastdr_forms()
  forms$x.form <- ~x1 + y
  expect_error(fastDR(forms, dat), "Outcome in y.form")
})

test_that("fastDR validates keys, treatment indicators, and weights", {
  forms <- make_fastdr_forms()

  dat <- make_fastdr_data()
  dat$id[1] <- NA
  expect_error(fastDR(forms, dat), "missing values in key")

  dat <- make_fastdr_data()
  dat$id[2] <- dat$id[1]
  expect_error(fastDR(forms, dat), "key has duplicates")

  dat <- make_fastdr_data()
  dat$treat[1] <- 2
  expect_error(fastDR(forms, dat), "does not only take values 0 or 1")

  forms <- make_fastdr_forms(weights = TRUE)
  dat <- make_fastdr_data()
  dat$sw[1] <- NA
  expect_warning(
    expect_error(fastDR(forms, dat), "missing values in weights"),
    "sample weighted data"
  )

  dat <- make_fastdr_data()
  dat$sw[1] <- -1
  expect_error(
    expect_warning(fastDR(forms, dat), "sample weighted data"),
    "negative"
  )
})
