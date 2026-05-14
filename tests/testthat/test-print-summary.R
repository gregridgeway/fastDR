test_that("print.fastDR formats model-specific outcome summaries", {
  x <- make_print_object("gaussian")

  expect_match(capture.output(print(x, model = "un")), "Unadjusted results", all = FALSE)
  expect_match(capture.output(print(x, model = "ps")), "Propensity score weighting results", all = FALSE)
  expect_match(capture.output(print(x, model = "dr")), "Doubly robust results", all = FALSE)
  expect_match(capture.output(print(x, model = "dr")), "effect \\(95% CI\\): 0.4", all = FALSE)
})

test_that("print.fastDR formats binary and count outcome scales", {
  expect_match(
    capture.output(print(make_print_object("quasibinomial"))),
    "percentage point difference",
    all = FALSE
  )
  expect_match(
    capture.output(print(make_print_object("quasipoisson"))),
    "rate difference",
    all = FALSE
  )
})

test_that("print.fastDR validates object, type, model, and missing effects", {
  expect_error(fastDR:::print.fastDR(list()), "object must be a fastDR object")
  expect_error(print(make_print_object(), type = "bad"), "type parameter")
  expect_error(print(make_print_object(), model = "bad"), "model parameter")
  expect_error(print(make_print_object(effects = FALSE)), "no estimated effects")
})

test_that("summary.fastDR prints effects and balance table", {
  output <- capture.output(summary(make_print_object()))

  expect_match(output, "Results", all = FALSE)
  expect_match(output, "Estimand:  ATT", all = FALSE)
  expect_match(output, "Balance table", all = FALSE)
  expect_match(output, "x1", all = FALSE)
})

test_that("summary.fastDR validates object and balance table", {
  expect_error(fastDR:::summary.fastDR(list()), "object must be a fastDR object")

  x <- make_print_object()
  x$balance.tab <- NULL
  expect_error(summary(x), "missing the balance table")
})
