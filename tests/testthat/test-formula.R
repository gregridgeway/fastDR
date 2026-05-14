test_that("make.fastDR.formula builds explicit formulas", {
  forms <- make.fastDR.formula(
    y.vars = c("y1", "y2"),
    t.var = "treat",
    x.vars = c("x1", "x2"),
    weights.var = "sw",
    key.var = "id"
  )

  expect_equal(deparse(forms$y.form), "~y1 + y2")
  expect_equal(deparse(forms$t.form), "~treat")
  expect_equal(deparse(forms$x.form), "~x1 + x2")
  expect_equal(deparse(forms$weights.form), "~sw")
  expect_equal(deparse(forms$key.form), "~id")
})

test_that("make.fastDR.formula expands dot covariates from data", {
  dat <- data.frame(
    id = 1:3,
    y = 1:3,
    y2 = 3:1,
    treat = c(0, 1, 0),
    sw = c(1, 2, 1),
    x1 = 4:6,
    x2 = letters[1:3]
  )

  forms <- make.fastDR.formula(
    y.vars = c("y", "y2"),
    t.var = "treat",
    x.vars = ".",
    weights.var = "sw",
    key.var = "id",
    data = dat
  )

  expect_equal(attr(terms(forms$x.form), "term.labels"), c("x1", "x2"))
})

test_that("make.fastDR.formula validates dot and warns on missing variables", {
  dat <- data.frame(id = 1:3, y = 1:3, treat = c(0, 1, 0), x = 4:6)

  expect_error(
    make.fastDR.formula("y", "treat", x.vars = ".", key.var = "id"),
    "data must be given"
  )
  expect_warning(
    make.fastDR.formula("missing_y", "treat", "x", key.var = "id", data = dat),
    "Not all y.vars appear in data"
  )
  expect_warning(
    make.fastDR.formula("y", "missing_t", "x", key.var = "id", data = dat),
    "t.var does not appear in data"
  )
  expect_warning(
    make.fastDR.formula("y", "treat", "x", weights.var = "missing_w", key.var = "id", data = dat),
    "weights.var does not appear in data"
  )
  expect_warning(
    make.fastDR.formula("y", "treat", "x", key.var = "missing_id", data = dat),
    "key.var does not appear in data"
  )
})
