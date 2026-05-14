make_fastdr_data <- function(n = 80) {
  set.seed(101)
  x1 <- rnorm(n)
  x2 <- factor(rep(c("a", "b"), length.out = n))
  treat <- rep(0:1, each = n / 2)

  data.frame(
    id = seq_len(n),
    y = 1 + 0.45 * treat + 0.25 * x1 + rnorm(n, sd = 0.2),
    y_bin = as.numeric(1 + 0.6 * treat + 0.2 * x1 + rnorm(n) > 1.2),
    treat = treat,
    x1 = x1,
    x2 = x2,
    sw = seq(1, 2, length.out = n)
  )
}

make_fastdr_forms <- function(weights = FALSE) {
  list(
    y.form = ~y,
    t.form = ~treat,
    x.form = ~x1 + x2,
    weights.form = if (weights) ~sw else NULL,
    key.form = ~id
  )
}

make_print_object <- function(y.dist = "gaussian", effects = TRUE) {
  x <- list(
    estimand = "ATT",
    y.dist = y.dist,
    balance.tab = matrix(
      c(0.2, 0.3, 0.1),
      nrow = 1,
      dimnames = list("x1", c("control", "treatment", "KS"))
    )
  )
  if (effects) {
    x$effects <- list(
      y = data.frame(
        E.y1 = c(1.2, 1.2, 1.2),
        E.y0 = c(1.0, 0.9, 0.8),
        se.y1 = c(0.1, 0.1, 0.1),
        se.y0 = c(0.1, 0.1, 0.1),
        TE = c(0.2, 0.3, 0.4),
        se.TE = c(0.05, 0.05, 0.05),
        p = c(0.04, 0.02, 0.01),
        row.names = c("un", "ps", "dr")
      )
    )
  }
  class(x) <- "fastDR"
  x
}
