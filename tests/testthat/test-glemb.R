# Integration tests for glemb()
# These tests require the 'mix' package and use a small synthetic dataset.

make_test_data <- function(n = 100, seed = 42) {
  set.seed(seed)
  df <- data.frame(
    cat1 = factor(sample(c("a","b"), n, replace = TRUE)),
    cat2 = factor(sample(c("x","y","z"), n, replace = TRUE)),
    cont1 = rnorm(n),
    cont2 = rnorm(n, mean = 5)
  )
  # Introduce ~15% missingness in cont1 and cont2
  df$cont1[sample(n, floor(n * 0.15))] <- NA
  df$cont2[sample(n, floor(n * 0.15))] <- NA
  df
}

test_that("glemb() returns a glemb object with correct structure", {
  skip_if_not_installed("mix")

  df  <- make_test_data()
  imp <- glemb(df, m = 3, seed = 1, p2s = 0)

  expect_s3_class(imp, "glemb")
  expect_equal(imp$m, 3L)
  expect_length(imp$imputations, 3L)
  expect_equal(imp$categorical, c("cat1", "cat2"))
  expect_equal(imp$continuous,  c("cont1", "cont2"))
  expect_equal(imp$cat.interact, 2L)
  expect_equal(imp$seed, 1L)
})

test_that("glemb() produces complete data (no missing values)", {
  skip_if_not_installed("mix")

  df  <- make_test_data()
  imp <- glemb(df, m = 3, seed = 2, p2s = 0)

  for (b in seq_len(imp$m)) {
    expect_false(anyNA(imp$imputations[[b]]),
                 label = paste("imputation", b, "has NA values"))
  }
})

test_that("glemb() imputed datasets have correct dimensions", {
  skip_if_not_installed("mix")

  df  <- make_test_data(n = 80)
  imp <- glemb(df, m = 4, seed = 3, p2s = 0)

  for (b in seq_len(imp$m)) {
    expect_equal(nrow(imp$imputations[[b]]), nrow(df))
    expect_equal(ncol(imp$imputations[[b]]), ncol(df))
    expect_equal(names(imp$imputations[[b]]), names(df))
  }
})

test_that("glemb() imputed factors have correct levels", {
  skip_if_not_installed("mix")

  df  <- make_test_data()
  imp <- glemb(df, m = 3, seed = 4, p2s = 0)

  for (b in seq_len(imp$m)) {
    expect_equal(levels(imp$imputations[[b]]$cat1), c("a", "b"))
    expect_equal(levels(imp$imputations[[b]]$cat2), c("x", "y", "z"))
  }
})

test_that("glemb() is reproducible with a seed", {
  skip_if_not_installed("mix")

  df   <- make_test_data()
  imp1 <- glemb(df, m = 3, seed = 99, p2s = 0)
  imp2 <- glemb(df, m = 3, seed = 99, p2s = 0)

  for (b in 1:3) {
    expect_equal(imp1$imputations[[b]], imp2$imputations[[b]])
  }
})

test_that("glemb() idvars are passed through unchanged", {
  skip_if_not_installed("mix")

  df     <- make_test_data()
  df$id  <- seq_len(nrow(df))
  imp    <- glemb(df, m = 3, idvars = "id", seed = 5, p2s = 0)

  # id column present and values are integers
  for (b in seq_len(imp$m)) {
    expect_true("id" %in% names(imp$imputations[[b]]))
  }
  # id not in categorical or continuous
  expect_false("id" %in% imp$categorical)
  expect_false("id" %in% imp$continuous)
})

test_that("glemb() cat.interact = 3 runs without error", {
  skip_if_not_installed("mix")

  df  <- make_test_data()
  expect_no_error(glemb(df, m = 2, cat.interact = 3, seed = 6, p2s = 0))
})

test_that("print.glemb and summary.glemb produce output", {
  skip_if_not_installed("mix")

  df  <- make_test_data()
  imp <- glemb(df, m = 2, seed = 7, p2s = 0)

  expect_output(print(imp),   "glemb")
  expect_output(summary(imp), "glemb")
})

test_that("output = 'mids' returns a mids object", {
  skip_if_not_installed("mix")
  skip_if_not_installed("mice")

  df  <- make_test_data()
  imp <- glemb(df, m = 3, seed = 8, output = "mids", p2s = 0)

  expect_s3_class(imp, "mids")
  expect_equal(imp$m, 3L)
})

test_that("output = 'mids' can be used directly with mice::pool()", {
  skip_if_not_installed("mix")
  skip_if_not_installed("mice")

  df  <- make_test_data()
  imp <- glemb(df, m = 3, seed = 9, output = "mids", p2s = 0)

  fit    <- with(imp, lm(cont1 ~ cat1 + cat2 + cont2))
  pooled <- mice::pool(fit)

  expect_s3_class(pooled, "mipo")
  # intercept + cat1 (1 contrast) + cat2 (2 contrasts) + cont2 (1) = 5 terms
  expect_equal(nrow(summary(pooled)), 5L)
})
