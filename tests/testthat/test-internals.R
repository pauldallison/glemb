test_that(".make_margins generates correct vectors", {

  # p=1: just the single variable
  expect_equal(.make_margins(1L, 2L), 1L)

  # p=2, 2-way: single pair, no separator needed
  expect_equal(.make_margins(2L, 2L), c(1L, 2L))

  # p=3, 2-way: all three pairs separated by zeros
  expect_equal(.make_margins(3L, 2L), c(1L, 2L, 0L, 1L, 3L, 0L, 2L, 3L))

  # p=4, 2-way: all six pairs
  expect_equal(.make_margins(4L, 2L),
               c(1L, 2L, 0L, 1L, 3L, 0L, 1L, 4L, 0L,
                 2L, 3L, 0L, 2L, 4L, 0L, 3L, 4L))

  # p=3, 3-way: single triple (implies all lower-order terms in mix)
  expect_equal(.make_margins(3L, 3L), c(1L, 2L, 3L))

  # p=2, cat.interact=3: falls back to 2-way (can't exceed p)
  expect_equal(.make_margins(2L, 3L), c(1L, 2L))
})

test_that(".make_design returns correct identity matrix", {
  fm <- list(
    a = list(name = "a", levels = c("x","y"),   nlevels = 2L),
    b = list(name = "b", levels = c("p","q","r"), nlevels = 3L)
  )
  D <- .make_design(fm)
  expect_true(is.matrix(D))
  expect_equal(dim(D), c(6L, 6L))
  expect_equal(D, diag(6L))
})

test_that(".check_inputs catches bad arguments", {
  df <- data.frame(x = factor(c("a","b","a")), y = rnorm(3))

  # Not a data frame
  expect_error(glemb(list()), "'data' must be a data frame")

  # Bad m
  expect_error(glemb(df, m = 0),    "'m' must be a positive integer")
  expect_error(glemb(df, m = 1.5),  "'m' must be a positive integer")

  # Bad cat.interact
  expect_error(glemb(df, cat.interact = 4), "'cat.interact' must be 2 or 3")

  # Bad idvars
  expect_error(glemb(df, idvars = "z"), "not columns in 'data'")

  # Bad cat.prior
  expect_error(glemb(df, cat.prior = -1), "'cat.prior' must be NULL")

  # Bad seed
  expect_error(glemb(df, seed = -5),  "'seed' must be a positive integer")
  expect_error(glemb(df, seed = 3e9), "'seed' must be a positive integer")

  # Bad output
  expect_error(glemb(df, output = "csv"), "'output' must be")
})

test_that(".preprocess_data classifies and recodes variables correctly", {
  df <- data.frame(
    id  = 1:6,
    cat = factor(c("a","b","a","b","a","b")),
    grp = c("x","y","x","y","x","y"),   # character -> factor
    val = c(1.1, 2.2, 3.3, NA, 5.5, 6.6)
  )

  expect_warning(
    prep <- .preprocess_data(df, idvars = "id"),
    "coerced to factor"
  )

  # Categorical variables identified correctly
  expect_equal(prep$cat_names, c("cat", "grp"))
  expect_equal(prep$cont_names, "val")
  expect_equal(prep$p, 2L)

  # Factors recoded to integers
  expect_true(is.integer(prep$data_model[["cat"]]))
  expect_true(is.integer(prep$data_model[["grp"]]))

  # Categorical columns come first
  expect_equal(names(prep$data_model)[1:2], c("cat", "grp"))

  # ID variable isolated correctly
  expect_equal(names(prep$data_ids), "id")

  # Factor metadata preserved
  expect_equal(prep$factor_meta$cat$levels, c("a", "b"))
})

test_that(".preprocess_data errors on too many factor levels", {
  df <- data.frame(
    x = factor(letters),     # 26 levels > 20
    y = rnorm(26)
  )
  expect_error(.preprocess_data(df, idvars = NULL), "26 levels")
})

test_that(".restore_factors recovers original factor levels", {
  factor_meta <- list(
    cat = list(name = "cat", levels = c("a","b","c"), nlevels = 3L)
  )
  df <- data.frame(cat = c(1L, 3L, 2L, 1L), val = 1:4)
  out <- .restore_factors(df, factor_meta)

  expect_s3_class(out$cat, "factor")
  expect_equal(levels(out$cat), c("a","b","c"))
  expect_equal(as.character(out$cat), c("a","c","b","a"))
})
