# S3 methods for the "glemb" class

# --------------------------------------------------------------------------- #
# print
# --------------------------------------------------------------------------- #

#' Print a glemb object
#'
#' @param x A `glemb` object returned by [glemb()].
#' @param digits Integer. Number of decimal places for numeric values.
#'   Default is `3`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.glemb <- function(x, digits = 3L, ...) {
  cat("glemb: General Location Model with EMB\n")
  cat(rep("-", 40L), "\n", sep = "")
  cat("Imputations       :", x$m, "\n")
  cat("Noms (categorical):", paste(x$categorical, collapse = ", "), "\n")
  cat("Continuous vars   :", paste(x$continuous,  collapse = ", "), "\n")
  if (length(x$idvars) > 0L)
    cat("ID vars           :", paste(x$idvars, collapse = ", "), "\n")
  cat("Interaction order :", x$cat.interact, "\n")
  cat("cat.prior         :", round(x$cat.prior, digits), "\n")
  cat("empri             :", round(x$empri,     digits), "\n")
  cat("maxits            :", x$maxits, "\n")
  if (!is.null(x$seed))
    cat("Seed              :", x$seed, "\n")
  invisible(x)
}

# --------------------------------------------------------------------------- #
# summary
# --------------------------------------------------------------------------- #

#' Summarise a glemb object
#'
#' Prints a summary of the imputation model and the missingness in the first
#' imputed dataset.
#'
#' @param object A `glemb` object returned by [glemb()].
#' @param digits Integer. Number of decimal places for numeric values.
#'   Default is `3`.
#' @param ... Ignored.
#' @return `object`, invisibly.
#' @export
summary.glemb <- function(object, digits = 3L, ...) {
  cat("glemb: General Location Model with EMB\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat("Imputations:", object$m, "\n\n")

  cat("Variables in imputation model:\n")
  cat("  Noms/categorical (", length(object$categorical), "): ",
      paste(object$categorical, collapse = ", "), "\n", sep = "")
  cat("  Continuous       (", length(object$continuous),  "): ",
      paste(object$continuous,  collapse = ", "), "\n", sep = "")
  if (length(object$idvars) > 0L)
    cat("  ID (excluded): ", paste(object$idvars, collapse = ", "), "\n",
        sep = "")

  cat("\nModel settings:\n")
  cat("  Interaction order :", object$cat.interact, "\n")
  cat("  cat.prior         :", round(object$cat.prior, digits), "\n")
  cat("  empri             :", round(object$empri,     digits), "\n")
  cat("  maxits            :", object$maxits, "\n")
  if (!is.null(object$seed))
    cat("  Seed              :", object$seed, "\n")

  cat("\nImputed dataset dimensions:",
      nrow(object$imputations[[1L]]), "rows x",
      ncol(object$imputations[[1L]]), "columns\n")

  invisible(object)
}

# --------------------------------------------------------------------------- #
# as.mids
# --------------------------------------------------------------------------- #

#' Convert a glemb object to a mids object
#'
#' Converts a `glemb` object to the `mids` format used by the `mice` package,
#' enabling use of [mice::pool()] for results pooling.
#'
#' @param x A `glemb` object returned by [glemb()].
#' @param data The original incomplete data frame that was passed to [glemb()].
#' @param ... Ignored.
#' @return A `mids` object.
#'
#' @details
#' The function stacks the original data (as imputation 0) and all `m`
#' imputed datasets in long format, adds `.imp` and `.id` columns, and
#' calls [mice::as.mids()] to construct the `mids` object.
#'
#' @seealso [mice::as.mids()], [mice::pool()]
#' @export
as.mids.glemb <- function(x, data, ...) {

  if (!requireNamespace("mice", quietly = TRUE))
    stop("Package 'mice' is required for as.mids.glemb(). ",
         "Install it with: install.packages(\"mice\")", call. = FALSE)

  # Build long-format data: imputation 0 = original, 1..m = completed datasets
  long_list <- lapply(seq_len(x$m), function(b) {
    df      <- x$imputations[[b]]
    df$.imp <- b
    df$.id  <- seq_len(nrow(df))
    df
  })

  orig      <- data
  orig$.imp <- 0L
  orig$.id  <- seq_len(nrow(data))

  long <- do.call(rbind, c(list(orig), long_list))

  mice::as.mids(long)
}
