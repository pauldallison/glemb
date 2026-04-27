# Internal preprocessing functions for glemb

# --------------------------------------------------------------------------- #
# Non-standard evaluation helper for noms
# --------------------------------------------------------------------------- #

# Parse an unevaluated noms expression (from substitute()) into a character
# vector of variable names.  Called only when evaluating noms directly fails
# (i.e. the user wrote unquoted names like c(pov, race, gender)).
#
.parse_noms_expr <- function(expr) {
  if (is.symbol(expr))
    return(as.character(expr))
  if (is.call(expr) && identical(as.character(expr[[1]]), "c")) {
    parts <- as.list(expr)[-1L]
    return(vapply(parts, function(p) {
      if (is.symbol(p))    return(as.character(p))
      if (is.character(p)) return(p)
      stop("'noms' contains an expression that is neither a name nor a ",
           "quoted string.", call. = FALSE)
    }, character(1L)))
  }
  stop("Invalid 'noms' specification. Use noms = c(var1, var2) or ",
       "noms = c(\"var1\", \"var2\").", call. = FALSE)
}

# --------------------------------------------------------------------------- #
# Input validation
# --------------------------------------------------------------------------- #

.check_inputs <- function(data, m, noms, idvars, cat.interact, cat.prior,
                          empri, maxits, seed, output, p2s) {

  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)

  if (!is.numeric(m) || length(m) != 1L || m < 1L || m != floor(m))
    stop("'m' must be a positive integer.", call. = FALSE)

  # --- noms -------------------------------------------------------------------
  if (is.null(noms) || length(noms) == 0L)
    stop("'noms' must specify at least one categorical variable.",
         call. = FALSE)
  if (!is.character(noms))
    stop("'noms' must be a character vector of column names.", call. = FALSE)
  bad_noms <- setdiff(noms, names(data))
  if (length(bad_noms) > 0L)
    stop("The following 'noms' variables are not columns in 'data': ",
         paste(bad_noms, collapse = ", "), call. = FALSE)

  # --- idvars -----------------------------------------------------------------
  if (!is.null(idvars)) {
    if (!is.character(idvars))
      stop("'idvars' must be a character vector of column names.",
           call. = FALSE)
    bad_id <- setdiff(idvars, names(data))
    if (length(bad_id) > 0L)
      stop("The following 'idvars' are not columns in 'data': ",
           paste(bad_id, collapse = ", "), call. = FALSE)
    overlap <- intersect(noms, idvars)
    if (length(overlap) > 0L)
      stop("The following variables appear in both 'noms' and 'idvars': ",
           paste(overlap, collapse = ", "), call. = FALSE)
  }

  if (!cat.interact %in% c(2L, 3L))
    stop("'cat.interact' must be 2 or 3.", call. = FALSE)

  if (!is.null(cat.prior)) {
    if (!is.numeric(cat.prior) || length(cat.prior) != 1L || cat.prior < 0)
      stop("'cat.prior' must be NULL, 0, or a positive number.",
           call. = FALSE)
  }

  if (!is.null(empri)) {
    if (!is.numeric(empri) || length(empri) != 1L || empri < 0)
      stop("'empri' must be NULL, 0, or a positive number.", call. = FALSE)
  }

  if (!is.numeric(maxits) || length(maxits) != 1L ||
      maxits < 1L || maxits != floor(maxits))
    stop("'maxits' must be a positive integer.", call. = FALSE)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || seed <= 0 ||
        seed != floor(seed) || seed > 2e9)
      stop("'seed' must be a positive integer no greater than 2,000,000,000.",
           call. = FALSE)
  }

  if (!output %in% c("list", "mids", "mitml"))
    stop("'output' must be \"list\", \"mids\", or \"mitml\".", call. = FALSE)

  if (output == "mids" && !requireNamespace("mice", quietly = TRUE))
    stop("output = \"mids\" requires the 'mice' package. ",
         "Install it with: install.packages(\"mice\")", call. = FALSE)

  if (output == "mitml" && !requireNamespace("mitml", quietly = TRUE))
    stop("output = \"mitml\" requires the 'mitml' package. ",
         "Install it with: install.packages(\"mitml\")", call. = FALSE)

  if (!p2s %in% c(0L, 1L))
    stop("'p2s' must be 0 or 1.", call. = FALSE)

  invisible(NULL)
}

# --------------------------------------------------------------------------- #
# Main preprocessing
# --------------------------------------------------------------------------- #

#' @keywords internal
.preprocess_data <- function(data, noms, idvars) {

  # --- Separate ID variables --------------------------------------------------
  if (!is.null(idvars)) {
    data_ids  <- data[, idvars,                       drop = FALSE]
    data_work <- data[, setdiff(names(data), idvars), drop = FALSE]
  } else {
    data_ids  <- NULL
    data_work <- data
  }

  # --- Check that there is something to impute --------------------------------
  if (!anyNA(data_work))
    stop("No missing values found in 'data'. glemb() is only needed when ",
         "data contain missing values.", call. = FALSE)

  # --- Coerce haven_labelled columns to numeric --------------------------------
  # haven_labelled is a common class from haven::read_sav() / read_dta().
  # Strip the labels so downstream type checks work correctly.
  haven_cols <- names(data_work)[
    vapply(data_work, function(x) inherits(x, "haven_labelled"), logical(1L))
  ]
  if (length(haven_cols) > 0L) {
    data_work[haven_cols] <- lapply(data_work[haven_cols],
                                    function(x) as.numeric(unclass(x)))
  }

  # --- Reject character variables ---------------------------------------------
  char_cols <- names(data_work)[vapply(data_work, is.character, logical(1L))]
  if (length(char_cols) > 0L)
    stop("The following variable(s) are character type: ",
         paste(char_cols, collapse = ", "), ". Convert to numeric before ",
         "calling glemb(), e.g. use model.matrix() or as.numeric().",
         call. = FALSE)

  # --- Reject factor variables not listed in noms ----------------------------
  factor_cols <- names(data_work)[vapply(data_work, is.factor, logical(1L))]
  unexpected  <- setdiff(factor_cols, noms)
  if (length(unexpected) > 0L)
    stop("The following factor variable(s) are not listed in 'noms': ",
         paste(unexpected, collapse = ", "), ". Either add them to 'noms' ",
         "or convert to numeric.", call. = FALSE)

  # --- Classify variables -----------------------------------------------------
  cat_names  <- noms
  cont_names <- setdiff(names(data_work), noms)

  if (length(cont_names) == 0L)
    stop("No continuous variables found. glemb() requires at least one ",
         "variable not listed in 'noms'.", call. = FALSE)

  # --- Validate continuous variables ------------------------------------------
  for (v in cont_names) {
    col  <- data_work[[v]]
    vals <- col[!is.na(col)]

    if (any(is.infinite(col)))
      stop("Continuous variable '", v, "' contains infinite values ",
           "(Inf or -Inf). Remove or replace them before calling glemb().",
           call. = FALSE)

    if (length(vals) == 0L)
      warning("Continuous variable '", v, "' is entirely missing. Imputed ",
              "values will be drawn from an unanchored distribution.",
              call. = FALSE)
    else if (length(unique(vals)) == 1L)
      stop("Continuous variable '", v, "' has no variance (all non-missing ",
           "values are identical). Remove it or recode it before calling ",
           "glemb().", call. = FALSE)
  }

  # --- Validate noms variables ------------------------------------------------
  factor_meta <- vector("list", length(cat_names))
  names(factor_meta) <- cat_names

  for (v in cat_names) {
    col <- data_work[[v]]

    if (is.factor(col)) {
      orig_type   <- "factor"
      orig_values <- levels(col)           # character vector of level labels
    } else {
      orig_type   <- "numeric"
      orig_values <- sort(unique(col[!is.na(col)]))   # sorted numeric values
    }

    nl <- length(orig_values)
    if (nl > 20L)
      stop("Variable '", v, "' has ", nl, " unique values. glemb() allows ",
           "at most 20 categories per noms variable.", call. = FALSE)
    if (nl < 2L)
      stop("Variable '", v, "' has fewer than 2 unique non-missing values.",
           call. = FALSE)

    factor_meta[[v]] <- list(
      name        = v,
      orig_type   = orig_type,
      orig_values = orig_values,
      nlevels     = nl
    )
  }

  # --- Save original column order (for output reconstruction) ----------------
  col_order <- names(data_work)

  # --- Reorder: categorical first, then continuous (required by mix) ----------
  data_work <- data_work[, c(cat_names, cont_names), drop = FALSE]

  # --- Recode noms variables to consecutive integers 1:k (double) -------------
  # mix::prelim.mix requires categorical variables coded 1, 2, ..., k.
  # We recode by mapping each original value to its rank in orig_values.
  for (v in cat_names) {
    col  <- data_work[[v]]
    meta <- factor_meta[[v]]
    if (meta$orig_type == "factor") {
      codes <- as.numeric(as.integer(col))          # factor → 1:k directly
    } else {
      # numeric: find position of each value in sorted unique values
      codes <- match(col, meta$orig_values)          # NA stays NA
      codes <- as.numeric(codes)
    }
    data_work[[v]] <- codes
  }

  # --- Coerce continuous columns to double ------------------------------------
  # Ensure no integer-typed columns reach prelim.mix; Fortran expects double.
  data_work[cont_names] <- lapply(data_work[cont_names], as.numeric)

  list(
    data_model  = data_work,
    data_ids    = data_ids,
    factor_meta = factor_meta,
    cat_names   = cat_names,
    cont_names  = cont_names,
    col_order   = col_order,
    p           = length(cat_names)
  )
}
