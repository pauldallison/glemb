# Functions for constructing mix arguments from data metadata

# --------------------------------------------------------------------------- #
# Margins vector
# --------------------------------------------------------------------------- #

# Build the margins vector for mix::ecm.mix.
#
# Encodes all combinations of size `cat.interact` among `p` categorical
# variables as a flat numeric vector with 0 separating each term, following
# Schafer's mix format.  For example:
#   p=3, cat.interact=2  ->  c(1,2, 0, 1,3, 0, 2,3)
#   p=3, cat.interact=3  ->  c(1,2,3)
#   p=2, cat.interact=3  ->  c(1,2)   (falls back: can't have 3-way with 2 vars)
#
# @param p           Number of categorical variables.
# @param cat.interact Maximum interaction order (2 or 3).
# @return Numeric vector in mix margins format.
#
.make_margins <- function(p, cat.interact) {

  if (p == 1L) return(1)

  order  <- min(as.integer(cat.interact), p)   # can't exceed number of vars
  combos <- utils::combn(p, order, simplify = FALSE)

  # Interleave with 0 separators; no trailing 0
  result <- numeric(0)
  for (i in seq_along(combos)) {
    result <- c(result, combos[[i]])
    if (i < length(combos)) result <- c(result, 0)
  }
  result
}

# --------------------------------------------------------------------------- #
# Design matrix
# --------------------------------------------------------------------------- #

# Build the design matrix for mix::ecm.mix.
#
# A diagonal identity matrix of size n_cells x n_cells specifies
# cell-specific intercepts (means) with a shared covariance matrix — the
# parameterisation used in glemb.
#
# @param factor_meta Named list of factor metadata (from .preprocess_data).
# @return Numeric identity matrix of dimension n_cells.
#
.make_design <- function(factor_meta) {
  n_cells <- prod(vapply(factor_meta, `[[`, integer(1L), "nlevels"))
  diag(rep(1, n_cells))
}

# --------------------------------------------------------------------------- #
# Dirichlet pseudocounts
# --------------------------------------------------------------------------- #

# Create pseudo-observations for the Dirichlet prior on cell probabilities.
#
# One complete (no-missing) observation is created for every categorical cell.
# Continuous variables are set to their observed column means. The rows are
# then replicated `n_pseudo` times (where n_pseudo = max(1, round(cat.prior))).
#
# These rows are appended to each bootstrap sample before calling
# mix::prelim.mix(), and removed from the imputed output afterwards.
#
# @param data_model  Working data frame (integer-coded factors, cat vars first).
# @param factor_meta Named list of factor metadata.
# @param cont_names  Character vector of continuous variable names.
# @param cat.prior   Pseudocount per cell (scalar > 0).
# @return Data frame of pseudo-observations in the same column order as
#   data_model, or NULL if cat.prior == 0.
#
.make_pseudo_obs <- function(data_model, factor_meta, cont_names, cat.prior) {

  if (cat.prior == 0) return(NULL)

  cat_names <- names(factor_meta)
  n_levels  <- vapply(factor_meta, `[[`, integer(1L), "nlevels")

  # All combinations of categorical levels (fully-crossed grid)
  level_grid        <- do.call(expand.grid, lapply(n_levels, seq_len))
  names(level_grid) <- cat_names

  # Continuous variables: column means (ignoring NA)
  cont_means <- vapply(data_model[cont_names], mean, numeric(1L), na.rm = TRUE)
  cont_df    <- as.data.frame(
    matrix(rep(cont_means, nrow(level_grid)),
           nrow  = nrow(level_grid),
           byrow = TRUE,
           dimnames = list(NULL, cont_names))
  )

  # Replicate rows by pseudocount (at least 1 row per cell)
  n_pseudo <- max(1L, round(cat.prior))
  pseudo   <- cbind(level_grid, cont_df)
  pseudo   <- pseudo[rep(seq_len(nrow(pseudo)), each = n_pseudo), ,
                     drop = FALSE]
  rownames(pseudo) <- NULL

  # Return in same column order as data_model
  pseudo[, names(data_model), drop = FALSE]
}

# --------------------------------------------------------------------------- #
# Prior resolution helpers
# --------------------------------------------------------------------------- #

# Automatic Dirichlet pseudocount: Jeffreys-like, 1/K per cell.
.resolve_cat_prior <- function(cat.prior, n, n_cells) {
  if (is.null(cat.prior)) return(1 / n_cells)
  cat.prior
}

# Automatic ridge prior: conservative, proportional to q/n.
.resolve_empri <- function(empri, n, q) {
  if (is.null(empri)) return(q / n)
  empri
}
