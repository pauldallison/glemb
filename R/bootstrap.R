# Bootstrap iteration and related helpers

# --------------------------------------------------------------------------- #
# Single bootstrap iteration
# --------------------------------------------------------------------------- #

# Run one complete bootstrap iteration: sample, EM, impute.
#
# The EMB algorithm has three steps per iteration:
#   1. Bootstrap the incomplete data.
#   2. Run EM on the bootstrap sample to estimate parameters.
#   3. Impute the ORIGINAL data using those parameters.
#
# Imputing the original data (not the bootstrap resample) is critical:
# if the bootstrap sample were imputed, duplicated rows with missing values
# would receive independently-drawn values, diluting the signal and
# attenuating regression coefficients.
#
# @param b           Iteration index (used for seeding).
# @param data_model  Working data frame (double-coded factors first, then
#                    continuous).  Used as the imputation target.
# @param data_mat    as.matrix(data_model) — passed in to avoid recomputing.
# @param orig_pre    prelim.mix(data_mat, p) computed once outside the loop.
# @param data_ids    ID variable columns, or NULL.
# @param n           Number of rows in data_model.
# @param p           Number of categorical variables.
# @param margins     Margins vector for mix::ecm.mix (from .make_margins).
# @param design      Design matrix for mix::ecm.mix (from .make_design).
# @param factor_meta Named list of factor metadata (from .preprocess_data).
# @param col_order   Original column order (for restoring output).
# @param cat.prior   Dirichlet pseudocount per cell.
# @param empri       Ridge prior strength for continuous covariance.
# @param maxits      Maximum ECM iterations (passed as `steps` to ecm.mix).
# @param seed        Integer seed, or NULL.
# @return A list with two elements:
#   \item{imp_df}{Completed data frame with the original n rows imputed.}
#   \item{converged}{Logical; FALSE if EM hit maxits without converging.}
#
.boot_iter <- function(b, data_model, data_mat, orig_pre, data_ids, n, p,
                       margins, design, factor_meta, col_order,
                       cat.prior, empri, maxits, seed) {

  cont_names <- setdiff(names(data_model), names(factor_meta))

  # ---- Seed R's RNG for reproducible bootstrap sampling ---------------------
  if (!is.null(seed)) set.seed(seed + b)

  # ---- Draw bootstrap sample (used only for EM parameter estimation) --------
  idx      <- sample(n, n, replace = TRUE)
  boot_dat <- data_model[idx, , drop = FALSE]
  rownames(boot_dat) <- NULL

  # ---- Add Dirichlet pseudo-observations to bootstrap (for EM only) ---------
  if (cat.prior > 0) {
    pseudo   <- .make_pseudo_obs(data_model, factor_meta, cont_names, cat.prior)
    boot_dat <- rbind(boot_dat, pseudo)
  }

  # ---- Run EM on bootstrap sample to get parameter estimates ----------------
  # maxits: maximum ECM iterations; we capture printed output (showits=TRUE) to
  #   detect whether the algorithm converged or hit the limit.
  # Note: ecm.mix does not expose a ridge prior for the continuous covariance
  #   (empri); that argument is accepted by glemb() for future use.
  boot_pre  <- mix::prelim.mix(as.matrix(boot_dat), p)
  em_output <- utils::capture.output(
    em_pars <- mix::ecm.mix(boot_pre,
                            margins = margins,
                            design  = design,
                            maxits  = maxits,
                            showits = TRUE)
  )

  # ---- Convergence check -----------------------------------------------------
  # ecm.mix prints the iteration number at convergence (or the last iteration
  # if it hits the limit).  Extract the last integer from the captured output
  # and flag as non-converged if it equals maxits.
  last_iter <- suppressWarnings(
    max(as.integer(
      unlist(regmatches(em_output,
                        gregexpr("[0-9]+", em_output)))
    ), na.rm = TRUE)
  )
  converged <- is.na(last_iter) || last_iter < maxits

  # ---- Seed mix's RNG for reproducible imputation draw ----------------------
  # mix::imp.mix() requires rngseed() to have been called before it can draw
  # imputed values. When seed is NULL we use a random integer from R's RNG so
  # that mix is always initialised, even on a fresh session.
  rng_seed <- if (!is.null(seed)) seed + b else sample.int(2e9L, 1L)
  mix::rngseed(rng_seed)

  # ---- Impute the ORIGINAL data using bootstrap parameters ------------------
  # orig_pre describes the missingness in the original n observations.
  # em_pars provides the bootstrap-derived parameters (proper MI uncertainty).
  imp_mat <- mix::imp.mix(orig_pre, em_pars, data_mat)

  # Convert to data frame (imp_mat has the same n rows as the original data)
  imp_df <- as.data.frame(imp_mat)

  # ---- Restore factor levels -------------------------------------------------
  imp_df <- .restore_factors(imp_df, factor_meta)

  # ---- Restore original column order ----------------------------------------
  imp_df <- imp_df[, col_order, drop = FALSE]
  rownames(imp_df) <- NULL

  # ---- Reattach ID variables (original row order) ---------------------------
  if (!is.null(data_ids)) {
    imp_df <- cbind(data_ids, imp_df)
  }

  list(imp_df = imp_df, converged = converged)
}

# --------------------------------------------------------------------------- #
# Factor restoration
# --------------------------------------------------------------------------- #

# Convert integer-coded categorical variables back to factors with original
# levels.  Imputed values are rounded and clamped to the valid level range
# before mapping, since mix may occasionally produce values slightly outside
# the integer grid.
#
# @param df          Data frame with integer-coded categorical columns.
# @param factor_meta Named list of factor metadata.
# @return Data frame with categorical columns restored as factors.
#
.restore_factors <- function(df, factor_meta) {
  for (v in names(factor_meta)) {
    meta  <- factor_meta[[v]]
    codes <- as.integer(round(df[[v]]))
    codes <- pmin(pmax(codes, 1L), meta$nlevels)   # clamp to valid range

    # Always restore noms variables as factors so that downstream lm()/glm()
    # calls treat them as categorical regardless of original storage type.
    # (A numeric noms variable coded 0/1 or 1/2/3 would otherwise be treated
    # as continuous in regression, giving incorrect coefficients.)
    df[[v]] <- factor(meta$orig_values[codes], levels = meta$orig_values)
  }
  df
}

# --------------------------------------------------------------------------- #
# Output builder
# --------------------------------------------------------------------------- #

# Assemble the final glemb return object.
#
.build_output <- function(imputations, m, call, cat_names, cont_names,
                          idvars, cat.interact, cat.prior, empri, maxits,
                          seed, output, data) {

  result <- structure(
    list(
      imputations  = imputations,
      m            = m,
      call         = call,
      categorical  = cat_names,
      continuous   = cont_names,
      idvars       = if (is.null(idvars)) character(0L) else idvars,
      cat.interact = cat.interact,
      cat.prior    = cat.prior,
      empri        = empri,
      maxits       = maxits,
      seed         = seed
    ),
    class = "glemb"
  )

  if (output == "mids") {
    return(as.mids.glemb(result, data))
  }

  if (output == "mitml") {
    return(mitml::as.mitml.list(imputations))
  }

  result
}
