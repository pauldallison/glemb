# Pure-R imputation step for the General Location Model.
#
# This replaces mix::imp.mix() to obtain full reproducibility. mix's Fortran
# RNG cannot be completely reset by rngseed() across successive calls; a
# persistent internal state causes alternating results when glemb() is called
# more than once in a session with odd m. Using R's own RNG (seeded via
# set.seed() in .boot_iter()) eliminates this dependency entirely.

# --------------------------------------------------------------------------- #
# Main imputation function
# --------------------------------------------------------------------------- #

# Draw one completed dataset from the posterior predictive distribution.
#
# For each row with any missing value:
#   1. Compute the posterior probability over cells given the observed data.
#   2. Draw a cell from that distribution.
#   3. Fill in missing categorical variables from the drawn cell.
#   4. Draw missing continuous variables from the conditional normal given the
#      drawn cell and observed continuous variables.
#
# All random draws use R's RNG; call set.seed() before .boot_iter() for
# reproducibility.
#
# Performance strategy
# --------------------
# Rows are grouped by their continuous missingness pattern (the set of
# observed continuous-variable indices).  Within each group the expensive
# linear-algebra operations are done once and shared:
#
#   * Cholesky R = chol(Sigma[obs_z, obs_z])           — one decomposition
#   * R^{-T} applied to all observed z-vectors         — one batched backsolve
#   * R^{-T} applied to all ncells mean vectors        — one backsolve
#   * Squared Mahalanobis distances for all (row, cell) — one matrix multiply
#       via  ||a - b||^2 = ||a||^2 - 2a'b + ||b||^2
#
# Per-row work (categorical masking, cell draw, conditional-normal draw) uses
# the precomputed Cholesky and coefficient matrix K = S_mo S_oo^{-1}.
#
# @param s      Output of mix::prelim.mix() on the original data.
# @param theta  Output of mix::ecm.mix(): list with $sigma, $mu, $pi.
# @param x      Original data matrix (n rows; first p cols = integer-coded
#               categorical variables, next q cols = continuous variables).
# @return       Completed numeric matrix with the same dimensions as x.
#
.glm_impute <- function(s, theta, x) {

  p      <- s$p
  q      <- s$q
  n      <- s$n
  ncells <- s$ncells
  jmp    <- s$jmp     # p-vector: jump sizes for mixed-radix cell indexing
  xbar   <- s$xbar    # q-vector: continuous column means (from prelim.mix)
  sdv    <- s$sdv     # q-vector: continuous column SDs  (from prelim.mix)

  # ---- Reconstruct full q x q covariance matrix ----------------------------
  Sigma <- matrix(0.0, q, q)
  if (q > 0L) {
    psi <- s$psi
    sv  <- theta$sigma
    for (i in seq_len(q))
      for (j in seq_len(q))
        Sigma[i, j] <- sv[psi[i, j]]
  }

  mu <- theta$mu    # q x ncells matrix of cell means (standardized scale)
  pi <- theta$pi    # ncells vector of cell probabilities

  result   <- x     # working copy; fill in missing values below
  cat_idx  <- seq_len(p)
  cont_idx <- p + seq_len(q)

  log_pi <- log(pmax(pi, .Machine$double.xmin))  # ncells; precomputed once

  # ---- Precompute cell category table --------------------------------------
  # cell_cats[cc, j] = category level of variable j in cell cc.
  # Decoding once here avoids repeated mixed-radix arithmetic in the row loop.
  cell_cats <- matrix(0L, nrow = ncells, ncol = p)
  for (cc in seq_len(ncells)) {
    rem <- cc - 1L
    for (j in rev(seq_len(p))) {
      lv               <- (rem %/% jmp[j]) + 1L
      rem              <- rem %% jmp[j]
      cell_cats[cc, j] <- lv
    }
  }

  # ---- Find rows that need imputation --------------------------------------
  need_imp <- which(apply(x, 1L, anyNA))
  if (length(need_imp) == 0L) return(result)

  # ---- Encode each row's continuous missingness pattern -------------------
  # Rows sharing the same obs_z can reuse the same Cholesky and batch their
  # cell-likelihood computations.
  miss_cont_mat <- is.na(x[need_imp, cont_idx, drop = FALSE])  # nrows x q

  # Build pattern key: comma-separated observed indices
  obs_z_list <- lapply(seq_len(nrow(miss_cont_mat)), function(i)
                  which(!miss_cont_mat[i, ]))
  pat_keys   <- vapply(obs_z_list, function(v) paste(v, collapse = ","),
                        character(1L))

  # ---- Process rows grouped by continuous missingness pattern --------------
  for (pat in unique(pat_keys)) {

    rows_in_pat <- need_imp[pat_keys == pat]
    obs_z       <- obs_z_list[[which(pat_keys == pat)[1L]]]
    miss_z_vec  <- setdiff(seq_len(q), obs_z)   # positions of missing conts

    has_cont_obs  <- length(obs_z)    > 0L
    has_cont_miss <- length(miss_z_vec) > 0L

    # -- Cholesky of Sigma[obs_z, obs_z] (once per pattern) -----------------
    if (has_cont_obs) {
      S_oo <- Sigma[obs_z, obs_z, drop = FALSE]
      R    <- chol(S_oo)          # upper-tri: S_oo = R' R
      ldet <- 2.0 * sum(log(diag(R)))

      # R^{-T} mu[obs_z, ] — transform all cell means at once
      RinvMu <- backsolve(R, mu[obs_z, , drop = FALSE], transpose = TRUE)
      sq_mu  <- colSums(RinvMu^2)   # ncells

      # Precompute K = S_mo S_oo^{-1} for conditional-normal means
      # (used below per row; reuse same R)
      if (has_cont_miss) {
        S_mo <- Sigma[miss_z_vec, obs_z, drop = FALSE]
        # K %*% dev = S_mo S_oo^{-1} dev = S_mo R^{-1} R^{-T} dev
        # Store S_mo R^{-1}: solve R x' = S_mo' -> x = R^{-T} S_mo'
        SmoRinv <- t(backsolve(R, t(S_mo), transpose = TRUE))  # |miss_z| x |obs_z|
        # Conditional variance (same for all rows in this pattern and cell)
        # cond_var = S_mm - S_mo S_oo^{-1} S_om = S_mm - SmoRinv %*% t(SmoRinv)
        S_mm     <- Sigma[miss_z_vec, miss_z_vec, drop = FALSE]
        cond_var <- S_mm - tcrossprod(SmoRinv)
        cond_var <- (cond_var + t(cond_var)) / 2.0   # enforce symmetry
        nm       <- length(miss_z_vec)
        if (nm == 1L) {
          cond_sd <- sqrt(max(cond_var[1L, 1L], 0.0))
        } else {
          R_c <- chol(cond_var)   # for multivariate draw; same for all rows
        }
      }

      # Standardise observed continuous values for all rows in pattern at once
      Z_obs_std <- sweep(x[rows_in_pat, cont_idx[obs_z], drop = FALSE],
                         2L, xbar[obs_z], `-`)
      Z_obs_std <- sweep(Z_obs_std, 2L, sdv[obs_z], `/`)
      # R^{-T} Z_obs_std': |obs_z| x n_pat
      RinvZ  <- backsolve(R, t(Z_obs_std), transpose = TRUE)
      sq_z   <- colSums(RinvZ^2)   # n_pat
      # Mahalanobis: maha[i, cc] = sq_z[i] - 2*(RinvZ'*RinvMu)[i,cc] + sq_mu[cc]
      maha   <- outer(sq_z, sq_mu, `+`) - 2.0 * crossprod(RinvZ, RinvMu)
                # n_pat x ncells
    }

    # -- Per-row: categorical masking, cell draw, conditional normal ---------
    for (ii in seq_along(rows_in_pat)) {
      i   <- rows_in_pat[ii]
      w_i <- as.integer(x[i, cat_idx])
      z_i <- x[i, cont_idx]

      miss_w <- which(is.na(w_i))
      obs_w  <- which(!is.na(w_i))
      # miss_z is the same for all rows in this pattern
      miss_z <- miss_z_vec

      # ---- Cell posterior -------------------------------------------------
      if (length(miss_w) == 0L) {
        # All categoricals observed: cell is deterministic
        c_i <- 1L + sum((w_i - 1L) * jmp)

      } else {
        log_post <- log_pi   # ncells; copy

        # Continuous likelihood: use precomputed maha row
        if (has_cont_obs) {
          log_post <- log_post - 0.5 * (ldet + maha[ii, ])
        }

        # Mask cells incompatible with observed categoricals
        if (length(obs_w) > 0L) {
          compat <- rep(TRUE, ncells)
          for (j in obs_w) compat <- compat & (cell_cats[, j] == w_i[j])
          log_post[!compat] <- -Inf
        }

        # Normalise and draw
        mx   <- max(log_post[is.finite(log_post)])
        post <- exp(log_post - mx)
        post[!is.finite(log_post)] <- 0.0
        post <- post / sum(post)
        c_i  <- sample.int(ncells, 1L, prob = post)

        # Fill imputed categorical values from drawn cell
        for (j in miss_w) result[i, j] <- cell_cats[c_i, j]
      }

      # ---- Conditional normal draw for missing continuous -----------------
      if (has_cont_miss && q > 0L) {
        mu_c <- mu[, c_i]   # q-vector, standardised

        if (!has_cont_obs) {
          # No observed continuous: draw from marginal N(mu_c, Sigma)
          R_S    <- chol(Sigma)
          z_draw <- mu_c + drop(t(R_S) %*% rnorm(q))
          y_imp  <- z_draw * sdv + xbar
          result[i, cont_idx[miss_z]] <- y_imp[miss_z]

        } else {
          # Conditional mean: mu_miss_c + SmoRinv %*% (R^{-T} dev)
          dev     <- Z_obs_std[ii, ] - mu_c[obs_z]
          Rdev    <- backsolve(R, dev, transpose = TRUE)   # R^{-T} dev
          cond_mu <- mu_c[miss_z] + SmoRinv %*% Rdev

          if (nm == 1L) {
            z_miss <- rnorm(1L, cond_mu, cond_sd)
          } else {
            z_miss <- drop(cond_mu + t(R_c) %*% rnorm(nm))
          }
          # Unstandardise
          y_imp <- z_miss * sdv[miss_z] + xbar[miss_z]
          result[i, cont_idx[miss_z]] <- y_imp
        }
      }
    }   # end per-row loop
  }   # end pattern loop

  result
}
