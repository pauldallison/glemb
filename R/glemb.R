#' Multiple Imputation Using the General Location Model with the EMB Algorithm
#'
#' @description
#' `glemb()` performs multiple imputation of mixed continuous and categorical
#' data using the General Location Model (GLM) combined with the
#' Expectation-Maximization with Bootstrapping (EMB) algorithm.
#'
#' @param data A data frame containing the data to be imputed. All variables
#'   not listed in `noms` or `idvars` are treated as continuous. There must
#'   be at least one continuous variable. Character variables are not accepted;
#'   convert to numeric first.
#' @param m Positive integer. Number of imputed datasets to generate.
#'   Default is `20`.
#' @param noms Variable names to treat as categorical (nominal) in the
#'   imputation model. Can be supplied quoted or unquoted:
#'   `noms = c("pov", "race")` or `noms = c(pov, race)`. These variables
#'   may be stored as numeric (e.g. 0/1 binary codes, 1/2/3 integer codes)
#'   or as `factor`. All `noms` variables are returned as factors after
#'   imputation, ensuring they are treated as categorical in downstream
#'   regression calls. At least one variable must be listed.
#' @param idvars Character vector of column names to be carried through to the
#'   output unchanged but excluded from the imputation model. Typically used
#'   for subject or case identifiers. Default is `NULL`.
#' @param cat.interact Integer, either `2` (default) or `3`. Maximum order of
#'   interactions among categorical variables in the log-linear model. `2`
#'   includes all pairwise interactions; `3` includes all three-way (and
#'   lower) interactions. Option `3` should be used with caution as it may
#'   greatly increase the number of parameters to be estimated.
#' @param cat.prior Numeric scalar or `NULL`. Pseudocount for the Dirichlet
#'   prior on categorical cell probabilities, used to stabilise estimation
#'   with sparse cells. `NULL` (default) selects an automatic value
#'   proportional to `1 / n_cells`. `0` applies no prior. A positive value
#'   specifies a fixed pseudocount per cell.
#' @param empri Numeric scalar or `NULL`. Reserved for a future ridge prior on
#'   the within-cell continuous covariance matrix, analogous to the `empri`
#'   argument in [Amelia::amelia()]. Currently accepted but not applied, as
#'   [mix::ecm.mix()] does not expose this parameter. Set to `0` or `NULL`.
#' @param maxits Positive integer. Maximum number of ECM iterations allowed
#'   per bootstrap sample. Default is `1000`. A warning is issued for any
#'   bootstrap iteration that reaches this limit without converging.
#' @param seed Positive integer or `NULL`. Random seed for reproducibility.
#'   When supplied, both R's random number generator (used for bootstrap
#'   sampling) and `mix`'s internal generator (used for imputation draws) are
#'   seeded deterministically as `seed + b` on iteration `b`. Must be a
#'   positive integer no greater than 2,000,000,000. Default is `NULL` (no
#'   seed set). For reproducibility, it is essential to set the seed within
#'   the `glemb()` function. Using `set.seed()` before calling `glemb()` will
#'   have no effect on the output.
#' @param output Character scalar. `"mids"` (default) returns a `mids` object
#'   compatible with [mice::pool()] (requires the `mice` package). `"mitml"`
#'   returns a `mitml.list` object compatible with [mitml::testEstimates()]
#'   (requires the `mitml` package). `"list"` returns a plain list of `m`
#'   completed data frames.
#' @param p2s Integer, `0` or `1`. Verbosity of console output. `0` = silent;
#'   `1` = show progress (default).
#'
#' @return An object of class `"glemb"`, which is a list with components:
#' \describe{
#'   \item{`imputations`}{A list of `m` completed data frames.}
#'   \item{`m`}{Number of imputations.}
#'   \item{`call`}{The matched call.}
#'   \item{`categorical`}{Names of variables treated as categorical.}
#'   \item{`continuous`}{Names of variables treated as continuous.}
#'   \item{`idvars`}{Names of ID variables excluded from the imputation model.}
#'   \item{`cat.interact`}{Interaction order used in the log-linear model.}
#'   \item{`cat.prior`}{Dirichlet pseudocount used (after automatic selection
#'     if applicable).}
#'   \item{`empri`}{Ridge prior strength used (after automatic selection if
#'     applicable).}
#'   \item{`seed`}{Seed used, or `NULL`.}
#' }
#'
#' @details
#' ## Statistical model
#'
#' The General Location Model (Olkin & Tate 1961; Little & Schluchter 1987)
#' factorises the joint distribution of continuous variables **y** and
#' categorical variables **x** as:
#' \deqn{p(\mathbf{y}, \mathbf{x}) =
#'       p(\mathbf{y} \mid \mathbf{x})\, p(\mathbf{x})}
#' The categorical margin \eqn{p(\mathbf{x})} is modelled with a log-linear
#' model restricted to at most `cat.interact`-way interactions among the
#' categorical variables. The continuous conditional
#' \eqn{p(\mathbf{y} \mid \mathbf{x} = c)} is multivariate normal with a
#' cell-specific mean vector and a covariance matrix **shared** across all
#' cells. This model implies that the conditional mean of every continuous
#' variable is a linear function of all the other variables (both continuous and
#' categorical). The conditional distribution of every categorical variable
#' is a logistic regression function of all the other variables.
#'
#' ## EMB algorithm
#'
#' Proper multiple imputations are generated via the EMB algorithm
#' (Honaker, King & Blackwell 2011):
#' 1. Draw a bootstrap sample (with replacement) from the incomplete data.
#' 2. Run the ECM algorithm ([mix::ecm.mix()]) on the bootstrap sample to
#'    obtain parameter estimates.
#' 3. Draw one set of imputed values from the posterior predictive
#'    distribution given those parameters ([mix::imp.mix()]).
#' 4. Repeat steps 1–3 `m` times to produce `m` completed datasets.
#'
#' ## Sparsity
#'
#' When the number of cells in the multiway table for the categorical variables
#' is large relative to sample size, some cells may contain few or no
#' observations. A Dirichlet prior (controlled by `cat.prior`) is applied by
#' adding `cat.prior` pseudo-observations per cell to the bootstrap sample
#' before estimation. `glemb()` warns when the number of cells exceeds
#' `n / 5`.
#'
#' @references
#' Olkin, I. & Tate, R.F. (1961). Multivariate correlation models with mixed
#' discrete and continuous variables. *Annals of Mathematical Statistics*,
#' 32, 448--465.
#'
#' Little, R.J.A. & Schluchter, M.D. (1987). Maximum likelihood estimation
#' for mixed continuous and categorical data with missing values.
#' *Biometrika*, 74, 165--173.
#'
#' Honaker, J., King, G. & Blackwell, M. (2011). Amelia II: A program for
#' missing data. *Journal of Statistical Software*, 45(7), 1--47.
#'
#' Schafer, J.L. (1997). *Analysis of Incomplete Multivariate Data*.
#' Chapman & Hall.
#'
#' @examples
#' \dontrun{
#' # Basic usage — specify categorical variables with noms (unquoted names)
#' imp <- glemb(mydata, m = 20, noms = c(race, gender), seed = 123)
#'
#' # Equivalent using quoted names
#' imp <- glemb(mydata, m = 20, noms = c("race", "gender"), seed = 123)
#'
#' # Access the third imputed dataset
#' imp$imputations[[3]]
#'
#' # Allow three-way interactions among categorical variables
#' imp <- glemb(mydata, m = 20, noms = c(race, gender),
#'              cat.interact = 3, seed = 123)
#'
#' # Exclude a subject ID from the imputation model
#' imp <- glemb(mydata, m = 20, noms = c(race, gender),
#'              idvars = "subject_id", seed = 123)
#' }
#'
#' @seealso [print.glemb()], [summary.glemb()], [as.mids.glemb()]
#'
#' @export
glemb <- function(data,
                  m            = 20L,
                  noms         = NULL,
                  idvars       = NULL,
                  cat.interact = 2L,
                  cat.prior    = NULL,
                  empri        = NULL,
                  maxits       = 1000L,
                  seed         = NULL,
                  output       = "mids",
                  p2s          = 1L) {

  cl <- match.call()

  # ---- Allow unquoted variable names in noms --------------------------------
  # Capture the expression before R evaluates it, so that
  # noms = c(pov, race) works as well as noms = c("pov", "race").
  noms_sub <- substitute(noms)
  noms <- tryCatch(
    noms,                                       # works if already quoted
    error = function(e) .parse_noms_expr(noms_sub)  # unquoted names
  )

  # ---- Input validation ------------------------------------------------------
  .check_inputs(data, m, noms, idvars, cat.interact, cat.prior, empri,
                maxits, seed, output, p2s)

  # ---- Preprocess data -------------------------------------------------------
  prep <- .preprocess_data(data, noms, idvars)
  # Returns: data_model, data_ids, factor_meta, cat_names, cont_names,
  #          col_order, p

  n       <- nrow(prep$data_model)
  p       <- prep$p
  q       <- length(prep$cont_names)
  margins <- .make_margins(p, cat.interact)
  design  <- .make_design(prep$factor_meta)
  n_cells <- ncol(design)

  # ---- Resolve prior strengths -----------------------------------------------
  cat.prior <- .resolve_cat_prior(cat.prior, n, n_cells)
  empri     <- .resolve_empri(empri, n, q)

  # ---- Sparsity warning ------------------------------------------------------
  if (n_cells > n / 5) {
    warning("The number of cells (", n_cells, ") is large relative to n (",
            n, "). Consider reducing the number of categorical variables or ",
            "their levels, or increase cat.prior.", call. = FALSE)
  }

  # ---- Guard against both priors being zero with many cells ------------------
  if (cat.prior == 0 && empri == 0 && n_cells > 20L) {
    warning("cat.prior = 0 and empri = 0 with ", n_cells, " cells may cause ",
            "numerical instability (singular covariance matrix). Consider ",
            "setting cat.prior > 0.", call. = FALSE)
  }

  if (p2s > 0L) {
    message("glemb: ", p, " categorical variable(s), ", q,
            " continuous variable(s), ", n_cells, " cell(s)")
    message("glemb: cat.prior = ", round(cat.prior, 4),
            ",  empri = ", round(empri, 4))
    message("glemb: running ", m, " bootstrap iterations...")
  }

  # ---- Pre-compute preliminary analysis on original data ---------------------
  # This is used by imp.mix in every iteration to impute the original
  # observations (not the bootstrap resample) — the correct EMB approach.
  data_mat <- as.matrix(prep$data_model)
  orig_pre <- mix::prelim.mix(data_mat, p)

  # ---- Bootstrap loop --------------------------------------------------------
  imputations   <- vector("list", m)
  nonconv_iters <- integer(0)   # iterations where EM may not have converged

  for (b in seq_len(m)) {
    if (p2s > 0L && (b == 1L || b %% 10L == 0L)) {
      message("  iteration ", b, " of ", m)
    }

    result <- .boot_iter(
      b           = b,
      data_model  = prep$data_model,
      data_mat    = data_mat,
      orig_pre    = orig_pre,
      data_ids    = prep$data_ids,
      n           = n,
      p           = p,
      margins     = margins,
      design      = design,
      factor_meta = prep$factor_meta,
      col_order   = prep$col_order,
      cat.prior   = cat.prior,
      empri       = empri,
      maxits      = maxits,
      seed        = seed
    )

    imputations[[b]] <- result$imp_df
    if (!result$converged) nonconv_iters <- c(nonconv_iters, b)
  }

  if (length(nonconv_iters) > 0L) {
    warning(length(nonconv_iters), " of ", m,
            " bootstrap iteration(s) may not have converged (hit maxits = ",
            maxits, "): iteration(s) ",
            paste(nonconv_iters, collapse = ", "), ". ",
            "Consider increasing maxits or checking for sparse cells.",
            call. = FALSE)
  }

  # ---- Check that all observed categories appear in at least one imputation --
  # If a category is so rare that bootstrap resampling never selects it,
  # it will be absent from all imputed values.
  for (v in prep$cat_names) {
    missing_rows <- which(is.na(prep$data_model[[v]]))
    if (length(missing_rows) == 0L) next

    imputed_vals <- unlist(lapply(imputations, function(df) {
      as.character(df[[v]][missing_rows])
    }))
    observed_cats <- as.character(prep$factor_meta[[v]]$orig_values)
    absent_cats   <- setdiff(observed_cats, unique(imputed_vals))

    if (length(absent_cats) > 0L) {
      warning("Variable '", v, "': the following categories never appeared ",
              "as imputed values across all ", m, " imputations: ",
              paste(absent_cats, collapse = ", "), ". This may indicate very ",
              "rare categories. Consider increasing m or cat.prior.",
              call. = FALSE)
    }
  }

  if (p2s > 0L) message("glemb: done.")

  # ---- Build and return output object ----------------------------------------
  .build_output(
    imputations  = imputations,
    m            = m,
    call         = cl,
    cat_names    = prep$cat_names,
    cont_names   = prep$cont_names,
    idvars       = idvars,
    cat.interact = cat.interact,
    cat.prior    = cat.prior,
    empri        = empri,
    maxits       = maxits,
    seed         = seed,
    output       = output,
    data         = data
  )
}
