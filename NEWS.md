# glemb 0.1.3

* Fixed a reproducibility bug where calling `glemb()` with an odd number of
  imputations (`m`) caused results to alternate between two distinct sets of
  imputations across successive calls with the same seed. The imputation step
  now uses R's own random number generator instead of the Fortran RNG in the
  `mix` package, whose internal state could not be fully reset between calls.

# glemb 0.1.2

* Added `max.cells` argument to `glemb()` (default `500000`). A pre-flight
  check now stops with an informative error before `glemb` allocates any memory
  when the nominal variables would produce more cells than this limit. Set
  `max.cells = Inf` to disable the check.
* Added `meanmodel` argument to `glemb()`. `"main"` (default) constrains the
  continuous variable means to depend only on the main effects of the
  categorical variables, reducing the number of mean parameters from
  `n_cells` to `1 + sum(K_j - 1)`. `"saturated"` restores the original
  cell-specific intercepts (unrestricted means).
* The `meanmodel` value is stored in the returned object and displayed by
  `print()` and `summary()`.

# glemb 0.1.1

* Added `output = "mitml"` option, returning a `mitml.list` object compatible
  with `mitml::testEstimates()`.
* Several bug fixes.

# glemb 0.1.0

Initial release.

* `glemb()`: multiple imputation for mixed continuous and categorical data
  using the General Location Model with the EMB algorithm.
* Factor variables automatically treated as categorical; all other variables
  treated as continuous.
* Log-linear model restricted to 2-way interactions by default; 3-way
  interactions available via `cat.interact = 3`.
* Dirichlet prior on cell probabilities (`cat.prior`) for sparse-cell
  stabilisation.
* Ridge prior on continuous covariance matrix (`empri`).
* Output as list of completed data frames (default) or `mids` object
  compatible with `mice::pool()`.
* `print()` and `summary()` methods for `glemb` objects.
