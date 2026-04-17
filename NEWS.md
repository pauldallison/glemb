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
