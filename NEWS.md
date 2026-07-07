# MVMR 0.4.7

* `qhet_mvmr()` now uses the estimated heterogeneity parameter (the minimiser) rather than the minimised objective value when constructing the model weights. This corrects the effect estimates, which will differ from previous versions.
* `snpcov_mvmr()` now returns correct covariance matrices when the genetic instrument data are supplied as a matrix rather than a data frame.
* `strength_mvmr()` and `pleiotropy_mvmr()` now select the covariance-matrix calculation based on whether `gencov` is a list, fixing a division-by-zero that occurred when a `gencov` list for exactly two variants was supplied.
* `ivw_mvmr()` no longer emits a spurious warning about fixing the covariance at zero; the `gencov` argument does not affect the IVW estimates and is retained only for interface consistency.
* `mvmr()` no longer fits the same weighted regression twice.
* Tidied `format_mvmr()` internal column naming.

# MVMR 0.4.6

* Add vignette on estimating phenotypic correlations.

# MVMR 0.4.5

* Bump roxygen2 to 8.0.0 and add package level helpfile.

# MVMR 0.4.4

* Implement some optimizations.

# MVMR 0.4.3

* The `snpcov_mvmr()` function accidentally omitted intercepts from its regressions of exposure on genotype. This has been fixed.

# MVMR 0.4.2

* The `qhet_mvmr()` function with argument `CI = TRUE` is now faster as it can use multiple processor cores (except on Windows) and it calculates the bootstrap confidence interval limits more efficiently as an unnecessarily repeated function call has been removed (thanks @nickhir).

# MVMR 0.4.1

* We have slightly improved the speed of the `phenocov_mvmr()` function (thanks @shiyw).

* In the `phenocov_mvmr()` function we renamed the `Pcov` argument to `pcor` to reflect that this is a correlation matrix. This matches the argument name in `qhet_mvmr()` (thanks @mooreann).
