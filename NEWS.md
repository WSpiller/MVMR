# MVMR (development version)

* `strength_mvmr()` and `strhet_mvmr()` now divide the conditional Q-statistic by L - (K - 1), where L is the number of variants and K the number of exposures, as in Equation 7 of Sanderson, Spiller and Bowden (2021, <https://doi.org/10.1002/sim.9133>). They previously divided by L, which slightly understated the conditional F-statistics. The `strhet_mvmr()` documentation now describes its estimator as iteratively reweighted least squares rather than exact Q-statistic minimisation.
* `pleiotropy_mvmr()` now computes its p-value on L - K degrees of freedom, as in Section 3.2 of Sanderson, Spiller and Bowden (2021). It previously used L - K - 1, which gave p-values that were too small. Its documentation now notes that the Q-statistic is evaluated at the IVW estimate rather than at the Q-minimising estimate used in the paper.
* `qhet_mvmr()` now chooses the heterogeneity parameter tau-squared so that the minimised Q-statistic equals L - K, as in Sanderson, Spiller and Bowden (2021). It previously targeted L - 2, which is only correct with two exposures. Tau-squared is now constrained to be non-negative (it was searched over -10 to 10, which could give negative weights and is poorly scaled for summary data), and is found by root finding with the Q-statistic minimised by BFGS from the IVW estimate. Effect estimates will differ from previous versions, including with two exposures.
* `qhet_mvmr()` gains a `CI_method` argument. `CI_method = "jackknife"` calculates leave-one-variant-out jackknife standard errors and normal-based 95% confidence intervals, as recommended by Sanderson, Spiller and Bowden (2021); the default `"bootstrap"` retains the BCa bootstrap intervals. The documentation no longer refers to a non-existent `se` argument.
* The `qhet_mvmr()` documentation now explains how `pcor` is used, and that the identity matrix should be supplied when the exposure GWAS were estimated in non-overlapping samples.

# MVMR 0.4.8

* `strhet_mvmr()` has been reimplemented. The previous version built a combinatorial grid of candidate coefficients with `utils::combn()` (which could exhaust memory with four or more exposures) and never actually minimised the Q-statistic, returning invalid conditional F-statistics. It now estimates the coefficients by iteratively reweighted least squares, giving distinct, well-identified conditional F-statistics for each exposure. Reported values will differ from previous versions.

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
