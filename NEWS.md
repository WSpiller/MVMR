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
