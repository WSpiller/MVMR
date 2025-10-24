# MVMR 0.4.2

* The `qhet_mvmr()` function with argument `CI = TRUE` is now faster as it can use multiple processor cores (except on Windows) and it calculates the bootstrap confidence interval limits more efficiently as an unnecessarily repeated function call has been removed (thanks @nickhir).

# MVMR 0.4.1

* We have slightly improved the speed of the `phenocov_mvmr()` function (thanks @shiyw).

* In the `phenocov_mvmr()` function we renamed the `Pcov` argument to `pcor` to reflect that this is a correlation matrix. This matches the argument name in `qhet_mvmr()` (thanks @mooreann).
