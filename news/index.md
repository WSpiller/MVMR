# Changelog

## MVMR 0.4.2

- The
  [`qhet_mvmr()`](https://wspiller.github.io/MVMR/reference/qhet_mvmr.md)
  function with argument `CI = TRUE` is now faster as it can use
  multiple processor cores (except on Windows) and it calculates the
  bootstrap confidence interval limits more efficiently as an
  unnecessarily repeated function call has been removed (thanks
  [@nickhir](https://github.com/nickhir)).

## MVMR 0.4.1

- We have slightly improved the speed of the
  [`phenocov_mvmr()`](https://wspiller.github.io/MVMR/reference/phenocov_mvmr.md)
  function (thanks [@shiyw](https://github.com/shiyw)).

- In the
  [`phenocov_mvmr()`](https://wspiller.github.io/MVMR/reference/phenocov_mvmr.md)
  function we renamed the `Pcov` argument to `pcor` to reflect that this
  is a correlation matrix. This matches the argument name in
  [`qhet_mvmr()`](https://wspiller.github.io/MVMR/reference/qhet_mvmr.md)
  (thanks [@mooreann](https://github.com/mooreann)).
