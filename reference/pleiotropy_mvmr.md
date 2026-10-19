# pleiotropy_mvmr

Calculates modified form of Cochran's Q statistic measuring
heterogeneity in causal effect estimates obtained using each genetic
variant. Observed heterogeneity is indicative of a violation of the
exclusion restriction assumption in MR (validity), which can result in
biased effect estimates. The function takes a formatted dataframe as an
input, obtained using the function
[`format_mvmr()`](https://wspiller.github.io/MVMR/reference/format_mvmr.md).
Additionally, covariance matrices for estimated effects of individual
genetic variants on each exposure can also be provided. These can be
estimated using external data by applying the
[`snpcov_mvmr()`](https://wspiller.github.io/MVMR/reference/snpcov_mvmr.md)
or
[`phenocov_mvmr()`](https://wspiller.github.io/MVMR/reference/phenocov_mvmr.md)
functions, are input manually. The function returns a dataframe
including the conditional Q-statistic for instrument validity, and a
corresponding P-value.

## Usage

``` r
pleiotropy_mvmr(r_input, gencov = 0, estimator = c("ivw", "qmin"))
```

## Arguments

- r_input:

  A formatted data frame using the
  [`format_mvmr()`](https://wspiller.github.io/MVMR/reference/format_mvmr.md)
  function or an object of class `MRMVInput` from
  [`MendelianRandomization::mr_mvinput()`](https://rdrr.io/pkg/MendelianRandomization/man/mr_mvinput.html)

- gencov:

  Calculating heterogeneity statistics requires the covariance between
  the effect of the genetic variants on each exposure to be known. This
  can either be estimated from individual level data, be assumed to be
  zero, or fixed at zero using non-overlapping samples of each exposure
  GWAS. A value of `0` is used by default.

- estimator:

  The causal effect estimates at which the Q-statistic is evaluated:
  `"ivw"` (the default) for the inverse variance weighted estimates, or
  `"qmin"` for the estimates that minimise the Q-statistic.

## Value

A Q-statistic for instrument validity and the corresponding p-value

## Details

By default the Q-statistic is evaluated at the inverse variance weighted
(IVW) estimates of the causal effects and compared with a chi-squared
distribution on L - K degrees of freedom, where L is the number of
genetic variants and K the number of exposures. This is the Q-statistic
for heterogeneity in MVMR described in Sanderson, Spiller and Bowden
(2021). It is not robust to weak instruments and will be inflated when
the instruments are conditionally weak, so instrument strength should
first be assessed with
[`strength_mvmr()`](https://wspiller.github.io/MVMR/reference/strength_mvmr.md);
if the instruments are weak the results should not be interpreted.
Setting `estimator = "qmin"` instead evaluates the Q-statistic at the
causal effect estimates that minimise it, which are robust to weak
instruments, as discussed in Section 3.2 of Sanderson, Spiller and
Bowden (2021).

## References

Sanderson, E., et al., An examination of multivariable Mendelian
randomization in the single-sample and two-sample summary data settings.
International Journal of Epidemiology, 2019, 48, 3, 713–727.
[doi:10.1093/ije/dyy262](https://doi.org/10.1093/ije/dyy262)

Sanderson, E., Spiller, W., and Bowden, J., Testing and correcting for
weak and pleiotropic instruments in two-sample multivariable Mendelian
randomization. Statistics in Medicine, 2021, 40, 25, 5434–5452.
[doi:10.1002/sim.9133](https://doi.org/10.1002/sim.9133)

## Author

Wes Spiller; Eleanor Sanderson; Jack Bowden.

## Examples

``` r
if (FALSE) { # \dontrun{
pleiotropy_mvmr(r_input, covariances)
} # }
```
