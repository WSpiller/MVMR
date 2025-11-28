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
pleiotropy_mvmr(r_input, gencov = 0)
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

## Value

A Q-statistic for instrument validity and the corresponding p-value

## References

Sanderson, E., et al., An examination of multivariable Mendelian
randomization in the single-sample and two-sample summary data settings.
International Journal of Epidemiology, 2019, 48, 3, 713–727.
[doi:10.1093/ije/dyy262](https://doi.org/10.1093/ije/dyy262)

## Author

Wes Spiller; Eleanor Sanderson; Jack Bowden.

## Examples

``` r
if (FALSE) { # \dontrun{
pleiotropy_mvmr(r_input, covariances)
} # }
```
