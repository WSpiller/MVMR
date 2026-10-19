# qhet_mvmr

Fits a multivariable Mendelian randomization model adjusting for weak
instruments. The functions requires a formatted dataframe using the
[`format_mvmr()`](https://wspiller.github.io/MVMR/reference/format_mvmr.md)
function, as well a phenotypic correlation matrix `pcor`. This should be
obtained from individual level phenotypic data, or constructed as a
correlation matrix where correlations have previously been reported. By
default, confidence intervals are not produced but can be calculated by
setting `CI = TRUE`. They are calculated using either a non-parametric
bootstrap over the genetic variants (`CI_method = "bootstrap"`, the
default), with the number of bootstrap iterations specified using the
`iterations` argument, or a leave-one-variant-out jackknife
(`CI_method = "jackknife"`). Sanderson, Spiller and Bowden (2021)
recommend the jackknife, as with a moderate number of variants bootstrap
resamples can give very weak instruments. Note that calculating
confidence intervals can take a substantial amount of time.

## Usage

``` r
qhet_mvmr(
  r_input,
  pcor,
  CI,
  iterations,
  ncores = parallelly::availableCores(omit = 1),
  CI_method = c("bootstrap", "jackknife")
)
```

## Arguments

- r_input:

  A formatted data frame using the
  [`format_mvmr()`](https://wspiller.github.io/MVMR/reference/format_mvmr.md)
  function or an object of class `MRMVInput` from
  [`MendelianRandomization::mr_mvinput()`](https://rdrr.io/pkg/MendelianRandomization/man/mr_mvinput.html)

- pcor:

  A phenotypic correlation matrix including the correlation between each
  exposure included in the MVMR analysis. It is used to construct the
  covariance between the estimated effects of each genetic variant on
  each exposure, as `pcor` multiplied by the outer product of the
  variant's standard errors. This is appropriate when the exposure GWAS
  were estimated in the same sample. If the exposure GWAS were estimated
  in non-overlapping samples these covariances are zero, and the
  identity matrix `diag(K)`, where K is the number of exposures, should
  be supplied instead.

- CI:

  Indicates whether 95 percent confidence intervals should be
  calculated.

- iterations:

  Specifies number of bootstrap iterations for calculating 95 percent
  confidence intervals. Not used when `CI_method = "jackknife"`.

- ncores:

  Number of cores to use for parallel processing in the bootstrap or
  jackknife. Default is `parallelly::availableCores(omit = 1)`. On
  Windows, this is automatically set to 1 regardless of user input. It
  is recommended to only set this to a maximum of
  `parallelly::availableCores(omit = 1)`.

- CI_method:

  Method used to calculate confidence intervals when `CI = TRUE`:
  `"bootstrap"` (the default) for bias-corrected and accelerated (BCa)
  non-parametric bootstrap intervals, or `"jackknife"` for normal-based
  intervals using the leave-one-variant-out jackknife standard error.

## Value

A dataframe containing effect estimates with respect to each exposure,
and when `CI = TRUE` their 95 percent confidence intervals. With
`CI_method = "jackknife"` the jackknife standard errors are also
returned.

## Details

Effects are estimated by minimising a Q-statistic that includes an
additive heterogeneity parameter tau-squared. Following Sanderson,
Spiller and Bowden (2021), tau-squared is chosen so that the minimised
Q-statistic equals L - K, where L is the number of genetic variants and
K the number of exposures. Tau-squared is constrained to be
non-negative, and is fixed at zero when the minimised Q-statistic
without heterogeneity is already no greater than L - K.

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
qhet_mvmr(r_input, pcor, CI = TRUE, iterations = 1000)
qhet_mvmr(r_input, pcor, CI = TRUE, CI_method = "jackknife")
} # }
```
