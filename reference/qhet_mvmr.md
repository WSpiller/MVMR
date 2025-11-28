# qhet_mvmr

Fits a multivariable Mendelian randomization model adjusting for weak
instruments. The functions requires a formatted dataframe using the
[`format_mvmr()`](https://wspiller.github.io/MVMR/reference/format_mvmr.md)
function, as well a phenotypic correlation matrix `pcor`. This should be
obtained from individual level phenotypic data, or constructed as a
correlation matrix where correlations have previously been reported.
Confidence intervals are calculated using a non-parametric bootstrap. By
default, standard errors are not produced but can be calculated by
setting `se = TRUE`. The number of bootstrap iterations is specified
using the `iterations` argument. Note that calculating confidence
intervals at present can take a substantial amount of time.

## Usage

``` r
qhet_mvmr(
  r_input,
  pcor,
  CI,
  iterations,
  ncores = parallelly::availableCores(omit = 1)
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
  exposure included in the MVMR analysis.

- CI:

  Indicates whether 95 percent confidence intervals should be calculated
  using a non-parametric bootstrap.

- iterations:

  Specifies number of bootstrap iterations for calculating 95 percent
  confidence intervals.

- ncores:

  Number of cores to use for parallel processing in bootstrap. Default
  is `parallelly::availableCores(omit = 1)`. On Windows, this is
  automatically set to 1 regardless of user input. It is recommended to
  only set this to a maximum of `parallelly::availableCores(omit = 1)`.

## Value

An dataframe containing effect estimates with respect to each exposure.

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
qhet_mvmr(r_input, pcor, CI = TRUE, iterations = 1000)
} # }
```
