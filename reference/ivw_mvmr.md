# ivw_mvmr

Fits an IVW multivariable Mendelian randomization model using first
order weights.

## Usage

``` r
ivw_mvmr(r_input, gencov = 0)
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

An dataframe containing MVMR results, including estimated coefficients,
their standard errors, t-statistics, and corresponding (two-sided)
p-values.

## References

Sanderson, E., et al., An examination of multivariable Mendelian
randomization in the single-sample and two-sample summary data settings.
International Journal of Epidemiology, 2019, 48, 3, 713–727.
[doi:10.1093/ije/dyy262](https://doi.org/10.1093/ije/dyy262)

## Author

Wes Spiller; Eleanor Sanderson; Jack Bowden.

## Examples

``` r
r_input <- format_mvmr(
  BXGs = rawdat_mvmr[, c("LDL_beta", "HDL_beta")],
  BYG = rawdat_mvmr$SBP_beta,
  seBXGs = rawdat_mvmr[, c("LDL_se", "HDL_se")],
  seBYG = rawdat_mvmr$SBP_se,
  RSID = rawdat_mvmr$SNP
)
ivw_mvmr(r_input)
#> Warning: Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.
#> 
#> Multivariable MR
#> 
#>               Estimate Std. Error    t value  Pr(>|t|)
#> exposure1 -0.031003996 0.01302925 -2.3795686 0.0186526
#> exposure2  0.006039167 0.01029181  0.5867933 0.5582678
#> 
#> Residual standard error: 2.209 on 143 degrees of freedom
#> 
#> 
#>               Estimate Std. Error    t value  Pr(>|t|)
#> exposure1 -0.031003996 0.01302925 -2.3795686 0.0186526
#> exposure2  0.006039167 0.01029181  0.5867933 0.5582678
```
