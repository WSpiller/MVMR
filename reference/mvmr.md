# mvmr (legacy)

Note: This function is from the old version of the MVMR package and will
be replaced in the future: The gencov argument should be set to zero
when using `mvmr()`.

## Usage

``` r
mvmr(r_input, gencov, weights)
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

- weights:

  A value specifying the inverse variance weights used to calculate IVW
  estimate and Cochran's Q statistic. Currently only first order weights
  are available (`1`).

## Value

An object of class `"MVMRIVW"` containing the following components:

- `summary`:

  A summary of the MVMR regression model, including estimated
  coefficients, standard errors, t-statistics, p-values, and
  heterogeneity statistics.

- `coef`:

  The estimated coefficients, their standard errors, t-statistics, and
  corresponding (two-sided) p-values.

- `Q_strength`:

  A data frame displaying modified Cochran's Q statistics for assessing
  instrument strength with respect to each exposure. The Q-statistic
  increases proportionally with instrument strength, and analogous to
  univariate MR analyses, a value equal to or greater than 10 can be
  used as a minimum threshold for instrument strength. Note that for
  these statistics it is not informative to evaluate p-values.

- `Q_valid`:

  A modified form of Cochran's Q statistic measuring heterogeneity in
  causal effect estimates obtained using each genetic variant. Observed
  heterogeneity is indicative of a violation of the exclusion
  restriction assumption in MR (validity), which can result in biased
  effect estimates.

- `p_valid`:

  A p-value corresponding to the heterogeneity measure for instrument
  validity (`Q_valid`)

## Details

Fits an IVW multivariable Mendelian randomization model using first
order weights. The function returns an object of class `"MVMRIVW"`,
containing regression estimates, estimated heterogeneity as a measure of
instrument strength (`Q_strength`), and estimated heterogeneity as a
measure of instrument validity (`Q_valid`).

## References

Sanderson, E., et al., An examination of multivariable Mendelian
randomization in the single-sample and two-sample summary data settings.
International Journal of Epidemiology, 2019, 48, 3, 713–727.
[doi:10.1093/ije/dyy262](https://doi.org/10.1093/ije/dyy262)

## Author

Wes Spiller; Eleanor Sanderson; Jack Bowden.

## Examples

``` r
# Example using format_mvmr formatted data
r_input <- format_mvmr(
  BXGs = rawdat_mvmr[, c("LDL_beta", "HDL_beta")],
  BYG = rawdat_mvmr$SBP_beta,
  seBXGs = rawdat_mvmr[, c("LDL_se", "HDL_se")],
  seBYG = rawdat_mvmr$SBP_se,
  RSID = rawdat_mvmr$SNP
)
mvmr(r_input, 0, 1)
#> 
#> Multivariable MR
#> 
#>               Estimate Std. Error    t value  Pr(>|t|)
#> exposure1 -0.031003996 0.01302925 -2.3795686 0.0186526
#> exposure2  0.006039167 0.01029181  0.5867933 0.5582678
#> 
#> Residual standard error: 2.209 on 143 degrees of freedom
#> 
#> F-statistic: 3.27 on 2 and 143 DF, p-value: 0.0408
#> 
#> ------------------------------
#> Q-Statistics for instrument strength:
#> 
#>   exposure1 exposure2
#> Q  9739.922  10401.49
#> 
#> ------------------------------
#> Q-Statistic for instrument validity:
#> 
#> 695.5924 on 142 DF , p-value: 7.338e-74

# Example using MRMVInput formatted data from the MendelianRandomization package
if (require("MendelianRandomization", quietly = TRUE)) {
  bx <- as.matrix(rawdat_mvmr[, c("LDL_beta", "HDL_beta")])
  bxse <- as.matrix(rawdat_mvmr[, c("LDL_se", "HDL_se")])
  dat <- MendelianRandomization::mr_mvinput(
    bx = bx,
    bxse = bxse,
    by = rawdat_mvmr$SBP_beta,
    byse = rawdat_mvmr$SBP_se,
    snps = rawdat_mvmr$SNP
  )
  mvmr(r_input = r_input, gencov = 0, weights = 1)
}
#> 
#> Multivariable MR
#> 
#>               Estimate Std. Error    t value  Pr(>|t|)
#> exposure1 -0.031003996 0.01302925 -2.3795686 0.0186526
#> exposure2  0.006039167 0.01029181  0.5867933 0.5582678
#> 
#> Residual standard error: 2.209 on 143 degrees of freedom
#> 
#> F-statistic: 3.27 on 2 and 143 DF, p-value: 0.0408
#> 
#> ------------------------------
#> Q-Statistics for instrument strength:
#> 
#>   exposure1 exposure2
#> Q  9739.922  10401.49
#> 
#> ------------------------------
#> Q-Statistic for instrument validity:
#> 
#> 695.5924 on 142 DF , p-value: 7.338e-74
```
