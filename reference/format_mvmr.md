# format_mvmr

Reads in summary data. Checks and organises columns for use in
calculating multivariable Mendelian Randomization analyses. Where
variant IDs are not provided, a vector is generated for variant
identification.

## Usage

``` r
format_mvmr(BXGs, BYG, seBXGs, seBYG, RSID)
```

## Arguments

- BXGs:

  A matrix containing beta-coefficient values for genetic associations
  with the each exposure. Columns should indicate exposure number, with
  rows representing estimates for a given genetic variant.

- BYG:

  A numeric vector of beta-coefficient values for genetic associations
  with the outcome.

- seBXGs:

  A matrix containing standard errors corresponding to the matrix of
  beta-coefficients `BXGs`.

- seBYG:

  A numeric vector of standard errors corresponding to the
  beta-coefficients `BYG`.

- RSID:

  A vector of names for genetic variants included in the analysis. If
  variant IDs are not provided a vector of ID numbers will be generated.

## Value

A formatted data frame of class `mvmr_format`.

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
names(r_input)
#> [1] "SNP"      "betaYG"   "sebetaYG" "betaX1"   "betaX2"   "sebetaX1" "sebetaX2"
class(r_input)
#> [1] "data.frame"  "mvmr_format"
```
