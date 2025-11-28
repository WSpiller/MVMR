# phenocov_mvmr

Uses an external phenotypic covariance matrix and summary data to
estimate covariance matrices for estimated effects of individual genetic
variants on each exposure. The phenotypic covariance matrix should be
constructed using standardised phenotype measures. The function returns
a number of covariance matrices equal to the number of SNPs, where SNP
and row numbers reference ordered exposures.

## Usage

``` r
phenocov_mvmr(pcor, seBXGs)
```

## Arguments

- pcor:

  A phenotypic correlation matrix using exposures, constructed using
  individual level exposure data. Columns should be ordered by exposure
  so as to match
  [`format_mvmr()`](https://wspiller.github.io/MVMR/reference/format_mvmr.md).

- seBXGs:

  A matrix containing standard errors corresponding in relation to the
  gene-exposure association for each SNP.

## Value

A list of covariance matrices with respect to each genetic variant,
retaining the ordering in `seBXGs`

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
phenocov_mvmr(pcor, summarydata[,c(3,4)])
} # }
```
