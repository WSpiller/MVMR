# snpcov_mvmr

Uses individual level genetic and exposure data to generate covariance
matrices for estimated effects of individual genetic variants on each
exposure. The function returns a number of covariance matrices equal to
the number of SNPs, where SNP and row numbers reference ordered
exposures.

## Usage

``` r
snpcov_mvmr(Gs, Xs)
```

## Arguments

- Gs:

  A matrix or dataframe containing genetic instrument measures. Columns
  should indicate genetic variant number, with rows representing an
  observed measure of the genetic variant.

- Xs:

  A matrix or dataframe containing exposure measures. Columns should
  indicate exposure number, with rows representing an observed measure
  for the given exposure.

## Value

A list of covariance matrices with respect to each genetic variant,
retaining the ordering in `Gs`

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
snpcov_mvmr(data[,1:10], data[,11:13])
} # }
```
