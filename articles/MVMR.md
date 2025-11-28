# Multivariable MR Tutorial

## Overview

Multivariable Mendelian Randomisation (MVMR) is a form of instrumental
variable analysis which estimates the direct effect of multiple
exposures on an outcome using genetic variants as instruments. The
`MVMR` R package facilitates estimation of causal effects using MVMR, as
well as including a range of sensitivity analyses evaluating the
underlying assumptions of the approach. The methods included in `MVMR`
originate from Sanderson, Spiller, and Bowden (2021), available
[here](https://doi.org/10.1002/sim.9133).

### Workflow

Fitting and interpreting MVMR models can be achieved by following the 5
steps given below:

1.  Obtain data
2.  Format data
3.  Assess instrument strength
4.  Assess horizontal pleiotropy
5.  Estimate causal effects

Each of these steps are shown in the figure below, highlighting the R
function used for each step.

![Figure showing the workflow of an MVMR analysis.](png/workflow.png)

The workflow of an MVMR analysis.

## Step 1: Obtain summary data

The following information is necessary to estimate causal effects using
MVMR:

1.  Gene-exposure associations for each variant selected as an
    instrument for any exposure.

2.  Corresponding standard errors for the gene-exposure associations.

3.  Gene-outcome associations for each instrument.

4.  Corresponding standard errors for the gene-outcome associations.

The data frame `rawdat_mvmr`, included in the `MVMR` package shows an
example of such data obtained from MRBase. When data is extracted from
MRBase or using the TwoSampleMR R package, the
[`mrmvinput_to_mvmr_format()`](https://wspiller.github.io/MVMR/reference/mrmvinput_to_mvmr_format.md)
function can be used to convert the data to the `MVMR` data format. This
specifically requires gene-exposure associations for all included SNPs,
gene-outcome associations, and corresponding standard errors. In this
case, low-density lipoprotein cholesterol (LDL-C), high-density
lipoprotein cholesterol (HDL-C), and triglycerides (Trg) have been
selected as exposures, while systolic blood pressure (SBP) is the
outcome of interest. Here the suffix `_beta` is used to denote
association estimates, while `_se` denotes standard errors. Please note
that the `MVMR` can take an arbirtary number of exposures (greater than
1), and that three exposures have been selected purely for illustration.

The first 6 rows of `rawdat_mvmr` are:

``` r
library(MVMR)
head(rawdat_mvmr)
#>   LDL_beta HDL_beta Trg_beta LDL_se HDL_se Trg_se    SBP_beta     SBP_se
#> 1  -0.0270   0.0182   0.0228 0.0046 0.0050 0.0045 -0.00426935 0.00280123
#> 2   0.1179   0.0020   0.0379 0.0038 0.0042 0.0039  0.00110389 0.00227690
#> 3  -0.0269   0.0079   0.0000 0.0039 0.0042 0.0038 -0.01317370 0.00222743
#> 4   0.0081  -0.0508   0.0085 0.0064 0.0069 0.0062 -0.00111303 0.00374264
#> 5   0.0191   0.0103  -0.0270 0.0034 0.0036 0.0033 -0.00854986 0.00207701
#> 6   0.0043   0.0320   0.0109 0.0038 0.0040 0.0037  0.00472509 0.00237133
#>          SNP
#> 1 rs10019888
#> 2 rs10468017
#> 3  rs1047891
#> 4 rs10490626
#> 5 rs10761762
#> 6 rs10832962
```

Note that the final column `SNP` contains the rsid numbers for each
genetic variant. These are not necessary for conducting MVMR, but assist
in follow-up analyses. Summary data for LDL-C, HDL-C, and Triglycerides
originate from [GLGC](https://www.nature.com/articles/ng.2797), while
SBP data was obtained using [UK
Biobank](https://www.nature.com/articles/ng.3768).

### Estimating pairwise covariances between SNP associations

The MVMR approach requires pairwise covariances between an instrument
and pairs of exposures to be known across all SNPs for testing and
sensitivity analyses. However, this is often not reported in published
GWAS analyses. Before continuing with MVMR it is therefore **necessary**
to select one of the following three solutions:

1.  Estimate the covariance terms using individual level data

    If individual level data is available from which the GWAS summary
    estimates were obtained, the
    [`snpcov_mvmr()`](https://wspiller.github.io/MVMR/reference/snpcov_mvmr.md)
    function can be used to calculate the necessary covariance terms.

2.  Estimate the covariance terms using phenotypic correlation between
    exposures. an estimate of the correlation between the (phenotypic)
    exposures is available, the
    [`phenocov_mvmr()`](https://wspiller.github.io/MVMR/reference/phenocov_mvmr.md)
    function can be used to provide an approximation for the necessary
    covariance terms. This function takes the phenotypic correlation
    between the exposures and the standard error of the SNP-exposure
    betas as inputs.

3.  Obtain gene-exposure associations from non-overlapping samples.

    If gene-exposure associations are estimated in seperate
    non-overlapping samples, then the covariances will be zero by
    design. It is therefore not necessary to calculate the set of
    covariances, although this approach can be difficult to apply due to
    a lack of suitable sources of data.

When the necessary data are provided to the
[`snpcov_mvmr()`](https://wspiller.github.io/MVMR/reference/snpcov_mvmr.md)
or
[`phenocov_mvmr()`](https://wspiller.github.io/MVMR/reference/phenocov_mvmr.md)
functions, a set of covariance matrices will be produced equal to the
number of SNPs used in estimation. By saving this output as an object,
it is possible to use this information in downstream senstivity analyses
and assumption testing. As the
[`phenocov_mvmr()`](https://wspiller.github.io/MVMR/reference/phenocov_mvmr.md)
function requires gene-exposure standard errors, it can be useful to
estimate the covariance matrices after initially formatiing the data. An
illustrative example is provided in step 6, creating an object `Xcovmat`
using the
[`phenocov_mvmr()`](https://wspiller.github.io/MVMR/reference/phenocov_mvmr.md)
function.

## Step 2: Format summary data

Downstream functions in the `MVMR` package rely upon prior formatting of
raw summary data using the
[`format_mvmr()`](https://wspiller.github.io/MVMR/reference/format_mvmr.md)
function. Specifically,
[`format_mvmr()`](https://wspiller.github.io/MVMR/reference/format_mvmr.md)
checks and organises summary data columns for use in MVMR analyses. The
[`format_mvmr()`](https://wspiller.github.io/MVMR/reference/format_mvmr.md)
function takes the following arguments:

- `BXGs`: A subset containing beta-coefficient values for genetic
  associations with each exposure. Columns should indicate exposure
  number, with rows representing estimates for a given genetic variant.
- `BYG`: A numeric vector of beta-coefficient values for genetic
  associations with the outcome.
- `seBXGs`: A subset containing standard errors corresponding to the
  subset of beta-coefficients `BXGs`.
- `seBYG`: A numeric vector of standard errors corresponding to the
  beta-coefficients `BYG`.
- `RSID`: A vector of names for genetic variants included in the
  analysis. If variant IDs are not provided (`RSID = "NULL"`), a vector
  of ID numbers will be generated.

Using the previous data `rawdat.mvmr`, we can format the data using the
following command:

``` r
F.data <- format_mvmr(
  BXGs = rawdat_mvmr[, c(1, 2, 3)],
  BYG = rawdat_mvmr[, 7],
  seBXGs = rawdat_mvmr[, c(4, 5, 6)],
  seBYG = rawdat_mvmr[, 8],
  RSID = rawdat_mvmr[, 9]
)
head(F.data)
#>          SNP      betaYG   sebetaYG  betaX1  betaX2  betaX3 sebetaX1 sebetaX2
#> 1 rs10019888 -0.00426935 0.00280123 -0.0270  0.0182  0.0228   0.0046   0.0050
#> 2 rs10468017  0.00110389 0.00227690  0.1179  0.0020  0.0379   0.0038   0.0042
#> 3  rs1047891 -0.01317370 0.00222743 -0.0269  0.0079  0.0000   0.0039   0.0042
#> 4 rs10490626 -0.00111303 0.00374264  0.0081 -0.0508  0.0085   0.0064   0.0069
#> 5 rs10761762 -0.00854986 0.00207701  0.0191  0.0103 -0.0270   0.0034   0.0036
#> 6 rs10832962  0.00472509 0.00237133  0.0043  0.0320  0.0109   0.0038   0.0040
#>   sebetaX3
#> 1   0.0045
#> 2   0.0039
#> 3   0.0038
#> 4   0.0062
#> 5   0.0033
#> 6   0.0037
```

In the above code we have provided the numbered columns for each
argument. For example, `BXGs = rawdat.mvmr[, c(1, 2, 3)]` indicates that
columns 1, 2, and 3 are the association estimates for exposures 1, 2,
and 3. It is important to note that standard error columns `seBXGs`
should be input in the same order as BXGs to ensure the correct matching
of association estimates with corresponding standard errors.

In subsequent steps, each exposure is numbered such that `X1`, `X2`, and
`X3` are the first, second, and third entries in the
`BXGs = rawdat.mvmr[, c(1, 2, 3)]` argument.

## Step 3: Test for weak instruments

In univariate two-sample summary MR, genetic variants selected as
instruments are required to be strongly associated with their
corresponding exposure. This is quantified by regressing the exposure
upon each instrument, and evaluating conditional dependence using the
F-statistic for the instrument. Conventionally, a F-statistic greater
than 10 is used as a threshold for sufficient instrument strength,
representing a 10% relative bias towards the null in the two-sample MR
setting.

Multivariable MR relies upon an extension of this assumption, requiring
instruments to be strongly associated with their corresponding exposure
conditioning on the remaining included exposures. Conditional instrument
strength is quantified by a conditional F-statistic which has the same
distribution as the univariate F-statistic. Consequently, the same
conventional instrument strength threshold of 10 can be used.

Further details are available [here](https://doi.org/10.1002/sim.9133).

The
[`strength_mvmr()`](https://wspiller.github.io/MVMR/reference/strength_mvmr.md)
function is used to evaluate instrument strength in the MVMR setting.
The function contains two arguments:

- `r_input`: A formatted data frame created using the
  [`format_mvmr()`](https://wspiller.github.io/MVMR/reference/format_mvmr.md)
  function or an object of class `MRMVInput` from the
  [`mr_mvinput()`](https://rdrr.io/pkg/MendelianRandomization/man/mr_mvinput.html)
  function in the MendelianRandomization package.
- `gencov`: A variance-covariance matrix for the effect of the genetic
  variants on each exposure. This is obtained from either
  [`snpcov_mvmr()`](https://wspiller.github.io/MVMR/reference/snpcov_mvmr.md),
  [`phenocov_mvmr()`](https://wspiller.github.io/MVMR/reference/phenocov_mvmr.md),
  or set to zero when omitted.

**Note**: The
[`strength_mvmr()`](https://wspiller.github.io/MVMR/reference/strength_mvmr.md)
function will output a warning if a variance-covariance matrix is not
provided. Please see Step 1 for further information.

Continuing with the previous example, we can evaluate the conditional
strength of the instruments for each exposure using the following
command

``` r
sres <- strength_mvmr(r_input = F.data, gencov = 0)
#> Warning in strength_mvmr(r_input = F.data, gencov = 0): Covariance between
#> effect of genetic variants on each exposure not specified. Fixing covariance at
#> 0.
#> 
#> Conditional F-statistics for instrument strength
#> 
#>             exposure1 exposure2 exposure3
#> F-statistic  46.33671  67.80463  38.80184
```

In this case the set of instruments is sufficiently strong for MVMR
estimation using the conventional F-statistic threshold of 10. However,
note that we have manually set `gencov` to zero, which would likely not
be appropriate given each SNP-exposure estimate was obtained from the
same sample. Using a random phenotypic correlation matrix, conditional
F-statistics can be calculated as

``` r
mvmrcormatrix <- matrix(c(1, -0.1, -0.05, -0.1, 1, 0.2, -0.05, 0.2, 1), nrow = 3, ncol = 3)
Xcovmat <- phenocov_mvmr(mvmrcormatrix, F.data[, 7:9])
sres2 <- strength_mvmr(r_input = F.data, gencov = Xcovmat)
#> 
#> Conditional F-statistics for instrument strength
#> 
#>             exposure1 exposure2 exposure3
#> F-statistic  48.20993  69.55193  39.77326
```

## Step 4: Test for horizontal pleiotropy using conventional Q-statistic estimation

Horizontal pleiotropy can be evaluated using a modified form of
Cochran’s Q statistic with respect to differences in MVMR estimates
across the set of instruments. In this case, observed heterogeneity is
indicative of a violation of the exclusion restriction assumption in MR
(validity), which can result in biased effect estimates.

Importantly, weak instruments can increase the false positive rate for
pleiotropy detection, as heterogeneity in effect estimates due to weak
instrument bias is conflated with heterogeneity as a result of
pleiotropic bias. As a correction it is possible to estimate
heterogeneity from pleiotropy through Q-statistic minimisation.

The function
[`pleiotropy_mvmr()`](https://wspiller.github.io/MVMR/reference/pleiotropy_mvmr.md)
can be used to test for heterogeneity, requiring the same arguments as
the
[`strength_mvmr()`](https://wspiller.github.io/MVMR/reference/strength_mvmr.md);
`r_input` and `gencov`.

``` r
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
#> Warning in pleiotropy_mvmr(r_input = F.data, gencov = 0): Covariance between
#> effect of genetic variants on each exposure not specified. Fixing covariance at
#> 0.
#> Q-Statistic for instrument validity:
#> 683.0807 on 141 DF , p-value: 4.880403e-72
```

And with the example covariance matrices from Step 3:

``` r
pres <- pleiotropy_mvmr(r_input = F.data, gencov = Xcovmat)
#> Q-Statistic for instrument validity:
#> 682.843 on 141 DF , p-value: 5.36533e-72
```

## Step 5: Estimate causal effects

Two MVMR estimation methods are provided in the `MVMR` package. The
first method fits an inverse variance weighted (IVW) MVMR model,
providing estimates of the direct effect of each exposure upon the
outcome. This is performed using the
[`ivw_mvmr()`](https://wspiller.github.io/MVMR/reference/ivw_mvmr.md)
function as shown below:

``` r
res <- ivw_mvmr(r_input = F.data)
#> Warning in ivw_mvmr(r_input = F.data): Covariance between effect of genetic
#> variants on each exposure not specified. Fixing covariance at 0.
#> 
#> Multivariable MR
#> 
#>               Estimate Std. Error   t value  Pr(>|t|)
#> exposure1 -0.021845061 0.01417255 -1.541364 0.1254538
#> exposure2  0.003735249 0.01033779  0.361320 0.7183973
#> exposure3  0.025572042 0.01601913  1.596344 0.1126351
#> 
#> Residual standard error: 2.197 on 142 degrees of freedom
```

In this case, the effect estimates are interpreted as the direct effects
of LDL-C (exposure 1), HDL-C (exposure 2), and Trg (exposure 3) on SBP.
Estimates are not robust to weak instruments of pleiotropic bias, and
therefore rely upon the underlying MVMR assumptions being satisfied.

## Step 6: Robust causal effect estimation.

Where the MVMR assumptions are potentially violated, specifically where
instruments are weak or exhibit pleiotropy, it is possible to obtain
more robust estimates through Q-statistic minimisation. This can be
performed using the
[`qhet_mvmr()`](https://wspiller.github.io/MVMR/reference/qhet_mvmr.md)
function.

``` r
res1 <- qhet_mvmr(F.data, mvmrcormatrix, CI = FALSE, iterations = 100)
#> Warning in qhet_mvmr(F.data, mvmrcormatrix, CI = FALSE, iterations = 100):
#> qhet_mvmr() is currently undergoing development.
res1
#>            Effect Estimates
#> Exposure 1    -0.0264865644
#> Exposure 2     0.0094372624
#> Exposure 3     0.0002575009
```

It is important to highlight that the phenotypic covariance matrix is
used as an input, and not the set of estimated covariance matrices which
previously formed the `gencov` argument. It should also be noted that as
the number of exposure and instruments increases, estimation using
[`qhet_mvmr()`](https://wspiller.github.io/MVMR/reference/qhet_mvmr.md)
may prove difficult, owing to the substantial amount of computing power
required. This can be initially relaxed by not computing 95% confidence
intervals as above.

## References

Sanderson, E., W. Spiller, and J. Bowden. 2021. “Testing and correcting
for weak and pleiotropic instruments in two-sample multivariable
Mendelian randomization.” *Statistics in Medicine* 40 (25): 5434–52.
<https://doi.org/10.1002/sim.9133>.
