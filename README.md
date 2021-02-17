# MVMR

## Installation

To install `MVMR` directly from the GitHub repository, first make sure you have the `remotes` package installed:

```r
install.packages("remotes")
```

Then the `MVMR` package can be installed using:
```r
library(remotes)
install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
```
To update the package just run the `remotes::install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)` command again.

## Description

We have written the `MVMR` R package to perform multivariable Mendelian randomization analyses, including heterogeneity
statistics for assessing instrument strength and validity. The package accommodates any number of exposures greater than 1,
and is currently includes a range of functions for estimating causal effects, as well as assessing conditional instrument strength and pleiotropic bias. For a detailed description regarding each function, please install the package and input `vignette("MVMR")`

## Citation

The corresponding paper has been published by the International Journal of Epidemiology, and can be accessed at:

[An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2018. [Internet]. 2018;dyy262.](https://dx.doi.org/10.1093/ije/dyy262)

## License

This project is licensed under GNU GPL v3.
