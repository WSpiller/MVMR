# MVMR

<!-- badges: start -->
[![R-CMD-check](https://github.com/WSpiller/MVMR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WSpiller/MVMR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Installation

MVMR can be installed from the [MRCIEU R-Universe](https://mrcieu.r-universe.dev/) with

```r
install.packages("MVMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
```

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

The corresponding paper has been published in Statistics in Medicine, and can be accessed at:

[Testing and correcting for weak and pleiotropic instruments in two-sample multivariable Mendelian randomization. Statistics in Medicine, 2021. [Internet]. 2021;sim.9133.]( https://doi.org/10.1002/sim.9133)

## License

This project is licensed under GNU GPL v3.
