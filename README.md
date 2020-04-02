# MVMR

## Installation

To install `MVMR` directly from the GitHub repository, first make sure you have the `remotes` package installed:

    install.packages("remotes")

Then the `MVMR` package can be installed using:

    library(remotes)
    install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
    
To update the package just run the `remotes::install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)` command again.

## Description

We have written the `MVMR` R package to perform multivariable Mendelian randomization analyses, including heterogeneity
statistics for assessing instrument strength and validity. The package accommodates any number of exposures greater than 1,
and is currently comprised of two functions:

1. The `format_mvmr` function is used to convert a data frame containing summary data into a set format for MVMR analyses.

2. The `mvmr` function fits an IVW multivariable Mendelian randomization model using first order weights. The function returns 
an object of class `MVMRIVW`, containing regression estimates, estimated heterogeneity as a measure of instrument strength
(`Q_strength`), and estimated heterogeneity as a measure of instrument validity (`Q_valid`). Estimation follows the method
outlined in [Sanderson et al, 2018](https://dx.doi.org/10.1093/ije/dyy262).

Multivariable Mendelian randomization is implemented in a number of R packages such as `TwoSampleMR`, and `MendelianRandomization`,
however, this package includes further sensitivity analyses leveraging information on causal effect heterogeneity across instruments.
To incorporate these features and future developments we will  continue to develop our own `MVMR` package to implement multivariable
Mendelian randomization.

## Citation

The corresponding paper has been published by the International Journal of Epidemiology, and can be accessed at:

[An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2018. [Internet]. 2018;dyy262.](https://dx.doi.org/10.1093/ije/dyy262)

## License

This project is licensed under GNU GPL v3.




