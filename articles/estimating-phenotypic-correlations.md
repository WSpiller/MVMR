# Estimating phenotypic correlations

## Covariance and Phenotypic correlations in MVMR

Summary data multivariable Mendelian randomization (MVMR) requires
estimation of the covariance between the effect of each instrument on
each exposure included in the MVMR model for effective strength testing,
heterogeneity estimation, and pleiotropy robust estimation (Sanderson et
al. (2021)).

This covariance is dependent on the level of sample overlap between the
GWAS used for each exposure and is zero if there is no sample overlap.
In that case `gencov = 0` can be used in the
[`pleiotropy_mvmr()`](https://wspiller.github.io/MVMR/reference/pleiotropy_mvmr.md)
and
[`strength_mvmr()`](https://wspiller.github.io/MVMR/reference/strength_mvmr.md)
functions.

If there is any level of sample overlap between the samples used for the
exposure GWAS then accurate estimation of the Q-statistic, conditional
F-statistics and robust estimation with weak or pleiotropic SNPs (under
the assumption of balanced pleiotropy) relies on an estimation of the
covariance between the errors in the GWAS results between the exposures
for each SNP. In a MVMR estimation with two exposures X_1 and X_2 the
estimated effects \hat{\beta}\_1 and \hat{\beta}\_2 are distributed as;

\left\[\begin{matrix} \hat{\beta}\_{X_1}\\ \hat{\beta}\_{X_2}
\end{matrix}\right\] \sim N\left(\left\[ \begin{matrix} \beta\_{X_1}\\
\beta\_{X_2} \end{matrix}\right\],\left\[ \begin{matrix} \sigma\_{x_1}^2
& \sigma\_{X_1,X_2}\\ \sigma\_{X_1,X_2} & \sigma\_{X_2}^2
\end{matrix}\right\]\right)

The estimate of the covariance \hat{\sigma}\_{X_1,X_2} is not available
from standard GWAS summary statistics and so needs to be obtained
elsewhere. This can be done in a couple of different ways:

1.  *Estimating the covariances from individual level data.*

    If there is complete sample overlap between the exposure samples,
    and the individual level data on which the GWAS summary statistics
    were generated is available, then the covariance required can be
    estimated directly from the individual level data by obtaining the
    covariance (\hat{\sigma}\_{12,i}) between the residuals
    \varepsilon\_{1i} and \varepsilon\_{2i} from the regressions;

    x_1 = \alpha_1 + \gamma\_{1i} G_i + \varepsilon\_{1i}

    and

    x_2 = \alpha_2 + \gamma\_{2i}G_i + \varepsilon\_{2i}

    These estimates can be obtained from the
    [`snpcov_mvmr()`](https://wspiller.github.io/MVMR/reference/snpcov_mvmr.md)
    function. Estimation in this way does not account for population
    structure controls included in the original GWAS which may affect
    the estimates of the covariance of the resulting residuals.

2.  ***RECOMMENDED:** Estimate* \lambda

    The distribution of the estimated effects can be rewritten as;

    \left\[\begin{matrix} \hat{\beta}\_{X_1}\\ \hat{\beta}\_{X_2}
    \end{matrix}\right\] \sim N\left(\left\[ \begin{matrix}
    \beta\_{X_1}\\ \beta\_{X_2} \end{matrix}\right\],\left\[
    \begin{matrix} \sigma\_{x_1}^2 & \lambda\sigma\_{X_1}\sigma\_{X_2}\\
    \lambda\sigma\_{X_1}\sigma\_{X_2} & \sigma\_{X_2}^2
    \end{matrix}\right\]\right)

    where;

    \lambda = \frac{n\_{\text{overlap}}\rho}{\sqrt{n\_{X_1}n\_{X_2}}}

    \rho is the correlation between X_1 and X_2, n\_{X_1} is the sample
    size for X_1, n\_{X_2} is the sample size for X_2 and
    n\_{\text{overlap}} is the sample size for the overlap between the
    samples for X_1 and X_2.

    Forde et al. (2026) provide a function to estimate this value of
    \lambda using the null SNPs for both X_1 and X_2 from the GWAS
    summary statistics. It therefore does not require access to any
    individual level data. This method works best if the full summary
    statistics are available. The function is `est_lambda()` in the
    [mr.simss package](https://amandaforde.github.io/mr.simss/).

    This value of \lambda should be used in the
    [`phenocov_mvmr()`](https://wspiller.github.io/MVMR/reference/phenocov_mvmr.md)
    function as the `pcor` argument to calculate the values
    \hat{\sigma}\_{X_1,X_2} for use in the required functions in the
    MVMR package.

3.  *Obtain* \rho *from phenotypic data.*

    If there is full sample overlap between the GWAS for X_1 and X_2
    then

    \lambda = \rho

    and so an estimate of \rho obtained from the phenotypic data that
    was used to generate the GWAS summary statistics can be used in
    place of \lambda. This requires access to the individual level data
    that was used to estimate the GWAS summary statistics.

## References

Forde, A., G. Hemani, and J. Ferguson. 2026. “Simulated sample splitting
approach to address biases due to instrument selection and participant
overlap in two-sample Mendelian Randomization studies.” *PLOS Genetics*
22 (5): e1011949. <https://doi.org/10.1371/journal.pgen.1011949>.

Sanderson, E., W. Spiller, and J. Bowden. 2021. “Testing and correcting
for weak and pleiotropic instruments in two-sample multivariable
Mendelian randomization.” *Statistics in Medicine* 40 (25): 5434–52.
<https://doi.org/10.1002/sim.9133>.
