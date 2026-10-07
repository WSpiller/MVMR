#' qhet_mvmr
#'
#' Fits a multivariable Mendelian randomization model adjusting for weak instruments. The functions requires a formatted dataframe using the [`format_mvmr()`] function, as well a phenotypic correlation matrix \code{pcor}. This should be obtained from individual level
#' phenotypic data, or constructed as a correlation matrix where correlations have previously been reported. Confidence intervals are calculated using a non-parametric bootstrap.
#' By default, standard errors are not produced but can be calculated by setting \code{se = TRUE}. The number of bootstrap iterations is specified using the \code{iterations} argument.
#' Note that calculating confidence intervals at present can take a substantial amount of time.
#'
#' Effects are estimated by minimising a Q-statistic that includes an additive heterogeneity parameter tau-squared. Following
#' Sanderson, Spiller and Bowden (2021), tau-squared is chosen so that the minimised Q-statistic equals L - K, where L is the
#' number of genetic variants and K the number of exposures. Tau-squared is constrained to be non-negative, and is fixed at zero
#' when the minimised Q-statistic without heterogeneity is already no greater than L - K.
#'
#' @param r_input A formatted data frame using the [`format_mvmr()`] function or an object of class `MRMVInput` from [`MendelianRandomization::mr_mvinput()`]
#' @param pcor A phenotypic correlation matrix including the correlation between each exposure included in the MVMR analysis.
#' @param CI Indicates whether 95 percent confidence intervals should be calculated using a non-parametric bootstrap.
#' @param iterations Specifies number of bootstrap iterations for calculating 95 percent confidence intervals.
#' @param ncores Number of cores to use for parallel processing in bootstrap. Default is `parallelly::availableCores(omit = 1)`. On Windows, this is automatically set to 1 regardless of user input. It is recommended to only set this to a maximum of `parallelly::availableCores(omit = 1)`.
#'
#' @return A dataframe containing effect estimates with respect to each exposure.
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#' @references Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2019, 48, 3, 713--727. \doi{10.1093/ije/dyy262}
#'
#' Sanderson, E., Spiller, W., and Bowden, J., Testing and correcting for weak and pleiotropic instruments in two-sample multivariable Mendelian randomization. Statistics in Medicine, 2021, 40, 25, 5434--5452. \doi{10.1002/sim.9133}
#' @export
#' @examples
#' \dontrun{
#' qhet_mvmr(r_input, pcor, CI = TRUE, iterations = 1000)
#' }

qhet_mvmr <- function(r_input, pcor, CI, iterations, ncores = parallelly::availableCores(omit = 1)) {
  # convert MRMVInput object to mvmr_format
  if ("MRMVInput" %in% class(r_input)) {
    r_input <- mrmvinput_to_mvmr_format(r_input)
  }

  # Perform check that r_input has been formatted using format_mvmr function
  if (
    !("mvmr_format" %in%
      class(r_input))
  ) {
    stop(
      'The class of the data object must be "mvmr_format", please resave the object with the output of format_mvmr().'
    )
  }

  warning("qhet_mvmr() is currently undergoing development.")

  if (missing(CI)) {
    CI <- FALSE
    warning("95 percent confidence interval not calculated")
  }

  if (missing(iterations)) {
    iterations <- 1000
    warning("Iterations for bootstrap not specified. Default = 1000")
  }

  # Check if on Windows and adjust ncores accordingly
  if (.Platform$OS.type == "windows" && ncores > 1) {
    warning("Multi-core processing is not supported on Windows. Setting ncores to 1.")
    ncores <- 1
  }

  if (ncores > parallelly::availableCores(omit = 1)) {
    stop('You have set the number of cores greater than the number available on the machine minus one. We recommend setting this to a maximum of parallelly::availableCores(omit = 1).')
  }

  exp.number <- length(names(r_input)[-c(1, 2, 3)]) / 2

  Qtemp <- function(r_input, pcor) {
    exp.number <- length(names(r_input)[-c(1, 2, 3)]) / 2
    stderr <- as.matrix(r_input[, (exp.number + 4):length(r_input)])
    correlation <- pcor
    covlist <- lapply(seq_len(nrow(r_input)), function(l) correlation * outer(stderr[l, ], stderr[l, ]))
    gammahat <- r_input$betaYG
    segamma <- r_input$sebetaYG
    pihat <- as.matrix(r_input[, c(4:(3 + exp.number))])

    nsnp <- nrow(r_input)

    # Q-statistic for causal effects b with additive heterogeneity tau2
    Qstat <- function(b, tau2) {
      w <- segamma^2 + vapply(covlist, function(cv) drop(t(b) %*% cv %*% b), numeric(1)) + tau2
      sum((1 / w) * ((gammahat - pihat %*% b)^2))
    }

    # Minimise the Q-statistic over b for a given tau2, starting from the IVW estimate
    bstart <- stats::lm.wfit(pihat, gammahat, 1 / segamma^2)$coefficients
    Qmin <- function(tau2) {
      stats::optim(bstart, Qstat, tau2 = tau2, method = "BFGS", control = list(reltol = 1e-12))
    }

    # Choose tau2 >= 0 so that the minimised Q-statistic equals its expectation
    # L - K under its chi-squared null distribution (Sanderson, Spiller and
    # Bowden, 2021). If the minimised Q-statistic is already no greater than
    # L - K there is no excess heterogeneity and tau2 is fixed at 0.
    target <- nsnp - exp.number
    tau_i <- 0
    if (Qmin(0)$value > target) {
      # Search upwards from the scale of the outcome variances for an upper bound
      upper <- stats::median(segamma^2)
      while (Qmin(upper)$value > target) {
        upper <- upper * 10
      }
      tau_i <- stats::uniroot(
        function(tau2) Qmin(tau2)$value - target,
        interval = c(0, upper),
        tol = upper * 1e-8
      )$root
    }

    limlhets <- Qmin(tau_i)$par

    Effects <- limlhets
    Effects <- data.frame(Effects)
    names(Effects) <- "Effect Estimates"
    for (i in 1:exp.number) {
      rownames(Effects)[i] <- paste("Exposure", i, sep = " ")
    }

    return(Effects)
  }

  if (!CI) {
    res <- Qtemp(r_input, pcor)
  }

  if (CI) {
    bootse <- function(data, indices) {
      bres <- Qtemp(data[indices, ], pcor)[, 1]

      return(bres)
    }

    if (ncores > 1) {
      b.results <- boot::boot(data = r_input, statistic = bootse, R = iterations, parallel = "multicore", ncpus = ncores)
    } else {
      b.results <- boot::boot(data = r_input, statistic = bootse, R = iterations)
    }

    lcb <- NULL
    ucb <- NULL
    ci <- NULL

    for (i in 1:exp.number) {
      boot_ci <- boot::boot.ci(b.results, type = "bca", index = i)
      lcb[i] <- round(boot_ci$bca[4], digits = 3)
      ucb[i] <- round(boot_ci$bca[5], digits = 3)

      ci[i] <- paste(lcb[i], ucb[i], sep = "-")
    }

    res <- data.frame(b.results$t0, ci)

    names(res) <- c("Effect Estimates", "95% CI")
    for (i in 1:exp.number) {
      rownames(res)[i] <- paste("Exposure", i, sep = " ")
    }
  }

  return(res)
}
