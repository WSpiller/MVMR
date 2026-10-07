#' pleiotropy_mvmr
#'
#' Calculates modified form of Cochran's Q statistic measuring heterogeneity in causal effect estimates obtained using each genetic variant. Observed heterogeneity is indicative of a violation of the exclusion restriction assumption in MR (validity), which can result in biased effect estimates.
#' The function takes a formatted dataframe as an input, obtained using the function [`format_mvmr()`]. Additionally, covariance matrices
#'  for estimated effects of individual genetic variants on each exposure can also be provided. These can be estimated using external data by
#'  applying the [`snpcov_mvmr()`] or [`phenocov_mvmr()`] functions, are input manually. The function returns a dataframe including the conditional
#'  Q-statistic for instrument validity, and a corresponding P-value.
#'
#'  By default the Q-statistic is evaluated at the causal effect estimates that minimise it and compared with a chi-squared distribution
#'  on L - K degrees of freedom, where L is the number of genetic variants and K the number of exposures, as in Section 3.2 of
#'  Sanderson, Spiller and Bowden (2021). This test does not over-reject when the instruments are weak, although it can then be
#'  conservative.
#'  Setting \code{estimator = "ivw"} instead evaluates the Q-statistic at the inverse variance weighted (IVW) estimate, as in earlier
#'  versions of the package. As the Q-statistic at the IVW estimate can be no smaller than its minimum, that test may over-reject,
#'  particularly when the instruments are conditionally weak.
#'
#'
#' @param r_input A formatted data frame using the [`format_mvmr()`] function or an object of class `MRMVInput` from [`MendelianRandomization::mr_mvinput()`]
#' @param gencov Calculating heterogeneity statistics requires the covariance between the effect of the genetic variants on each exposure to be known. This can either be estimated from individual level data, be assumed to be zero, or fixed at zero using non-overlapping samples of each exposure GWAS. A value of \code{0} is used by default.
#' @param estimator The causal effect estimates at which the Q-statistic is evaluated: \code{"qmin"} (the default) for the estimates
#' that minimise the Q-statistic, or \code{"ivw"} for the inverse variance weighted estimates.
#'
#' @return A Q-statistic for instrument validity and the corresponding p-value
#'
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#' @references Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2019, 48, 3, 713--727. \doi{10.1093/ije/dyy262}
#'
#' Sanderson, E., Spiller, W., and Bowden, J., Testing and correcting for weak and pleiotropic instruments in two-sample multivariable Mendelian randomization. Statistics in Medicine, 2021, 40, 25, 5434--5452. \doi{10.1002/sim.9133}
#' @export
#' @examples
#' \dontrun{
#' pleiotropy_mvmr(r_input, covariances)
#' }
#'

pleiotropy_mvmr <- function(r_input, gencov = 0, estimator = c("qmin", "ivw")) {
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

  #gencov is the covariance between the effect of the genetic variants on each exposure.
  #By default it is set to 0.

  if (!is.list(gencov) && gencov == 0) {
    warning(
      "Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0."
    )
  }

  estimator <- match.arg(estimator)

  #Determine the number of exposures included in the model

  exp.number <- length(names(r_input)[-c(1, 2, 3)]) / 2
  nsnp <- nrow(r_input)

  betas <- as.matrix(r_input[, 4:(3 + exp.number)])
  sebetas <- as.matrix(r_input[, (exp.number + 4):length(r_input)])

  # Per-SNP covariance matrices of the exposure effect estimates. When gencov is
  # a scalar these are diagonal with the squared standard errors; when gencov is
  # a list, one covariance matrix is supplied per SNP.
  if (is.list(gencov)) {
    covlist <- gencov
  } else {
    covlist <- lapply(seq_len(nsnp), function(l) diag(sebetas[l, ]^2, exp.number))
  }

  ########################
  ## Instrument Validity #
  ########################

  # Q-statistic for instrument validity at causal effects b, with
  # sigma^2_A = se(betaYG)^2 + t(b) %*% Sigma_l %*% b
  Qstat <- function(b) {
    sigma2A <- r_input[, 3]^2 +
      vapply(covlist, function(S) drop(t(b) %*% S %*% b), numeric(1))
    sum((1 / sigma2A) * (r_input[, 2] - betas %*% b)^2)
  }

  # Fit the IVW MVMR model
  bivw <- stats::lm.wfit(betas, r_input[, 2], 1 / r_input[, 3]^2)$coefficients

  if (estimator == "qmin") {
    # Minimise the Q-statistic over the causal effects, starting from the IVW
    # estimate, as in Section 3.2 of Sanderson, Spiller and Bowden (2021)
    Q_valid <- stats::optim(bivw, Qstat, method = "BFGS", control = list(reltol = 1e-12))$value
  } else {
    Q_valid <- Qstat(bivw)
  }

  #Calculates p_value for instrument validity on L - K degrees of freedom
  Q_chiValid <- stats::pchisq(Q_valid, nsnp - exp.number, lower.tail = FALSE)

  ##########
  # Output #
  ##########

  cat("Q-Statistic for instrument validity:")

  cat("\n")

  cat(
    Q_valid,
    "on",
    nsnp - exp.number,
    "DF",
    ",",
    "p-value:",
    Q_chiValid
  )

  cat("\n")

  multi_return <- function() {
    Out_list <- list("Qstat" = Q_valid, "Qpval" = Q_chiValid)

    #Defines class of output object

    return(Out_list)
  }
  OUT <- multi_return()
}
