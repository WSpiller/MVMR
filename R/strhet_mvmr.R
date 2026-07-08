#' strhet_mvmr
#'
#' Calculates the conditional F-statistic for assessing instrument strength in two sample summary multivariable Mendelian randomization through minimisation of Q-statistics.
#' The function takes a formatted dataframe as an input, obtained using the function [`format_mvmr()`]. Additionally, covariance matrices
#' for estimated effects of individual genetic variants on each exposure can also be provided. These can be estimated using external data by
#' applying the [`snpcov_mvmr()`] or [`phenocov_mvmr()`] functions, or are input manually. The function returns a dataframe including the conditional
#' F-statistic with respect to each exposure. A conventional F-statistic threshold of 10 is used in basic assessments of instrument strength.
#'
#' @param r_input A formatted data frame using the [`format_mvmr()`] function or an object of class `MRMVInput` from [`MendelianRandomization::mr_mvinput()`]
#' @param gencov Calculating heterogeneity statistics requires the covariance between the effect of the genetic variants on each exposure to be known. This can either be estimated from individual level data, be assumed to be zero, or fixed at zero using non-overlapping samples of each exposure GWAS. A value of \code{0} is used by default.
#'
#' @return A dataframe showing the conditional F-statistic for each exposure.
#'
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#' @references Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2019, 48, 3, 713--727. \doi{10.1093/ije/dyy262}
#' @export
#' @examples
#' \dontrun{
#' strhet_mvmr(r_input, covariances)
#' }

strhet_mvmr <- function(r_input, gencov) {
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

  if (missing(gencov)) {
    gencov <- as.numeric(0)
    warning(
      "Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0."
    )
  }

  #Determine the number of exposures included in the model

  exp.number <- length(names(r_input)[-c(1, 2, 3)]) / 2

  # Full matrices of exposure effect estimates and their standard errors.
  betas_all <- as.matrix(r_input[, 4:(3 + exp.number)])
  sebetas_all <- as.matrix(r_input[, (exp.number + 4):length(r_input)])
  nsnp <- nrow(r_input)

  # Per-SNP covariance matrices of the exposure effect estimates. When gencov is
  # a scalar (typically 0) these are diagonal with the squared standard errors;
  # when gencov is a list, one covariance matrix is supplied per SNP.
  if (is.list(gencov)) {
    covlist <- gencov
  } else {
    covlist <- lapply(
      seq_len(nsnp),
      function(l) diag(sebetas_all[l, ]^2, exp.number)
    )
  }

  #############################################
  # Generalised instrument strength het.stats #
  #############################################

  # For each exposure the delta coefficients (regressing that exposure's
  # associations on the others) are estimated by minimising the Q-statistic. The
  # variance weights sigma^2 are held fixed at the current delta within each
  # iteration and updated until convergence (iteratively reweighted least
  # squares), which keeps the per-exposure statistics well identified. Writing
  # the residual as a contrast v with -1 in the exposure's own position and delta
  # elsewhere, sigma^2_l = t(v) %*% Sigma_l %*% v.
  qminvec <- numeric(exp.number)
  for (m in seq_len(exp.number)) {
    y <- betas_all[, m]
    X <- betas_all[, -m, drop = FALSE]
    d <- rep(0, exp.number - 1)

    for (iter in seq_len(100)) {
      v <- numeric(exp.number)
      v[m] <- -1
      v[-m] <- d
      sig <- vapply(covlist, function(S) drop(t(v) %*% S %*% v), numeric(1))
      w <- 1 / sig
      d_new <- as.vector(solve(crossprod(X, w * X), crossprod(X, w * y)))
      if (max(abs(d_new - d)) < 1e-10) {
        d <- d_new
        break
      }
      d <- d_new
    }

    v <- numeric(exp.number)
    v[m] <- -1
    v[-m] <- d
    sig <- vapply(covlist, function(S) drop(t(v) %*% S %*% v), numeric(1))
    resid <- as.vector(y - X %*% d)
    qminvec[m] <- sum(resid^2 / sig) / nsnp
  }

  Q_strength <- data.frame(t(qminvec))
  names(Q_strength) <- paste0("exposure", seq_len(exp.number))
  rownames(Q_strength) <- "F-statistic"

  ##########
  # Output #
  ##########

  # Print a few summary elements that are common to both lm and plm model summary objects
  cat("\n")

  cat("Conditional F-statistics for instrument strength\n")

  cat("\n")

  print(Q_strength)

  return(Q_strength)
}
