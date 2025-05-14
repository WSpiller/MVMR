#' ivw_mvmr
#'
#' Fits an IVW multivariable Mendelian randomization model using first order weights.
#'
#' @param r_input A formatted data frame using the [`format_mvmr()`] function or an object of class `MRMVInput` from [`MendelianRandomization::mr_mvinput()`]
#' @param gencov Calculating heterogeneity statistics requires the covariance between the effect of the genetic variants on each exposure to be known. This can either be estimated from individual level data, be assumed to be zero, or fixed at zero using non-overlapping samples of each exposure GWAS. A value of \code{0} is used by default.
#'
#' @return An dataframe containing MVMR results, including estimated coefficients, their standard errors, t-statistics, and corresponding (two-sided) p-values.
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#' @references Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2019, 48, 3, 713--727. \doi{10.1093/ije/dyy262}
#' @export
#' @examples
#' r_input <- format_mvmr(
#'   BXGs = rawdat_mvmr[, c("LDL_beta", "HDL_beta")],
#'   BYG = rawdat_mvmr$SBP_beta,
#'   seBXGs = rawdat_mvmr[, c("LDL_se", "HDL_se")],
#'   seBYG = rawdat_mvmr$SBP_se,
#'   RSID = rawdat_mvmr$SNP
#' )
#' ivw_mvmr(r_input)
# Define IVW Multivariable MR function: This takes the formatted dataframe from
# the format_MVMR function as an input, as well as the covariance between the effect of the
# genetic variants on each exposure.

ivw_mvmr<-function(r_input,gencov=0){

  # convert MRMVInput object to mvmr_format
  if ("MRMVInput" %in% class(r_input)) {
    r_input <- mrmvinput_to_mvmr_format(r_input)
  }

  # Perform check that r_input has been formatted using format_mvmr function
  if(!("mvmr_format" %in%
       class(r_input))) {
    stop('The class of the data object must be "mvmr_format", please resave the object with the output of format_mvmr().')
  }

  if(!is.list(gencov) && gencov == 0) {
    warning("Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.")
  }

  #If weights is missing, first order weights are used by default.


  # Inverse variance weighting is used.

    Wj<-1/r_input[,3]^2

  #Determine the number of exposures included in the model

  exp.number<-length(names(r_input)[-c(1,2,3)])/2

  #Fit the IVW MVMR model

  A_sum<-summary(stats::lm(stats::as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,weights=Wj,data=r_input))

  A<-summary(stats::lm(stats::as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,weights=Wj,data=r_input))$coef

  #Rename the regressors for ease of interpretation
  for(i in 1:exp.number){
    dimnames(A)[[1]][i]<- paste0("exposure",i,collapse="")
  }


  ##########
  # Output #
  ##########

  # Print a few summary elements that are common to both lm and plm model summary objects
  cat("\n")

  cat("Multivariable MR\n")

  cat("\n")

  print(A)

  cat("\nResidual standard error:", round(A_sum$sigma,3), "on", A_sum$df[2], "degrees of freedom")

  cat("\n")

  cat("\n")

  cat("\n")


  return(A)

}
