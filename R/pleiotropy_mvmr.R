#' pleiotropy_mvmr
#'
#' Calculates modified form of Cochran's Q statistic measuring heterogeneity in causal effect estimates obtained using each genetic variant. Observed heterogeneity is indicative of a violation of the exclusion restriction assumption in MR (validity), which can result in biased effect estimates.
#' The function takes a formatted dataframe as an input, obtained using the function [`format_mvmr`]. Additionally, covariance matrices
#'  for estimated effects of individual genetic variants on each exposure can also be provided. These can be estimated using external data by
#'  applying the [`snpcov_mvmr`] or [`phenocov_mvmr`] functions, are input manually. The function returns a dataframe including the conditional
#'  Q-statistic for instrument validity, and a corresponding P-value.
#'
#'
#' @param r_input A formatted data frame using the [`format_mvmr`] function or an object of class `MRMVInput` from [`MendelianRandomization::mr_mvinput`]
#' @param gencov Calculating heterogeneity statistics requires the covariance between the effect of the genetic variants on each exposure to be known. This can either be estimated from individual level data, be assumed to be zero, or fixed at zero using non-overlapping samples of each exposure GWAS. A value of \code{0} is used by default.
#'
#' @return A Q-statistic for instrument validity and the corresponding p-value
#'
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#' @references Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2019, 48, 3, 713-727. \doi{10.1093/ije/dyy262}
#' @export
#' @examples
#' \dontrun{
#' pleiotropy_mvmr(r_input, covariances)
#' }
#'

pleiotropy_mvmr<-function(r_input,gencov=0){

  # convert MRMVInput object to mvmr_format
  if ("MRMVInput" %in% class(r_input)) {
    r_input <- mrmvinput_to_mvmr_format(r_input)
  }

  # Perform check that r_input has been formatted using format_mvmr function
  if(!("mvmr_format" %in%
       class(r_input))) {
    stop('The class of the data object must be "mvmr_format", please resave the object with the output of format_mvmr().')
  }

  #gencov is the covariance between the effect of the genetic variants on each exposure.
  #By default it is set to 0.

  if(!is.list(gencov) && gencov == 0) {
    warning("Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.")
  }

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


  #Create a subset containing only standard errors for exposure effect estimates
  sebetas<-r_input[,(exp.number + 4):length(r_input)]


  ########################
  ## Instrument Validity #
  ########################

  if(length(gencov) <2){

    # Generate Sigma^2_A values
    sigma2A<-r_input[,3]^2
    for(i in 1:exp.number){
      sigma2A<-sigma2A + (A[i]^2 * sebetas[,i]^2)
    }

    #Create a subset of exposure effect estimates
    betas<-r_input[,c(4:(3+exp.number))]

    #Generates the component of the Q statistic to be subtracted from the outcome estimates
    temp.sub2<-0
    for(i in 1:exp.number){
      temp.sub2<-temp.sub2 + (betas[,i] * A[i])
    }

    #Calculates Q statistic for instrument validity
    Q_valid<- sum((1/sigma2A)*(r_input[,2]-temp.sub2)^2)

    #Calculates p_value for instrument validity
    Q_chiValid<-stats::pchisq(Q_valid,length(r_input[,2])-exp.number-1,lower.tail = FALSE)


  }

  if(length(gencov) >2){

    # Generate Sigma^2_A values
    sigma2A<-r_input[,3]^2
    for(i in seq_along(r_input[,3])){
      sigma2A[i]<-sigma2A[i] + (t(as.matrix(A[,1])) %*% gencov[[i]]%*% as.matrix(A[,1]))
    }

    #Create a subset of exposure effect estimates
    betas<-r_input[,c(4:(3+exp.number))]

    #Generates the component of the Q statistic to be subtracted from the outcome estimates
    temp.sub2<-0
    for(i in 1:exp.number){
      temp.sub2<-temp.sub2 + (betas[,i] * A[i])
    }

    #Calculates Q statistic for instrument validity
    Q_valid<- sum((1/sigma2A)*(r_input[,2]-temp.sub2)^2)

    #Calculates p_value for instrument validity
    Q_chiValid<-stats::pchisq(Q_valid,length(r_input[,2])-exp.number-1,lower.tail = FALSE)


  }


  ##########
  # Output #
  ##########

  cat("Q-Statistic for instrument validity:")

  cat("\n")

  cat(Q_valid, "on", length(r_input[,2])-exp.number-1, "DF",",", "p-value:" , Q_chiValid)

  cat("\n")

  multi_return <- function() {
    Out_list <- list("Qstat" = Q_valid, "Qpval"=Q_chiValid)

    #Defines class of output object

    return(Out_list)
  }
  OUT<-multi_return()
}
