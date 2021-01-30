#' mvmr (legacy)
#'
#' Note: This function is from the old version of the MVMR package and will be replaced in the future: The gencov argument should be set to zero when using \code{mvmr()}.
#'
#' Fits an IVW multivariable Mendelian randomization model using first order weights. The function returns an object of class \code{"MVMRIVW"}, containing regression estimates, estimated heterogeneity as a measure of instrument strength (\code{Q_strength}), and estimated heterogeneity as a measure of instrument validity (\code{Q_valid}).
#'
#' @param r_input A formatted data frame using the \code{format_mvmr} function or an object of class `MRMVInput` from [`MendelianRandomization::mr_mvinput`]
#' @param gencov Calculating heterogeneity statistics requires the covariance between the effect of the genetic variants on each exposure to be known. This can either be estimated from individual level data, be assumed to be zero, or fixed at zero using non-overlapping samples of each exposure GWAS. A value of \code{0} is used by default.
#' @param weights A value specifying the inverse variance weights used to calculate IVW estimate and Cochran's Q statistic. Currently only first order weights are available (\code{1}).
#'
#' @return An object of class \code{"MVMRIVW"} containing the following components:\describe{
#' \item{\code{summary}}{A summary of the MVMR regression model, including estimated coefficients, standard errors, t-statistics, p-values, and heterogeneity statistics.}
#' \item{\code{coef}}{The estimated coefficients, their standard errors, t-statistics, and corresponding (two-sided) p-values.}
#' \item{\code{Q_strength}}{A data frame displaying modified Cochran's Q statistics for assessing instrument strength with respect to each exposure. The Q-statistic increases proportionally with instrument strength, and analogous to univariate MR analyses, a value equal to or greater than 10 can be used as a minimum threshold for instrument strength. Note that for these statistics it is not informative to evaluate p-values.}
#' \item{\code{Q_valid}}{A modified form of Cochran's Q statistic measuring heterogeneity in causal effect estimates obtained using each genetic variant. Observed heterogeneity is indicative of a violation of the exclusion restriction assumption in MR (validity), which can result in biased effect estimates.}
#' \item{\code{p_valid}}{A p-value corresponding to the heterogeneity measure for instrument validity (\code{Q_valid})}
#'}
#'@author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#'@references Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2019, 48, 3, 713-727. <https://dx.doi.org/10.1093/ije/dyy262>
#' @importFrom stats lm as.formula pchisq pf
#' @export
#' @examples
#' # Example using format_mvmr formatted data
#' r_input <- format_mvmr(
#'     BXGs = rawdat_mvmr[,c("LDL_beta","HDL_beta")],
#'     BYG = rawdat_mvmr$SBP_beta,
#'     seBXGs = rawdat_mvmr[,c("LDL_se","HDL_se")],
#'     seBYG = rawdat_mvmr$SBP_se,
#'     RSID = rawdat_mvmr$SNP)
#' mvmr(r_input, 0, 1)
#'
#' # Example using MRMVInput formatted data from the MendelianRandomization package
#' bx <- as.matrix(rawdat_mvmr[,c("LDL_beta", "HDL_beta")])
#' bxse <- as.matrix(rawdat_mvmr[,c("LDL_se", "HDL_se")])
#' dat <- MendelianRandomization::mr_mvinput(bx = bx,
#'                                           bxse = bxse,
#'                                           by = rawdat_mvmr$SBP_beta,
#'                                           byse = rawdat_mvmr$SBP_se,
#'                                           snps = rawdat_mvmr$SNP)
#' mvmr(r_input = r_input, gencov = 0, weights = 1)

# Define IVW Multivariable MR function: This takes the formatted dataframe from
# the format_MVMR function, as an input, the covariance between the effect of the
# genetic variants on each exposure, and a value indicating the weights to
# used in the analysis.

mvmr<-function(r_input,gencov,weights){

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

  if(missing(gencov)) {
    gencov<-0
    warning("Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.")
  }

  #If weights is missing, first order weights are used by default.

  if(missing(weights)) {
    weights<-1
    warning("Weights not specified: Adopting first-order weights")
  }

  # A value of 1 for weights indicates inverse variance weighting is to be used.

  if(weights==1){
    Wj<-1/r_input[,3]^2
  }

  #Determine the number of exposures included in the model

  exp.number<-length(names(r_input)[-c(1,2,3)])/2

  #Fit the IVW MVMR model

  A_sum<-summary(lm(as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,weights=Wj,data=r_input))

  A<-summary(lm(as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,weights=Wj,data=r_input))$coef

  #Rename the regressors for ease of interpretation
  for(i in 1:exp.number){
    dimnames(A)[[1]][i]<- paste0("exposure",i,collapse="")
  }

  #############################################
  # Generalised instrument strength het.stats #
  #############################################

  #Create an empty matrix for delta value (coefficients regressing each set of exposure effects upon other
  #exposure effects)

  delta_mat<-matrix(0,ncol=exp.number,nrow=exp.number-1)

  #Obtain delta values fitting regression models for each set of exposure effects upon other exposure effects
  for(i in 1:exp.number){
    regressand<-names(r_input[3 + i])
    regressors<-names(r_input)[-c(1,2,3,
                                  4+exp.number:length(names(r_input)))]
    C<-paste(regressand, "~", "-1 +", paste(regressors[-i], collapse="+"))
    D.reg<-lm(C,data=r_input)
    delta_mat[,i]<-D.reg$coefficients
  }

  #Create an empty matrix for sigma2xj values
  sigma2xj_dat<-matrix(ncol=exp.number,nrow=length(r_input[,1]),0)

  #Create a subset containing only standard errors for exposure effect estimates
  sebetas<-r_input[,(exp.number + 4):length(r_input)]

  #Generates the sigma2xj values for each exposure
  for(i in 1:exp.number){
    se.temp<-as.matrix(sebetas[,-i])
    for(j in 1:(exp.number-1)){
      sigma2xj_dat[,i]<- sigma2xj_dat[,i] + (se.temp[,j]^2 * delta_mat[j,i]^2)
    }
    sigma2xj_dat[,i]<- sigma2xj_dat[,i] + sebetas[,i]^2 - gencov

  }

  #Create an empty matrix for instrument strength Q statistics
  Q_strength<-matrix(ncol=exp.number,nrow=1,0)

  #Generates the component of the Q statistic to be subtracted from the exposure estimates
  for(i in 1:exp.number){
    betas<-r_input[,c(4:(3+exp.number))]
    betas<-data.frame(betas[,-i])
    temp.sub <- 0
    for(j in 1:(exp.number-1)){
      temp.sub<-temp.sub + (delta_mat[j,i] * betas[,j])
    }

    #Populates matrix of Q statistics with respect to instrument strength
    Q_strength[i]<- sum( (1/sigma2xj_dat[,i]) * ((r_input[,3+i] - temp.sub)^2) )
  }

  Q_strength<-data.frame(Q_strength)
  names(Q_strength)<-dimnames(A)[[1]]
  rownames(Q_strength)<-"Q"

  ########################
  ## Instrument Validity #
  ########################

  # Generate Sigma^2_A values
  sigma2A<-r_input[,3]^2
  for(i in 1:exp.number){
    sigma2A<-sigma2A + (A[i]^2 * sebetas[,i]^2)
  }

  if(gencov != 0){

    warning("Validity statistics do not currently account for covariance between the effect of the genetic variants on each exposure")

  }

  #Create a subset of exposure effect estimates
  betas<-r_input[,c(4:(3+exp.number))]

  #Generates the component of the Q statistic to be subtracted from the outcome estimates
  temp.sub2<-0
  for(i in 1:exp.number){
    temp.sub2<-temp.sub2 + (betas[,i] * A[i])
  }

  #Calculates Q statistic for instrument validity
  Q_valid<- sum ((1/sigma2A)*(r_input[,2]-temp.sub2)^2)

  #Calculates p_value for instrument validity
  Q_chiValid<-pchisq(Q_valid,length(r_input[,2])-exp.number-1,lower.tail = FALSE)

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

  cat(paste(c("\nF-statistic:", " on"," and"), round(A_sum$fstatistic,2), collapse=""),
      "DF, p-value:",
      format.pval(pf(A_sum$fstatistic[1L], A_sum$fstatistic[2L], A_sum$fstatistic[3L],
                     lower.tail = FALSE), digits=3))

  cat("\n")

  cat("\n")

  cat("------------------------------")

  cat("\n")

  cat("Q-Statistics for instrument strength:")

  cat("\n")

  cat("\n")

  print(Q_strength)

  cat("\n")

  cat("------------------------------")

  cat("\n")

  cat("Q-Statistic for instrument validity:")

  cat("\n")

  cat("\n")

  cat(Q_valid, "on", length(r_input[,2])-exp.number-1, "DF",",", "p-value:" , Q_chiValid)

  cat("\n")


  #A function allowing multiple objects to be accessed from the output object
  multi_return <- function() {
    Out_list <- list("coef" = A, "Q_strength"=Q_strength, "Q_valid"=Q_valid, "p_valid"=Q_chiValid)

    #Defines class of output object
    class(Out_list)<-"MVMRIVW"
    return(Out_list)
  }
  OUT<-multi_return()
}
