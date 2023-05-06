#' strength_mvmr
#'
#' Calculates the conditional F-statistic for assessing instrument strength in two sample summary multivariable Mendelian randomization.
#' The function takes a formatted dataframe as an input, obtained using the function [`format_mvmr`]. Additionally, covariance matrices
#' for estimated effects of individual genetic variants on each exposure can also be provided. These can be estimated using external data by
#' applying the [`snpcov_mvmr`] or [`phenocov_mvmr`] functions, are input manually. The function returns a dataframe including the conditional
#' F-statistic with respect to each exposure. A conventional F-statistic threshold of 10 is used in basic assessments of instrument strength.
#'
#' @param r_input r_input A formatted data frame using the [`format_mvmr`] function or an object of class `MRMVInput` from [`MendelianRandomization::mr_mvinput`]
#' @param gencov Calculating heterogeneity statistics requires the covariance between the effect of the genetic variants on each exposure to be known. This can either be estimated from individual level data, be assumed to be zero, or fixed at zero using non-overlapping samples of each exposure GWAS. A value of \code{0} is used by default.
#'
#' @return A dataframe showing the conditional F-statistic for each exposure.
#'
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#' @references Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2018, 48, 3, 713-727. Available from: \doi{10.1093/ije/dyy262}
#' @export
#' @examples
#' \dontrun{
#' strength_mvmr(data, covariances)
#' }

strength_mvmr<-function(r_input,gencov){

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
    gencov<-as.numeric(0)
    warning("Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.")
  }


  # Inverse variance weighting is used.

  Wj<-1/r_input[,3]^2

  #Determine the number of exposures included in the model

  exp.number<-length(names(r_input)[-c(1,2,3)])/2

  A<-summary(stats::lm(stats::as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,data=r_input))$coef

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
    D.reg<-stats::lm(C,data=r_input)
    delta_mat[,i]<-D.reg$coefficients
  }

  sigma2xj_dat<-matrix(ncol=exp.number,nrow=length(r_input[,1]),0)

  if(length(gencov) < 2){

    #Create a subset containing only standard errors for exposure effect estimates
    sebetas<-r_input[,(exp.number + 4):length(r_input)]

    #Generates the sigma2xj values for each exposure
    for(i in 1:exp.number){
      se.temp<-as.matrix(sebetas[,-i])
      for(j in 1:(exp.number-1)){
        sigma2xj_dat[,i]<- sigma2xj_dat[,i] + (se.temp[,j]^2 * delta_mat[j,i]^2)
      }
      sigma2xj_dat[,i] <- sigma2xj_dat[,i] + sebetas[,i]^2

    }


  }

  if(length(gencov) > 2){
    sigma2xj_dat<-matrix(ncol=exp.number,nrow=length(r_input[,1]),0)
    delta.temp <- matrix(0,ncol=exp.number,nrow=exp.number)


    #Generates the sigma2xj values for each exposure
    for(i in 1:exp.number){

      if(i == 1) {
        delta.temp[,i] <- c(-1, delta_mat[,i])
      }

      if(i>1 & i<exp.number){
        delta.temp[,i] <- c(delta_mat[1:(i-1),i],-1,delta_mat[i:(exp.number-1),i])
      }

      if(i == exp.number){
        delta.temp[,i] <- c(delta_mat[,i],-1)
      }


      for(l in 1:length(r_input[,1])){

        sigma2xj_dat[l,i]<- sigma2xj_dat[l,i] +  t(delta.temp[,i])%*%gencov[[l]]%*%delta.temp[,i]

      }


    }
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
    Q_strength[i]<-Q_strength[i]/nrow(r_input)

  }

  Q_strength<-data.frame(Q_strength)
  names(Q_strength)<-dimnames(A)[[1]]
  rownames(Q_strength)<-"F-statistic"


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
