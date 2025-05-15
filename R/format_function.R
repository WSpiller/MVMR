#' format_mvmr
#'
#' Reads in summary data. Checks and organises columns for use in calculating multivariable Mendelian Randomization analyses. Where variant IDs are not provided, a vector is generated for variant identification.
#'
#' @param BXGs A matrix containing beta-coefficient values for genetic associations with the each exposure. Columns should indicate exposure number, with rows representing estimates for a given genetic variant.
#' @param BYG A numeric vector of beta-coefficient values for genetic associations with the outcome.
#' @param seBXGs A matrix containing standard errors corresponding to the matrix of beta-coefficients \code{BXGs}.
#' @param seBYG A numeric vector of standard errors corresponding to the beta-coefficients \code{BYG}.
#' @param RSID A vector of names for genetic variants included in the analysis. If variant IDs are not provided a vector of ID numbers will be generated.
#' @return A formatted data frame of class `mvmr_format`.
#'
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
#' names(r_input)
#' class(r_input)

#Function for formatting MVMR data frame

#Define formating function: This takes as inputs a matrix of instrument-exposure
#associations, a vector of instrument-outcome associations, a matrix of
#corresponding instrument-exposure standard errors, and a vector of identification
#numbers for the instruments.

format_mvmr<-function(BXGs,BYG,seBXGs,seBYG,RSID){

  #If no instrument-identification vector is provided, a set of placeholder values
  #is produced. A warning is also given to indicate no values were provided

  if(missing(RSID)) {
    RSID<-seq(from=1,to=length(BYG),by=1)
    warning("Missing SNP IDs; Generating placeholders")
  }

  #This loop names each column of instrument-exposure associations in the order
  #they appear in the provided matrix, labeling each exposure with an X and
  #subsequent index number. The first exposure provided is labelled betaX1.

  BXGs<-data.frame(BXGs)
  seBXGs<-data.frame(seBXGs)
  BYG<-data.frame(BYG)
  seBYG<-data.frame(seBYG)
  RSID<-data.frame(RSID)



  for(i in seq_len(ncol(BXGs))){
    names(BXGs)[i]<-paste0("betaX",i,collapse=",")
  }

  #This loop names each column of instrument-exposure standard errors in the order
  #they appear in the provided matrix, labeling each exposure with an X and
  #subsequent index number. The standard error for the first exposure provided
  #is labelled sebetaX1.

  for(i in seq_len(ncol(seBXGs))){
    names(seBXGs)[i]<-paste0("sebetaX",i,collapse=",")
  }

  # A dataframe containing all the necessary information for performing
  #multivariable MR is created, placing the variables in a specific order.

  dat<-cbind(RSID,BYG,seBYG,BXGs,seBXGs)

  # The columns of the dataframe are renamed so as to be interpretable in
  #subsequent functions.

  for(i in seq_len(ncol(dat))){

    if(i>1){
      dat[,i]<-as.numeric(dat[,i])

    }

  }

  names(dat) <- c("SNP", "betaYG", "sebetaYG",
                  names(BXGs), names(seBXGs))

  class(dat) <- append(class(dat), "mvmr_format")

  # The function returns the formatted dataframe

  return(dat)
}
