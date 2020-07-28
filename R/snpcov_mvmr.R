#' snpcov_mvmr
#'
#' Uses individual level genetic and exposure data to generate covariance matrices for estimated effects of individual genetic
#' variants on each exposure. The function returns a number of covariance matrices equal to the number of SNPs, where SNP and
#' row numbers reference ordered exposures.
#'
#' @param Gs A matrix or dataframe containing genetic instrument measures. Columns should indicate genetic variant number, with rows representing an observed measure of the genetic variant.
#' @param Xs A matrix or dataframe containing exposure measures. Columns should indicate exposure number, with rows representing an observed measure for the given exposure.
#'
#' @return A list of covariance matrices with respect to each genetic variant, retaining the ordering in \code{Gs}
#' 
#'@author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#'@references Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2018. [Internet]. 2018;dyy262. Available from: https://dx.doi.org/10.1093/ije/dyy262
#'@export
#'@examples
#'
#' snpcov_mvmr(data[,1:10],data[,11:13])
#' 
#' @importFrom stats lm

snpcov_mvmr<-function(Gs,Xs){
  
  betas<-matrix(0,ncol=length(Xs[1,]),nrow=length(Gs[1,]))
  
  resmat<-data.frame(rep(0,length(Gs[,1])))
  
  for(j in 1:length(Xs[1,])){
    
    for(i in 1:length(Gs[1,])){
      
      betas[i,j]<-lm(Xs[,j]~-1 + Gs[,i])$coefficients
      
      resids<-data.frame(lm(Xs[,j]~-1 + Gs[,i])$residuals)
      
      resmat<-cbind(resmat,resids)
      
    }
    
  }
  
  resmat<-resmat[,-1]
  
  sigmalist <- vector("list", length(Gs[1,]))
  
  for(i in 1:length(Gs[1,])){
    
    sigma_mattemp<-matrix(((t(Gs[,i]) %*% Gs[,i])^-1)/length(Gs[,i]),ncol=length(Xs[1,])
                          ,nrow=length(Xs[1,]))
    
    for(j in 1:length(Xs[1,])){
      
      for(k in 1:length(Xs[1,])){
        
        sigma_mattemp[k,j]<-sigma_mattemp[k,j] * sum(resmat[,i+((k-1)*length(Gs))] * resmat[,i+((j-1)*length(Gs))])
        
        
      }
    }
    
    sigmalist[[i]] <-sigma_mattemp
    
  }
 
  return(sigmalist)
   
}