#' phenocov_mvmr
#'
#' Uses an external phenotypic covariance matrix and summary data to estimate covariance matrices for estimated effects of individual genetic
#' variants on each exposure. The phenotypic covariance matrix should be constructed using standardised phenotype measures. The function returns a number of covariance matrices equal to the number of SNPs, where SNP and
#' row numbers reference ordered exposures.
#'
#' @param Pcov A phenotypic matrix using exposures, constructed using individual level exposure data. Columns should be ordered by exposure so as to match \code{format_mvmr}.
#' @param seBXGs A matrix containing standard errors corresponding in relation to the gene-exposure association for each SNP.
#'
#' @return A list of covariance matrices with respect to each genetic variant, retaining the ordering in \code{seBXGs}
#' 
#'@author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#'@references Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2018. [Internet]. 2018;dyy262. Available from: https://dx.doi.org/10.1093/ije/dyy262
#'@export
#'@examples
#'
#' phenocov_mvmr(Pcov,summarydata[,c(3,4)])
#' 

phenocov_mvmr<-function(Pcov,seBXGs){
  
  sigmalist <- vector("list", length(seBXGs[1,]))
  
  for(i in 1:length(seBXGs[,1])){
    
    sigma_mattemp<-Pcov
    
    for(j in 1:length(seBXGs[1,])){
      
      for(k in 1:length(seBXGs[1,])){
        
        sigma_mattemp[j,k]<-sigma_mattemp[j,k] * seBXGs[i,j] * seBXGs[i,k]
        
      }
      
    }
    
    sigmalist[[i]] <-sigma_mattemp
    
  }
  
  return(sigmalist)
  
}
  
  
  
  
  