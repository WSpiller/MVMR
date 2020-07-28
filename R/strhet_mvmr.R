#' strhet_mvmr
#'
#' Calculates the conditional F-statistic for assessing instrument strength in two sample summary multivariable Mendelian randomization through minimisation of Q-statistics.
#' The function takes a formatted dataframe as an input, obtained using the function \code{format_mvmr}. Additionally, covariance matrices
#'  for estimated effects of individual genetic variants on each exposure can also be provided. These can be estimated using external data by
#'  applying the \code{snpcov_mvmr} or \code{phenocov_mvmr} functions, are input manually. The function returns a dataframe including the conditional
#'  F-statistic with respect to each exposure. A conventional F-statistic threshold of 10 is used in basic assessments of instrument strength.
#'
#'
#' @param r_input A formatted data frame using the \code{format_mvmr} function.
#' @param gencov Calculating heterogeneity statistics requires the covariance between the effect of the genetic variants on each exposure to be known. This can either be estimated from individual level data, be assumed to be zero, or fixed at zero using non-overlapping samples of each exposure GWAS. A value of \code{0} is used by default.
#'
#' @return A dataframe showing the conditional F-statistic for each exposure.
#' 
#'@author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#'@references Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2018. [Internet]. 2018;dyy262. Available from: https://dx.doi.org/10.1093/ije/dyy262
#'@export
#'@examples
#'
#' strhet_mvmr(data,covariances)
#' 

strhet_mvmr<-function(r_input,gencov){
  
  # Inverse variance weighting is used. 
  
  Wj<-1/r_input[,3]^2
  
  if(missing(gencov)) {
    gencov<-as.numeric(0)
    warning("Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.")
  }
  
  #Determine the number of exposures included in the model
  
  exp.number<-length(names(r_input)[-c(1,2,3)])/2
  
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
  
  Q_strength<-matrix(ncol=exp.number,nrow=1,0)
  
  qminvec<-NULL
  
  valvec<-combn(-100:100,c(exp.number-1))
  
  for(m in 1:exp.number){
    
    qmin<-function(a){
      
      delta_mat[,m]<-valvec[a]
      
      sigma2xj_dat<-matrix(ncol=exp.number,nrow=length(r_input[,1]),0)
      
      if(length(gencov) > 0){
        
        #Create a subset containing only standard errors for exposure effect estimates
        sebetas<-r_input[,(exp.number + 4):length(r_input)]
        
        #Generates the sigma2xj values for each exposure
        for(i in 1:exp.number){
          se.temp<-as.matrix(sebetas[,-i])
          for(j in 1:(exp.number-1)){
            sigma2xj_dat[,i]<- sigma2xj_dat[,i] + (se.temp[,j]^2 * delta_mat[j,i]^2)
          }
          sigma2xj_dat[,i]<- sigma2xj_dat[,i] + sebetas[,i]^2
          
        }
        
        
      }
      
      if(length(gencov) > 2){
        
        #Create a subset containing only standard errors for exposure effect estimates
        sebetas<-r_input[,(exp.number + 4):length(r_input)]
        
        #Generates the sigma2xj values for each exposure
        for(i in 1:exp.number){
          se.temp<-as.matrix(sebetas[,-i])
          for(j in 1:(exp.number-1)){
            sigma2xj_dat[,i]<- sigma2xj_dat[,i] + (se.temp[,j]^2 * delta_mat[j,i]^2)
          }
          sigma2xj_dat[,i]<- sigma2xj_dat[,i] + sebetas[,i]^2 - 0
          
          for(k in 1:(exp.number-1)){
            
            temp<-seq(1:exp.number)
            
            temp<-temp[-i]
            
            for(l in 1:length(r_input[,1])){
              
              sigma2xj_dat[l,i]<- sigma2xj_dat[l,i] - delta_mat[j,i]*2*gencov[[l]][i,temp[k]]
              
            }
            
          }
          
        }
        
      }
      
      
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
      
      return(Q_strength[m])
      
    }
    
    qminvec[m] <- optimize(qmin, interval = c(1,length(valvec)))$objective
    
  }
  
  Q_strength<-data.frame(qminvec)
  Q_strength<-data.frame(t(Q_strength))
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