#' qhet_mvmr
#'
#' Fits a multivariable Mendelian randomization model adjusting for weak instruments. The functions requires a formatted dataframe using the \code{format_mvmr} function, as well a phenotypic correlation matrix \code{pcor}. This should be obtained from individual level
#' phenotypic data, or constructed as a correlation matrix where correlations have previously been reported. Confidence intervals are calculated using a non-parametric bootstrap.
#' By default, standard errors are not produced but can be calculated by setting \code{se = TRUE}. The number of bootstrap iterations is specified using the \code{iterations} argument.
#' Note that calculating confidence intervals at present can take a substantial amount of time.
#' 
#'
#'
#' @param r_input A formatted data frame using the \code{format_mvmr} function.
#' @param pcor A phenotypic correlation matrix including the correlation between each exposure included in the MVMR analysis.
#' @param CI Indicates whether 95 percent confidence intervals should be calculated using a non-parametric bootstrap.
#' @param iterations Specifies number of bootstrap iterations for calculating 95 percent confidence intervals.
#'
#' @return An dataframe containing effect estimates with respect to each exposure.
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#' @references Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2019, 48, 3, 713-727. <https://dx.doi.org/10.1093/ije/dyy262>
#' @importFrom stats optim optimize
#' @export
#' @examples
#'
#' \dontrun{
#' qhet_mvmr(r_input, pcor, CI = TRUE, iterations = 1000)
#' }
#'

qhet_mvmr<-function(r_input,pcor,CI,iterations){

  requireNamespace("boot", quietly = TRUE)

  warning("qhet_mvmr() is currently undergoing development.")
  
  if(missing(CI)) {
    CI<-F
    warning("95 percent confidence interval not calculated")
  }
  
  if(missing(iterations)) {
    iterations<-1000
    warning("Iterations for bootstrap not specified. Default = 1000")
  }
  
  exp.number<-length(names(r_input)[-c(1,2,3)])/2
  
  Qtemp<-function(r_input,pcor){
    
    exp.number<-length(names(r_input)[-c(1,2,3)])/2
    stderr<-as.matrix(r_input[,(exp.number + 4):length(r_input)])
    correlation<-pcor
    gammahat<-r_input$betaYG
    segamma<-r_input$sebetaYG
    pihat<-as.matrix(r_input[,c(4:(3+exp.number))])
    
    #estimation with the heterogeneity statistic
    
    PL_MVMR = function(a){
      tau2   = a[1]
      
      PL2_MVMR = function(ab){
        b<-ab
        
        cov = matrix(nrow = exp.number, ncol = exp.number)
        w=NULL
        for(l in 1:nrow(r_input)){
          for(pp in 1:exp.number){
            for(p2 in 1:exp.number){
              cov[pp,p2] <- correlation[pp,p2]*stderr[l,pp]*stderr[l,p2]
            }
          }
          
          segamma<-r_input$sebetaYG
          
          w[l] <- segamma[l]^2+t(b)%*%cov%*%b + tau2
        }
        
        q =  sum((1/w)*((gammahat - pihat%*%b)^2))
        
        return(q)
        
      }
      
      st_PL2  =  rep(0,exp.number)
      
      bc    = optim(st_PL2, PL2_MVMR)
      
      bcresults <- bc$par
      
      cov = matrix(nrow =exp.number, ncol = exp.number)
      w=NULL
      for(l in 1:nrow(r_input)){
        for(pp in 1:exp.number){
          for(p2 in 1:exp.number){
            cov[pp,p2] <- correlation[pp,p2]*stderr[l,pp]*stderr[l,p2]
          }
        }
        
        w[l] <- segamma[l]^2+t(bcresults)%*%cov%*%bcresults + tau2
      }
      
      q = (sum((1/w)*((gammahat - pihat%*%bcresults)^2))-(nrow(r_input)-2))^2
      
    }
    
    PL2_MVMR = function(ab){
      b=ab
      
      w=NULL
      cov = matrix(nrow = exp.number, ncol = exp.number)
      for(l in 1:nrow(r_input)){
        for(pp in 1:exp.number){
          for(p2 in 1:exp.number){
            cov[pp,p2] <- correlation[pp,p2]*stderr[l,pp]*stderr[l,p2]
          }
        }
        
        w[l] <- segamma[l]^2+t(b)%*%cov%*%b + tau_i
      }
      
      q =  sum((1/w)*((gammahat - pihat%*%b)^2))
      
    }
    
    limltauest  = optimize(PL_MVMR,interval=c(-10,10))
    tau_i = limltauest$objective
    
    tau = tau_i
    
    liml_het2<- optim(rep(0,exp.number), PL2_MVMR)
    limlhets <- liml_het2$par
    Qexact_het <- liml_het2$value
    
    Effects<-limlhets
    Effects<-data.frame(Effects)
    names(Effects)<-"Effect Estimates"
    for(i in 1:exp.number){
      rownames(Effects)[i]<-paste("Exposure",i,sep=" ")
      
    }
    
    return(Effects)
    
    
  }
  
  if(CI==F){
    
    res<-Qtemp(r_input,pcor)
    
  }
  
  if(CI==T){
    
    bootse<-function(data, indices){
      
      bres<-Qtemp(data[indices,],pcor)[,1]
      
      return(bres)
      
    }
    
    b.results<-boot(data=r_input,statistic=bootse,R=iterations)
    
    lcb<-NULL
    ucb<-NULL
    ci<-NULL
    
    for(i in 1:exp.number){
      
      lcb[i]<-round(boot.ci(b.results, type="bca",index=i)$bca[4],digits=3)
      ucb[i]<-round(boot.ci(b.results, type="bca",index=i)$bca[5],digits=3)
      
      ci[i]<-paste(lcb[i],ucb[i],sep="-")
      
    }
    
    res<-data.frame(b.results$t0,ci)
    
    names(res)<-c("Effect Estimates","95% CI")
    for(i in 1:exp.number){
      rownames(res)[i]<-paste("Exposure",i,sep=" ")
      
    }
    
  }
  
  return(res)
 
 
}



