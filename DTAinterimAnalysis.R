# Methods for running the exact group sequential method for DTA data
# S Fleming and L Mwandigha, 2022 

library(dplyr)

# this function implements the exact group sequential method as described in
# Fleming (1982) and Zhau (2007)
# k is the number of stages
# beta is nominal type II error
# alpha is the one sided nominal type I error
# p0 is the proportion such that H0: p <= p0
# nk are the sample size in each stage (not cumulated across the stages)
# Ek are the number of successes / events in each stage (not cumulated across the stages)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

discreteInterimFleming <- function(k,
                                    alpha, 
                                    p0, 
                                    nk, 
                                    Ek) 
{
  ### Prerequisites

  
  if(length(nk)!=k) {stop("Every stage (k) specified should be assigned a corresponding sample size (nk)")}
  if(length(Ek)> k) {stop("Number of endpoint events (E) exceeds number of stages (k)")}
  
  # Fleming TR, page 146 (last paragraph before Evaluation section)
  # commented out as there is no justification for k<5
  #if(k > 5) {stop(" The stages specified should be utmost 5")}
  
  # Validity of nominal type I/ II error 
  if(alpha <= 0 | alpha >= 1 ) {stop("Nominal type I error rate must be between 0 and 1")}


  ### Fleming's (1982) critical boundaries 

  
  z.alpha <- qnorm(1-alpha)
  Nk <- sum(nk)
  pa <- (sqrt(sum(nk)*p0) + z.alpha*sqrt(1-p0))^2/(sum(nk) + (z.alpha)^2)
  
  rg <- NULL
  ag <- NULL
  for(g in seq(length(nk)))
  {
    ### Fleming TR, Equation 5 (recall sum(ai*b)=b*sum(ai))
    rg[g] <- round(p0*sum(nk[seq(g)]) + z.alpha*sqrt(Nk*p0*(1-p0)))+1 
    
    ###Fleming TR, Equation 6
    ### up to k-1
    if(g<k)
    {
      ag[g] <- round(pa*sum(nk[seq(g)]) - z.alpha*sqrt(Nk*pa*(1-pa)))
    }
    ### up to k
    if(g==k)
    {
      ag[g] <- rg[g]-1
    }
  }
  

  ### Conclusion based on critical boundaries

  
  Sk <- cumsum(Ek)
  
  terminationStage <- 0
  Conclude <- NULL
  for(g in seq(length(Sk)))
  {
    
    ### up to k-1
    if(g<k)
    {
      ###
      if(Sk[g] <= ag[g]){
        Conclude[g] <- paste0("Stop sampling and accept H0:p<=", p0)
        if(terminationStage==0) { terminationStage <- g}
      }
      if(Sk[g] >= rg[g]){
        Conclude[g] <- paste0("Stop sampling and accept H1:p>=", p0)
        if(terminationStage==0){terminationStage <- g}
      }
      if(Sk[g] > ag[g] & Sk[g] < rg[g]){
        Conclude[g] <- paste0("Continue to stage ", g+1)}
    }
    ### up to k
    if(g==k)
    {
      if(Sk[g] <= ag[g]){
        Conclude[g] <- paste0("Accept H0:p<=", p0," at study end")
        if(terminationStage==0){ terminationStage <- g}
      }
      if(Sk[g] >= rg[g]){
        Conclude[g] <- paste0("Accept H1:p>=", p0, " at study end")
        if(terminationStage==0){ terminationStage <- g}
      }
    }
  }
  
  ##### Output results
  
  output <- list("terminationStage"= terminationStage,
                 "terminationType" =Conclude[terminationStage],
                 "p" = p0,
                 "details"=cbind.data.frame("Stage"=seq(nk),
                                                          "Cumulative_N"=cumsum(nk),
                                                          "Cumulative_E"=cumsum(Ek),
                                                          "Crit_ag"=ag,
                                                          "Crit_rg"=rg,
                                                          "Decision"=Conclude))
  
  return(output)
  
}

# this applies Flemings method to DTA data
# data: data in the form of the test datasets created by createTestData.R
# analysispoints: vector of points at which (interim) analyses will be made
# pSe: threshold p0 for sensitivity data
# pSp: threshold p0 for specificity data
# alpha: one sided nominal type I error (default 0.05)
# simpleOutput: chooses between simplified (default), and detailed output
DTAdiscreteInterimAnalysis <- function(data,analysispoints,pSe,pSp,alpha=0.05, simpleOutput=TRUE){

  # first create a counter variable so that I can cut on it
  data$n <- 1:nrow(data)
  data$gp <- cut(data$n,c(0,analysispoints))
  # remove any rows with NAs (fall outside analysis)
  data <- data[!is.na(data$gp),]
  # create useful variables
  # we actually are analysing on 1-Se and 1-Sp, so Ek are adjusted for this
  interimvars <- data %>% group_by(gp) %>% 
    summarise(nkSe=sum(reference==TRUE),
              nkSp=sum(reference==FALSE),
              EkSe=nkSe - sum(TP==TRUE),  
              EkSp=nkSp - sum(TN==TRUE))
  k <- length(analysispoints) 
  # actually 1-Se and 1-Sp, so p0 needs to be 1-p
  SeInterim <- discreteInterimFleming(k=k,
                                       alpha=alpha,
                                       p0=1-pSe,   
                                       nk=interimvars$nkSe,
                                       Ek=interimvars$EkSe)
  SpInterim <- discreteInterimFleming(k=k,
                                       alpha=alpha,
                                       p0=1-pSp,
                                       nk=interimvars$nkSp,
                                       Ek=interimvars$EkSp)
  res <- list("Sensitivity" = SeInterim, "Specificity"= SpInterim)
  
  if(simpleOutput){
    ressimp <- list("Sensitivity" = simplifiedDiscreteInterimOutput(res$Sensitivity),
        "Specificity" = simplifiedDiscreteInterimOutput(res$Specificity))
    return(ressimp)
  }
  else{ return(res)
  }
  
  
}

# helper function for DTAdiscreteInterimAnalysis() 
## interpret output of main interim analysis code into a data frame suitable for DTA
simplifiedDiscreteInterimOutput <- function(detailedResults) {
  termpoint <- detailedResults$terminationStage
  termtype <- detailedResults$terminationType
  p0 <- 1-detailedResults$p
  endstop <- termpoint==length(detailedResults$details$Decision)
  futility <- length(grep("H1",termtype))>0 # true if H1 present
  
  termstring <- ""
  if(endstop) { 
    if(futility) termstring <- paste0(termstring, "Accept H1 p<=", p0, " at study end")
    else termstring <- paste0(termstring, "Accept H0 p>=", p0, " at study end")
  }
  else { termstring <- paste0(termstring, "Stop for ")
    if(futility) { termstring <- paste0(termstring, "futility with p<=", p0, " at k=", termpoint)}
    else { termstring <- paste0(termstring, "efficacy with p>=", p0 , " at k=", termpoint)}
  }
   
  res <- data.frame(terminateAt = termpoint,
                    futility = futility,
                    termstring = termstring)
  
  return(res)  
  
}

# model termination lines for Fleming method
# termination by accepting H0
# when Sk <= ag
# therefore termination line is where Sk=ag
# note that at end, ag would be set to just below rg
# termination by accepting H1
# when Sk >= rg
# therefore termination line is where Sk=rg
# Sk is cumulative number of events
# therefore proportion at given point (i) for termination is either rg[i]/n[i] or ag[i]/n[i]
# p0 is proportion to be tested
# n is number of points (default 200)
modelFlemingTerminationThresholds <- function(p0, n=200, alpha=0.05) {
  # set up variables
  nk = rep(1,n) # vector of additional n up to N=n
  rg <- ag <- H0lim <- H1lim <- NULL
  
  # copied code from discreteInterimFleming
  z.alpha <- qnorm(1-alpha)
  Nk <- sum(nk)
  pa <- (sqrt(sum(nk)*p0) + z.alpha*sqrt(1-p0))^2/(sum(nk) + (z.alpha)^2)
  
  # calculate rg and ag at every point
  for(g in seq(n)) {
    # copied code from discreteInterimFleming
    ### Fleming TR, Equation 5 (recall sum(ai*b)=b*sum(ai))
    rg[g] <- round(p0*sum(nk[seq(g)]) + z.alpha*sqrt(Nk*p0*(1-p0)))+1 
    ###Fleming TR, Equation 6
    ag[g] <- round(pa*sum(nk[seq(g)]) - z.alpha*sqrt(Nk*pa*(1-pa)))
    
    H0lim[g] <- ag[g]/sum(nk[seq(g)])
    H1lim[g] <- rg[g]/sum(nk[seq(g)])
  }
  
  #values of H0lim and H1lim greater than 1 or less than 0 are meaningless
  H0lim[H0lim<0] <- NA
  H1lim[H1lim>1] <- NA
  
  # add inverse limits of 1- for use in DTA
  res <- data.frame(rg=rg,
                    ag=ag,
                    H0limit=H0lim,
                    invH0limit = 1-H0lim,
                    H1limit=H1lim,
                    invH1limit = 1-H1lim)
  return(res)
}
  








