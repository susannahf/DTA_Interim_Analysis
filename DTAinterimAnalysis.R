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

  # sanity checking  
  if(length(nk)!=k) {stop("discreteInterimFleming: Every stage (k) should be assigned a corresponding sample size (nk)")}
  if(length(Ek)!=k) {stop("discreteInterimFleming: Every stage (k) should be assigned a corresponding number of events (Ek)")}
  if(p0 <= 0 | p0 >=1) {stop("discreteInterimFleming: p0 must be between 0 and 1")}
  if(any(Ek>nk)) {stop("discreteInterimFleming: there should not be more events than datapoints at any timepoint")}
  
  # Fleming TR, page 146 (last paragraph before Evaluation section)
  # commented out as there is no justification for k<5
  #if(k > 5) {stop(" The stages specified should be utmost 5")}
  
  # Validity of nominal type I/ II error 
  if(alpha <= 0 | alpha >= 1 ) {stop("discreteInterimFleming: Nominal type I error rate must be between 0 and 1")}


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
                                                          "Estimated proportion"=cumsum(Ek)/cumsum(nk),
                                                          "Decision"=Conclude))
  
  return(output)
  
}


# run Fleming discrete interim analysis for one or more interim analyses and cumulative data
# does not do the final adjustment of rg at k
# ns are the numbers of data points at each interim analysis (cumulative)
# events are the numbers of events at each interim analysis (cumulative)
# ns and events can be a vector or scalar but must have same length
# finaln is the planned final sample size
# p0is the proportion such that H0: p <= p0
# alpha is the one sided nominal type I error (default 0.05)
cumulDiscreteInterimFleming <- function(ns, events, finaln, p0, alpha=0.05) 
{
  # check inputs
  if(length(ns)!=length(events)) {stop("cumulDiscreteInterimFleming: there should be an equal number of ns and events")}
  if(max(ns)>finaln) {stop("cumulDiscreteInterimFleming: finaln should not be smaller than an interim n")}
  if(p0 <= 0 | p0 >=1) {stop("cumulDiscreteInterimFleming: p0 must be between 0 and 1")}
  if(alpha <= 0 | alpha >= 1 ) {stop("cumulDiscreteInterimFleming: Nominal type I error rate must be between 0 and 1")}
  if(any(events>ns)) {stop("cumulDiscreteInterimFleming: there should not be more events than datapoints at any timepoint")}
  
  ### Fleming's (1982) critical boundaries 
  z.alpha <- qnorm(1-alpha)
  Nk <- finaln
  pa <- (sqrt(Nk*p0) + z.alpha*sqrt(1-p0))^2/(Nk + (z.alpha)^2)
  
  rg <- NULL
  ag <- NULL
  for(g in seq(length(ns)))
  {
    ### Fleming TR, Equation 5 (recall sum(ai*b)=b*sum(ai))
    rg[g] <- round(p0*ns[g] + z.alpha*sqrt(Nk*p0*(1-p0)))+1 
    ag[g] <- round(pa*ns[g] - z.alpha*sqrt(Nk*pa*(1-pa)))
  }
  
  
  ### Conclusion based on critical boundaries
  Sk <- events
  
  terminationStage <- 0
  Conclude <- NULL
  for(g in seq(length(Sk)))
  {
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
  
  ##### Output results
  
  output <- list("terminationStage"= terminationStage,
                 "terminationType" =Conclude[terminationStage],
                 "p" = p0,
                 "details"=cbind.data.frame("Stage"=seq(length(Sk)),
                                            "Cumulative_N"=ns,
                                            "Cumulative_E"=events,
                                            "Crit_ag"=ag,
                                            "Crit_rg"=rg,
                                            "Estimated proportion"=events/ns,
                                            "Decision"=Conclude))
  
  return(output)
  
}





# this applies Flemings method to DTA data
# data: data in the form of the test datasets created by createTestData.R
# ie. a data frame with binary reference, TP and TN values
# analysispoints: vector of points at which (interim) analyses will be made
# pSe: threshold p0 for sensitivity data
# pSp: threshold p0 for specificity data
# prevalence: predicted / known prevalence for study
# N / positiveN: planned / actual sample sizes.  N is total sample size. positiveN is sample size
#   of positive cases.  One of these must be provided.  If both are provided, positiveN is used.
# alpha: one sided nominal type I error (default 0.05)
# simpleOutput: chooses between simplified (default), and detailed output
DTAdiscreteInterimAnalysis <- function(data,analysispoints,pSe,pSp, prevalence, positiveN=NULL, N=NULL, alpha=0.05, simpleOutput=TRUE){

  # sanity checking
  if (!is.data.frame(data)) {stop("DTAdiscreteInterimAnalysis: data should be a data frame")}
  needdatacols <- c("reference", "TP", "TN")
  if (!all(needdatacols %in% colnames(data))) {stop("DTAdiscreteInterimAnalysis: data should have reference TP and TN columns")}
  if (!is.logical(data$reference) | !is.logical(data$TP) | !is.logical(data$TN))
    {stop("DTAdiscreteInterimAnalysis: data should have binary reference TP and TN columns")}
  if(pSe <= 0 | pSe >=1 | pSp <= 0 | pSp >=1) {stop("DTAdiscreteInterimAnalysis: p0 thresholds must be between 0 and 1")}
  if(prevalence <= 0 | prevalence >=1) {stop("DTAdiscreteInterimAnalysis: prevalence must be between 0 and 1")}
  if(alpha <= 0 | alpha >=1) {stop("DTAdiscreteInterimAnalysis: alpha must be between 0 and 1")}
  
  # warn if both positiveN and N are provided
  if(!is.null(positiveN) & !is.null(N)) {
    warning("both N and positive N have been provided.  Using positive N only.")
  }
  # calculate NSe and NSp
  if(!is.null(positiveN)){
    NSe <- positiveN
    NSp <- NSe/prevalence - NSe
  } else {
    NSe <- N * prevalence
    NSp <- N - NSe
  }
  # check that Ns don't need inflating
  if(sum(data$reference==TRUE) > NSe) {
    NSe <- sum(data$reference==TRUE)
    warning("NSe inflated as more positives than expected")
  }
  if(sum(data$reference==FALSE) > NSp) {
    NSp <- sum(data$reference==FALSE)
    warning("NSp inflated as more negatives than expected")
  }
  
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
  
  # actually 1-Se and 1-Sp, so p0 needs to be 1-p
  SeInterim <- cumulDiscreteInterimFleming(ns=cumsum(interimvars$nkSe), 
                                           events = cumsum(interimvars$EkSe),
                                           finaln = NSe,
                                           p0=1-pSe,   
                                           alpha=alpha)
  SpInterim <- cumulDiscreteInterimFleming(ns=cumsum(interimvars$nkSp), 
                                           events = cumsum(interimvars$EkSp),
                                           finaln = NSp,
                                           p0=1-pSp,   
                                           alpha=alpha)
  
  res <- list("Sensitivity" = SeInterim, "Specificity"= SpInterim)
  
  if(simpleOutput){
    ressimp <- list("Sensitivity" = simplifiedDiscreteInterimOutput(res$Sensitivity, N),
        "Specificity" = simplifiedDiscreteInterimOutput(res$Specificity, N))
    return(ressimp)
  }
  else{ return(res)
  }
  
  
}

# helper function for DTAdiscreteInterimAnalysis() 
## interpret output of main interim analysis code into a data frame suitable for DTA
simplifiedDiscreteInterimOutput <- function(detailedResults, N) {
  termpoint <- detailedResults$terminationStage
  termtype <- detailedResults$terminationType
  p0 <- 1-detailedResults$p
  # termination happens? 
  term <- termpoint>0
  futility <- length(grep("H1",termtype))>0 # true if H1 present
  # no termination at study end
  endstop = (!term & length(detailedResults$details$Decision)==N)
  
  # if no termination, set termpoint and futility to NA
  if(!term) {
    termpoint=NA
    futility=NA
  }
  
  # set termstring
  termstring <- ""
  if(!term) { 
    if(endstop) termstring <- paste0(termstring, "Reached study end without termination")
    else termstring <- paste0(termstring, "No termination by final interim analysis")
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
  
# run DTA interim analysis using only cumulative data
# data should be of the form of a data frame with the following columns (each corresponding to the data up to an interim analysis)
# - N number of data points
# - RefT number of positive reference cases
# - TP number of true positives
# - TN number of true negatives
DTAcumulativeInterimAnalysis <- function(data, pSe, pSp, prevalence, positiveN=NULL, N=NULL, alpha=0.05, simpleOutput = TRUE){
  
  # sanity checking
  if (!is.data.frame(data)) {stop("DTAcumulativeInterimAnalysis: data should be a data frame")}
  needdatacols <- c("N", "RefT", "TP", "TN")
  if (!all(needdatacols %in% colnames(data))) {stop("DTAcumulativeInterimAnalysis: data should have N RefT TP and TN columns")}
  if (!is.numeric(data$N) | !is.numeric(data$RefT) | !is.numeric(data$TP) | !is.numeric(data$TN))
  {stop("DTAcumulativeInterimAnalysis: data should have numeric N RefT TP and TN columns")}
  if(pSe <= 0 | pSe >=1 | pSp <= 0 | pSp >=1) {stop("DTAcumulativeInterimAnalysis: p0 thresholds must be between 0 and 1")}
  if(prevalence <= 0 | prevalence >=1) {stop("DTAcumulativeInterimAnalysis: prevalence must be between 0 and 1")}
  if(alpha <= 0 | alpha >=1) {stop("DTAcumulativeInterimAnalysis: alpha must be between 0 and 1")}
  if(any(data$RefT>data$N) | any(data$TP>data$N) | any(data$TN>data$N)) 
    {stop("DTAcumulativeInterimAnalysis: RefT TP and TN must be less than N at each analysis point")}
  
  # warn if both positiveN and N are provided
  if(!is.null(positiveN) & !is.null(N)) {
    warning("both N and positive N have been provided.  Using positive N only.")
  }
  
  # calculate NSe and NSp
  if(!is.null(positiveN)){
    NSe <- positiveN
    NSp <- NSe/prevalence - NSe
  } else {
    NSe <- N * prevalence
    NSp <- N - NSe
  }
  # check that Ns don't need inflating
  if(max(data$RefT) > NSe) {
    NSe <- max(data$RefT)
    warning("NSe inflated as more positives than expected")
  }
  if(max(data$N-data$RefT) > NSp) {
    NSp <- max(data$N-data$RefT)
    warning("NSp inflated as more negatives than expected")
  }
  
    # create useful variables
  # we actually are analysing on 1-Se and 1-Sp, so Ek are adjusted for this
  nkSe <- data$RefT
  nkSp <- data$N-data$RefT
  EkSe <- nkSe - data$TP
  EkSp <- nkSp - data$TN
  
  # actually 1-Se and 1-Sp, so p0 needs to be 1-p
  SeInterim <- cumulDiscreteInterimFleming(ns=nkSe, 
                                           events = EkSe,
                                           finaln = NSe,
                                           p0=1-pSe,   
                                           alpha=alpha)
  SpInterim <- cumulDiscreteInterimFleming(ns=nkSp, 
                                           events = EkSp,
                                           finaln = NSp,
                                           p0=1-pSp,   
                                           alpha=alpha)
  
  res <- list("Sensitivity" = SeInterim, "Specificity"= SpInterim)
  
  if(simpleOutput){
    ressimp <- list("Sensitivity" = simplifiedDiscreteInterimOutput(res$Sensitivity, N),
                    "Specificity" = simplifiedDiscreteInterimOutput(res$Specificity, N))
    return(ressimp)
  }
  else{ return(res)
  }
  
}


# plot sensitivity and specificity against Fleming thresholds, allowing for which points contribute to each.
# data should contain reference, Se and Sp columns
# assumes that N is perfectly known for the data (i.e. this is the whole data and the planned sample size and prevalence)
plotSeSpFlemingThresholds <- function(data, p0Se, p0Sp) {
  
  # sanity checking
  if (!is.data.frame(data)) {stop("plotSeSpFlemingThresholds: data should be a data frame")}
  needdatacols <- c("reference", "Se", "Sp")
  if (!all(needdatacols %in% colnames(data))) {stop("plotSeSpFlemingThresholds: data should have reference Se and Sp columns")}
  if (!is.logical(data$reference) | !is.numeric(data$Se) | !is.numeric(data$Se))
  {stop("plotSeSpFlemingThresholds: data should have logical reference and numeric Se and Sp columns")}
  if(pSe <= 0 | pSe >=1 | pSp <= 0 | pSp >=1) {stop("plotSeSpFlemingThresholds: p0 thresholds must be between 0 and 1")}
  
  # calculate Ns and Fleming thresholds
  
  
  # plot 
  
  
  
}




