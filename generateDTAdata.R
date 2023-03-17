# generate DTA data for interim analysis simulations
# S Fleming 2022

# need to generate 2 binary variables: actual outcome and test outcome 
# want to be able to specify: 
# - n
# - prevalence
# - sensitivity
# - specificity

# creates a 2x2 table from n, p (prevalence), se (sensitivity), sp (specificity)
# 2x2 table shape is as follows:
# TP  FN
# FP  TN
# these are defaulted as follows: n=100, p=0.25, se=0.7, sp=0.9 
# if verbose=T, will print actual values after rounding.
create2x2 <- function(n=100,p=0.25,se=0.7,sp=0.9, digits=3,verbose=F) {
  
  # Preliminary checks 
  if(p <= 0 | p >= 1 ) {stop("The specified prevalence must be between 0 and 1")}
  if(se <= 0 | se >= 1 ) {stop("The specified sensitivity must be between 0 and 1")}
  if(sp <= 0 | sp >= 1 ) {stop("The specified specificity must be between 0 and 1")}
  
  
  d = n*(se-p)/(1+se/sp-1/sp)
  a = se*(n-d/sp)
  b = d*(1/sp-1)
  c = a*(1/se-1)
  # round all values
  d=round(d)
  a=round(a)
  b=round(b)
  c=round(c)
  #create 2x2 for output
  twobytwo <- matrix(data=c(a,b,c,d),2,2)
  if(verbose){
    dimnames(twobytwo) <- list(c("row1","row2"),
                               c("column1","column2"))
    print(twobytwo)
    
    cat("Total sample size: ",a+b+c+d,sep="", "\n")
    cat("Output rounded to ", digits, " significant figures",sep="", "\n")
    cat("Prevalence: ",signif((a+b)/(a+b+c+d),digits=digits),sep="", "\n")
    cat("Sensitivity: ",signif(a/(a+c),digits=digits),sep="", "\n")
    cat("Specificity: ",signif(d/(d+b),digits=digits),sep="", "\n")
    
  }
  if(any(twobytwo<0)) stop("some cells require negative numbers to fulfil desired specifications")
  return(twobytwo)
}


# given a 2x2 table, creates a data set corresponding to that table
# 2x2 table shape is as follows:
# TP  FN
# FP  TN
# only full ordering (values in order TP FN FP TN) or random ordering is available from this function
DTAdataFrom2x2 <- function(twobytwo, randomise=T) {
  # check input 
  if(!is.numeric(twobytwo)) stop("2x2 table must be numeric")
  if(nrow(twobytwo)!=2 | ncol(twobytwo)!=2) stop("2x2 table must have dimensions of 2 by 2")
  # create basic dataset
  # create in order TP FN FP TN
  n = sum(twobytwo)
  refresult <- logical(n) # false by default
  indexresult <- logical(n) # false by default
  # ref is T in TP and FN
  endofFN <- sum(twobytwo[1,])
  refresult[1:endofFN] <- TRUE
  # index is T in TP and FP
  indexresult[1:twobytwo[1,1]] <- TRUE #TP
  indexresult[(endofFN+1): (endofFN+ twobytwo[2,1])] <- TRUE # FP
  # combine
  res <- data.frame(
    reference = refresult,
    index = indexresult)
  # shuffle dataset if necessary
  if(randomise){
    k<-sample(1:nrow(res)) # random permutation of indices
    res <- res[k,]
  }
  return(res)
}


# take DTA data (with reference and index data)
# calculate sensitivity and specificity at every point
continuousSeSp <- function(DTAdata, ref="reference", index="index", Se="Se", Sp="Sp") {
  # check that reference and index data exist and are logical
  if(!is.logical(DTAdata[,ref])) stop("reference data is not logical")
  if(!is.logical(DTAdata[,index])) stop("index data is not logical")
  
  # at each point, calculate se and sp up to that point
  ref <- DTAdata[,ref]
  index <- DTAdata[,index]
  TP <- ref & index
  TN <- !ref & !index
  sp <- se <- numeric(length(ref))
  # se = TP/P  , Sp = TN/N
  for(i in 1:length(ref)) {
    se[i] <- sum(TP[1:i])/sum(ref[1:i])
    sp[i] <- sum(TN[1:i])/sum(!ref[1:i])
  }
  
  # add the TP, TN, se and sp to the input data and output
  DTAdata$TP <- TP
  DTAdata$TN <- TN
  DTAdata[,Se] <- se
  DTAdata[,Sp] <- sp
  
  library(tidyverse)
  
  DTAdata_output <- DTAdata%>%
                     dplyr::mutate(Se=dplyr::case_when(is.nan(Se)~NA_real_,
                                                       !is.nan(Se)~Se),
                                   Sp=dplyr::case_when(is.nan(Sp)~NA_real_,
                                                       !is.nan(Sp)~Sp))
  return(DTAdata_output)
}







