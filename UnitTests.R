# unit tests
# not everything is unit tested, but any new code from 24 Oct 2023 should be.

# clear environment first
rm(list=ls())

source("DTAinterimAnalysis.R")
source("generateDTAdata.R")

# unit tests for cumulDiscreteInterimFleming
# create a few random data sets 
# run both DiscreteInterimFleming and cumulDiscreteInterimFleming on the data, not including k=n
# should have the same rg, ag, and so on in details of output

# run 10 times
for(i in seq(10)) {
  # create random non-cumulative event data for n=100
  rdata <- rbinom(100,1,0.4)
  totaln <- 100
  k=5
  p0 <- 0.4
  # k random interim analysis points
  anaPoints <- c(sort(sample(10:(totaln-10), k-1)), totaln)
  # cumulative
  cdata <- cumsum(rdata)[anaPoints]
  # non-cumulative
  nk <- c(anaPoints[1], diff(anaPoints))
  Ek <- c(cdata[1], diff(cdata))
  
  discdata <- discreteInterimFleming(k=k, alpha=0.05, p0=p0, nk=nk, Ek=Ek)
  cumdata <- cumulDiscreteInterimFleming(ns=anaPoints, events=cdata, finaln=totaln, p0=p0, alpha=0.05)
  
  # tests - these should be identical
  disctst <- discdata$details
  cumtst <- cumdata$details
  if(any(disctst!=cumtst)) stop("Test failed: different results from discreteInterimFleming and cumulDiscreteInterimFleming")
}
# if not stopped by here, test passed
print("Test passed: discreteInterimFleming and cumulDiscreteInterimFleming agree")

# unit tests for DTAdiscreteInterimAnalysis
# compare DTAdiscreteInterimAnalysis with cumulDiscreteInterimFleming

# run 10 times
for(i in seq(10)) {
  p0Se <- 0.65
  p0Sp <- 0.85
  prevalence <- 0.35
  N <- 200

  # create data
  tbt <- create2x2(n=10*N, p=prevalence, se=p0Se, sp=p0Sp, digits=4, verbose=FALSE)
  dtadata <- DTAdataFrom2x2(tbt, randomise = T)
  dtadata <- continuousSeSp(dtadata[1:N, ])
  
  # analysis points
  k <- 5
  anaPoints <- sort(sample(10:150, k))
  
  # pull out data for cumulDiscreteInterimFleming and run it
  cumframe <- data.frame(N=numeric(),
                         RefT=numeric(), 
                         TP=numeric(),
                         TN=numeric())
  for(i in seq(length(anaPoints))) {
    k <- anaPoints[i]
    ifrm <- data.frame(N = k, 
                       RefT = sum(dtadata$reference[1:k]), 
                       TP = sum(dtadata$TP[1:k]),
                       TN = sum(dtadata$TN[1:k]))
    cumframe <- rbind(cumframe,ifrm)
  }
  nsSe <- cumframe$RefT
  nsSp <- cumframe$N - cumframe$RefT
  ekSe <- nsSe-cumframe$TP
  ekSp <- nsSp-cumframe$TN
  NSe <- sum(dtadata$reference)
  NSp <- N-NSe
  
  
  SeFleming <- cumulDiscreteInterimFleming(ns=nsSe, events=ekSe, finaln=NSe, p0=1-p0Se, alpha=0.05)
  SpFleming <- cumulDiscreteInterimFleming(ns=nsSp, events=ekSp, finaln=NSp, p0=1-p0Sp, alpha=0.05)
  
  # pull out data for DTAdiscreteInterimAnalysis and run it
  actprev <- NSe/N
  
  dtares <- DTAdiscreteInterimAnalysis(dtadata,analysispoints = anaPoints, pSe=p0Se, pSp=p0Sp, prevalence = actprev, N=N, simpleOutput = FALSE)
  # compare results
  SeFlemingTest <- SeFleming$details
  SpFlemingTest <- SpFleming$details
  dtaSeTest <- dtares$Sensitivity$details
  dtaSpTest <- dtares$Specificity$details
  
  if(any(SeFlemingTest!=dtaSeTest, na.rm=T)) stop("Test failed: different results for Se from DTAdiscreteInterimAnalsysis and cumulDiscreteInterimFleming")
  if(any(SpFlemingTest!=dtaSpTest, na.rm=T)) stop("Test failed: different results for Sp from DTAdiscreteInterimAnalsysis and cumulDiscreteInterimFleming")
}
# if not stopped by here, test passed
print("Test passed: DTAdiscreteInterimAnalysis and cumulDiscreteInterimFleming agree")


# unit tests for simplifiedDiscreteInterimOutput
# check that the output matches what's expected
# use known data
# load test data
testdata <- readRDS("testData1.rds")
# these should nominally have Se= 65 and Sp = 85, prev=0.4 but won't 
# calculate continuous Se and Sp
N <- 100
prev <- 0.4
contshort <- continuousSeSp(testdata[1:N, ])


#suppress warnings because of final prevalence
suppressWarnings({
# full outputs
# early termination for high Se at k=37 and Sp at k=30
lowterm <- DTAdiscreteInterimAnalysis(contshort,c(30, 40, 50), pSe=0.55, pSp=0.7, prevalence = prev, N=N, simpleOutput = FALSE)
# early termination for low Se at k=63  and Sp at k=32
highterm <- DTAdiscreteInterimAnalysis(contshort,c(30, 50, 70), pSe=0.8, pSp=0.95, prevalence = prev, N=N, simpleOutput = FALSE)
# no termination by end point for either
noterm <- DTAdiscreteInterimAnalysis(contshort,c(15, 25), pSe=0.85, pSp=0.9, prevalence = prev, N=25, simpleOutput = FALSE)
# pure interim, no termination for either
interim <- DTAdiscreteInterimAnalysis(contshort,c(15, 25), pSe=0.85, pSp=0.9, prevalence = prev, N=N, simpleOutput = FALSE)


# simple outputs
# early termination for low Se at k=37 and Sp at k=30 
lowterms <- DTAdiscreteInterimAnalysis(contshort,c(30, 40, 50), pSe=0.55, pSp=0.7, prevalence = prev, N=N, simpleOutput = TRUE)
# early termination for high Se at k=63  and Sp at k=32
highterms <- DTAdiscreteInterimAnalysis(contshort,c(30, 50, 70), pSe=0.8, pSp=0.95, prevalence = prev, N=N, simpleOutput = TRUE)
# no termination by end point for either
noterms <- DTAdiscreteInterimAnalysis(contshort,c(15, 25), pSe=0.85, pSp=0.9, prevalence = prev, N=25, simpleOutput = TRUE)
# pure interim, no termination for either
interims <- DTAdiscreteInterimAnalysis(contshort,c(15, 25), pSe=0.85, pSp=0.9, prevalence = prev, N=N, simpleOutput = TRUE)
})

# low termination
#print("Testing termination for low Se/Sp")
# check terminateAt matches
if(lowterm$Sensitivity$terminationStage != lowterms$Sensitivity$terminateAt) 
  stop("Test failed: different termination stages in full and simple outputs")
if(lowterm$Specificity$terminationStage != lowterms$Specificity$terminateAt) 
  stop("Test failed: different termination stages in full and simple outputs")
# check futility == FALSE
if(lowterms$Sensitivity$futility!=FALSE) stop("Test failed: wrong futility in simple output")
if(lowterms$Specificity$futility!=FALSE) stop("Test failed: wrong futility in simple output")
# check string has correct number
if(length(grep(lowterm$Sensitivity$terminationStage, lowterms$Sensitivity$termstring))==0)
  stop("Test failed: string has incorrect number")
if(length(grep(lowterm$Specificity$terminationStage, lowterms$Specificity$termstring))==0)
  stop("Test failed: string has incorrect number")
# check string contains "efficacy"
if(length(grep("efficacy", lowterms$Sensitivity$termstring))==0) stop("Test failed: string conclusion is incorrect for efficacy")
if(length(grep("efficacy", lowterms$Specificity$termstring))==0) stop("Test failed: string conclusion is incorrect for efficacy")

# high termination
#print("Testing termination for high Se/Sp")
# check terminateAt matches
if(highterm$Sensitivity$terminationStage != highterms$Sensitivity$terminateAt) 
  stop("Test failed: different termination stages in full and simple outputs")
if(highterm$Specificity$terminationStage != highterms$Specificity$terminateAt) 
  stop("Test failed: different termination stages in full and simple outputs")
# check futility == TRUE
if(highterms$Sensitivity$futility!=TRUE) stop("Test failed: wrong futility in simple output")
if(highterms$Specificity$futility!=TRUE) stop("Test failed: wrong futility in simple output")
# check string has correct number
if(length(grep(highterm$Sensitivity$terminationStage, highterms$Sensitivity$termstring))==0)
  stop("Test failed: string has incorrect number")
if(length(grep(highterm$Specificity$terminationStage, highterms$Specificity$termstring))==0)
  stop("Test failed: string has incorrect number")
# check string contains "futility"
if(length(grep("futility", highterms$Sensitivity$termstring))==0) stop("Test failed: string conclusion is incorrect for futility")
if(length(grep("futility", highterms$Specificity$termstring))==0) stop("Test failed: string conclusion is incorrect for futillity")

# no termination
#print("Testing non termination")
# check terminateAt = NA
if(!is.na(noterms$Sensitivity$terminateAt))  stop("Test failed: termination in simple output is not NA")
if(!is.na(noterms$Specificity$terminateAt))  stop("Test failed: termination in simple output is not NA")
# check futility = NA
if(!is.na(noterms$Sensitivity$futility)) stop("Test failed: wrong futility in simple output")
if(!is.na(noterms$Specificity$futility)) stop("Test failed: wrong futility in simple output")
# check string does not contain futility or efficacy
if(length(grep("efficacy|futility", noterms$Sensitivity$termstring))>0) stop("Test failed: string conclusion is incorrect for not having futility or efficacy")
if(length(grep("efficacy|futility", noterms$Specificity$termstring))>0) stop("Test failed: string conclusion is incorrect for not having futility or efficacy")
# check string contains end
if(length(grep("final", noterms$Sensitivity$termstring))==0) stop("Test failed: string conclusion is incorrect at end")
if(length(grep("final", noterms$Specificity$termstring))==0) stop("Test failed: string conclusion is incorrect at end")

# pure interim
#print("Testing pure interim analysis")
# check terminateAt = NA
if(!is.na(interims$Sensitivity$terminateAt))  stop("Test failed: termination in simple output is not NA")
if(!is.na(interims$Specificity$terminateAt))  stop("Test failed: termination in simple output is not NA")
# check futility = NA
if(!is.na(interims$Sensitivity$futility)) stop("Test failed: wrong futility in simple output")
if(!is.na(interims$Specificity$futility)) stop("Test failed: wrong futility in simple output")
# check string does not contain futility, efficacy, Se or Sp
if(length(grep("efficacy|futility|Se|Sp", interims$Sensitivity$termstring))>0) stop("Test failed: string conclusion is incorrect")
if(length(grep("efficacy|futility|Se|Sp", interims$Specificity$termstring))>0) stop("Test failed: string conclusion is incorrect")
# check string contains interim
if(length(grep("interim", interims$Sensitivity$termstring))==0) stop("Test failed: string conclusion is incorrect")
if(length(grep("interim", interims$Specificity$termstring))==0) stop("Test failed: string conclusion is incorrect")

# if not stopped by here, test passed
print("Test passed: simplifiedDiscreteInterimOutput gives expected outputs")


# check that DTAdiscreteInterimAnalysis and DTAcumulativeInterimAnalysis give the same results
# run 10 times
for(i in seq(10)) {
  # create random data for Se and Sp with p of 0.8 and 0.6
  p0Se <- 0.65
  p0Sp <- 0.85
  prevalence <- 0.35
  N <- 200
  k <- 5
  analysispoints <- sort(sample(10:150, k))
  
  
  # create data
  tbt <- create2x2(n=10*N, p=prevalence, se=p0Se, sp=p0Sp, digits=4, verbose=FALSE)
  dtadata <- DTAdataFrom2x2(tbt, randomise = T)
  dtadata <- continuousSeSp(dtadata[1:N, ])
  
  actualprev <- sum(dtadata$reference)/N
  
  # run on interim analysis
  suppressWarnings(
    dtares <- DTAdiscreteInterimAnalysis(dtadata, analysispoints, pSe=p0Se, pSp=p0Sp, prevalence = actualprev, N=N, simpleOutput = FALSE))
  
  # pull out data for DTAcumulativeInterimAnalysis
  # want a data frame with the following:
  # N, RefT, TP, TN
  cumframe <- data.frame(N=numeric(),
                         RefT=numeric(), 
                         TP=numeric(),
                         TN=numeric())
  for(i in seq(length(analysispoints))) {
    k <- analysispoints[i]
    ifrm <- data.frame(N = k, 
                       RefT = sum(dtadata$reference[1:k]), 
                       TP = sum(dtadata$TP[1:k]),
                       TN = sum(dtadata$TN[1:k]))
    cumframe <- rbind(cumframe,ifrm)
  }
  
  # run on cumulative analysis
  suppressWarnings(
    cumres <- DTAcumulativeInterimAnalysis(cumframe, pSe=p0Se, pSp=p0Sp, prevalence = actualprev, N=N, simpleOutput = FALSE))
  
  # compare results
  if(any(dtares$Sensitivity$details!=cumres$Sensitivity$details, na.rm=T)) 
    stop("Test failed: different results from DTAdiscreteInterimAnalsysis and DTAcumulativeInterimAnalysis Se")
  if(any(dtares$Specificity$details!=cumres$Specificity$details, na.rm=T)) 
    stop("Test failed: different results from DTAdiscreteInterimAnalsysis and DTAcumulativeInterimAnalysis Sp")

}
# if not stopped by here, test passed
print("Test passed: DTAdiscreteInterimAnalysis and DTAcumulativeInterimAnalysis agree")


# testing plotting
# does it run?
plotSeSpFlemingThresholds(dtadata, 0.6, 0.9)
# if not stopped by here, test passed
print("Test passed: plotSeSpFlemingThresholds runs")
