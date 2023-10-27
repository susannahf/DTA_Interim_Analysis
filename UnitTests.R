# unit tests
# not everything is unit tested, but any new code from 24 Oct 2023 should be.
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
  # cumulative
  cdata <- cumsum(rdata)
  totaln <- 100
  testn <- seq(100)
  noncumn <- rep(1,100)
  k <- 100
  p0 <- 0.4
  
  discdata <- discreteInterimFleming(k=k, alpha=0.05, p0=p0, nk=noncumn, Ek=rdata)
  cumdata <- cumulDiscreteInterimFleming(ns=testn, events=cdata, finaln=totaln, p0=p0, alpha=0.05)
  
  # tests - these should be identical
  disctst <- discdata$details[1:99, ]
  cumtst <- cumdata$details[1:99, ]
  if(any(disctst!=cumtst)) stop("Test failed: different results from discreteInterimFleming and cumulDiscreteInterimFleming")
}
# if not stopped by here, test passed
print("Test passed: discreteInterimFleming and cumulDiscreteInterimFleming agree")

# unit tests for DTAdiscreteInterimAnalysis
# compare DTAdiscreteInterimAnalysis with cumulDiscreteInterimFleming

# run 10 times
for(i in seq(10)) {
  # create random data for Se and Sp with p of 0.8 and 0.6
  p0Se <- 0.8
  p0Sp <- 0.6
  # prevalence = 0.2, N=100, so 20 Se data points and 80 Sp data points
  prevalence <- 0.2
  N <- 100
  NSe <- prevalence*N
  NSp <- N-NSe
  TP <- rbinom(NSe,1,p0Se)
  TN <- rbinom(NSp,1, p0Sp)
  # pull out data for cumulDiscreteInterimFleming and run it
  SeNs <- c(seq(NSe), rep(NSe,NSp))
  SpNs <- c(rep(0,NSe), seq(NSp))
  SeEvents <- c(cumsum(1-TP), rep(sum(1-TP),NSp))
  SpEvents <- c(rep(0,NSe), cumsum(1-TN))
  SeFleming <- cumulDiscreteInterimFleming(ns=SeNs, events=SeEvents, finaln=NSe, p0=1-p0Se, alpha=0.05)
  SpFleming <- cumulDiscreteInterimFleming(ns=SpNs, events=SpEvents, finaln=NSp, p0=1-p0Sp, alpha=0.05)
  
  # pull out data for DTAdiscreteInterimAnalysis and run it
  dtadata <- data.frame(reference=c(rep(TRUE, NSe), rep(FALSE, NSp)),
                        index=c(as.logical(TP), !as.logical(TN)))
  dtadata <- continuousSeSp(dtadata)
  dtares <- DTAdiscreteInterimAnalysis(dtadata,seq(N), pSe=p0Se, pSp=p0Sp, prevalence = prevalence, N=N, simpleOutput = FALSE)
  # compare results
  SeFlemingTest <- SeFleming$details
  SpFlemingTest <- SpFleming$details
  dtaSeTest <- dtares$Sensitivity$details
  dtaSpTest <- dtares$Specificity$details
  
  if(any(SeFlemingTest!=dtaSeTest)) stop("Test failed: different results for Se from DTAdiscreteInterimAnalsysis and cumulDiscreteInterimFleming")
  if(any(SpFlemingTest!=dtaSpTest)) stop("Test failed: different results for Sp from DTAdiscreteInterimAnalsysis and cumulDiscreteInterimFleming")
}
# if not stopped by here, test passed
print("Test passed: DTAdiscreteInterimAnalysis and cumulDiscreteInterimFleming agree")
  