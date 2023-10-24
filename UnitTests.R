# unit tests
# not everything is unit tested, but any new code from 24 Oct 2023 should be.
source("DTAinterimAnalysis.R")

# unit tests for cumulDiscreteInterimFleming
# TODO
# run both DiscreteInterimFleming and cumulDiscreteInterimFleming on the data, not including k=n
# should have the same rg, ag, and so on in details of output

# create random non-cumulative event data for n=100
rdata <- rbinom(100,1,0.4)
# cumulative
cdata <- cumsum(rdata)
totaln <- 100
testn <- seq(99)
noncumn <- rep(1,99)
k <- 99
p0 <- 0.4

discdata <- discreteInterimFleming(k=k, alpha=0.05, p0=p0, nk=noncumn, Ek=rdata)
cumdata <- cumulDiscreteInterimFleming(ns=testn, events=cdata, finaln=totaln, p0=p0, alpha=0.05)
# doesn't work now