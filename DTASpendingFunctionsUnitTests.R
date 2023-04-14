# unit testing for DTAspendingFunctions.R
# S Fleming 2023

# because re-implementing a huge C code base without testing is just stupid
# and I finally realised that

# Test 1: will it even source?
test1 <- try(source("DTAspendingFunctions.R"))
if(!inherits(test1, "try-error")) {
  print("Test1 passed: can source DTAspendingFunctions.R")
}

# reasonable test values
nmax=smax=9999
p0=0.003
p1=0.006 # we may not need this
power=0.8 # we may not need this
alpha=0.05

# Test 2: does the code even run?
test2 <- try(uberfunction(p0, p1, alpha, power, nmax, smax))
if(!inherits(test2, "try-error")) {
  print("Test2 passed: uberfunction runs.")
}

