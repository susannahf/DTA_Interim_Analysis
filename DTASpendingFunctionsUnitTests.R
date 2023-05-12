# unit testing for DTAspendingFunctions.R
# S Fleming 2023

# because re-implementing a huge C code base without testing is just stupid
# and I finally realised that

print("")
print("Testing DTASpendingFunctions.R")

# Test 1: will it even source?
test1 <- try(source("DTAspendingFunctions.R"))
if(!inherits(test1, "try-error")) {
  print("Test1 passed: can source DTAspendingFunctions.R")
} else { print("Test 2 failed.")}

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
} else { print("Test 2 failed.")}

# now let's worry about the output
#testout <- uberfunction(p0, p1, alpha, power, nmax, smax)

#############################################

print("")
print("Testing SpendingFunctions.R")

# also test the spending functions themselves
# Test 1: will it even source?
test1 <- try(source("SpendingFunctions.R"))
if(!inherits(test1, "try-error")) {
  print("Test 1: can source SpendingFunctions.R")
} else { print("Test 1 failed.")}

# reasonable test values
alpha=0.05
t <- seq(0,1,0.01)

# Test 2: Kim and DeMets alpha_u
test2 <- try(KimDeMets_alpha_u(alpha, t))
if(!inherits(test1, "try-error")) {
  print("Test 2: can calculate Kim and DeMets alpha_u")
  plot(t,test2,'l',main="Test 2: Kim/Demets alpha_u")
} else { print("Test 2 failed.")}

if(test2[t==0] == 0) {
  print("Test 2: alpha_u(0)=0")
} else  print("Test 2 failed.")

if(test2[t==1] == alpha) {
  print("Test 2: alpha_u(1)=alpha")
} else  print("Test 2 failed.")

# Test 3: Kim and DeMets alpha_l
test3 <- try(KimDeMets_alpha_l(alpha, t))
if(!inherits(test1, "try-error")) {
  print("Test 3: can calculate Kim and DeMets alpha_l")
  plot(t,test3,'l',main="Test 3: Kim/Demets alpha_l")
} else { print("Test 3 failed.")}

if(test3[t==0] == 0) {
  print("Test 3: alpha_l(0)=0")
} else  print("Test 3 failed.")

if(test3[t==1] == (1-alpha)) {
  print("Test 3: alpha_l(1)= 1-alpha")
} else  print("Test 3 failed.")


