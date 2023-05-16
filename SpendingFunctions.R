# implement various spending function forms
# S Fleming 2023

# methods from Stallard/Todd 2000 
# all mistakes are my own

# simple kim and Demets approach with rho=1

# this implements equation 5 in Stallard/Todd 2000
KimDeMets_alpha_u <- function(alpha, t){
  # check that alpha and t are both between 0 and 1
  if(any(alpha<0 | alpha>1 | t<0 | t>1)) {
    error("alpha and t must both be between 0 and 1")
  }
  return(alpha*t)
}

# this implements equation 6 in Stallard/Todd 2000
KimDeMets_alpha_l <- function(alpha, t){
  # check that alpha and t are both between 0 and 1
  if(any(alpha<0 | alpha>1 | t<0 | t>1)) {
    error("alpha and t must both be between 0 and 1")
  }
  return((1-alpha)*t)
}

# the horrendously complicated Stallard/Todd approach

# helper function for Mills Ratio (R(x) in the paper)
MillsRatio <- function(x) {
  numer <- 1-pnorm(x)
  denom <- dnorm(x)
  return(numer/denom)
}

# this implements equation 3 in Stallard/Todd 2000
# note that this function can't cope with vector t
TriangularTest_alpha_u <- function(alpha, t){
  # check that alpha and t are both between 0 and 1
  if(alpha<0 | alpha>1 | t<0 | t>1) {
    error("alpha and t must both be between 0 and 1")
  }
  # check they're both length 1
  if(length(alpha)>1 | length(t)>1)
    error("alpha and t must both be single numbers")

  # if t > 0.97, return alpha
  #if(t>0.9999) return(alpha)
  # not needed for this implementation
  
  # constant
  c <- sqrt(-1 * log(2*alpha) / 2)
  #print(paste0("c=",c))
  
  # first term
  term1arg <- c * (1- 3*t) / sqrt(t)
  term1 <- dnorm(term1arg) # lower case phi is std normal density
  #print(paste0("term1=",term1))
  
  # handle where term1 is zero already...
  if(term1 == 0) return(0)
  
  # helper function to calculate the thing that's summed
  sumarg <- function(i, c, t) {
    # first term
    term1 <- ((-1)^i) * exp(2 * i * (i-1) * (c^2) * (1- (1/t)) )
    # second term
    term2 <- MillsRatio( c * (2*i + 1 - 3*t) / sqrt(t) )
    # third term
    term3 <- MillsRatio( c * (2*i + 1 + 3*t) / sqrt(t) )
    # whole thing
    return(term1 * (term2 + term3))
  }
    
  # sum term
  maxi = 500 # changing maxi doesn't seem to affect the oddness 
  is <- 0:maxi 
  tolerance <- 1e-11 # nor does changing the tolerance
  sumterms <- rep(0, maxi+1)
  for(i in is) {
    sumterms[i+1] <- sumarg(i, c, t)
    if(abs(sumterms[i+1]) <tolerance) break # stop the loop when convergence
  }
  #print("sumterms:")
  #print(sumterms)
  term2 <- sum(sumterms)
  #print(paste0("term2=",term2))
  
  
  return(term1 * term2)
}