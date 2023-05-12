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
  if(t>0.97) return(alpha)
  
  # constant
  c <- sqrt(-1 * log(2*alpha) / 2)
  
  # first term
  term1arg <- c * (1- 3*t) / sqrt(t)
  term1 <- dnorm(term1arg) # lower case phi is std normal density
  
  # helper function to calculate the thing that's summed
  sumarg <- function(i, c, t) {
    # first term
    term1 <- ((-1)^i) * exp(2 * i * (i-1) * (c^2) * (1- (1/t)) )
    # second term
    term2 <- MillsRatio( c * (2*i + 1 - t) / sqrt(t) )
    # third term
    term3 <- MillsRatio( c * (2*i + 1 + t) / sqrt(t) )
    # whole thing
    return(term1 * (term2 + term3))
  }
    
  # sum term
  maxi = 500
  is <- 0:maxi
  sumterms <- sapply(is, sumarg, c=c, t=t)
  term2 <- sum(sumterms)
  
  return(term1 * term2)
}