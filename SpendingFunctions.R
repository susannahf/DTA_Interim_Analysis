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