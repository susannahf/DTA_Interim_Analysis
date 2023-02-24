# code for running spending function methods for DTA interim analysis
# S Fleming 2023

# methods from Stallard/Todd 2000 and Lan/DeMets1983

# try both triangular and simpler Kim/DeMets spending functions
# this is re-implemented from Stallard and Todd's original c code
# All mistakes are mine.



# REMEMBER, c vectors start at 0, and R vectors start at 1!

# equivalent of Stallard and Todd's get_boundary function
getBoundary <- function(nmax, smax){
  # initialise lower and upper boundaries
  lowerbound <- upperbound <- rep(NA, nmax)
  # for n=1, there boundaries are at their extremes
  lowerbound[1] = -1
  upperbound[1] = 1
  
  for(i in 2:(nmax+1)) {
    # update_probs(smax)
    # get_bounds_i(i, nmax, smax)
  }
  
  # get_power_and_n(nmax, smax, p, power, 0, output) <- probably not necessary for me
}



















# internal testing (to remove later)

nmax=smax=9999
p0=0.003
p1=0.006 # we may not need this
power=0.8 # we may not need this

