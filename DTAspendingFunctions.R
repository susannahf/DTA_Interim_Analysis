# code for running spending function methods for DTA interim analysis
# S Fleming 2023

# methods from Stallard/Todd 2000 and Lan/DeMets1983

# try both triangular and simpler Kim/DeMets spending functions
# this is re-implemented from Stallard and Todd's original c code
# All mistakes are mine.



# REMEMBER, c vectors start at 0, and R vectors start at 1!

# output from c program is:
# i, lower[i], upper[i], 1-pcross_low2, pcross_low, alpha_u(i,nmax), alpha_l(i, nmax)

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

# equivalent of Stallard and Todd's get_bounds_i()
# searches for upper and lower boundaries at inspection i to satisfy spending function requirements
# needs get_prob, alpha_l, alpha_u


# equivalent of Stallard and Todd's get_prob()
# finds probability of getting to (si,i) with P(success)=p

# equivalent of Stallard and Todd's alpha_l()
# spending function for triangular test
# needs set_std_para_grp, prob_cross_lwr

# equivalent of Stallard and Todd's alpha_u()
# spending function for triangular test
# needs set_std_para_grp, prob_cross_lwr

# equivalent of Stallard and Todd's set_std_para_grp()
# sets some constants afaict

# equivalent of Stallard and Todd's prob_cross_lwr()
# calculates the probability of crossing the lower boundary (p_L?)
# needs integrate, mills_ratio

# equivalent of Stallard and Todd's mills_ratio()
# calculates R(x) for the triangular test

# equivalent of Stallard and Todd's integrate()
# integrates... I guess. May be able to use an internal R function, like, um, integrate()?
# needs gauleg

# equivalent of Stallard and Todd's gauleg()
# some sort of integration helper, probably don't need

















# internal testing (to remove later)

nmax=smax=9999
p0=0.003
p1=0.006 # we may not need this
power=0.8 # we may not need this

