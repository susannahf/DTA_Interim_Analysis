# code for running spending function methods for DTA interim analysis
# S Fleming 2023

# methods from Stallard/Todd 2000 and Lan/DeMets1983

# try both triangular and simpler Kim/DeMets spending functions
# this is re-implemented from Stallard and Todd's original c code
# All mistakes are mine.



# REMEMBER, c vectors start at 0, and R vectors start at 1!
# so any initial starting points are coded here as var0 rather than var[0]
# this will need a bit of care, but I think is safer than doing everything as i+1

# there are a lot of global variables being set, which I don't like, so I'm going to have an uber-function
# which will allow them to be in lexical scope
# they will then be able to be accessed as normal, and set using <<-

# output from c program is: (starting with i=1) is from get_power_and_n
# i, lower[i], upper[i], 1-pcross_low2, pcross_low, alpha_u(i,nmax), alpha_l(i, nmax)

uberfunction <- function(){
  # things that need to be globally accessible happen here
  # this is kind of the equivalent of main() in Stallard and Todd
  
  # define globally accessible variables
  prob <- prob_last <- numeric()
  
  
  
  # do some stuff that ultimately calls getBoundary()
  
  
  
  # equivalent of Stallard and Todd's get_boundary function
  getBoundary <- function(nmax, smax){
    # initialise lower and upper boundaries
    lowerbound <- upperbound <- rep(NA, nmax)
    # for n=0, their boundaries are at their extremes
    lowerbound0 = -1 # this makes sure i matches
    upperbound0 = 1
    
    # init_probs(smax)
    # equivalent of Stallard and Todd's init_probs(smax)
    # initialises probability vector
    # prob is a vector of size smax+1 with all elements 0 except prob[0]=1
    prob <<- rep(0,smax)
    prob0 = 1
    
    for(i in 1:nmax) {
      # update_prob(smax)
      # equivalent of Stallard and Todd's update_prob(smax)
      # copy prob to prob_last
      prob_last <<- prob
      # get upper and lower boundaries at inspection i
      getBoundsAt_i(i, nmax, smax)
    }
    
    # get_power_and_n(nmax, smax, p, power, 0, output) <- probably not necessary for me
  }
  
  # equivalent of Stallard and Todd's get_bounds_i()
  # searches for upper and lower boundaries at inspection i to satisfy spending function requirements
  # /* searches for upper and lower boundaries at inspection i to satisfy 
  # spending function requirements 
  # 
  # The calculations are based on the distribution of the number of successes
  # at each stage given the stopping boundary thus far.  These probabilities
  # (given by prob[si]) are found recursively by get_prob.  By running over the loop
  # for si, the full distribution for this inspection is obtained.  The function update_prob 
  # then stores this as last_prob so that the recursion can continue
  # */
  # needs get_prob, alpha_l, alpha_u
  getBoundsAt_i <- function(i, nmax, smax) {
    
    # if i=1, then the crossing probabilities must be 0
    if(i==1) {
      pcross_low=0
      pcross_low2=0 # not yet sure what these two are
    }
    
    # loop from lower[i-1] to i+1, stopping if stop=1
      # if si<0, probsi=0
      # if si>=0
        ####get_prob(i,si,p0)
        ## probsi = prob[si]
    
    ### more here
    
    
  }
  
  # equivalent of Stallard and Todd's get_prob(i, si, p)
  # finds probability of getting to (si,i) with P(success)=p
  # I think this is p_n(s) in the paper, which uses "a simple path counting method"
  # set prob[si] = 0
  # if si>0 and si-1>lower[i-1] and si-1<upper[i-1] # si-1 did not cross bdy at i-1
  # prob[si] = prob[si] + prob_last[si-1]*p
  # if si>lower[i-1] and si<upper[i-1] # si did not cross bdy at i-1
  # prob[si] = prob[si] + prob_last[si]*(1-p)
  
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
  
  
  
  
}












# internal testing (to remove later)

nmax=smax=9999
p0=0.003
p1=0.006 # we may not need this
power=0.8 # we may not need this

