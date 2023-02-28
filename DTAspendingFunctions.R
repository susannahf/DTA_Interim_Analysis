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

uberfunction <- function(p0, p1, alpha, power, nmax, smax){
  # things that need to be globally accessible happen here
  # this is kind of the equivalent of main() in Stallard and Todd
  
  # define globally accessible variables
  prob <- prob_last <- numeric()
  
  
  
  # call getBoundary() in a searching loop to find nmax such that nmax<NMAX+1 
  # and power is closest to specified power (but a bit less)
  # this is needed to get the defined power because the lower boundary is defined by (1-alpha) and not beta
  # need to determine what nmax actually is...
  # /* search starts looking every GINC and is refined to GINC/10  ... 1 */
  #   for (g_inc=GINC; g_inc>=1; g_inc/=10)
  #   {
  #     if (g_inc==GINC) g_start=GSTART;
  #     else g_start=(int)(nmax/GROUP)-19*g_inc;
  #     power1[0]=0.;
  #     for (nmax=g_start*GROUP; nmax<NMAX_PLUS_1 && power1[0]<POWER;
  #          nmax+=g_inc*GROUP)
  #     {
  #       get_boundary(nmax, smax, power1, output);
  #       printf("nmax=%i: power at p=%f = %f\n", nmax,P1,power1[0]);
  #     }
  #   }
  # 
  # nmax-=GROUP;
  
  #get_power_and_n(nmax, smax, p0, power0, 1, output); #may not need this other than for output
  # this is the point where the output is given
  
  
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
    
    # get_power_and_n(nmax, smax, p, power, 0, output) <- needed for the power search
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
    
    
      # this is the  c code:
      # need to know what these stand for...
    # p=P0;
    # stop=0;
    # upset=0; # is upper boundary set for this i
    # lowset=0;  # is lower boundary set for this i
     
      # if i=1, then the crossing probabilities must be 0
    if(i==1) {
      pcross_low=0
      pcross_low2=0 # not yet sure what these two are
    }
    
      # loop si from lower[i-1] to i+1, incrementing, stop if stop==1
      
    # for (si=lower[i-1]; stop!=1 && si<=i+1; si++)
    # {
      # calculate probability of getting to (si,i) with P(success)=p
      # I think probsi is p_n(s) in the paper, in which case get_prob uses "a simple path counting method"
    #   /* find prob of getting to si */
    #     if (si>=0) 
    #     {
    #       get_prob(i, si, p);
    #       probsi=prob[si];
    #     }
      # prob is 0 if si is <0
    #   else
    #     probsi=0.;
    #   
    #   I think this bit below (only checking every group) is not relevant for our case
    #   /* boundaries for binomial data with GROUP>1 are given by calculations
    #   as for Bernoulli case except that stopping is only made possible
    #   every GROUP observations - when the code below is executed.
    #   For all other i, the lower boundary remains the same and the
    #   upper boundary is increased by 1 to make stopping impossible
    #   */
    #     if (i%GROUP==0) 
    #     {
    #       /* check to see if P(Si>si) is <= 1-alpha_u */
    #         if (si==lower[i-1]) pcross_low2=pcross_low; # if si is on the lower boundary, set pcross_low2 ???
    #         if (upset==0)  # if upper boundary not set
    #         {
    #           pcross_low2+=probsi; # increment pcrosslow2 by probsci 
      # the test below is eq 9 from paper, iff pcross_low2 is sum(p_L, 1 to n-1) and probsi is sum(p_n, 0 to u_n -1)
      # if si is u_n -1 above, then the below makes sense, because u_n is the upper boundary.
      # if the test is met, then upper bdy is set to si+1 (=u_n) and upset is set
    #           if ( pcross_low2 >= 1.-alpha_u(i,nmax) )  
    #           {
    #             upper[i]=si+1;
    #             upset=1;
    #           }
    #         }
    #         
      # this should implement eq 10 from paper.
    #         /* and to see if P(Si<si) is > alpha_l */
    #             if (lowset==0) # if lower boundary not set
    #             {
    #               pcross_low+=probsi; # increment pscrosslow by probsci <- this creates the next value of p_l
      # the test below is eq 10 from paper, iff pcross_low is sum(p_L, 1 to n-1) and probsi is sum(p_n, 0 to l_n)
      # if si is ln above, then the below makes sense, because l_n is the lower boundary
      # and we are testing for > rather than <= so we want the boundary to be si-1
      # if the test is met, then lower bdy is set to si-1 and lowset is set
      # pcross_low is set to the probability -probsi, which is the p_l for l_n
    #               if ( pcross_low > alpha_l(i,nmax) )
    #               {
    #                 lower[i]=si-1;
    #                 pcross_low-=probsi;
    #                 lowset=1;
    #               }
    #             }
    #         if (lowset==1 && upset==1) stop=1; # stop if both set
    #     }
    #   else # this deals with situations where this i isn't divisible by group
    #   {
    #     upper[i]=upper[i-1]+1;
    #     lower[i]=lower[i-1];
    #   }   
    # }
    
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

