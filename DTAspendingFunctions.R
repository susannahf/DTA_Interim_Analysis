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
# i, lower[i], upper[i], 1-pcross_up, pcross_low, alpha_u(i,nmax), alpha_l(i, nmax)

# load spending functions
source("SpendingFunctions.R")

STbounds <- function(p0, p1, alpha, power, nmax, smax){
  # things that need to be globally accessible happen here
  # this is kind of the equivalent of main() in Stallard and Todd
  
  # define globally accessible variables
  prob <- prob_last <- lowerbound <- upperbound <- rep(NA, nmax)
  prob0 <- prob_last0 <- lowerbound0 <- upperbound0 <- NA
  pcross_up <- pcross_low <- NA
  
  # spending functions
  spendfunc <- "Simple"

  # now we need to define all the subfunctions.
  # the actual function that runs the code is main(), which is called right at the end

  
  main <- function(p0, p1, alpha, power, nmax, smax, spend = "Simple") {
    spendfunc <<- spend
    
    # testing only
    getBoundary(nmax, smax)
    
    
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
    # 
    
    #get_power_and_n(nmax, smax, p0, power0, 1, output); #may not need this other than for output
    # this is the point where the output is given
  }
  
  ## this is the end of the main calculations
  ## from here on are the sub-functions
  
  # calculate upper and lower boundaries
  # equivalent of Stallard and Todd's get_boundary function
  getBoundary <- function(nmax, smax){
    # for n=0, boundaries are at their extremes
    lowerbound0 <<- -1 
    upperbound0 <<- 1
    
    # initialise prob[1:smax] and prob0
    initProbs(smax)

    for(i in 1:nmax) {
      # copy prob to prob_last for 0:smax
      updateProb(smax)
      # get upper and lower boundaries at inspection i
      getBoundsAt_i(i, nmax, smax)
    }
    
    # get_power_and_n(nmax, smax, p1, power, 0, output) <- needed for the power search
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

    #print(paste("setting bounds at",i))
    # p=P0; <- lets just use p0 to make this clear....
    stopflag=FALSE; # should the search stop
    upset=FALSE; # is upper boundary set for this i
    lowset=FALSE;  # is lower boundary set for this i
    probsi = 0 # I think probsi is p_n(s) in the paper
    
    # if i=1, then the crossing probabilities must be 0
    if(i==1) {
      pcross_low<<-0
      pcross_up<<-0 
      loopStart<-lowerbound0 # loop from lowerbound[i-1], but if i=1, this is lowerbound0
    } else {
      loopStart<-lowerbound[i-1] # loop from this
    }
    
    
    # loop si from lowerbound[i-1] to i+1
    # note that loopStart could be -1 (initial setting of lowerbound0)
    for(si in loopStart:(i+1)) {  
      # stop if necessary (this is if everything is set)
      if(stopflag) break
      # calculate probability of getting to (si,i) with P(success)=p
      # I think probsi is p_n(s) in the paper, in which case get_prob uses "a simple path counting method"
      #/* find prob of getting to si */
      if(si>0) {
        getProb(i, si, p0);
        probsi<-prob[si];
      } else if(si==0) {
        getProb(i, si, p0)
        probsi<- prob0
      } else { probsi<-0 } # prob is 0 if si is <0
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
      # check to see if P(Si>si) is <= 1-alpha_u */
      if (si==loopStart) pcross_up <<- pcross_low; # if si is on the lower boundary, set pcross_up ???
     
      if (!upset)  # if upper boundary not set
      {
        #print(paste0("i:", i," si:", si, " pcross up:", pcross_up, " probsi:", probsi))
        pcross_up <<- pcross_up + probsi; # increment pcrosslow2 by probsci 
        
        # the test below is eq 9 from paper, iff pcross_up is sum(p_L, 1 to n-1) and probsi is sum(p_n, 0 to u_n -1)
        # if si is u_n -1 above, then the below makes sense, because u_n is the upper boundary.
        # if the test is met, then upper bdy is set to si+1 (=u_n) and upset is set
        #print(paste0("si:", si, " pcross up:", pcross_up, ", 1-alphau:", 1-alpha_u(i, nmax)))
        if ( pcross_up >= 1-alpha_u(i,nmax) )  
        {
          upperbound[i]<<-si+1;
          upset=T;
        }
      }
      #         
      # this should implement eq 10 from paper.
      #        /* and to see if P(Si<si) is > alpha_l */
      if (!lowset) # if lower boundary not set
      {
        pcross_low <<- pcross_low + probsi; # increment pscrosslow by probsci <- this creates the next value of p_l
        # the test below is eq 10 from paper, iff pcross_low is sum(p_L, 1 to n-1) and probsi is sum(p_n, 0 to l_n)
        # if si is ln above, then the below makes sense, because l_n is the lower boundary
        # and we are testing for > rather than <= so we want the boundary to be si-1
        # if the test is met, then lower bdy is set to si-1 and lowset is set
        # pcross_low is set to the probability -probsi, which is the p_l for l_n
        #print(paste0("i", i, "si:", si, " pcross low:", pcross_low, ", alphal:", alpha_l(i, nmax)))
        if ( pcross_low > alpha_l(i,nmax) )
        {
          lowerbound[i] <<-si-1;
          pcross_low <<- pcross_low -probsi;
          lowset=T;
        }
      }
      if (lowset && upset) stopflag=T; # stop if both set
    }
    # if lowset or upset isn't set...
    if(!lowset) lowerbound[i] <<- lowerbound[i-1]
    if(!upset) upperbound[i] <<- upperbound[i-1]+1
    #   else # this deals with situations where this i isn't divisible by group
    #   {
    #     upper[i]=upper[i-1]+1;
    #     lower[i]=lower[i-1];
    #   }   
  }
  

  
  # initialise prob[1:smax] and prob0
  # implements Stallard and Todd's init_probs(smax)
  initProbs <- function(smax) {
    # initialises probability vector up to smax
    # all element of prob up to smax = 0, prob0 = 1
    prob[1:smax] <<- 0
    prob0 <<- 1
  }
  
  # copy prob to prob_last for 0 to smax
  # implements Stallard and Todd's update_prob(smax)
  updateProb <- function(smax){
    # copy prob to prob_last up to smax
    prob_last[1:smax] <<- prob[1:smax]
    prob_last0 <<- prob0
  }
  
  # finds probability of getting to (si,i) with P(success)=p
  # I think this is p_n(s) in the paper, which uses "a simple path counting method"
  # implements Stallard and Todd's get_prob(i, si, p)
  getProb <- function(i, si, p) {
    #initialise prob[si]
    probsi <<- 0
    
    # handle boundaries
    if(i==1) {
      lastlowerbound <- lowerbound0
      lastupperbound <- upperbound0
      lastproblast <- prob_last0
    } else {
      lastlowerbound <-lowerbound[i-1]
      lastupperbound <- upperbound[i-1]
      lastproblast <- prob_last[i-1]
    }
    if(si==0) problastsi = prob_last0
    else problastsi = prob_last[si]
    
    if(si>0 && (si-1)>lastlowerbound && (si-1)<lastupperbound) { # si-1 did not cross bdy at i-1
      probsi <- probsi + lastproblast*p }
    if(si>lastlowerbound && si<lastupperbound) { # si did not cross bdy at i-1
      probsi <- probsi + problastsi*(1-p) }
    
    if(si==0)
      prob0 <<- probsi
    else
      prob[si] <<-probsi
    
  }
  
  # implements Stallard and Todd's get_power_and_n(nmax, smax, p, *results, printall, output)
  # void get_power_and_n(int nmax, int smax, double p, double *results, int printall, FILE* output)
  # /* calls get_prob with p to get power for boundary obtained previously */
  #   {
  #     int i, si;
  #     double power, n, probsi;
  #     double pcross_low, pcross_low2, pcross_low_old, pcross_low2_old;
  #     
  #     power=0.;
  #     n=0.;
  #     pcross_low=0.;
  #     pcross_low2=1.;
  #     
  #     init_probs(smax);
  #     
  #     for (i=1; i<=nmax; i++)
  #     {
  #       update_prob(smax);
  #       pcross_low_old=pcross_low;
  #       pcross_low2_old=pcross_low2;
  #       for (si=lower[i-1]; si<upper[i]; si++)
  #       {
  #         if (si>=0) 
  #         {
  #           get_prob(i, si, p);
  #           probsi=prob[si];
  #         }
  #         else
  #           probsi=0.;
  #         
  #         if (si==lower[i-1])
  #         {
  #           pcross_low2=pcross_low;
  #         }
  #         pcross_low2+=probsi;
  #         if (si<=lower[i]) pcross_low+=probsi;
  #       }
  #       power+=(pcross_low2_old-pcross_low2);
  #       n+=i*(pcross_low2_old-pcross_low2+pcross_low-pcross_low_old);
  #       if (printall==1)
  #       {
  #         fprintf(output,"%i %i %i  %g %g",
  #                 i,lower[i],upper[i],1.-pcross_low2,pcross_low);
  #         fprintf(output,"           %g %g\n",
  #                 alpha_u(i,nmax),alpha_l(i,nmax));
  #       }
  #     }
  #     results[0]=power;
  #     results[1]=n;
  #   }
  
  # alpha_l calculates the spending function for the lower bound
  # using functions from SpendingFunctions.R
  alpha_l <- function(i, nmax) {
    t <- i/nmax
    if(spendfunc=="Triangular")
      return(TriangularTest_alpha_l(alpha,t))
    else if(spendfunc=="Simple")
      return(KimDeMets_alpha_l(alpha,t))
    else
      stop(paste("unrecognised spending function: ",spendfunc))
  }
  
  # alpha_u calculates the spending function for the upper bound
  # using functions from SpendingFunctions.R
  alpha_u <- function(i, nmax) {
    t <- i/nmax
    if(spendfunc=="Triangular")
      return(TriangularTest_alpha_u(alpha,t))
    else if(spendfunc=="Simple")
      return(KimDeMets_alpha_u(alpha,t))
    else
      stop(paste("unrecognised spending function: ",spendfunc))
  }
  
  # now that we've defined all the functions, we can call main()
  main(p0, p1, alpha, power, nmax, smax)

  # output all global variables 
  # (initially for testing, but probably useful in the end too)
  # define globally accessible variables
  outlist <- list()
  outlist$lowerbound <- c(lowerbound0, lowerbound)
  outlist$upperbound <- c(upperbound0, upperbound)
  outlist$prob <- c(prob0, prob)
  outlist$prob_last <- c(prob_last0, prob_last)
  
  return(outlist)
  
}














