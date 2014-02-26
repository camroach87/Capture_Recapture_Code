#### Preamble ####

# PROGRAM: Closed and open population simulation CODER: Cameron Roach DATE:
# 2/4/2013 DESCRIPTION: Simulates what is actually happening in closed and open 
# populations and then compares estimators for each. Assumes sampling occasions 
# are equally spaced. OpenSim uses a super population that operates as a set of
# fish that may or may not enter/exit the current population between sampling
# occasions. While the size of the superpopulation does affect the number of
# fish born/passed away, so too does the input for pDeath and pBirth which can
# be changed to get suitable numbers for fish born/deceased.
# 
# NOTES: 1. Currently, does not take into consideration the marked population
# size approaching the total population size. May cause errors to be produced
# once marked population exceeds total population.
# 2. Superpopulation is 10 times the size of actual population in closed model.
# If we have significantly more births than deaths then this will be an issue if
# population goes above superpopulation limit.
# 
# FIX: Open and closed sims not using the same Fish prob matrix
# ERROR: birth and death matrices not considering fish that have already been born/died.


#### Console commands ####
# These are not part of the program but are useful

#for (i in 1:100) {bla<-bla+ClosedSim(N_0, t, pCapture)}; bla/100
#for (i in 1:100) {bla<-bla+OpenSim(N_0, t, pCapture)}; bla/100
Initialise <- function() {
  N_0 <-1000
  t <- 10
  pBirth <- 0.03   
  pDeath <- 0.02
  pCapture <- 0.05
  N_final <- OpenSim(N_0, t, pCapture, pBirth, pDeath)
  ClosedSim(N_final, t, pCapture)
}

#### Code #####

ClosedSim <- function(N_0, t, pCapture) {
  
  # Initialise variables
  n <- rep(NA, t)
  m <- rep(NA, t)
  M <- rep(NA, t)
  
  
  
  # Setup fish matrix
  # Assume capture probabilities unifiormly distributed on each sampling occasion
  Fish <- matrix(runif(N_0*t), ncol=t)
  FishCaught <- ifelse(Fish<=pCapture,1,0)
  fVector <- cbind(apply(FishCaught,1,sum))
  
  # Chao estimator for time variation model M_t
  f1 <- sum(ifelse(fVector==1,1,0))
  f2 <- sum(ifelse(fVector==2,1,0))
  S <- sum(ifelse(fVector>=1,1,0))
  
  
  Z <- FishCaught[fVector==1,]
  Z <- colSums(Z)
  
  
  N = S + (f1^2-sum(Z^2))/(2*f2)
  
  return(N)
  
}


OpenSim <- function(N_0, t, pCapture, pBirth, pDeath) {
  # Simulates open population. N_0 is initial population size.
  # Assumes each animal only capable of one offspring per sampling interim.
  # Assumes death equally likely at any point in time of fish's lifespan
  
  # Initialise variables
  N <- rep(NA, t)
  n <- rep(NA, t)
  m <- rep(NA, t)
  M <- rep(NA, t)
  nbirths <- rep(0, t)
  ndeaths <- rep(0, t)
  Nsuper <- 10*N_0
  FishExists <- matrix(rep(0,Nsuper*t), ncol=t)
  FishGivesBirth <- matrix(runif(Nsuper*t), ncol=t)
  FishDies <- matrix(runif(Nsuper*t), ncol=t)
  
  N[1] <- N_0
  NewFish <- N_0
  FishExists[1:N_0,1] <- 1
  for (i in 2:t) {
    
    # Use 1-pBirth as probability since non-existing fish have 0 used to identify them
    # i.e. there will be lots of zeros when taking inner product of FishExists and FishGivesBirth
    births <- ifelse(FishExists[1:NewFish,i-1]*FishGivesBirth[1:NewFish,i-1]>=(1-pBirth),1,0)
    nbirths[i] <- sum(births)
    FishExists[(NewFish+1):(NewFish+nbirths[i]),i] <- 1
    
    deaths <- ifelse(FishExists[1:NewFish,i-1]*FishDies[1:NewFish,i-1]>=(1-pBirth),1,0)
    ndeaths[i] <- sum(deaths)
    FishExists[1:NewFish,i] <- FishExists[1:NewFish,i-1]*(1-deaths)
    
    N[i] <- N[i-1] + nbirths[i] - ndeaths[i]
    
    NewFish <- NewFish + nbirths[i]
  }

  plot(1:t,N)
  
  
  return(N[t])
  
  
  
  
  
  
  
#   # Simulate birth death process at each time step to get actual population.
#   N[1] <- N_0
#   for (i in 2:t) {
#     # Simulate births
#     Birth <- t(runif(N[i-1]))
#     FishBirth <- ifelse(Birth<=pBirth,1,0)
#     YearlyBirths[i] <- sum(FishBirth)
#     
#     # Simulate deaths
#     Death <- t(runif(N[i-1]))
#     FishDeath <- ifelse(Death<=pDeath,1,0)
#     YearlyDeaths[i] <- sum(FishDeath)
#     
#     N[i] <- N[i-1] + YearlyBirths[i] - YearlyDeaths[i]
#     
#     
#     # Simulate captures
#     Fish <- runif(N[i])
#     FishCaught[i] <- sum(ifelse(Fish<=pCapture,1,0))
#     
#     
#   }
#   
# 
#   plot(1:t,N, 1:t, n)
#   
#   
#   
#   # Simulate captures
#   Fish <- runif(max(N)*t, ncol=t)
#   
#   for (i in 1:t) {
#     #This is to stop non-existent fish being captured
#     Fish[N[i]:max(N)] <- 1
#   }
#   
#   FishCaught <- ifelse(Fish<=pCapture,1,0)
#   fVector <- cbind(apply(FishCaught,1,sum))
#   
  
  
#   # Setup fish matrix
#   # Assume capture probabilities uniformly distributed on each sampling occasion
#   Fish <- matrix(runif(superN*t), ncol=t)
#   FishCaught <- ifelse(Fish<=pCapture,1,0)
#   fVector <- cbind(apply(FishCaught,1,sum))
#   
#   # Fish population must be equal to N_0 (close population size) at time t
#   # Work backwards for births and deaths. i.e. start at t then check t-1 to see
#   # if any births/deaths occurred. If death occurred between t-1 and t add 1 to
#   # population size at t-1. If birth occurred between t-1 and t subtract 1 from
#   # population at t-1.
# 
#   # Work out number of births for each year leading up to t
#   
#   # ERROR: Not quite right - need to check if fish is already born. Can't have a
#   # fish be born twice or more, which is what can happen here.
#   
#   Birth <- matrix(runif(superN*t), ncol=t)
#   FishBirth <- ifelse(Birth<=pBirth,1,0)
#   YearlyBirths <- rbind(apply(FishBirth,2,sum))
#   
#   
#   # Work out number of deaths for each year leading up to t
#   
#   # ERROR: Deaths need to be examined differently - can't have non existing fish die which is what is happening here.
#   # Need to check if a fish is alive before it can die.
#   
#   Death <- matrix(runif(superN*t), ncol=t)
#   FishDeath <- ifelse(Death<=pDeath,1,0)
#   YearlyDeaths <- rbind(apply(FishDeath,2,sum))
#   # We know there are N_0 fish at time=t
#   N[t] <- N_0
#   for (i in t-1:1) {
#     N[i] <- N[i+1]+YearlyDeaths[i+1]-YearlyBirths[i+1]
#   }
  
  return(N)
  
}

