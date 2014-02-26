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
# number of new fish goes above 9*N (super population - initial population)
# limit.
# 
# Fix: actual Mt being used instead of estimated version - use Jolly Seber
# approximation!






#### Console commands ####
# These are not part of the program but are useful

#for (i in 1:100) {bla<-bla+ClosedSim(N_0, t, pCapture)}; bla/100
#for (i in 1:100) {bla<-bla+OpenSim(N_0, t, pCapture)}; bla/100
Initialise <- function() {
  # put these in a function so they didn't run every time I saved/sourced
  N_0 <-1000
  t <- 10
  pBirth <- 0.03   
  pDeath <- 0.02
  pCapture <- 0.05
  #list[Nt_open, Mt] <- OpenSim(N_0, t, pCapture, pBirth, pDeath)
  
  iterations <- 500
  finalpop <- rep(NA,iterations)
  openEstimate <- rep(NA,iterations)
  ChaoEstimate <- rep(NA,iterations)
  
  for (i in 1:iterations) {

    # Open test
    openOut <- OpenSim(N_0, t, pCapture, pBirth, pDeath)
    finalpop[i] <- openOut[1]
    ChaoEstimate[i] <- openOut[2]
    openEstimate[i] <- openOut[3] 
    
    
#     # Closed test
#     closedOut <- ClosedSim(N_0, t, pCapture)
#     ChaoEstimate[i] <- closedOut[1]
    
    
    
  }
  
  #Variance test for open
  mean(finalpop)
  mean(ChaoEstimate)
  mean(openEstimate)
  var(finalpop-ChaoEstimate)
  var(finalpop-openEstimate)
  
#   # Variance test for closed
#   mean(ChaoEstimate)
#   var(N-ChaoEstimate)
  
  
}






#### Code #####

ClosedSim <- function(N_0, t, pCapture) {
  
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
  
  
  N_Chao = S + (f1^2-sum(Z^2))/(2*f2)
  
  return(c(N_Chao, S))
  
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
  
  # Simulate births and deaths starting from initial population N_0
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

  #plot(1:t,N)
  
  
  # Simulate captures
  FishProbs <- FishExists*matrix(runif(Nsuper*t), ncol=t)
  FishCaught <- ifelse(FishProbs>=(1-pCapture),1,0)
  fVector <- cbind(apply(FishCaught,1,sum))
  n <- colSums(FishCaught)
  # These can be worked out later if need be
  #   m <- 
  #   M <-
  
  Mt <- sum(ifelse(FishExists[,t]*fVector>=1,1,0))
  
  # This is the estimate that would be produced fitting Chao's estimator for Mt
  # model to an open population.
  f1 <- sum(ifelse(fVector==1,1,0))
  f2 <- sum(ifelse(fVector==2,1,0))
  S <- sum(ifelse(fVector>=1,1,0))
  Z <- FishCaught[fVector==1,]
  Z <- colSums(Z)
  N_Chao = S + (f1^2-sum(Z^2))/(2*f2)
  
  
  # My modified estimate for population size would be.
  # FIX: shouldn't be using Mt (as we don't really know this in practice) -
  # should use Jolly - Seber estimate: Mt=mt+(Rt*Zt)/rt
  N_Cam = Mt/S*N_Chao
  
  
  
  return(c(N[t], N_Chao, N_Cam))
  #return(c(N[t],Mt))
  
}