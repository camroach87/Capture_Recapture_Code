#### Preamble ####

# PROGRAM: Closed and open population simulation
# CODER: Cameron Roach
# DATE: 2/4/2013
# DESCRIPTION: Simulates what is actually happening in closed and open 
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
# 3. Assumption: In open pop. model have assumed all captured fish released,
# i.e., no fish die during capture.
# 4. Calcr assumes all animals released after capture, i.e., R[i]=n[i]. If this
# is not the case, then Calcr will need to be modified. Will need to input a
# release matrix? Maybe call Calcr() from inside a new function CalcR() that
# calculates both r and R?
# 
# 
# Fix: actual Mt being used instead of estimated version - use Jolly Seber
# approximation!






#### Console commands ####
# These are not part of the program but are useful

#for (i in 1:100) {bla<-bla+ClosedSim(N.0, t, pCapture)}; bla/100
#for (i in 1:100) {bla<-bla+OpenSim(N.0, t, pCapture)}; bla/100
Initialise <- function() {
  # put these in a function so they didn't run every time I saved/sourced
  N.0 <-20
  t <- 10
  pBirth <- 0.03   
  pDeath <- 0.02
  pCapture <- 0.3
  #list[Nt_open, Mt] <- OpenSim(N.0, t, pCapture, pBirth, pDeath)
  
  iterations <- 500
  finalpop <- rep(NA,iterations)
  openEstimate <- rep(NA,iterations)
  ChaoEstimate <- rep(NA,iterations)
  
  for (i in 1:iterations) {

    # Open test
    openOut <- OpenSim(N.0, t, pCapture, pBirth, pDeath)
    finalpop[i] <- openOut[1]
    ChaoEstimate[i] <- openOut[2]
    openEstimate[i] <- openOut[3] 
    
    
#     # Closed test
#     closedOut <- ClosedSim(N.0, t, pCapture)
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

ClosedSim <- function(N.0, t, pCapture) {
  
  # Setup fish matrix
  # Assume capture probabilities unifiormly distributed on each sampling occasion
  Fish <- matrix(runif(N.0*t), ncol=t)
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




OpenSim <- function(N.0, t, pCapture, pBirth, pDeath) {
  # Simulates open population. N.0 is initial population size.
  # Assumes each animal only capable of one offspring per sampling interim.
  # Assumes death equally likely at any point in time of fish's lifespan
  
  # Initialise variables
  N <- rep(NA, t)
  n <- rep(NA, t)
  m <- rep(NA, t)
  M <- rep(NA, t)
  nbirths <- rep(0, t)
  ndeaths <- rep(0, t)
  Nsuper <- 10*N.0
  FishExists <- matrix(rep(0,Nsuper*t), ncol=t)
  FishGivesBirth <- matrix(runif(Nsuper*t), ncol=t)
  FishDies <- matrix(runif(Nsuper*t), ncol=t)
  
  # Simulate births and deaths starting from initial population N.0
  N[1] <- N.0
  NewFish <- N.0
  FishExists[1:N.0,1] <- 1
  for (i in 2:t) {
    
    # Use 1-pBirth as probability since non-existing fish have 0 used to identify them
    # i.e. there will be lots of zeros when taking inner product of FishExists and FishGivesBirth
    births <- ifelse(FishExists[1:NewFish,i-1]*FishGivesBirth[1:NewFish,i-1]>=(1-pBirth),1,0)
    nbirths[i] <- sum(births)
    FishExists[(NewFish+1):(NewFish+nbirths[i]),i] <- 1
    
    # ERROR: Should be 1-pDeath here
    #deaths <- ifelse(FishExists[1:NewFish,i-1]*FishDies[1:NewFish,i-1]>=(1-pBirth),1,0)
    deaths <- ifelse(FishExists[1:NewFish,i-1]*FishDies[1:NewFish,i-1]>=(1-pDeath),1,0)
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
  m <- Calcm(FishCaught, t, N)
  R <- n  # Note 3: Assumption: All fish captured are released
  Z <- CalcZ(FishCaught, t, N)
  r <- Calcr(FishCaught, t, N)
  
  # remember, values for 2:(t-1) are only valid values
  # Using Seber unbiased estimator for M (pg 41 Handbook)
  # old M <- m + R*Z/r
  M <- m + (R+1)*Z/(r+1)
    
  M_actual <- sum(ifelse(FishExists[,t]*fVector>=1,1,0))
  
  # This is the estimate that would be produced fitting Chao's estimator for M_actual
  # model to an open population.
  f1 <- sum(ifelse(fVector==1,1,0))
  f2 <- sum(ifelse(fVector==2,1,0))
  S <- sum(ifelse(fVector>=1,1,0))
  Z_Chao <- FishCaught[fVector==1,]
  Z_Chao <- colSums(Z_Chao)
  N_Chao = S + (f1^2-sum(Z_Chao^2))/(2*f2)
  
  
  # My modified estimate for population size would be.
  # FIX: shouldn't be using M_actual (as we don't really know this in practice) -
  # should use Jolly - Seber estimate: Mt=mt+(Rt*Zt)/rt
  #N_Cam = M_actual/S*N_Chao
  N_Cam = M[t-1]/S*N_Chao
  
  
  return(c(N[t], N_Chao, N_Cam))
  #return(c(N[t],M_actual))
  
}





Calcm <- function(CaptureMatrix, t, N) {
  # calculates number of marked fish captured, m_j, at time j
  
  m <- rep(0,t)  
  
  
  # start at i=2 because m[1]=0
  for (i in 2:t) {
    # Checks each fish to see if previously captured
    for (j in 1:N[i]) {
      # Checks if fish j caught at time i and if it was previously marked at times 1:(i-1)
      if (CaptureMatrix[j, i] == 1 && sum(CaptureMatrix[j,1:(i-1)])>=1) {
        m[i] = m[i] + 1
      }
    }
  }
  return(m) 
}


CalcZ <- function(CaptureMatrix, t, N) {
  # Calculates z_j, the number of members of the marked population not captured
  # at j (Mj-mj) captured again later.
  
  z <- rep(0,t)
  
  # Start at i=2 and at t-1 because z[1]=z[t] = 0
  for (i in 2:(t-1)) {
    for (j in 1:N[i]) {
      # Checks each fish to see if fish j not captured at time i but captured
      # later on.
      if (CaptureMatrix[j, i] == 0
          && sum(CaptureMatrix[j,1:(i-1)])>=1
          && sum(CaptureMatrix[j,(i+1):t])>=1) {
        z[i] = z[i] + 1
      }
    }
  }
  return(z)
}


Calcr <- function(CaptureMatrix, t, N) {
  # Calculates r_j, the number of members of R_j captured again later
  # Note 4: if R_j is not simply all animals captured, then this function will need to change.
  # This can be fixed by inputting a 'ReleaseMatrix' rather than a capture
  # matrix. Or, if we know that only p percent of animals survive, can probably
  # just scale the current r vector calculated.
  
  r <- rep(0,t)  
  
  # end at t-1 because r[t] = 0
  for (i in 1:(t-1)) {
    for (j in 1:N[i]) {
      # Checks each fish to see if fish j released at time i and captured
      # again later on.
      if (CaptureMatrix[j, i] == 1 && sum(CaptureMatrix[j,(i+1):t])>=1) {
        r[i] = r[i] + 1
      }
    }
  }
  return(r)
}