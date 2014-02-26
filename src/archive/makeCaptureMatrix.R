############################### PREAMBLE ################################
# PROGRAM: Make capture matrices
# CODER: Cameron Roach
# DATE: 9/7/2013
# DESCRIPTION: Constructs various capture matrices using either a datafile of
# capture occasions or by simulation.






#### Functions ####
mkCptrMtrx_years <- function(datafile = "TC_Data_Charles.csv") {
  # Constructs a capture matrix using sourcefile in format of TC_Data_Charles.csv.
  data <- read.csv(file.path(dataDir,datafile));
  
  D <- as.Date(data$surveydate, format="%d/%m/%Y")
  tabD <- table(D)
  
  # add year to data
  data[10:11, "year"] <- NA
  data$year <- as.numeric(format(D, "%Y"))
  
  
  # useful
  mtrxCapt <- table(data$idfish, data$year)
  
  return(mtrxCapt)
}

mkCptrMtrx <- function(datafile = "TC_Data_Charles.csv") {
  # Constructs a capture matrix using sourcefile in format of TC_Data_Charles.csv.
  data <- read.csv(file.path(dataDir,datafile));
  
  D <- as.Date(data$surveydate, format="%d/%m/%Y")
  tabD <- table(D)
  
  # add day to data
  data[10:11, "day"] <- NA
  data$day <- as.numeric(D)
  
  
  # useful
  mtrxCapt <- table(data$idfish, data$day)
  
  return(mtrxCapt)
}





mkCloseSimMtrx <- function(N.0, t, pCapture) {
  # Simulates closed population. N.0 is initial population size.
  # Assumes capture probabilities are uniformly distributed.
  mtrxCaptProb <- matrix(runif(N.0*t),N.0,t)
  mtrxCapt <- ifelse(mtrxCaptProb<=pCapture,1,0)
  
  return(mtrxCapt)  
}
  



mkOpenSimMtrx <- function(N.0, t, p) {
  beta <- 0.01;
  
  R <- 0;

  #Fix +1 and R
  Pop <- Pop.Mat(N.0,t+1,R,beta);
  
  mtrxP <- Pop*matrix(runif(N.0*(t-1)),N.0,t)
  
  mtrxCapt <- ifelse(mtrxP>=(1-p),1,0) 
  
  return(list(mtrxCapt, colSums(Pop)))
  
}



mkOpenSimMtrx1 <- function(N.0, t, p, phi, pBirth, pImmigration) {
  # Simulates open population. N.0 is initial population size.
  # phi: Probability of survival
  # p: Capture probability
  # Population calculed from N_t = X_t+Z_t+Y_t 
  # X_t =Number of individuals surviving from t-1 to t
  # Z_t = Number of offspring of individuals in population at t-1
  # Y_t = Number of immigrants in t-1,t
  
  
  # Generate population at each sampling occasion
  N <- rep(NA, t)
  X <- rep(NA, t)
  Z <- rep(NA, t)
  Y <- rep(NA, t)
  N[1] <- N.0
  X[1] <- 0
  Z[1] <- 0
  Y[1] <- 0
  
  for (i in 2:t) {
    X[i] <- rbinom(1,N[i-1],phi)
    Z[i] <- rpois(1, N[i-1]*pBirth)
    Y[i] <- rpois(1, pImmigration)
    N[i] <- X[i]+Z[i]+Y[i]
  }
  
  #cat(X, "\n")
  
  # Create capture probability matrix for each animal
  NupperB <- max(N)+sum(Z)+sum(Y)
  mtrxP <- matrix(runif(NupperB*t),NupperB,t)
  deaths <- 0
  for (i in 1:t) {
    # zero out deaths from beginning of matrix (oldest animals die first)
    if(i>=2) {
      deaths <- deaths + N[i-1] - X[i]
      #cat(deaths, "\n")
      if (deaths>0) {mtrxP[1:deaths,i] <- 0}
    }
    
    # Only want animals that have exist to have capture probability
    if ((N[i]+1+deaths)<=NupperB) {
      for(j in (N[i]+1+deaths):NupperB){
        mtrxP[j,i]<-0
      }
    }
  }
  
  mtrxCapt <- ifelse(mtrxP>=(1-p),1,0) 
  
  return(list(mtrxCapt, N))
}


#### old ####
# mkOpenSimMtrx_old1 <- function(N.0, t, pCapture, pBirth, pDeath) {
#   # Simulates open population. N.0 is initial population size.
#   # Assumes each animal only capable of one offspring per sampling interim.
#   # Assumes death equally likely at any point in time of fish's lifespan
#   
#   # NOTE: Fish death and birth probability matrices are not uniformly 
#   # distributed in real life. e.g. deaths are more likely to be normally 
#   # distributed around the mean lifespan of animal, births could be uniformly 
#   # distributed once animal reaches sexual maturity. Something to consider 
#   # fixing up down the line....
#   
#   # Initialise variables
#   N <- rep(NA, t)
#   n <- rep(NA, t)
#   m <- rep(NA, t)
#   M <- rep(NA, t)
#   nbirths <- rep(0, t)
#   ndeaths <- rep(0, t)
#   Nsuper <- 10*N.0
#   FishExists <- matrix(rep(0,Nsuper*t), ncol=t)
#   FishGivesBirth <- matrix(runif(Nsuper*t), ncol=t)
#   FishDies <- matrix(runif(Nsuper*t), ncol=t)
#   
#   # Simulate births and deaths starting from initial population N.0
#   N[1] <- N.0
#   NewFishMarker <- N.0
#   FishExists[1:N.0,1] <- 1
#   for (i in 2:t) {
#     
#     # Use 1-pBirth as probability since non-existing fish have 0 used to identify them
#     # i.e. there will be lots of zeros when taking inner product of FishExists and FishGivesBirth
#     
#     #deaths happen first so that a new fish does not immediately die
#     deaths <- ifelse(FishExists[1:NewFishMarker,i-1]*FishDies[1:NewFishMarker,i-1]>=(1-pDeath),1,0)
#     ndeaths[i] <- sum(deaths)
#     FishExists[1:NewFishMarker,i] <- FishExists[1:NewFishMarker,i-1]*(1-deaths)
#     
#     births <- ifelse(FishExists[1:NewFishMarker,i-1]*FishGivesBirth[1:NewFishMarker,i-1]>=(1-pBirth),1,0)
#     nbirths[i] <- sum(births)
#     FishExists[(NewFishMarker+1):(NewFishMarker+nbirths[i]),i] <- 1
#     
#     N[i] <- N[i-1] + nbirths[i] - ndeaths[i]
#     
#     NewFishMarker <- NewFishMarker + nbirths[i]
#   }
#   
#   #plot(1:t,N)
#   
#   
#   # Simulate captures
#   FishProbs <- FishExists*matrix(runif(Nsuper*t), ncol=t)
#   # Use 1-pCapture so that fish that don't exist aren't counted.
#   mtrxCapt <- ifelse(FishProbs>=(1-pCapture),1,0)
# 
#   return(list(mtrxCapt, FishExists))
# }
# mkOpenSimMtrx_old2 <- function(N.0, t, p, phi, pBirth, pImmigration) {
#   # Simulates open population. N.0 is initial population size.
#   # phi: Probability of survival
#   # p: Capture probability
#   # Population calculed from N_t = X_t+Z_t+Y_t 
#   # X_t =Number of individuals surviving from t-1 to t
#   # Z_t = Number of offspring of individuals in population at t-1
#   # Y_t = Number of immigrants in t-1,t
#   
#   
#   # Note: Still not quite right - first few animals will live forever. When 
#   # calculating the population at next sampling stage need to decide which
#   # animals have died (set probability of capture to zero for dead animal in
#   # mtrxP)
#   
#   # Generate population at each sampling occasion
#   N <- rep(NA, t)
#   N[1] <- N.0
#   for (i in 2:t) {
#     X <- rbinom(1,N[i-1],phi)
#     Z <- rpois(1, N[i-1]*pBirth)
#     Y <- rpois(1, pImmigration)
#     N[i] <- X+Z+Y
#   }
#   
#   # Create capture probability matrix for each fish
#   maxN <- max(N)
#   mtrxP <- matrix(runif(maxN*t),maxN,t)
#   for (i in 1:t) {
#     if ((N[i]+1)<=maxN) {
#       for(j in (N[i]+1):maxN){
#         mtrxP[j,i]<-0
#       }
#     }
#   }
#   
#   mtrxCapt <- ifelse(mtrxP>=(1-p),1,0) 
#   
#   return(list(mtrxCapt, N))
# }