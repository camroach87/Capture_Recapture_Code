#### Preamble ####

# PROGRAM: Get estimators
# CODER: Cameron Roach
# DATE: 24/2/2014
# DESCRIPTION: Various functions to produce capture-recapture estimators






######################### Estimators #####################

CR_RobustDesign <- function(mtrxCapt, window) {
  # DESCRIPTION: Estimates population by splitting sampling timeline into 
  # overlapping sections that are assumed to be closed. Calculate population in
  # each of these windows using closed methods.
  #
  #
  # 2*window+1 should be the number of sampling occasions in a row where each
  # animal is only captured at most twice starting from any time. If we take x0
  # as the centre of the window, we take "window" number of values left of x0
  # and "window" number of values right of x0. Bit crappy towards end points,
  # use duplicate values of end points so that we can get "window" number of
  # values.
  
  
  T <- ncol(mtrxCapt)
  
  if(T==2*window) {
    stop("Window is too big! Must be less than half the number of sampling occasions.")
  }
  
  # Calculate Chao's sparse data estimator
  #aCaptWindow <- array(0,c(nrow(mtrxCapt),2*window+1,T-2*window))
  aChao <- array(0,c(T-2*window,2))
  
  for (i in 1:(T-2*window)) {
    #aCaptWindow[,1:(1+2*window),i] <- mtrxCapt[,i:(i+2*window)]
    #aChao[i,] <- unlist(calcChaoMt (aCaptWindow[,,i]))
    aChao[i,] <- unlist(calcChaoMt (mtrxCapt[,i:(i+2*window)]))
  }
  
  # Copy values for end points
  aChao_endfix <- array(0,c(T,2))
  aChao_endfix[1:window,] <- t(matrix(aChao[1,],2,window))
  aChao_endfix[(window+1):(T-window),] <- aChao[,]
  aChao_endfix[(T-window+1):T,] <- t(matrix(aChao[dim(aChao)[1],],2,window))
  # Kernel smoothing
  sChao <- ksmooth(1:T,aChao_endfix[,1],bandwidth=window)
  
  
  # Calculate variance
  
  
  
  
  return(sChao)
}



calcChaoMt <- function(mtrxCapt) {
  # DESCRIPTION: Chao's sparse data estimator for closed populations
  
  fVector <- cbind(apply(mtrxCapt,1,sum))
  
  # Chao estimator for time variation model M_t
  f1 <- sum(ifelse(fVector==1,1,0))
  #Count f3, f4, etc in f2 just in case they are not equal to zero.
  f2 <- sum(ifelse(fVector>=2,1,0))
  S <- sum(ifelse(fVector>=1,1,0))
  
  
  Z <- mtrxCapt[fVector==1,]
  Z <- colSums(Z)
  
  
  N_Chao = S + (f1^2-sum(Z^2))/(2*f2)
  
  output <- list(N_Chao, S)
  
  return(output)
  
}



calcJS <- function(mtrxCapt) {
    #Jolly Seber
    n <- colSums(mtrxCapt)
    
    m <- NA
    R <- NA
    Z <- NA
    r <- NA
    
    for (i in 2:ncol(mtrxCapt)){
      m[i] <- calcMarked(mtrxCapt[,1:i])
      R <- m # assuming all released
      #       Z <- calcZ(mtrxCapt)
      #       r <- calcRecapt(mtrxCapt)
    }
    

    M <- m + R*Z/r
    N <- n*M[1:length(n)]/m
    
    plot(N ~ days, type="l", col="green")
}