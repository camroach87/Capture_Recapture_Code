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
  sChao <- ksmooth(1:T,aChao_endfix[,1],kernel="normal",
                   bandwidth=window,x.points=1:T)
  sChao <- sChao[[2]]
    
  return(sChao)
}


CR.bs <- function(mtrxCapt, window, indices) {
  mtrxCapt <- mtrxCapt[indices,]
  N.bs <- CR_RobustDesign(mtrxCapt, window)
  return(N.bs)
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
    k <- ncol(mtrxCapt)
    
    m <- rep(NA,k)
    R <- rep(NA,k)
    Z <- rep(NA,k)
    r <- rep(NA,k)
    
    for (i in 1:k){
      #cat(i, "\n")
      # gets marked animals that are not captured at time i
      curCapt <- as.matrix(mtrxCapt[mtrxCapt[,i]>=1,])
      curNotCapt <- as.matrix(mtrxCapt[mtrxCapt[,i]==0,])
      
      # R converts curCapt to an array when only one animal is captured (or not
      # captured) which causes trhe resulting matrix to be around the wrong way.
      # Need to transpose to fix.
      if (sum(mtrxCapt[,i]>=1)==1) {curCapt <- t(curCapt)}
      if (sum(mtrxCapt[,i]==0)==1) {curNotCapt <- t(curNotCapt)}
      
      # gets statistics
      if (i>1) {
        m[i] <- checkCapt(curCapt[,1:(i-1)])
        R[i] <- n[i] # assuming all released
        Z[i] <- calcZ(curNotCapt, i)  
      }
      
      if (i<k) {
        r[i] <- checkCapt(curCapt[,(i+1):k]) #assumes R=m        
      }
    }
    
#     M <- m + R*Z/r
#     N <- n*M/m
        
    # bias adjusted
    M <- m + (R+1)*Z/(r+1)
    N <- (n+1)*M/(m+1)
    
    return(N)
}