#### Preamble ####

# PROGRAM: Get estimators
# CODER: Cameron Roach
# DATE: 24/2/2014
# DESCRIPTION: Various functions to produce capture-recapture estimators






######################### Estimators #####################

calcCR <- function(mtrxCapt, window, xType="Occasion", timeBandwidth=180, dates.occ) {
  # DESCRIPTION: Estimates population by splitting sampling timeline into 
  # overlapping periods that are assumed to be closed. Calculate population in
  # each of these windows using closed methods.
  #
  #   
  # 2*window+1 is the period we assume the population is closed. 2*window+1
  # should be the number of sampling occasions in a row where each animal is
  # only captured at most twice starting from any time. If we take x0 as the
  # centre of the window, we take "window" number of values left of x0 and
  # "window" number of values right of x0. Bit crappy towards end points, use
  # duplicate values of end points so that we can get "window" number of 
  # values.
  #
  # Use direct plug-in method, dpill(), to calculate bandwidth for kernel 
  # regression when dealing with occasions. For times, bandwidth of 180 has been
  # selected as default as it appears to give the best tradeoff between bias and
  # variance for the trout cod data. Obviously, if different data is used this
  # value may need to be updated.
  #
  # xType specifies if we are basing the estimates on capture occasion or
  # capture time. Takes "Occasion" and "Time" as inputs.
  
  
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
  if (xType == "Occasion") {
    h <- dpill(1:T,aChao_endfix[,1])
    sChao <- locpoly(1:T,aChao_endfix[,1], degree=1, bandwidth=h, gridsize=T)
    
    output <- sChao$y
  } else if (xType == "Time") {
    h <- timeBandwidth # works well for the TC data
    sChao <- locpoly(dates.occ$Day,aChao_endfix[,1], degree=1, bandwidth=h, gridsize=T)
    
    # converts the x values from locpoly to dates
    sChao$x <- min(dates.occ$Date) + sChao$x*60*60*24
    
    output <- sChao
  }
  
  
  return(output)
}


CR.bs <- function(mtrxCapt, window, indices) {
  mtrxCapt <- mtrxCapt[indices,]
  N.bs <- calcCR(mtrxCapt, window)
  return(N.bs)
}

CR.bs.time <- function(mtrxCapt, window, dates.occ, indices) {
  mtrxCapt <- mtrxCapt[indices,]
  N.bs <- calcCR(mtrxCapt, window, xType="Time", dates.occ=dates.occ)
  return(N.bs$y)
}

calcChaoMt <- function(mtrxCapt, bc=TRUE) {
  # DESCRIPTION: Chao's sparse data estimator for closed populations
  
  fVector <- cbind(apply(mtrxCapt,1,sum))
  
  # Chao estimator for time variation model M_t
  f1 <- sum(ifelse(fVector==1,1,0))
  #Count f3, f4, etc in f2 just in case they are not equal to zero.
  f2 <- sum(ifelse(fVector>=2,1,0))
  S <- sum(ifelse(fVector>=1,1,0))
  
  
  Z <- mtrxCapt[fVector==1,]
  Z <- colSums(Z)
  
  if (bc==TRUE) {
    N_Chao = S + (f1^2-sum(Z^2))/(2*(f2+1))
  } else {
    N_Chao = S + (f1^2-sum(Z^2))/(2*f2)
  }
  
  
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
      # captured) which causes the resulting matrix to be around the wrong way. 
      # Need to transpose to fix. UPDATE: could have just used drop=FALSE when 
      # subsetting....
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
        
    # bias corrected
    M <- m + (R+1)*Z/(r+1)
    N <- (n+1)*M/(m+1)
    
    N.ci <- calcManlyCi(N,n,M,m,R,r)

    return(list("N"=N, "N.ci"=N.ci))
}


calcManlyCi <- function(N,n,M,m,R,r) {
  # Calculates 95% confidence intervals for JS using Manly's method
  
  p <- n/N
  
  T_1.N <- log(N) + log(0.5*(1-0.5*p+(1-p)^2))
  T_1.N.var <- (M-m+R+1)/(M+1)*(1/(r+1)-1/(R+1)) + 1/(m+1) - 1/(n+1)
  
  T_1l <- T_1.N - 1.6*sqrt(T_1.N.var)
  T_1u <- T_1.N + 2.4*sqrt(T_1.N.var)
  
  N.ci.l <- (4*exp(T_1l) + n)^2/(16*exp(T_1l))
  N.ci.u <- (4*exp(T_1u) + n)^2/(16*exp(T_1u))
  
  output <- data.frame("ci.l" = N.ci.l,
                       "ci.u" = N.ci.u)
  
  return(output)  
}