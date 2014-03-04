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
    n <- calcTotal(mtrxCapt)
    cat("n complete \n")
    m <- calcMarked(mtrxCapt)
    cat("m complete \n")
    R <- calcReleased(mtrxCapt)
    cat("R complete \n")
    Z <- calcZ(mtrxCapt)
    cat("Z complete \n")
    r <- calcRecapt(mtrxCapt)
    cat("r complete \n")
    
    M <- m + R*Z/r
    N <- n*M[1:length(n)]/m
    
    plot(N ~ days, type="l", col="green")
}





##################### Assorted functions #################


# Calculates number of animals captured exactly j times, j=1,...,k
calcExact <- function(data) {
  uniqueFish <- unique(data$idfish)
  counts <- c()
  f <- c()
  days <- sort(unique(data$day))
  k <- max(days)-min(days) + 1
  
  #First calculate number of times each fish caught
  
  for (i in 1:length(uniqueFish)) counts[i] <- sum(uniqueFish[i] == data$idfish)
  
  #Then calculate number of animals caught exactly j times
  for (j in 1:k) f[j] <- sum(counts == j)
  
  return(f)
  
}



calcTotal <- function(data) {
  
  captFish <- NA
  n <- NA
  M <- NA
  
  
  
  for (j in 1:k) {
    n[j] <- length(data$idfish[data$day == j + minDay - 1])
  }
  
  return(n)
}




# Calculates number of distinct animals caught prior to jth capture occasion.
# By definition M_0 = 0 and M_k+1 = total distinct animals captured
# This assumes a closed population
calcDistinct <- function(data) {
  
  u <- calcUnmarked(data)
  
  cat(u, "\n")
  
  M <- c()
  M[1] <- 0
  
  for (j in 1:length(u)) M[j+1] <- sum(u[1:j])
  
  return(M)
  
}






# calculates number of unmarked fish, u_j
calcUnmarked <- function(data) {
  u <- rep(0,k)
  markedFish <- c()
  
  
  for (j in 1:k) {
    
    # calculate number of unmarked animals for capture occasion j
    if (j != 1) {
      
      for (i in 1:length(fishDay[[j]])) {
        
        # checks if fish i in day j was captured in previous days (or already caught this day)
        
        if (fishDay[[j]][i] %in% markedFish) {
          fishNotCaught = 0
          #cat(fishDay[[j]][i],"\n")
          
        } else {
          fishNotCaught = 1
        }
        
        markedFish <- c(markedFish, fishDay[[j]][i])
        
        if (fishNotCaught == 1) u[j] = u[j] + 1
        
      }
      
    } else {
      u[j] = length(fishDay[[j]])
      markedFish <- c(fishDay[[j]])
    }
    
    
  }
  
  return(u)
  
}





# calculates number of marked fish, m_j
calcMarked <- function(data) {
  m <- rep(0,k)
  markedFish <- c(fishDay[[1]])
  
  
  
  # start at j=2 because m_1=0
  for (j in 2:k) {
    # DEBUG cat("\n j = ", j)
    
    # checks if any fish were caught on day j
    if (length(fishDay[[j]]) != 0) {
      # calculate number of unmarked animals for capture occasion j
      for (i in 1:length(fishDay[[j]])) {
        # DEBUG cat(" i = ", i)
        # checks if fish i in day j was captured in previous days (or already captured this day)
        if (fishDay[[j]][i] %in% markedFish) {
          m[j] = m[j] + 1
        }
        
        markedFish <- c(markedFish, fishDay[[j]][i])    
      }
    }  
  }
  
  return(m)
  
}







# calculates R_j, total number of animals capture at sampling occasion j that are released
# Note: if this changes see note for calcRecapt() as it will probably need to change as well (see definition in book)
calcReleased <- function(data) {
  # Until I hear otherwise, I will assume that all animals that are captured are relased
  R <- calcTotal(data)
  
  return(R)  
}




# Calculates r_j, the number of members of R_j captured again later
# Note: doesn't include fish from sampling occasion j captured again during sampling occasion j
# Note: if R_j is not simply all animals captured, then this function will need to change.
calcRecapt <- function(data) {
  r <- c()
  
  
  # Check day j for later captures in days i
  for (j in 1:(k-1)) {
    
    # refreshes recaptured matrix each day, j
    recaptured <- matrix(0)
    
    for (i in (j+1):k) {
      recaptured <- recaptured + unique(fishDay[[j]]) %in% fishDay[[i]]
    }
    
    # don't care if fish recaptured twice or more - only need to know they were captured
    recaptured[recaptured>=2] <- 1
    
    r[j] <- sum(recaptured)
    
  }
  
  r[k] <- 0
  
  return(r)
  
}







# Calculates z_j, the number of members of the marked population not captured at j (Mj-mj) captured again later
calcZ <- function(data) {
  z <- c()
  
  # Check day j for later captures in days i
  for (j in 1:(k-1)) {
    
    # refreshes markedPop and recaptured matrix each day, j
    markedPop <- c(NA)
    recaptured <- matrix(0)
    
    cat("\n j = ", j)
    
    # gets marked population excluding those fish captured in current day j
    if (j == 1) {
      recaptured <- 0
    } else {
      for (i in 1:(j-1)) {
        markedPop <- c(markedPop, fishDay[[i]][!(fishDay[[i]] %in% fishDay[[j]])])
        markedPop <- unique(markedPop)
        
        # DEBUG cat(" i=", i)
      }
      
      for (i in (j+1):k) {
        recaptured <- recaptured + markedPop %in% fishDay[[i]]
        # DEBUG cat(" i=", i)
      }
      
      # don't care if fish recaptured twice or more - only need to know they were captured
      recaptured[recaptured>=2] <- 1
    }
    
    z[j] <- sum(recaptured)
    
  }
  
  z[k] <- 0
  
  return(z)
  
}

# Calculates z_j, the number of members of the marked population not captured at j (Mj-mj) captured again later
calcZ_fast <- function(data) {
  
  t1 <- Sys.time()
  
  z <- rep(NA, k)     # CHANGE c()
  
  # Check day j for later captures in days i
  for (j in 1:300) {        #(k-1)) {
    
    # refreshes markedPop and recaptured matrix each day, j
    markedPop <- rep(NA,5000)    # CHANGE can't remember - check above function
    recaptured <- 0   # CHANGE matrix(0)
    
    # DEBUG cat("\n j = ", j)
    
    # gets marked population excluding those fish captured in current day j
    if (j == 1) {
      recaptured <- 0
    } else {
      for (i in 1:(j-1)) {
        markedPop <- c(markedPop, fishDay[[i]][!(fishDay[[i]] %in% fishDay[[j]])])
        markedPop <- unique(markedPop)
        
        # DEBUG cat(" i=", i)
      }
      
      for (i in (j+1):k) {
        recaptured <- recaptured + markedPop %in% fishDay[[i]]
        # DEBUG cat(" i=", i)
      }
      
      # don't care if fish recaptured twice or more - only need to know they were captured
      recaptured[recaptured>=2] <- 1
    }
    
    z[j] <- sum(recaptured)
    
  }
  
  z[k] <- 0
  
  (T1 <- Sys.time() - t1)
  
  return(z)
  
}