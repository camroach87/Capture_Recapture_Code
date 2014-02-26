#### Preamble ####

# PROGRAM: Robust Design
# CODER: Cameron Roach
# DATE: 7/7/2013
# DESCRIPTION: Estimates population by splitting sampling timeline into 
# overlapping sections that are assumed to be closed. Calculate population in
# each of these windows using closed methods.
# 
# aCaptWindow(animal, year in window, window) is a 3d array showing capture
# history in window. Each window has a length given by the 'window'
# variable.




#### Header ####
# Make sure the correct directory is chosen
# switch(Sys.info()["sysname"],
#   Darwin =  {setwd("/Volumes/Data/Users/chroach/Dropbox/Uni/Research/")},
#   Windows = {dirBootcamp <- "C:/Users/chroach/Dropbox/Uni/Research"
#              dirWindows <- "C:/Users/Cameron/Dropbox/Uni/Research"
#              if (file.info(dirBootcamp)$isdir) {setwd(dirBootcamp)}
#              else {setwd(dirWindows)}}
# )

source('Programs/MakeCaptureMatrix.R')
library(RMark)



#### Main #####
CR_RobustDesign <- function(window, sim="N", N.0=1000, nsampOcc=14, p=0.05, phi=0.9, pBirth=0.1, pImmigration=0.1) {
  # 2*window should be the number of sampling occasions in a row (minus one)
  # where each animal is only captured at most twice starting from any time. If
  # we take x0 as the centre of the window, we take "window" number of values
  # left of x0 and "window" number of values right of x0. Bit crappy towards end
  # points, use duplicate values of end points so that we can get "window"
  # number of values.
  # 
  # sim is "Y" or "N" variable that designates if this is a simulation
  # 
  # TO DO
  # - Would be nice if capture history matrix could be input for non 
  # simulations.
  # - A function that automatically determines the optimal window (longest
  # stretch of time where fish is caught at most twice for any window start
  # point).
  
  
  
  if (sim=="N")  {
    mtrxCapt <- mkCptrMtrx()
    sampOcc <- as.numeric(colnames(mtrxCapt))
    nsampOcc <- length(sampOcc)
    #rownames(aCaptWindow) <- rownames(mtrxCapt)
  } else if (sim=="Y") {
    sampOcc <- 1:nsampOcc
    output <- mkOpenSimMtrx(N.0, nsampOcc, p, phi, pBirth, pImmigration)
    mtrxCapt <- output[[1]]
    N.t <- output[[2]]
  } else {
    return("Didn't input Y or N for sim variable. Try again.")
  }
  
  
  # Calculate Chao's sparse data estimator
  aCaptWindow <- array(0,c(nrow(mtrxCapt),2*window+1,nsampOcc))
  aChao <- array(0,c(nsampOcc,2))
  
  
  
  
  
  for (i in 1:nsampOcc)) {
    if (i<=(window)) {
      # start points: create capture window matrix with first i columns a copy
      # of the (i+1)th column (which corresponds to the first column in the
      # capture history matrix).
      
      aCaptWindow[,1:(window-i+1),i] <- mtrxCapt[,1]
      aCaptWindow[,(window-i+2):(1+2*window),i] <- mtrxCapt[,1:(i+window)]
      aChao[i,] <- unlist(ChaoEst(aCaptWindow[,,i]))
      
    } else if(i>=(nsampOcc-window) {
      # end points: create capture window matrix with last i columns a copy of
      # the (nsampOcc-(i+1))th column (which corresponds to the first column in
      # the capture history matrix).
      
      
      
      
      
    } else {
      # center points: take window number of points left and right of ith
      # sampling occasion
      aCaptWindow[,1:(1+2*window),i] <- mtrxCapt[,(i-window):(i+window)]
      aChao[i,] <- unlist(ChaoEst(aCaptWindow[,,i]))
    }
  }
  
  # Kernel smoothing
  sChao <- ksmooth(sampOcc[1:nsampOcc],aChao[,1])
  
  
  # Get CJS estimates from RMark
  p_rmark <- TestRMark(mtrxCapt)
  n.t <- colSums(mtrxCapt)
  N_rmark <- n.t/p_rmark
  
  #Plot estimates
  plot(N_rmark ~ sampOcc, type="l", col="red", ylim=c(0,max(sChao)))
  lines(sampOcc[1:nsampOcc], sChao,col=colors()[25+6*window])
  #And actual if a simulation
  if (sim=="Y") {lines(sampOcc,N.t,col="green")}
  
  
  rm(output)
  if (sim=="Y") {
    output <- list(sChao, N.t)
  } else {
    output <- list(sChao)
  }
  
  return(output)
  
}





#### Functions ####
ChaoEst <- function(mtrxCapt) {
  
  fVector <- cbind(apply(mtrxCapt,1,sum))
  
  # Chao estimator for time variation model M_t
  f1 <- sum(ifelse(fVector==1,1,0))
  #f2 <- sum(ifelse(fVector==2,1,0))
  #Count f3, f4, etc in f2 just in case they are not equal to zero.
  f2 <- sum(ifelse(fVector>=2,1,0))
  S <- sum(ifelse(fVector>=1,1,0))
  
  
  Z <- mtrxCapt[fVector==1,]
  Z <- colSums(Z)
  
  
  N_Chao = S + (f1^2-sum(Z^2))/(2*f2)
  
  output <- list(N_Chao, S)
  
  return(output)
  
}




TestRMark <- function(mtrxCapt) {
  #RMark doesn't like values different to one or zero. Also doesn't like rows with only zeros.
  mtrxCapt[mtrxCapt>1]<-1
  mtrxCapt<-mtrxCapt[rowSums(mtrxCapt)!=0,]
  mtrxCapt_f <- apply(format(mtrxCapt), 1, paste, collapse="")
  mtrxCapt_f <- paste(mtrxCapt_f, " 1;")
  write.table(as.matrix(mtrxCapt_f),"Programs/Data/MarkTest1.inp",col.names = F, row.names = F, sep="", quote=FALSE)
  mtrxCapt_f <- convert.inp("Programs/Data/MarkTest1.inp")
  output <- mark(mtrxCapt_f)
  p_est <- output$results$real$estimate[2]  
  
  return(p_est)
}
