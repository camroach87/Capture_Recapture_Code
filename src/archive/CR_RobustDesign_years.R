#### Preamble ####

# PROGRAM: Robust Design
# CODER: Cameron Roach
# DATE: 7/7/2013
# DESCRIPTION: Estimates population by splitting sampling timeline into 
# overlapping sections that are assumed to be closed. Calculate population in
# each of these "Year Groups" (YG) using closed methods.
# 
# aCaptYG(animal, year in YG, YG) is a 3d array showing capture
# history in year groups. Each year group has a length given by the 'length'
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
CR_RobustDesign_years <- function(length, sim="N", N.0=1000, t=14, p=0.05, phi=0.9, pBirth=0.1, pImmigration=0.1) {
  # length is length of time (plus one) where each animal is only captured at
  # most twice starting from any time
  # sim is "Y" or "N" variable that designates if this is a simulation
  #
  # TO DO: Would be nice if capture history matrix could be input for non simulations
  
  
  
  if (sim=="N")  {
    mtrxCapt <- mkCptrMtrx_years()
    years <- as.numeric(colnames(mtrxCapt))
    nyears <- length(years)
    #rownames(aCaptYG) <- rownames(mtrxCapt)
  } else if (sim=="Y") {
    nyears <- t
    years <- 1:nyears
    output <- mkOpenSimMtrx(N.0, nyears, p, phi, pBirth, pImmigration)
    mtrxCapt <- output[[1]]
    N.t <- output[[2]]
  } else {
    return("Didn't input Y or N for sim variable. Try again.")
  }
  
  
  # Calculate Chao's sparse data estimator
  aCaptYG <- array(0,c(nrow(mtrxCapt),length+1,nyears-length))
  aChao <- array(0,c(nyears-length,2))
  
  for (i in 1:(nyears-length)) {
    aCaptYG[,1:(1+length),i] <- mtrxCapt[,i:(i+length)]
    
    aChao[i,] <- unlist(ChaoEst(aCaptYG[,,i]))
  }
  
  # My basic smoothing algorithm
  sChao <- SmoothChao(aChao[,1],length,t)
  
  
  # Get CJS estimates from RMark
  p_rmark <- TestRMark(mtrxCapt)
  n.t <- colSums(mtrxCapt)
  N_rmark <- n.t/p_rmark
  
  #Plot estimates
  plot(N_rmark ~ years, type="l", col="red", ylim=c(0,max(sChao)))
  lines(years[1:nyears], sChao,col=colors()[25+6*length])
  #And actual if a simulation
  if (sim=="Y") {lines(years,N.t,col="green")}
  
  
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


SmoothChao <- function(aChao, length, t) {
  mSmooth <- matrix(NA, length+1, t)
  aCounts <- rep(NA,t)
  
  if(length(aChao)!=(t-length)) {
    cat("aChao array too small")
    return()
  }
  
  for(i in 1:(length+1)) {
    mSmooth[i,i:(i+(t-length-1))]<-aChao
    
#     for(j in i:(i+length-1)) {
#       mSmooth[i,j] <- aChao[j-i+1]
#     }
  }
  

  for(i in 1:t) {
    aCounts[i] <- sum(mSmooth[,i]!="NA", na.rm=TRUE)
  }
  aSmooth <- colSums(mSmooth, na.rm=TRUE)
  aSmooth <- aSmooth/aCounts
  
  return(aSmooth)
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
