############################### PREAMBLE ################################
# PROGRAM: main
# CODER: Cameron Roach
# DATE: 24/2/2014
# DESCRIPTION: Sets directories and calls programs

rm(list=ls())

require(lubridate)
require(ggplot2)
require(plyr)
require(reshape2)


#### Set directories ####
dataDir <- "./data"
srcDir <- "./src"
outputDir <- "./output"



##### Source functions ####
source(file.path(srcDir,"mkPopMatrix.R"))     # population matrix functions
source(file.path(srcDir,"mkCaptureMatrix.R")) # capture matrix functions
source(file.path(srcDir,"getStats.R"))        # calculates various statistics
source(file.path(srcDir,"getEstimators.R"))   # estimator functions



testOpenSim  <- function() {
  #### Initialise parameters ####
  N.0   <- 5000
  t     <- 20
  p     <- 0.02
  beta  <- 0.003
  
  #### Generate population and capture matrices ####
  #Uncomment scenario to run
  
  ##### Open population sim ####
  
  # simulate population matrix
  burnIn  <- 500 #burn in for population stability
  mtrxPop <- Pop.Mat(N.0,t,burnIn,beta)
  actN    <- apply(mtrxPop,2,sum)
  nCaptSims   <- 100
  
  
  estN     <- list("CR" = matrix(NA,nCaptSims,t),
                   "JS" = matrix(NA,nCaptSims,t))
  estMSE   <- list("CR" = rep(NA,nCaptSims),
                   "JS" = rep(NA,nCaptSims))
  estVar   <- list("CR" = rep(NA,nCaptSims),
                   "JS" = rep(NA,nCaptSims))
  estNmean <- list("CR" = rep(NA,t),
                   "JS" = rep(NA,t))
  
  
  
  
  for (iS in 1:nCaptSims) {
    #simulate a different capture matrix on each loop
    mtrxCapt <- mkOpenCaptMtrx(mtrxPop, t, p, beta)
    
    # Calculate estimators
    estN[["CR"]][iS,] <- CR_RobustDesign(mtrxCapt,4)
    estN[["JS"]][iS,] <- calcJS(mtrxCapt)
    
    estMSE[["CR"]][iS] <- sum((estN[["CR"]][iS,]-actN)^2)/t
    estMSE[["JS"]][iS] <- sum((estN[["JS"]][iS,]-actN)^2)/t
    
    estVar[["CR"]][iS] <- sum((estN[["CR"]][iS,]-mean(estN[["CR"]][iS,]))^2)/t # pretty sure it is /t and not /(t-1) since we are looking at population variance, not sample variance (all t included)
    estVar[["JS"]][iS] <- sum((estN[["JS"]][iS,]-mean(estN[["JS"]][iS,]))^2)/t
    
  }
  
  estNmean[["CR"]] <- apply(estN[["CR"]],2,mean)
  estNmean[["JS"]] <- apply(estN[["JS"]],2,mean)
  
  #estBias <- how do I calculate bias here?
  
  
  #par(mfrow=c(1,1))
  yMin <- min(unlist(estNmean), actN, na.rm=T)
  yMax <- max(unlist(estNmean), actN, na.rm=T)
  plot(estNmean[["CR"]], type="p", col="red", 
       ylim=c(yMin,yMax))
  points(actN, col="blue")
  points(estNmean[["JS"]], col="green")
  legend("bottomright",legend=c("CR estimate", "JS estimate", "Actual population"),
         pch=c(1,1,1),col=c("red","green","blue"),cex=0.5)
}




# Trout Cod data analysis
testTroutCod <- function() {
  data <- read.csv(file.path(dataDir,"TC_Data_Charles.csv"))
  data$surveydate <- dmy_hm(data$surveydate)
  
  # add year to data
  data$year <- year(data$surveydate)
  data$month <- month(data$surveydate)
  
  mtrxCaptY <- table(data$idfish, data$year)
  mtrxCaptD <- table(data$idfish, data$surveydate)
  
  cat(length(data$idfish), "captures made.")
  cat(length(unique(data$idfish)), "unique fish.")
  
  # table of yearly and monthly captures
  acast(data,year~month,length)
  
  # table of f_i
  table(cbind(apply(mtrxCaptY,1,sum)))
  
  # table of next recapture year
  mtrxRecaptY <- mtrxCaptY
  for (i in 2:ncol(mtrxRecaptY)) {
    mtrxRecaptY[,i] <- mtrxRecaptY[,i]+mtrxRecaptY[,i-1]
  }
  mtrxRecaptY <- mtrxRecaptY - 1
  
  
  
  # Jolly Seber population estimates
  N_Y <- calcJS(mtrxCaptY)
  N_D <- calcJS(mtrxCaptD)
  
  plot(N_Y, type="l")
  plot(N_D, type="l")
  
  # apply Chao's sparse data estimator
  N_ChaoY <- calcChaoMt(mtrxCaptY)
  N_ChaoD <- calcChaoMt(mtrxCaptD)
}








#### Results ####

#Plot estimators
# #plot(N_rmark ~ sampOcc, type="l", col="red", ylim=c(0,max(sChao)))
# #lines(sampOcc[1:nsampOcc], sChao,col=colors()[25+6*window])
# if (sim=="Y") {
#   y.up <- max(N_rmark, sChao$y, N.t)
#   y.down <- min(N_rmark, sChao$y, N.t)
#   plot(sChao$y ~ c(1:nsampOcc), type="l", col="blue", ylim=c(y.down,y.up))
#   lines(N_rmark ~ c(1:nsampOcc), type="l", col="pink")
#   bla<-loess(N_rmark~c(1:nsampOcc))
#   lines(predict(bla),type="l", col="red")
#   lines(sampOcc,N.t,col="green")
# } else {
#   y.up <- max(N_rmark, sChao$y)
#   y.down <- min(N_rmark, sChao$y)
#   plot(sChao$y ~ c(1:nsampOcc), type="l", col="blue", ylim=c(y.down,y.up))
#   lines(N_rmark ~ c(1:nsampOcc), type="l", col="pink")
#   bla<-loess(N_rmark~c(1:nsampOcc))
#   lines(predict(bla),type="l", col="red")
# }