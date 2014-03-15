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





#### Initialise parameters ####
N.0   <- 5000
t     <- 20
p     <- 0.02
beta  <- 0.002

#### Generate population and capture matrices ####
#Uncomment scenario to run

##### Open population sim ####

# simulate population matrix
burnIn  <- 500 #burn in for population stability
mtrxPop <- Pop.Mat(N.0,t,burnIn,beta)
actN    <- apply(mtrxPop,2,sum)


nSims   <- 100
estNcr  <- matrix(NA,nSims,t)
estMSE  <- rep(NA,nSims)
estVar  <- rep(NA,nSims)
estBias <- rep(NA,nSims)

for (iS in 1:nSims) {
  #simulate a different capture matrix on each loop
  mtrxCapt <- mkOpenSimMtrx(mtrxPop, t, p, beta)
  
  # Calculate estimators
  estNcr[iS,] <- CR_RobustDesign(mtrxCapt,4)
  estMSE[iS] <- sum((estNcr[iS,]-actN)^2)/t
  estVar[iS] <- sum((estNcr[iS,]-mean(estNcr[iS,]))^2)/t # pretty sure it is /t and not /(t-1) since we are looking at population variance, not sample variance (all t included)
  #estBias[iS] <- 
}

estNmean <- apply(estNcr,2,mean)
estBias <- estN

#par(mfrow=c(1,1))
plot(estNmean, type="p", col="red", 
     ylim=c(min(actN,estNmean),max(actN,estNmean)))
points(actN, col="blue")
legend("bottomright",legend=c("Estimated population","Actual population"),
       lty=c(1,1),col=c("red","blue"),cex=0.5)


# Closed population sim
# To do (maybe not required)



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
  
  N_Y <- calcJS(mtrxCaptY)
  N_D <- calcJS(mtrxCaptD)
  
  plot(N_Y, type="l")
  plot(N_D, type="l")
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