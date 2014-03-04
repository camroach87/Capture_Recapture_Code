############################### PREAMBLE ################################
# PROGRAM: main
# CODER: Cameron Roach
# DATE: 24/2/2014
# DESCRIPTION: Sets directories and calls programs

rm(list=ls())

require(lubridate)
require(ggplot2)
require(plyr)


#### Set directories ####
dataDir <- "./data"
srcDir <- "./src"
outputDir <- "./output"



##### Source functions ####
source(file.path(srcDir,"mkPopMatrix.R"))     # population matrix functions
source(file.path(srcDir,"mkCaptureMatrix.R")) # capture matrix functions
source(file.path(srcDir,"getEstimators.R"))   # estimator functions





#### Initialise parameters ####
N.0   <- 5000
t     <- 20
p     <- 0.03
beta  <- 0.002

#### Generate population and capture matrices ####
#Uncomment scenario to run

# Open population sim
R <- 500 #burn in for population stability
mtrxPop <- Pop.Mat(N.0,t,R,beta)
mtrxCapt <- mkOpenSimMtrx(mtrxPop, t, p, beta)



# # Closed population sim




# # Trout Cod data analysis

data <- read.csv(file.path(dataDir,"TC_Data_Charles.csv"))
data$surveydate <- dmy_hm(data$surveydate)

# add year to data
data$year <- format(data$surveydate, "%Y")

mtrxCapt <- table(data$idfish, data$year)






# Calculate estimators
estCrRobust <- CR_RobustDesign(mtrxCapt[[1]],6)

#par(mfrow=c(1,1))
plot(estCrRobust, type="l", col="red", 
     ylim=c(min(mtrxCapt[[2]],estCrRobust[[2]]),max(mtrxCapt[[2]],estCrRobust[[2]])))
points(mtrxCapt[[2]], col="blue")



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