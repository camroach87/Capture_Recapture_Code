############################### PREAMBLE ################################
# PROGRAM: main
# CODER: Cameron Roach
# DATE: 24/2/2014
# DESCRIPTION: Sets directories and calls programs


#### Set directories ####
dataDir <- file.path("./data")
srcDir <- file.path("./src")
outputDir <- file.path("./output")



##### Source functions ####
source(file.path(srcDir,"mkPopMatrix.R"))     # population matrix functions
source(file.path(srcDir,"mkCaptureMatrix.R")) # capture matrix functions
source(file.path(srcDir,"getEstimators.R"))   # estimator functions





#### Usage ####

# Parameters
N.0   <- 5000
t     <- 20
p     <- 0.03
beta  <- 0.002

# Generate capture matrix
mtrxCapt <- mkOpenSimMtrx(N.0, t, p, beta)

# Calculate estimators
popEst <- CR_RobustDesign(mtrxCapt[[1]],6)

#par(mfrow=c(1,1))
plot(popEst, type="l", col="red", 
     ylim=c(min(mtrxCapt[[2]],popEst[[2]]),max(mtrxCapt[[2]],popEst[[2]])))
points(mtrxCapt[[2]], col="blue")



#### Results ####

# #Plot estimates
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