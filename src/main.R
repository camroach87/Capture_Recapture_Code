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
  # Open population sim
  timer1 <- Sys.time()
  
  # Initialise parameters
  set.seed(1234)
  N.0   <- 5000
  t     <- 20
  p     <- 0.02
  beta  <- 0.003
  window <- 8
  
  
  # simulate population matrix
  burnIn  <- 500 #burn in for population stability
  mtrxPop <- Pop.Mat(N.0,t,burnIn,beta)
  actN    <- apply(mtrxPop,2,sum)
  nCaptSims   <- 100
  
  
  # initialise lists
  estN     <- list("CR" = matrix(NA,nCaptSims,t),
                   "JS" = matrix(NA,nCaptSims,t))
  estMSE   <- list("CR" = rep(NA,nCaptSims),
                   "JS" = rep(NA,nCaptSims))
  estVar   <- list("CR" = rep(NA,nCaptSims),
                   "JS" = rep(NA,nCaptSims))
  estN.mean <- list("CR" = rep(NA,t),
                    "JS" = rep(NA,t))
  estN.bs.sd <- matrix(NA,nCaptSims,t)
  
  
  
  for (iS in 1:nCaptSims) {
    cat("Running capture simulation",iS,"of",nCaptSims,"...\n")
    
    #simulate a different capture matrix on each loop
    mtrxCapt <- mkOpenCaptMtrx(mtrxPop, t, p, beta)
    
    # Calculate estimators
    estN[["CR"]][iS,] <- CR_RobustDesign(mtrxCapt,window)
    estN[["JS"]][iS,] <- calcJS(mtrxCapt)
    
    #     # we are ignoring the effects of correlation on MSE and Var - needs further investigation
    #     estMSE[["CR"]][iS] <- sum((estN[["CR"]][iS,]-actN)^2)/t
    #     estMSE[["JS"]][iS] <- sum((estN[["JS"]][iS,]-actN)^2)/t
    #     
    #     estVar[["CR"]][iS] <- sum((estN[["CR"]][iS,]-mean(estN[["CR"]][iS,]))^2)/t # pretty sure it is /t and not /(t-1) since we are looking at population variance, not sample variance (all t included)
    #     estVar[["JS"]][iS] <- sum((estN[["JS"]][iS,]-mean(estN[["JS"]][iS,]))^2)/t
    
    # Get bootstrap samples
    nB <- 999
    estN.bs <- matrix(NA,nB,t)
    for (iB in 1:nB) {
      idx <- sample(nrow(mtrxCapt),nrow(mtrxCapt),replace=TRUE)
      mtrxCapt.bs <- mtrxCapt[idx,]
      estN.bs[iB,] <- CR_RobustDesign(mtrxCapt.bs,window)
    }
    
    # Calculate bootstrap estimates (first step)
    estN.bs.mean <- apply(estN.bs,2,mean)
    tmp <- NULL
    for (iB in 1:nB) {
      tmp <- rbind(tmp,(estN.bs.mean-estN.bs[iB,])^2)
    }
    estN.bs.sd[iS,] <- sqrt(1/(nB-1)*apply(tmp,2,sum))
    
  }
  
  # Calculate loop time and save workspace
  timer2 <- Sys.time()
  print(difftime(timer2,timer1,units="mins"))
  fId <- file.path(outputDir,paste0("afterSimLoop_seed_1234_window_",window,".RData"))
  save.image(file = fId)
  
  # calculate means from repeated simulations
  estN.mean[["Actual"]]  <- actN
  estN.mean[["CR"]] <- apply(estN[["CR"]],2,mean)
  estN.mean[["JS"]] <- apply(estN[["JS"]],2,mean)
  
  # Effectively, these bootstrapped se values have been averaged across many
  # different capture matrices.
  estN.bs.sd.mean <- apply(estN.bs.sd,2,mean)
  
  
  # Calculate actual sd and mse
  estN.CR.sd  <- apply(estN[["CR"]],2,sd)
  estN.JS.sd  <- apply(estN[["JS"]],2,sd)
  
  estN.CR.mse  <- sapply(1:t, function(i) mse.f(estN[["CR"]][,i], actN[i]))
  estN.JS.mse  <- sapply(1:t, function(i) mse.f(estN[["JS"]][,i], actN[i]))
  
  # Convert to dataframe
  tmp1 <- data.frame("Method" = "CR",
                     "Period" = 1:t,
                     "sd" = estN.CR.sd,
                     "sd.bs" = estN.bs.sd.mean,
                     "mse" = estN.CR.mse)
  tmp2 <- data.frame("Method" = "JS",
                     "Period" = 1:t,
                     "sd" = estN.JS.sd,
                     "sd.bs" = NA,
                     "mse" = estN.JS.mse)
  
  error.df <- rbind(tmp1,tmp2)
  rm(list=ls(pattern="tmp"))
  
  # calculate two sided 95% CI
  z.val <- qnorm(1-.05/2)
  error.df$ci.95 <- z.val*error.df$sd/sqrt(nCaptSims)
  error.df$ci.bs.95 <- z.val*error.df$sd.bs/sqrt(nCaptSims)
  
  #estBias <- how do I calculate bias here?
  
  # Munge data and then plot
  estN.df <- as.data.frame(estN.mean)
  estN.df$Period <- c(1:t)
  estN.tidy <- melt(estN.df, "Period", variable.name = "Method", value.name = "N")
  estN.tidy$Method <- factor(estN.tidy$Method, levels=c("Actual","CR","JS"))
  
  #add errors to N estimates
  estN.tidy <- merge(estN.tidy,error.df,by=c("Method","Period"),all.x=TRUE)
  
  # plots
  pd <- position_dodge(0.1)
  
  # plot of bootstrapped confidence intervals for CR
  plot1 <- ggplot(estN.tidy[estN.tidy$Method!="JS",], aes(x=Period, y=N, colour=Method)) + 
    geom_errorbar(aes(ymin=N-ci.bs.95, ymax=N+ci.bs.95), width=.5, alpha=0.4) +
    geom_line() + 
    geom_point(size=3,shape=21,fill="white") +
    ggtitle("Abundance estimates of simulated open population. Bootstrap CI.") +
    theme_bw()
  
  # plot of actual confidence intervals for CR
  plot2 <- ggplot(estN.tidy[estN.tidy$Method!="JS",], aes(x=Period, y=N, colour=Method)) + 
    geom_errorbar(aes(ymin=N-ci.95, ymax=N+ci.95), width=.5, alpha=0.4) +
    geom_line() + 
    geom_point(size=3,shape=21,fill="white") +
    ggtitle("Abundance estimates of simulated open population. Actual CI") +
    theme_bw()
  
  # Plot of JS and CR with CIs
  plot3 <- ggplot(estN.tidy, aes(x=Period, y=N, colour=Method)) +
    geom_errorbar(aes(ymin=N-ci.95, ymax=N+ci.95), position=pd, width=.5, alpha=0.4) +
    geom_line() + 
    geom_point(position=pd, size=3,shape=21,fill="white") +
    ggtitle("Abundance estimates of simulated open population") +
    theme_bw()
  
  
  print(plot1)
  print(plot2)
  print(plot3)
  
}




# Trout Cod data analysis
testTroutCod <- function() {
  data <- read.csv(file.path(dataDir,"TC_Data_Charles.csv"))
  data$surveydate <- dmy_hm(data$surveydate)
  data <- data[order(data$surveydate),]
  
  # add year to data
  data$year <- year(data$surveydate)
  data$month <- month(data$surveydate)
  
  mtrxCaptY <- table(data$idfish, data$year)
  mtrxCaptD <- table(data$idfish, data$surveydate)
  
  cat(length(data$idfish), "captures made.\n")
  cat(length(unique(data$idfish)), "unique fish.\n\n")
  
  # table of yearly and monthly captures
  cat("Table of yearly and monthly captures...\n")
  print(acast(data,year~month,length, value.var="month"))
  cat("\n")
  
  # table of f_i
  cat("Table of f_i...")
  print(table(cbind(apply(mtrxCaptY,1,sum))))
  cat("\n")
  
  
  
  #   # table of next recapture year
  #   mtrxRecaptY <- mtrxCaptY
  #   for (i in 2:ncol(mtrxRecaptY)) {
  #     mtrxRecaptY[,i] <- mtrxRecaptY[,i]+mtrxRecaptY[,i-1]
  #   }
  #   mtrxRecaptY <- mtrxRecaptY - 1
  
  
  
  
  # Jolly Seber population estimates
  N_Y <- calcJS(mtrxCaptY)
  N_Y.df <- as.data.frame(N_Y)
  #Taking 1st June to be estimate time for years as most sampling happens
  #around/before this time.
  N_Y.df$Date <- ymd(paste0(rownames(N_Y.df),"-05-01"))
  N_Y.df$Occasion <- c(1:(dim(N_Y.df)[1]))
  N_Y.df$Method <- "JS on yearly grouped data"
  rownames(N_Y.df)  <- NULL
  colnames(N_Y.df)[1] <- "N"
  
  N_D <- calcJS(mtrxCaptD)
  N_D.df <- as.data.frame(N_D)
  N_D.df$Date <- ymd(rownames(N_D.df))
  N_D.df$Occasion <- c(1:(dim(N_D.df)[1]))
  N_D.df$Method <- "JS on daily grouped data"
  rownames(N_D.df)  <- NULL
  colnames(N_D.df)[1] <- "N"
  
  
  # Combine all data frames
  estN.tidy <- rbind(N_D.df, N_Y.df)
  
  # Plots
  tmp <- estN.tidy[estN.tidy$Method=="JS on yearly grouped data",]
  plot1 <- ggplot(tmp, aes(x=Date, y=N, colour=Method)) + geom_line() + 
    ggtitle("JS estimate of abundance for yearly grouping of TC capture data.") +
    theme_bw() +
    theme(legend.position="bottom")
  
  tmp <- estN.tidy[estN.tidy$Method=="JS on daily grouped data",]
  plot2 <- ggplot(tmp, aes(x=Date, y=N, colour=Method)) + geom_line() + 
    ggtitle("JS estimate of abundance for daily grouping of TC capture data.") +
    theme_bw() +
    theme(legend.position="bottom")
  
  plot3 <- ggplot(estN.tidy, aes(x=Date, y=N, colour=Method)) + geom_line() + 
    ggtitle("Comparison of JS estimate of abundance for daily and yearly grouping of TC capture data.") +
    theme_bw() +
    theme(legend.position="bottom")
  
  print(plot1)
  print(plot2)
  print(plot3)
  ggsave(plot1,file=file.path(outputDir,"plots/TC_JS_y.png"),width=14,height=8)
  ggsave(plot2,file=file.path(outputDir,"plots/TC_JS_d.png"),width=14,height=8)
  ggsave(plot3,file=file.path(outputDir,"plots/TC_JS_y_d.png"),width=14,height=8)
  rm(tmp)
  
  
  
  
  #   # Chao's sparse data estimator for closed populations
  #   N_ChaoY <- calcChaoMt(mtrxCaptY)
  #   N_ChaoD <- calcChaoMt(mtrxCaptD)
  
  # Get dates and occasion from N_D.df. Bit hacky...
  dates.occ <- N_D.df[,c("Date","Occasion")]
  #dates.diff <- max(dates.occ$Date)-min(dates.occ$Date)+1
  #dates.all <- min(dates.occ$Date) + days(0:dates.diff)
  #dates.all <- data.frame("Date" = dates.all)
  #dates.all <- merge(dates.all,dates.occ,all.x=TRUE)
  
  
  
  
  
  # Apply CR estimator to daily trout cod data
  # Investigate impact of different window sizes
  estN.CR <- data.frame()
  CR.levels <- NULL
  for (iW in seq(10,100,10)) {
    cat("Calculating CR abundance estimate for window size",iW,"...\n")
    tmp <- CR_RobustDesign(mtrxCaptD, iW)
    tmp <- as.data.frame(tmp)
    tmp$Occasion <- c(1:(dim(tmp)[1]))
    tmp <- merge(tmp, dates.occ)
    tmp$Method <- paste("CR window size",iW)
    colnames(tmp)[2] <- "N"
    estN.CR <- rbind(estN.CR,tmp)
    CR.levels <- rbind(CR.levels,paste("CR window size",iW))
  }
  
  # Plot window size impact for CR estimate of TC population
  estN.CR$Method <- factor(estN.CR$Method, levels=CR.levels)
  plot1 <- ggplot(estN.CR, aes(x=Occasion, y=N, colour=Method)) + geom_line() + 
    ggtitle("CR estimates of abundance for TC capture data.") +
    theme_bw()
  tmp <- estN.CR[estN.CR$Method=="CR window size 20",]
  plot2 <- ggplot(tmp, aes(x=Occasion, y=N, colour=Method)) + geom_line() + 
    ggtitle("CR estimates of abundance for TC capture data.") +
    theme_bw()
  
  print(plot1)
  print(plot2)
}


