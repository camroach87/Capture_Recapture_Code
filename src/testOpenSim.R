# Open population sim
timer1 <- Sys.time()

# Initialise parameters
set.seed(1234)
N.0   <- 5000
t     <- 20
p     <- 0.02
beta  <- 0.003
window.val <- 8


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
estN.bs <- NULL
estN.bs.ci <- NULL



for (iS in 1:nCaptSims) {
  cat("Running capture simulation",iS,"of",nCaptSims,"...\n")
  
  #simulate a different capture matrix on each loop
  mtrxCapt <- mkOpenCaptMtrx(mtrxPop, t, p, beta)
  
  # Calculate estimators
  estN[["CR"]][iS,] <- calcCR(mtrxCapt,window.val)
  estN[["JS"]][iS,] <- calcJS(mtrxCapt)
  
  
  # Get bootstrap estimates
  nB <- 5000
  if (.Platform$OS.type == "windows") {
    # For windows use snow. Need to work out clustering to use more than one cpu - to do.
    estN.bs[[iS]] <- boot(data=mtrxCapt,statistic=CR.bs,R=nB,parallel="snow",ncpus=1,window=window.val)
  } else if (.Platform$OS.type == "unix") {
    # For linux use multicore. Don't need to worry about cluster stuff.
    estN.bs[[iS]] <- boot(data=mtrxCapt,statistic=CR.bs,R=nB,parallel="multicore",ncpus=3,window=window.val)
  }
  
}


# Calculate loop time and save workspace
timer2 <- Sys.time()
print(difftime(timer2,timer1,units="mins"))
fId <- file.path(outputDir,paste0("afterSimLoop_seed_1234_window_",window.val,".RData"))
save.image(file = fId)



## Bootstrap checks. Remember: index is the variable of interest, i.e. the 
## sampling occasion number. Interesting to look at edges vs middle.
# boot.ci(estN.bs[[iS]], index=10)
# plot(estN.bs[[iS]], index=1) 

## check if evidence to support normal
#   shapiro.test(estN.bs[[iS]]$t[,10])


# Average bootstrap confidence intervals for different capture matrices
estN.bs.se <- laply(estN.bs, function(x) apply(x$t, 2, sd))
estN.bs.se.mean <- apply(estN.bs.se, 2, mean)


# calculate means from repeated simulations
estN.mean[["Actual"]]  <- actN
estN.mean[["CR"]] <- apply(estN[["CR"]],2,mean)
estN.mean[["JS"]] <- apply(estN[["JS"]],2,mean)

# Calculate actual se and mse
estN.CR.se  <- apply(estN[["CR"]],2,sd)
estN.JS.se  <- apply(estN[["JS"]],2,sd)
estN.CR.mse  <- sapply(1:t, function(i) mse.f(estN[["CR"]][,i], actN[i]))
estN.JS.mse  <- sapply(1:t, function(i) mse.f(estN[["JS"]][,i], actN[i]))



# Convert to dataframe
tmp1 <- data.frame("Method" = "CR",
                   "Period" = 1:t,
                   "se" = estN.CR.se,
                   "se.bs" = estN.bs.se.mean,
                   "mse" = estN.CR.mse)
tmp2 <- data.frame("Method" = "JS",
                   "Period" = 1:t,
                   "se" = estN.JS.se,
                   "se.bs" = NA,
                   "mse" = estN.JS.mse)

error.df <- rbind(tmp1,tmp2)
rm(list=ls(pattern="tmp"))

# calculate two sided 95% CI
z.val <- qnorm(1-.05/2)
error.df$ci.95 <- z.val*error.df$se
error.df$ci.bs.95 <- z.val*error.df$se.bs

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