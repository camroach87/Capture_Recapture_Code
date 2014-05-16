source("./src/init.R")

# Open population sim
timer1 <- Sys.time()

# Initialise parameters
set.seed(1234)
N.0   <- 5000
t     <- 20
p     <- 0.02
beta  <- 0.003
window.val <- 4


# simulate population matrix
burnIn  <- 500 #burn in for population stability
mtrxPop <- Pop.Mat(N.0,t,burnIn,beta)
actN    <- apply(mtrxPop,2,sum)
nCaptSims   <- 1000


# initialise lists
estN     <- list("CR" = matrix(NA,nCaptSims,t),
                 "JS" = matrix(NA,nCaptSims,t))
estMSE   <- list("CR" = rep(NA,nCaptSims),
                 "JS" = rep(NA,nCaptSims))
estVar   <- list("CR" = rep(NA,nCaptSims),
                 "JS" = rep(NA,nCaptSims))
estN.mean <- list("CR" = rep(NA,t),
                  "JS" = rep(NA,t))
estN.JS.ci.l <- matrix(NA,nCaptSims,t)
estN.JS.ci.u <- matrix(NA,nCaptSims,t)
estN.bs <- NULL
estN.bs.ci <- NULL



for (iS in 1:nCaptSims) {
  cat("Running capture simulation",iS,"of",nCaptSims,"...\n")
  
  #simulate a different capture matrix on each loop
  mtrxCapt <- mkOpenCaptMtrx(mtrxPop, t, p, beta)
  
  # Calculate JS pop est. and Manly C.I.
  tmp <- calcJS(mtrxCapt)
  estN[["JS"]][iS,] <- tmp[[1]]
  estN.JS.ci.l[iS,] <- tmp[[2]]$ci.l
  estN.JS.ci.u[iS,] <- tmp[[2]]$ci.u
  
  
  # Calculate CR population estimate and then bootstrap
  estN[["CR"]][iS,] <- calcCR(mtrxCapt,window.val)
  
  nB <- 5000
  if (.Platform$OS.type == "windows") {
    # For windows use snow an clustering.
    cl <- makeCluster(4, type = "SOCK")
    registerDoSNOW(cl)
    clusterEvalQ(cl, library(boot))
    clusterEvalQ(cl, library(KernSmooth))
    clusterExport(cl, c("calcCR","calcChaoMt"))
    
    estN.bs[[iS]] <- boot(data=mtrxCapt,statistic=CR.bs,R=nB,parallel="snow",ncpus=4,window=window.val,cl=cl)
    
    stopCluster(cl)
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
estN.bs.se.med <- apply(estN.bs.se, 2, median)

# calculate means from repeated simulations
estN.mean[["Actual"]]  <- actN
estN.mean[["CR"]] <- apply(estN[["CR"]],2,mean)
estN.mean[["JS"]] <- apply(estN[["JS"]],2,mean)

# Calculate actual se and mse
estN.CR.se  <- apply(estN[["CR"]],2,sd)
estN.JS.se  <- apply(estN[["JS"]],2,sd)
estN.CR.mse  <- sapply(1:t, function(i) mse.f(estN[["CR"]][,i], actN[i]))
estN.JS.mse  <- sapply(1:t, function(i) mse.f(estN[["JS"]][,i], actN[i]))

# Get confidence intervals
estN.JS.ci.l.mean <- apply(estN.JS.ci.l,2,mean)
estN.JS.ci.u.mean <- apply(estN.JS.ci.u,2,mean)

CR.ci.l <- estN.mean[["CR"]] - qnorm(1-.05/2)*estN.CR.se
CR.ci.u <- estN.mean[["CR"]] + qnorm(1-.05/2)*estN.CR.se


#Calculate BCa confidence intervals. Note: Clustering hasn't been tested on
#linux systems. Recommend breaking foreach loop into chunks and joining lists.
cl <- makeCluster(4, type = "SOCK")
registerDoSNOW(cl)
clusterEvalQ(cl, library(boot))
clusterExport(cl, c("estN.bs"))
writeLines(c(""), "./output/log.txt")

bs.bca.ci <- foreach (iS=1:nCaptSims) %dopar% {
  sink("./output/log.txt", append=TRUE)
  ci.bca <- sapply(1:t, 
                   function(i) {
                     cat(paste0("Running sim ", iS," index ", i, "...\n"))
                     ci.tmp <- boot.ci(estN.bs[[iS]], index = i, type="bca")
                     ci.tmp <- ci.tmp$bca[c(4,5)]
                   })
  sink()
  output <- data.frame("bs.ci.l" = ci.bca[1,],
                       "bs.ci.u" = ci.bca[2,],
                       "Occasion" =1:t)
}
stopCluster(cl)

fId <- file.path(outputDir,paste0("bs_bca_ci_window_",window.val,".RData"))
save(bs.bca.ci, file=fId)


# calculate bootstrap percentile confidence intervals
bs.perc.ci <- foreach (iS=1:nCaptSims) %do% {
  ci.perc <- sapply(1:t, 
                    function(i) {
                      cat(paste0("Running sim ", iS," index ", i, "...\n"))
                      ci.tmp <- boot.ci(estN.bs[[iS]], index = i, type="perc")
                      ci.tmp <- ci.tmp$perc[c(4,5)]
                    })
  output <- data.frame("bs.ci.l" = ci.perc[1,],
                       "bs.ci.u" = ci.perc[2,],
                       "Occasion" =1:t)
}


# calculate coverage probabilities for BCa and percentile methods
calcCoverage(actN, bs.bca.ci)
calcCoverage(actN, bs.perc.ci)


# calculate mean confidence intervals for bootstraps
CR.bs.ci <- do.call("rbind", bs.perc.ci)
CR.bs.ci.mean <- ddply(CR.bs.ci,
                       .(Occasion),
                       function(x) {
                         bs.ci.l <- mean(x$bs.ci.l)
                         bs.ci.u <- mean(x$bs.ci.u)
                         data.frame("bs.ci.l"=bs.ci.l,
                                    "bs.ci.u"=bs.ci.u)
                       })




# Convert to dataframe
CR.df <- data.frame("Method" = "CR",
                   "Occasion" = 1:t,
                   "se" = estN.CR.se,
                   "se.bs" = estN.bs.se.mean,
                   "mse" = estN.CR.mse,
                   "ci.l" = CR.ci.l,
                   "ci.u" = CR.ci.u,
                   "ci.l.bs" = CR.bs.ci.mean$bs.ci.l,
                   "ci.u.bs" = CR.bs.ci.mean$bs.ci.u)
JS.df <- data.frame("Method" = "JS",
                   "Occasion" = 1:t,
                   "se" = estN.JS.se,
                   "se.bs" = NA,
                   "mse" = estN.JS.mse,
                   "ci.l" = estN.JS.ci.l.mean,
                   "ci.u" = estN.JS.ci.u.mean,
                   "ci.l.bs" = NA,
                   "ci.u.bs" = NA)

error.df <- rbind(CR.df,JS.df)
rm(CR.df,JS.df)





# Munge data and then plot
estN.df <- as.data.frame(estN.mean)
estN.df$Occasion <- c(1:t)
estN.tidy <- melt(estN.df, "Occasion", variable.name = "Method", value.name = "N")
estN.tidy$Method <- factor(estN.tidy$Method, levels=c("Actual","CR","JS"))

#add errors to N estimates
estN.tidy <- merge(estN.tidy,error.df,by=c("Method","Occasion"),all.x=TRUE)

# plots
pd <- position_dodge(0.1)

# plot of bootstrapped confidence intervals for CR
idx <- estN.tidy$Method!="JS"
plot1 <- ggplot(estN.tidy[idx,], aes(x=Occasion, y=N, linetype=Method)) + 
  geom_errorbar(aes(ymin=ci.l.bs, ymax=ci.u.bs), width=.5, alpha=0.4) +
  geom_line() + 
  geom_point(size=3,shape=21,fill="white") +
  ylim(1500,4000) +
  ggtitle("Abundance estimates of simulated open population. Bootstrap CI.") +
  theme_bw()

# plot of actual confidence intervals for CR
idx <- estN.tidy$Method!="JS"
plot2 <- ggplot(estN.tidy[idx,], aes(x=Occasion, y=N, linetype=Method)) + 
  geom_errorbar(aes(ymin=ci.l, ymax=ci.u), width=.5, alpha=0.4) +
  geom_line() + 
  geom_point(size=3,shape=21,fill="white") +
  ylim(1500,4000) +
  ggtitle("Abundance estimates of simulated open population. Actual CI.") +
  theme_bw()

# Plot of JS and CR with CIs
idx <- estN.tidy$Method!="CR"
plot3 <- ggplot(estN.tidy[idx,], aes(x=Occasion, y=N, linetype=Method)) +
  geom_errorbar(aes(ymin=ci.l, ymax=ci.u), position=pd, width=.5, alpha=0.4) +
  geom_line() + 
  geom_point(position=pd, size=3,shape=21,fill="white") +
  ggtitle("Abundance estimates of simulated open population") +
  theme_bw()

plot4 <- qplot(x=estN.CR.se,y=estN.bs.se.med) + 
  #geom_smooth(method="lm", colour="black") +
  ggtitle("Median bootstrap standard error vs standard deviation for N for 20 capture occasions.") +
  xlab("Median bootstrap standard error") +
  ylab("Standard deviation of N") +
  geom_abline(intercept=0, slope=1, linetype=2) +
  theme_bw()
  


print(plot1)
print(plot2)
print(plot3)
print(plot4)
ggsave(plot1,file=file.path(plotDir,"OpenSim_CR_bs.png"),width=14,height=8)
ggsave(plot2,file=file.path(plotDir,"OpenSim_CR.png"),width=14,height=8)
ggsave(plot3,file=file.path(plotDir,"OpenSim_JS.png"),width=14,height=8)
ggsave(plot4,file=file.path(plotDir,"BS_se_vs_N_sd.png"),width=14,height=8)
