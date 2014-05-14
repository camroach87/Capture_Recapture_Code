source("./src/init.R")

# input
# choose from young, adult or all
dataFilter <- "adult"

# load data
data <- read.csv(file.path(dataDir,"TC_Data_Charles.csv"))
data$surveydate <- dmy_hm(data$surveydate)
data <- data[order(data$surveydate),]

# add year to data
data$year <- year(data$surveydate)
data$month <- month(data$surveydate)

# select data filter
if (dataFilter == "young") {
  # 276 is median length. Used to split young and adults.
  data <- data[data$totallength<=276,]
} else if (dataFilter == "adult") {
  data <- data[data$totallength>276,]
}

  
  

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
tmp <- calcJS(mtrxCaptY)
N_Y <- tmp[[1]]
N_Y.ci <- tmp[[2]]
N_Y.df <- as.data.frame(N_Y)
#Taking 1st June to be estimate time for years as most sampling happens
#around/before this time.
N_Y.df$Date <- ymd(paste0(rownames(N_Y.df),"-05-01"))
N_Y.df$Occasion <- c(1:(dim(N_Y.df)[1]))
N_Y.df$Method <- "JS on yearly grouped data"
rownames(N_Y.df)  <- NULL
colnames(N_Y.df)[1] <- "N"

tmp <- calcJS(mtrxCaptD)
N_D <- tmp[[1]]
N_D.ci <- tmp[[2]]
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
plot1 <- ggplot(tmp, aes(x=Date, y=N)) + geom_line() + 
  ggtitle("JS estimate of abundance for yearly grouping of TC capture data.") +
  theme_bw()

tmp <- estN.tidy[estN.tidy$Method=="JS on daily grouped data",]
plot2 <- ggplot(tmp, aes(x=Date, y=N)) + geom_line() + 
  ggtitle("JS estimate of abundance for daily grouping of TC capture data.") +
  theme_bw()

plot3 <- ggplot(estN.tidy, aes(x=Date, y=N, linetype=Method)) + geom_line() + 
  ggtitle("Comparison of JS estimate of abundance for daily and yearly grouping of TC capture data.") +
  theme_bw() +
  theme(legend.position="bottom")

print(plot1)
print(plot2)
print(plot3)
ggsave(plot1,file=file.path(plotDir,paste0("TC_JS_y_",dataFilter,".png")),width=14,height=8)
ggsave(plot2,file=file.path(plotDir,paste0("TC_JS_d_",dataFilter,".png")),width=14,height=8)
ggsave(plot3,file=file.path(plotDir,paste0("TC_JS_y_d_",dataFilter,".png")),width=14,height=8)
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
  tmp <- calcCR(mtrxCaptD, iW)
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
plot1 <- ggplot(estN.CR, aes(x=Occasion, y=N)) + 
  geom_line() + 
  facet_wrap(~Method,ncol=2) +
  ggtitle("CR estimates of abundance for TC capture data.") +
  theme_bw()
estN.CR.20 <- estN.CR[estN.CR$Method=="CR window size 20",]
plot2 <- ggplot(estN.CR.20, aes(x=Occasion, y=N)) + 
  geom_line() + 
  ggtitle("CR estimates of abundance for TC capture data.") +
  theme_bw()

print(plot1)
print(plot2)
ggsave(plot1,file=file.path(plotDir,paste0("TC_CR_window_all_",dataFilter,".png")),width=8.3,height=11.7)
ggsave(plot2,file=file.path(plotDir,paste0("TC_CR_window_20_",dataFilter,".png")),width=14,height=8)



# Produce same plots but with BCa confidence intervals

# First bootstrap
nB <- 5000
window.val <- 20
timer1 <- Sys.time()
if (.Platform$OS.type == "windows") {
  cl <- makeCluster(4, type = "SOCK")
  registerDoSNOW(cl)
  clusterEvalQ(cl, library(boot))
  clusterEvalQ(cl, library(KernSmooth))
  clusterExport(cl, c("calcCR","calcChaoMt"))
  
  estN.bs <- boot(data=mtrxCaptD,statistic=CR.bs,R=nB,parallel="snow",ncpus=4,window=window.val,cl=cl)
  
  stopCluster(cl)
} else if (.Platform$OS.type == "unix") {
  estN.bs <- boot(data=mtrxCaptD,statistic=CR.bs,R=nB,parallel="multicore",ncpus=4,window=window.val)
}

timer2 <- Sys.time()
print(difftime(timer2,timer1,units="mins"))
fId <- file.path(outputDir,paste0("bs_TC_window_",window.val,"_",dataFilter,".RData"))
save(estN.bs, file=fId)


# Then calculate bootstrap percentile CIs
t <- ncol(mtrxCaptD)
ci.perc <- sapply(1:t, 
                  function(i) {
                    cat(paste0("Running index ", i, "...\n"))
                    ci.tmp <- boot.ci(estN.bs, index = i, type="perc")
                    ci.tmp <- ci.tmp$perc[c(4,5)]
                  })
bs.perc.ci <- data.frame("bs.ci.l" = ci.perc[1,],
                     "bs.ci.u" = ci.perc[2,],
                     "Occasion" =1:t)

# plot trout cod estimate and percentile confidence intervals
estN.CR.20 <- merge(estN.CR.20,bs.perc.ci)

p1 <- ggplot(estN.CR.20, aes(x=Occasion, y=N)) + 
  #geom_errorbar(aes(ymin=bs.ci.l, ymax=bs.ci.u), width=.5, alpha=0.4) +
  geom_ribbon(aes(ymin=bs.ci.l, ymax=bs.ci.u), alpha=0.2) +
  geom_line() + 
  ggtitle("CR estimates of abundance for TC capture data.") +
  theme_bw()
print(p1)
ggsave(file=file.path(plotDir,paste0("TC_CR_window_20_CI_",dataFilter,".png")),width=14,height=8)
