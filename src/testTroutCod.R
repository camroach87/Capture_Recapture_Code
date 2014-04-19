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
ggsave(plot1,file=file.path(plotDir,"TC_JS_y.png"),width=14,height=8)
ggsave(plot2,file=file.path(plotDir,"TC_JS_d.png"),width=14,height=8)
ggsave(plot3,file=file.path(plotDir,"TC_JS_y_d.png"),width=14,height=8)
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