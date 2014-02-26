#### 23/9/2013 ####
CR_RobustDesign(10,sim="Y", N.0=4000, p=0.02, nsampOcc=200, pBirth=0.03)

#### 22/7/2013 ####
mkOpenSimMtrx(N.0=5000, t=14, p=0.05, phi=0.9,pBirth=0.25, pImmigration=0.25)
table(rowSums(output[[1]][,1:8]))
CR_RobustDesign(1, sim="Y", N.0=5000, p=0.02, phi=0.9,pBirth=0.25, pImmigration=0.25)


#### 19/7/2013 ####
#Comparing Chao estimator to CJS
mtrxCapt <- mkCloseSimMtrx(1000,15,0.01)
p <- TestRMark(mtrxCapt)
colSums(mtrxCapt)/p
mean(colSums(mtrxCapt)/p)
ChaoEst(mtrxCapt)






#### 16/7/2013 ####
# output capture matrix to a text file
# Need to source CR_RobustDesign.R
mtrxCapt[mtrxCapt>1]<-1
write.table(as.matrix(mtrxCapt),"Programs/Data/MarkTest1.inp",col.names = F, row.names = F, sep="", quote=FALSE)
#Format for RMARK
library(RMark)
mtrxCapt_f <- apply(format(mtrxCapt), 1, paste, collapse="")
mtrxCapt_f <- paste(mtrxCapt_f, " 1;")
write.table(as.matrix(mtrxCapt_f),"Programs/Data/MarkTest1.inp",col.names = F, row.names = F, sep="", quote=FALSE)
# Read by RMark
mtrxCapt_f <- convert.inp("Programs/Data/MarkTest1.inp")
mark(mtrxCapt_f)








#### 9/7/2013 ####
# Testing CR_RobustDesign.R
bla <- CR_RobustDesign(3,"Y", pCapture=0.03, pBirth=0.2)
plot(1:11, bla[[1]][,1], type="l", col="green", ylim=c(0,5000))
lines(1:14, bla[[2]])







############################### 24//2012 ###############################

#adult fish
data1 <- data[data$totallength>=250,]
JollySeber(data1)

data2 <- data[data$totallength<250,]
JollySeber(data2)







#########################################################################
############################### 24/9/2012 ###############################
#########################################################################


#Jolly Seber
n <- CalcTotal(data)
m <- CalcMarked(data)
R <- CalcReleased(data)
Z <- CalcZ(data)
r <- CalcRecapt(data)

M <- m + R*Z/r
N <- n*M[1:14]/m

plot(N ~ years, type="l", col="green")

#Peterson Lincoln estimator
Npl <-c()
for (i in 2:13) Npl[i] <- n[i]n[i-1]/m[i]






#########################################################################
############################### 23/9/2012 ###############################
#########################################################################

# 23/9/2012 old version of CalcDistinct - not functioning correctly
CalcDistinct <- function(data) {
  years <- sort(unique(data$year))
  minYear <- min(years)
  maxYear <- max(years)
  curYear <- minYear
  k <- maxYear - minYear
  
  captFish <- NA
  fishYear <- list(NA)
  M <- NA
  
  
  for (j in 1:(k+1)) {
    fishYear[[j]] <- data$idfish[data$year == j + minYear - 1]
  }
  
  
  # calculate distinct fish captured each year
  for (j in 1:(k+1)) {
    
    if (j != 1) {
      
      ### ERROR HERE
      
      M[j] <- sum(!(fishYear[[j]] %in% captFish))
      
      ###not handling first distinct occasion correctly
      
    } else {
      M[1] <- 0 
    }
    
    if (j == 1) {
      captFish <- fishYear[[1]]
    } else {
      captFish <- append(captFish, fishYear[[j]])
    }
    
  }
  
  M[k+1] <- sum(M)
  
  return(M)
  
}



#########################################################################
############################### OLD #####################################
#########################################################################


#plot(table(D)) date values are spaced sequentially rather than in their correct time position.


# thought - maybe just add a new column to existing data frame with the year each fish was caught - will then be able to easily sort
# the below program would be useful if captures did not occur yearly, but since they do there is little reason to do this (although it would be fun)
GetYearCounts <- function(tabD) {
  dStart <- as.Date("1999/1/1")
  dEnd <- as.Date("2012/1/1")
  
  dCurrent <- startD
  dNext <- seq(dCurrent, by = "1 year", length = 2)[2]
  
  
  
}
