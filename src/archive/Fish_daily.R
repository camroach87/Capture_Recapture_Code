############################### PREAMBLE ################################

# Make sure the correct directory is chosen
if (Sys.info() ["sysname"] == "Darwin") {
  # MAC
  setwd("/Volumes/Data/Users/chroach/Dropbox/Uni/Research/")
} else {
  # PC
  setwd("C:/Users/Cameron/Dropbox/Uni/Research");
}

data <- read.csv("TC_Data_Charles.csv");

D <- as.Date(data$surveydate, format="%d/%m/%Y")
minDate <- min(D)
tabD <- table(D)

# add day to data
data[10:11, "day"] <- NA
data$day <- as.numeric(D-minDate)

# useful
tabCapt <- table(data$idfish, data$day)




days <- sort(unique(data$day))
minDay <- min(days)
maxDay <- max(days)
curDay <- minDay
k <- maxDay - minDay + 1


#fishDay <- list(days)
fishDay <- list(NA)

for (j in 1:k) {
  fishDay[[j]] <- data$idfish[data$day == j + minDay - 1]
}




################################### MAIN ##################################

JollySeber <- function(inData) {
  #Jolly Seber
  n <- CalcTotal(inData)
  cat("n complete \n")
  m <- CalcMarked(inData)
  cat("m complete \n")
  R <- CalcReleased(inData)
  cat("R complete \n")
  Z <- CalcZ(inData)
  cat("Z complete \n")
  r <- CalcRecapt(inData)
  cat("r complete \n")
  
  M <- m + R*Z/r
  N <- n*M[1:length(n)]/m
  
  plot(N ~ days, type="l", col="green")
}
  

CJS <- function() {
  
  
  
}





############################### FUNCTIONS ###############################





# Calculates number of animals captured exactly j times, j=1,...,k
CalcExact <- function(data) {
  uniqueFish <- unique(data$idfish)
  counts <- c()
  f <- c()
  days <- sort(unique(data$day))
  k <- max(days)-min(days) + 1
  
  #First calculate number of times each fish caught
  
  for (i in 1:length(uniqueFish)) counts[i] <- sum(uniqueFish[i] == data$idfish)
  
  #Then calculate number of animals caught exactly j times
  for (j in 1:k) f[j] <- sum(counts == j)
  
  return(f)
  
}



CalcTotal <- function(data) {
  
  captFish <- NA
  n <- NA
  M <- NA
  
  
  
  for (j in 1:k) {
    n[j] <- length(data$idfish[data$day == j + minDay - 1])
  }
  
  return(n)
}




# Calculates number of distinct animals caught prior to jth capture occasion.
# By definition M_0 = 0 and M_k+1 = total distinct animals captured
# This assumes a closed population
CalcDistinct <- function(data) {
  
  u <- CalcUnmarked(data)
  
  cat(u, "\n")
  
  M <- c()
  M[1] <- 0
  
  for (j in 1:length(u)) M[j+1] <- sum(u[1:j])
    
  return(M)
  
}






# calculates number of unmarked fish, u_j
CalcUnmarked <- function(data) {
  u <- rep(0,k)
  markedFish <- c()
  

  for (j in 1:k) {
    
    # calculate number of unmarked animals for capture occasion j
    if (j != 1) {
      
      for (i in 1:length(fishDay[[j]])) {
        
        # checks if fish i in day j was captured in previous days (or already caught this day)
        
        if (fishDay[[j]][i] %in% markedFish) {
          fishNotCaught = 0
          #cat(fishDay[[j]][i],"\n")
          
        } else {
          fishNotCaught = 1
        }
        
        markedFish <- c(markedFish, fishDay[[j]][i])
        
        if (fishNotCaught == 1) u[j] = u[j] + 1
        
      }
      
    } else {
      u[j] = length(fishDay[[j]])
      markedFish <- c(fishDay[[j]])
    }
    
    
  }
  
  return(u)
  
}





# calculates number of marked fish, m_j
CalcMarked <- function(data) {
  m <- rep(0,k)
  markedFish <- c(fishDay[[1]])
  
  
  
  # start at j=2 because m_1=0
  for (j in 2:k) {
    # DEBUG cat("\n j = ", j)
    
    # checks if any fish were caught on day j
    if (length(fishDay[[j]]) != 0) {
      # calculate number of unmarked animals for capture occasion j
      for (i in 1:length(fishDay[[j]])) {
        # DEBUG cat(" i = ", i)
        # checks if fish i in day j was captured in previous days (or already captured this day)
        if (fishDay[[j]][i] %in% markedFish) {
          m[j] = m[j] + 1
        }
        
        markedFish <- c(markedFish, fishDay[[j]][i])    
      }
    }  
  }
  
  return(m)
  
}







# calculates R_j, total number of animals capture at sampling occasion j that are released
# Note: if this changes see note for CalcRecapt() as it will probably need to change as well (see definition in book)
CalcReleased <- function(data) {
  # Until I hear otherwise, I will assume that all animals that are captured are relased
  R <- CalcTotal(data)
  
  return(R)  
}




# Calculates r_j, the number of members of R_j captured again later
# Note: doesn't include fish from sampling occasion j captured again during sampling occasion j
# Note: if R_j is not simply all animals captured, then this function will need to change.
CalcRecapt <- function(data) {
  r <- c()
  
  
  # Check day j for later captures in days i
  for (j in 1:(k-1)) {
    
    # refreshes recaptured matrix each day, j
    recaptured <- matrix(0)
    
    for (i in (j+1):k) {
      recaptured <- recaptured + unique(fishDay[[j]]) %in% fishDay[[i]]
    }
    
    # don't care if fish recaptured twice or more - only need to know they were captured
    recaptured[recaptured>=2] <- 1
    
    r[j] <- sum(recaptured)
    
  }
  
  r[k] <- 0
  
  return(r)
  
}







# Calculates z_j, the number of members of the marked population not captured at j (Mj-mj) captured again later
CalcZ <- function(data) {
  z <- c()
  
  # Check day j for later captures in days i
  for (j in 1:(k-1)) {
    
    # refreshes markedPop and recaptured matrix each day, j
    markedPop <- c(NA)
    recaptured <- matrix(0)
    
    cat("\n j = ", j)
    
    # gets marked population excluding those fish captured in current day j
    if (j == 1) {
      recaptured <- 0
    } else {
      for (i in 1:(j-1)) {
        markedPop <- c(markedPop, fishDay[[i]][!(fishDay[[i]] %in% fishDay[[j]])])
        markedPop <- unique(markedPop)
        
        # DEBUG cat(" i=", i)
      }
      
      for (i in (j+1):k) {
        recaptured <- recaptured + markedPop %in% fishDay[[i]]
        # DEBUG cat(" i=", i)
      }
      
      # don't care if fish recaptured twice or more - only need to know they were captured
      recaptured[recaptured>=2] <- 1
    }
    
    z[j] <- sum(recaptured)
    
  }
  
  z[k] <- 0
  
  return(z)
  
}

# Calculates z_j, the number of members of the marked population not captured at j (Mj-mj) captured again later
CalcZ_fast <- function(data) {
  
  t1 <- Sys.time()
  
  z <- rep(NA, k)     # CHANGE c()
  
  # Check day j for later captures in days i
  for (j in 1:300) {        #(k-1)) {
    
    # refreshes markedPop and recaptured matrix each day, j
    markedPop <- rep(NA,5000)    # CHANGE can't remember - check above function
    recaptured <- 0   # CHANGE matrix(0)
    
    # DEBUG cat("\n j = ", j)
    
    # gets marked population excluding those fish captured in current day j
    if (j == 1) {
      recaptured <- 0
    } else {
      for (i in 1:(j-1)) {
        markedPop <- c(markedPop, fishDay[[i]][!(fishDay[[i]] %in% fishDay[[j]])])
        markedPop <- unique(markedPop)
        
        # DEBUG cat(" i=", i)
      }
      
      for (i in (j+1):k) {
        recaptured <- recaptured + markedPop %in% fishDay[[i]]
        # DEBUG cat(" i=", i)
      }
      
      # don't care if fish recaptured twice or more - only need to know they were captured
      recaptured[recaptured>=2] <- 1
    }
    
    z[j] <- sum(recaptured)
    
  }
  
  z[k] <- 0
  
  (T1 <- Sys.time() - t1)
  
  return(z)
  
}