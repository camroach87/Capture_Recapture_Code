############################### PREAMBLE ################################

data <- read.csv("TC_Data_Charles.csv");

D <- as.Date(data$surveydate, format="%d/%m/%Y")
tabD <- table(D)

# add year to data
data[10:11, "year"] <- NA
data$year <- as.numeric(format(D, "%Y"))

# useful
tabCapt <- table(data$idfish, data$year)




years <- sort(unique(data$year))
minYear <- min(years)
maxYear <- max(years)
curYear <- minYear
k <- maxYear - minYear + 1


#fishYear <- list(years)
fishYear <- list(NA)

for (j in 1:k) {
  fishYear[[j]] <- data$idfish[data$year == j + minYear - 1]
}




################################### MAIN ##################################

JollySeber <- function(inData) {
  #Jolly Seber. inData is just "data" variable as defined above.
  n <- CalcTotal(inData)
  m <- CalcMarked(inData)
  R <- CalcReleased(inData)
  Z <- CalcZ(inData)
  r <- CalcRecapt(inData)
  
  M <- m + R*Z/r
  N <- n*M[1:14]/m
  
  plot(N ~ years, type="l", col="green")#, ylim=c(0,30000))
}
  

CJS <- function() {
  
  
  
}





############################### FUNCTIONS ###############################





# Calculates number of animals captured exactly j times, j=1,...,k
CalcExact <- function(data) {
  uniqueFish <- unique(data$idfish)
  counts <- c()
  f <- c()
  years <- sort(unique(data$year))
  k <- max(years)-min(years) + 1
  
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
    n[j] <- length(data$idfish[data$year == j + minYear - 1])
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
      
      for (i in 1:length(fishYear[[j]])) {
        
        # checks if fish i in year j was captured in previous years (or already caught this year)
        
        if (fishYear[[j]][i] %in% markedFish) {
          fishNotCaught = 0
          #cat(fishYear[[j]][i],"\n")
          
        } else {
          fishNotCaught = 1
        }
        
        markedFish <- c(markedFish, fishYear[[j]][i])
        
        if (fishNotCaught == 1) u[j] = u[j] + 1
        
      }
      
    } else {
      u[j] = length(fishYear[[j]])
      markedFish <- c(fishYear[[j]])
    }
    
    
  }
  
  return(u)
  
}





# calculates number of marked fish, m_j
CalcMarked <- function(data) {
  m <- rep(0,k)
  markedFish <- c(fishYear[[1]])
  
  
  
  # start at j=2 because m_1=0
  for (j in 2:k) {
    
    # calculate number of unmarked animals for capture occasion j
    for (i in 1:length(fishYear[[j]])) {
      
      # checks if fish i in year j was captured in previous years (or already captured this year)
      if (fishYear[[j]][i] %in% markedFish) {
        m[j] = m[j] + 1
      }
      
      markedFish <- c(markedFish, fishYear[[j]][i])    
          
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
  
  
  # Check year j for later captures in years i
  for (j in 1:(k-1)) {
    
    # refreshes recaptured matrix each year, j
    recaptured <- matrix(0)
    
    for (i in (j+1):k) {
      recaptured <- recaptured + unique(fishYear[[j]]) %in% fishYear[[i]]
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
  
  # Check year j for later captures in years i
  for (j in 1:(k-1)) {
    
    # refreshes markedPop and recaptured matrix each year, j
    markedPop <- c(NA)
    recaptured <- matrix(0)
    
    # gets marked population excluding those fish captured in current year j
    if (j == 1) {
      recaptured <- 0
    } else {
      for (i in 1:(j-1)) {
        markedPop <- c(markedPop, fishYear[[i]][!(fishYear[[i]] %in% fishYear[[j]])])
        markedPop <- unique(markedPop)
      }
      
      for (i in (j+1):k) {
        recaptured <- recaptured + markedPop %in% fishYear[[i]]
      }
      
      # don't care if fish recaptured twice or more - only need to know they were captured
      recaptured[recaptured>=2] <- 1
    }
    
    z[j] <- sum(recaptured)
    
  }
  
  z[k] <- 0
  
  return(z)
  
}