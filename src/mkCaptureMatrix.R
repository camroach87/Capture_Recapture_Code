############################### PREAMBLE ################################
# PROGRAM: Make capture matrices
# CODER: Cameron Roach
# DATE: 9/7/2013
# DESCRIPTION: Constructs various capture matrices using either a datafile of
# capture occasions or by simulation.






#### Functions ####

mkCaptureMtrx_years <- function(datafile = "TC_Data_Charles.csv") {
  # Constructs a capture matrix using sourcefile in format of TC_Data_Charles.csv.
  data <- read.csv(file.path(dataDir,datafile));
  
  D <- as.Date(data$surveydate, format="%d/%m/%Y")
  tabD <- table(D)
  
  # add year to data
  data[10:11, "year"] <- NA
  data$year <- as.numeric(format(D, "%Y"))
  
  
  # useful
  mtrxCapt <- table(data$idfish, data$year)
  
  return(mtrxCapt)
}





mkCaptureMtrx <- function(datafile = "TC_Data_Charles.csv") {
  # Constructs a capture matrix using sourcefile in format of TC_Data_Charles.csv.
  data <- read.csv(file.path(dataDir,datafile));
  
  D <- as.Date(data$surveydate, format="%d/%m/%Y")
  tabD <- table(D)
  
  # add day to data
  data[10:11, "day"] <- NA
  data$day <- as.numeric(D)
  
  
  # useful
  mtrxCapt <- table(data$idfish, data$day)
  
  return(mtrxCapt)
}





mkCloseSimMtrx <- function(N.0, t, pCapture) {
  # Simulates closed population. N.0 is initial population size.
  # Assumes capture probabilities are uniformly distributed.
  mtrxCaptProb <- matrix(runif(N.0*t),N.0,t)
  mtrxCapt <- ifelse(mtrxCaptProb<=pCapture,1,0)
  
  return(mtrxCapt)  
}
  



mkOpenSimMtrx <- function(Pop, t, p, beta) {  

  mtrxP <- Pop*matrix(runif(nrow(Pop)*t),nrow(Pop),t)
  
  mtrxCapt <- ifelse(mtrxP>=(1-p),1,0) 
  
  return(list(mtrxCapt, colSums(Pop)))
  
}



