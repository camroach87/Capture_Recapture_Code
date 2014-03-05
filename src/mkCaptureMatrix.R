############################### PREAMBLE ################################
# PROGRAM: Make capture matrices
# CODER: Cameron Roach
# DATE: 9/7/2013
# DESCRIPTION: Constructs various capture matrices using either a datafile of
# capture occasions or by simulation.




################################## Functions ##################################

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
  
  return(mtrxCapt)
  
}