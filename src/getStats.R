#### Preamble ####

# PROGRAM: Get capture recapture stats
# CODER: Cameron Roach
# DATE: 5/3/2014
# DESCRIPTION: Functions to calculate capture-recapture statistics from capture
# matrix.


################################## Functions ##################################

calcMarked <- function(mtrxCapt) {
  # Checks if animal has been captured more than once. Only can be used to
  # determine if animal is captured prior to the latest sampling occasion.
  nCaptures <- apply(mtrxCapt,1,sum)
  
  m <- length(nCaptures[nCaptures>1])
  
  return(m)
}