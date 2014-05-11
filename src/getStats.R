#### Preamble ####

# PROGRAM: Get capture recapture stats
# CODER: Cameron Roach
# DATE: 5/3/2014
# DESCRIPTION: Functions to calculate capture-recapture statistics from capture
# matrix.


################################## Functions ##################################

checkCapt <- function(mtrxCapt) {
  # calculates how many animals have been captured in capture matrix
  
  mtrxCapt <- as.matrix(mtrxCapt)
  
  nCaptures <- apply(mtrxCapt,1,sum)
  
  m <- length(nCaptures[nCaptures>=1])
  
  return(m)
}



calcZ <- function(curNotCapt, i) {
  # Takes capture matrix as input where no animals caught at time i. Returns
  # number of members of marked population not captured at sampling occasion i
  # that are captured again later.
  
  k  <- ncol(curNotCapt)
  
  if (i<k & i>1) {
        
    preCapt <- as.matrix(curNotCapt[,1:(i-1)])
    preCapt <- apply(preCapt,1,sum)>=1
    
    postCapt <- as.matrix(curNotCapt[,(i+1):k])
    postCapt <- as.matrix(postCapt[preCapt,])
    
    Z <- sum(apply(postCapt,1,sum)>=1)
  } else {
    Z <- NA
  }
  
  return(Z)
}



mse.f <- function(x, x.act) {
  sum((x - x.act) ^ 2) / length(x)
}



calcCoverage <- function(actN, bs.ci) {
  # calculate coverage probability
  nCaptSims <- length(bs.ci)
  
  coverage <- NULL
  for (i in 1:nCaptSims) {
    x <- bs.ci[[i]]
    tmp <- actN<=x$bs.ci.u & actN>=x$bs.ci.l
    coverage.tmp <- data.frame("withinCi"=tmp,"Period"=1:20)
    coverage <- rbind(coverage,coverage.tmp)
  }
  
  
  coverage.prob <- aggregate(withinCi~Period, coverage, sum)
  coverage.prob$probWithinCi <- coverage.prob$withinCi/nCaptSims
  
  return(coverage.prob)
}