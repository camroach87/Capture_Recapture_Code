#### Preamble ####

# CODER: Cameron Roach
# DATE: 2/4/2013
# DESCRIPTION: Simulates what is actually happening in closed and open
# populations and then compares estimators for each. Assumes sampling occasions
# are equally spaced.
#
# NOTE: Currently, does not take into consideration the marked population size
# approaching the total population size. Can cause NaN's to be produced once
# marked population exceeds total population.


# Constants
N_1 <-1000
t <- 20
birthRate <- 10   
deathRate <-5
captureRate <- 30   #rate of captures

# N_1 <- as.numeric(readline("Initial population? "))
# t <- as.numeric(readline("Number of sampling occasions? "))
# birthRate <- as.numeric(readline("Rate of migrations (births/deaths)? "))
# captureRate <- as.numeric(readline("Rate of captures? "))

# Initialise variables
#N <- rep(NA, t)
n <- rep(NA, t)
m <- rep(NA, t)
M <- rep(NA, t)


##### Closed population ####

# Calculate population
# OPEN: No change (migrations will only affect rate - so lambda1 might change based on migration rate)
N <- N_1



# Calculate number of animals captured at each sampling occasion
# OPEN: some of these will die
n <- rpois(t,captureRate)


# Calculate number of captured animals that are marked for each sampling occasion
m[1] <- n[1]
M[1] <- 0
for (i in 2:t) {
  
  #ERROR HERE  
    M[i] <- M[i-1] + (n[i-1] - m[i-1])
  
    p <- M[i]/N
  m[i] <- rbinom(1,n[i],p)
}
cat("n ", n,"\n M ",M,"\n m ",m,"\n")


# This is the estimate that would have been derived based on the capture recapture data
N_est <- M*n/m
cat(N_est)

# Chao's estimator for sparse capture recapture data
# Geometric distribution? Calculate p based on number of marked fish?
pGeom <- N/(N+sum(n[1:t-1]))
NoTimesFishCapt <- rgeom(N,pGeom)

# calc f1 and f2 then do estimator













#### Open population ####

# Calculate population
# N[1] <- N_1
# for (i in 2:t) {
#   N[i] <- N[i-1] + rpois(1, birthRate) - rpois(1,deathRate)
# }

