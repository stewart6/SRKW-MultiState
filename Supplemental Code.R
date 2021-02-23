# Supplemental code for Stewart et al., 2020
# Survival of the Fattest: Linking body condition to prey availability and survivorship of killer whales
# This code will run body condition multi-state models for Southern Resident killer whale L Pod.

# Required packages:
library(R2jags)

# Required Data Files:
# SRKW Body Condition Classes.RData
# Chinook Salmon Indices.RData
# SRKW_MultiState.jags

#######################
# 1) LOAD DATA FILES
#######################

# Load the Body Condition Matrices for L pod
# (make sure this .RData file is in your working directory)
load("SRKW Body Condition Classes.RData")

# Notes: 2007 is the 'initialization' year. All whales start in condition class 3
# and we disregard the transition probabilities from 2007-2008 and don't include 
# salmon data for 2008 so the covariate fits are not influenced by initialization.
# Body condition classes 1-5, mortality is logged at 6, unmeasured whales logged as NA.
# No measurements taken 2009-2012, 2014, but known deaths are included for all years 2008-2019.


# Load the Chinook Salmon abundance indices
load("Chinook Salmon Indices.RData")

# Notes: 
# Salmon abundance data are from the Fishery Regulation Assessment Model (FRAM - refer to Methods)
# Fraser, Columbia, and Puget are aggregate abundances of all Chinook stocks returning to
# the Fraser River, Columbia River, and Puget Sound, respectively. 
# NOF, SWVI, OR, and Salish are regional indices that include all Chinook salmon from ANY stock
# that are present in the North of Cape Falcon (Washington coast), Southwest Vancouver Island,
# Oregon coast, and Salish Sea regions, respectively. 
# Named vectors are Z-scored abundance indices, created by subtracting mean abundance from annual abundance and dividing by the standard deviation, within each region
# The Chinook_Indices object contains all of the Chinook abundance and Z-scored data,
# while each abundance index also has its own named vector object for inclusion in the model loop below.
# The leading NA value in each named abundance vector is to account for the initialization year (see above)



# Load the Age Matrices and Sex matrices for L pod:
load("SRKW AgeClasses.RData")

# Notes:
# 1 = Not yet born (in all models, mortality probability is forced to 0 for age/sex class 1)
# 2 = Calf
# 3 = Juvenile
# 4 = Young Female
# 5 = Old Female
# 6 = Young Male
# 7 = Old Male




#############################
# 2) RUN THE NULL MODEL
#############################

# Choose which pod to run
BC <- as.matrix(LBC)
AgeSex <- as.matrix(LAgeClass)


jags.data <- list(n.occasions = dim(BC)[2], #number of time steps
                  n.ind = dim(BC)[1], #number of animals
                  n.bc = 5, #number of condition classes
                  BC = BC, #body condition data
                  Params = 3, #total number of transition parameters (G, S, D)
                  AgeSex = AgeSex) #vector of which whales are male (1 / 0)

parameters <- c( "G", "S", "D", "M", "Base.M", "Mean.M")

# Set some initial values for mortality probability
inits <- function(){
  Base.M = runif(7,-10,-5)
  M = runif(5,-10,-5)
  list(Base.M=Base.M,
       M=M)
}
nc= 3 #number of chains
ni = 100000 #number of iterations
nb = 50000 #burn-in length
nt = 50 #thinning

SRKW_LPod_Null <- jags(jags.data, inits=inits, parameters, "SRKW_MultiState_Null.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

attach.jags(SRKW_LPod_Null)


#############################
# 3) RUN THE TIME-ONLY MODEL
#############################

# Choose which pod to run
BC <- as.matrix(LBC)
AgeSex <- as.matrix(LAgeClass)


jags.data <- list(n.occasions = dim(BC)[2], #number of time steps
                  n.ind = dim(BC)[1], #number of animals
                  n.bc = 5, #number of condition classes
                  BC = BC, #body condition data
                  Params = 3, #total number of transition parameters (G, S, D)
                  AgeSex = AgeSex) #vector of which whales are male (1 / 0)

parameters <- c( "G", "S", "D", "M", "Base.M", "Mean.M")

# Set some initial values for mortality probability
inits <- function(){
  Base.M = runif(7,-10,-5)
  M = runif(5,-10,-5)
  list(Base.M=Base.M,
       M=M)
}
nc= 3 #number of chains
ni = 100000 #number of iterations
nb = 50000 #burn-in length
nt = 50 #thinning

SRKW_LPod_Time <- jags(jags.data, inits=inits, parameters, "SRKW_MultiState_Time.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

attach.jags(SRKW_LPod_Time)



#############################
# 4) RUN THE COVARIATE MODEL
#############################

BC <- as.matrix(LBC)
Cov <- Puget
AgeSex <- LAgeClass



jags.data <- list(n.occasions = dim(BC)[2], #number of time steps
                  n.ind = dim(BC)[1], #number of animals
                  n.bc = 5, #number of condition classes
                  BC = BC, #body condition data
                  cov = Cov, #covariate data
                  Params = 3, #total number of transition parameters (EG, G, S, G, ED, M)
                  AgeSex = AgeSex) #vector of which whales are male (1 / 0)

parameters <- c( "G", "S", "D", "M", "slope", "intercept", "Base.M", "Mean.M")

# Set some initial values for mortality probability
inits <- function(){
  Base.M = runif(7,-10,-5)
  M = runif(5,-10,-5)
  list(Base.M=Base.M,
       M=M)
}

nc= 3 #number of chains
ni = 100000 #number of iterations
nb = 50000 #burn-in length
nt = 50 #thinning

SRKW_LPod_Puget <- jags(jags.data, inits=inits, parameters, "SRKW_MultiState_Cov.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

attach.jags(SRKW_LPod_Puget)


