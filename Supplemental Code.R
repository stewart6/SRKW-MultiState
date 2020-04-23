# Supplemental code for Stewart et al., 2020
# Survival of the Fattest: Using aerial photogrammetry to monitor nutritional stress in killer whales
# This code will run body condition multi-state models for all Southern Resident killer whale pods,
# and will run the full model selection loop with all candidate Chinook salmon covariates.

# Required packages:
library(R2jags)
library(loo)

# Required Data Files:
# SRKW Body Condition Classes.RData
# Chinook Salmon Indices.RData
# SRKW_MultiState.jags

#######################
# 1) LOAD DATA FILES
#######################

# Load the Body Condition Matrices for J, K and L pods 
# (make sure this .RData file is in your working directory)
load("SRKW Body Condition Classes.RData")

# Notes: 2007 is the 'initialization' year. All whales start in condition class 3
# and we disregard the transition probabilities from 2007-2008 and don't include 
# salmon data for 2008 so the covariate fits are not influenced by initialization.
# Body condition classes 1-5, mortality is logged at 6, unmeasured whales logged as NA.
# No measurements taken 2009-2012, 2014, but deaths are included for those years.

# Load the Chinook Salmon abundance indices
load("Chinook Salmon Indices.RData")

# Notes: 
# Salmon abundance data are from the Fishery Regulation Assessment Model (FRAM - refer to Methods)
# Fraser, Columbia, and Puget are aggregate abundances of all Chinook stocks returning to
# the Fraser River, Columbia River, and Puget Sound, respectively. 
# NOF, SWVI, OR, and Salish are regional indices that include all Chinook salmon from ANY stock
# that are present in the North of Cape Falcon (Washington coast), Southwest Vancouver Island,
# Oregon coast, and Salish Sea regions, respectively. 
# Abundance indices were created by dividing the annual abundance by the mean abundance, within each region
# The Chinook_Indices object contains all of the Chinook abundance and abundance index data,
# while each abundance index also has its own object for inclusion in the model loop below.
# The leading NA value in each abundance index is to account for the initialization year (see above)


#######################
# 2) RUN A SINGLE MODEL
#######################

# Choose which pod & salmon index to run
BC <- JBC
Cov <- Fraser

# Format the input data
jags.data <- list(n.occasions = dim(BC)[2], #number of time steps
                  n.ind = dim(BC)[1], #number of animals
                  n.bc = 5, #number of condition classes
                  BC = BC, #body condition data
                  cov = Cov, #covariate data
                  Params = 6) #total number of transition parameters (EG, G, S, G, ED, M)


# Parameters monitored
parameters <- c("EG", "G", "S", "D", "ED", "M", "slope", "intercept", "mean.M", "M.intercept", "log_lik")

# Model setup
nc= 3 #number of chains
ni = 400000 #number of iterations
nb = 200000 #burn-in length
nt = 300 #thinning

# Run the model:
SRKW_Condition <- jags(jags.data, inits=NULL, parameters, "SRKW_MultiState.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Explore results
print(SRKW_Condition)

attach.jags(SRKW_Condition) #plot/inspect transition probabilities, mortality probabilities, covariate slopes, etc.


#######################
# 3) RUN A MODEL LOOP
# FOR MODEL SELECTION
#######################

# List of covariates to run
Stock_List <- c("Fraser","Columbia","Puget","NOF","OR","Salish","SWVI")

# Make a matrix of pod x covariate combinations to loop through,
# with empty cells to fill with information criteria
Combos<-as.matrix(data.frame(Pod=c(rep("J",length(Stock_List)),
                                   rep("K",length(Stock_List)),
                                   rep("L",length(Stock_List))),
                             Stock=c(as.character(rep(Stock_List,3))),
                             LOOIC=NA,
                             LOOIC_SE=NA,
                             WAIC=NA,
                             WAIC_SE=NA,
                             DIC=NA))

# Run model loop
# NOTE: at 400,000 iterations, each model takes several hours to run (strap in)

for(i in 1:(dim(Combos)[1])){
  
  
  # Set up data files for each run:
  
  BC <- get(paste0(Combos[i,1],"BC"))
  Cov <- get(Combos[i,2])
  
  jags.data <- list(n.occasions = dim(BC)[2],
                    n.ind = dim(BC)[1],
                    n.bc = 5,
                    BC = BC,
                    cov = Cov,
                    Params = 6)
  
  
  
  # Parameters monitored
  parameters <- c("EG", "G", "S", "D", "ED","M", "slope","intercept","mean.M","M.intercept","log_lik")
  
  # Assign the JAGS object and run the model
  assign(paste0("SRKW_",Combos[i,1],"Pod_",Combos[i,2]), jags(jags.data, inits=NULL, parameters, "SRKW_MultiState.jags", n.chains = 3, n.thin = 300, n.iter = 400000, n.burnin = 200000, working.directory = getwd()))
  
  # Attach the model output
  attach.jags(get(paste0("SRKW_",Combos[i,1],"Pod_",Combos[i,2]))) #attach the model for LOOIC 
  
  # Calculate LOOIC / WAIC
  ModelLOOIC <- loo(log_lik) #Calculate LOOIC
  ModelWAIC <- waic(log_lik) #Calculate WAIC
  
  # Fill combination matrix with information criteria
  Combos[i,3] <- ModelLOOIC$estimates[3,1] #LOOIC
  Combos[i,4] <- ModelLOOIC$estimates[3,2] #LOOIC SE
  Combos[i,5] <- ModelWAIC$estimates[3,1] #WAIC
  Combos[i,6] <- ModelWAIC$estimates[3,2] #WAIC SE
  Combos[i,7] <- get(paste0("SRKW_",Combos[i,1],"Pod_",Combos[i,2]))$BUGSoutput$DIC #DIC
  
  # Save the workspace as an .RData file so you can inspect the results of each model run as desired
  save.image(paste0("SRKW_",Combos[i,1],"Pod_",Combos[i,2],"_",Sys.Date(),".RData")) #
  
  detach.jags()
  
  # Remove the model output from the workspace after the model selection matrix is populated
  # and the model file is saved. Otherwise you'll have a humungous workspace object to save
  # after each model run, and it will get bigger and bigger with duplicates as the loop goes on
  rm(list=paste0("SRKW_",Combos[i,1],"Pod_",Combos[i,2]))
  
  # Write out the model selection table as a .csv so it's updated and saved separately after each model run
  write.csv(Combos,"SRKW Loop LOOICs.csv",row.names = F)
}


