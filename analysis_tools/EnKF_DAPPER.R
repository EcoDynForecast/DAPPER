#rm(list = ls())
if (!"mvtnorm" %in% installed.packages()) install.packages("mvtnorm")
if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4")
if (!"glmtools" %in% installed.packages()) install.packages('glmtools', repos=c('http://cran.rstudio.com', 'http://owi.usgs.gov/R'))
library(mvtnorm)
library(glmtools)
library(ncdf4)

#### set working directory
pathDAPPER = '/Users/laurapuckett/Documents/Research/GIT/DAPPER/' 
EnKF_directory = '/Users/laurapuckett/Documents/Spring_2018/Research_Spring_2018/'
setwd(EnKF_directory)

#### define variables needed for DAPPER functions
working_directory = pathDAPPER # this is references in prepare_update_states(), and prepare_obs()
input_directory = '/Users/laurapuckett/Documents/Research/GIT/DAPPER_inputdata/'
output_directory = '/Users/laurapuckett/Documents/Research/GIT/DAPPER_projects/'
plotlist = c(40001) #THIS IS Duke Forest
focal_plotID = plotlist#

#### source required functions
source('prepare_update_states.R')
source('update_states.R')

#---------------------------------------------------------
###### MAIN EnKF CODE ######

#See table 1 in Rastetter et al. 2010
nmembers = 50
NO_UNCERT = FALSE
ADD_NOISE_TO_OBS = FALSE
USE_SYNTHETIC_DATA = FALSE
USE_OBS_COORD = FALSE
USE_OBS_CONTRAINT = TRUE


start_forecast_step = 1


# Define variables needed for prepare_update_states or update_states
plotnum = 1
all_studies = c('/Duke/TIER4_Duke')
nsteps = 1


# Retreieve variables from prepare_update_states() needed for update_states()
prep_update_states <- prepare_update_states(plotnum)
output_dim = prep_update_states$output_dim
pars = prep_update_states$pars
site = prep_update_states$site
nopars = prep_update_states$nopars
nosite = prep_update_states$nosite
met = prep_update_states$met
thin_event = prep_update_states$thin_event
exclude_hardwoods = prep_update_states$exclude_hardwoods
median_pars = prep_update_states$median_pars
ASW_max = prep_update_states$ASW_max


#USE FIRST OBSERVATION AS THE INITIAL CONDITIONS
if(!USE_SYNTHETIC_DATA){
  realdata = read.csv(file = "/Users/laurapuckett/Documents/Spring_2018/Research_Spring_2018/fake_LAI_data.csv", header = TRUE, sep = ",")
}else{
  realdata = read.csv(file = "/Users/laurapuckett/Documents/Spring_2018/Research_Spring_2018/fake_LAI_data.csv", header = TRUE, sep = ",")
}
observedLAI = realdata[1,2] # (one obs for now) change this when there is actual data


#NUMBER OF STATE SIMULATED = SPECIFIED DEPTHS
nstates <- 6 #nlayers_init
# LAI, WS, WCR, WFR, Stem_Count, ASW

#Observations for each observed state at each time step
#an observation with at least 1 observation but without an observation in a time-step gets assigned an NA
z <- t(matrix(rep(NA,length(observedLAI)), nrow = length(observedLAI), ncol = nsteps))

for(i in 1:nsteps){
  #index = which(as.POSIXct(realdata$datetime) == full_time[i])
  #if(length(index) > 0){
    index = 1
    z[i,] = unlist(realdata[1,2]) # fix - reformat later for correct file format
    #if(ADD_NOISE_TO_OBS & i > 1){
    #  z[i,] = rnorm(length( z[i,]), z[i,],0.1)
    #}
  #}
}

if(!USE_OBS_CONTRAINT){
  z_obs = z
  z[,] = NA
}

#FIGURE OUT WHICH LOCATIONS HAVE OBSERVATIONS
obs_index = rep(NA,length(observedLAI))
#for(i in 1:length(observedLAI)){
 # obs_index[i] = which.min(abs(the_depths_init - observedLAI[i])) # how does this work?
#}

#A matrix for knowing which state the observation corresponds to
z_states <- t(matrix(obs_index, nrow = length(obs_index), ncol = nsteps))

#Process error 
# fix - this needs a lot of work. Supposed to be sigma = gamma + rho*value

# Define variables for Qt matrix
# need to incorporate rho into sigma later
sigma = array(data = NA, dim = nstates)
sigma[1] = median_pars[52] # gamma_LAI
sigma[2] = median_pars[53] +  median_pars[64] # WS, gamma_WS+ rho_WS * last stem measurement
sigma[3] = median_pars[54] # WCR
sigma[4] = median_pars[55] ### check this one ## WFR
sigma[5] = .1 # gamma_Stem_Count - look into this, should have a gamma
sigma[6] = .1 # gamma_ASW

rho_LAI = 0 #median_pars[63]
rho_WS = 0  #median_pars[64]
rho_WCR =0 # median_pars[65]
rho_WFR = 0 #median_pars[66]
rho_Stem_Count = 0
rho_ASW = 0

Qt = diag(sigma) # fix this part


#Measurement error 
psi = rep(0.0001,length(obs_index)) # update this later

#Initial states
init_LAI = 1 # WFi # change these to starting values for Duke forest (at age 11ish) 
init_WS = 2 # WSi
init_WCR = 0.4 # WCRi
init_WR = 0.16 # WRi
init_Stem_Count = 1200 #"StemNoi"
init_ASW = ASW_max 
init_states <- array(data = c(init_LAI, init_WS, init_WCR, init_WR, init_Stem_Count, init_ASW), dim = c(nstates,1))
x <- array(NA,dim=c(nsteps,nmembers,nstates))
x[1,,] <- rmvnorm(n=nmembers, mean=init_states, sigma=as.matrix(Qt)) 


#Matrix to store ensemble specific deviations and innovations
dit = array(NA,dim=c(nmembers,nstates))
#dit_star = array(NA,dim=c(nmembers,nstates)) #Adaptive noise estimation



#loop through time steps
for(i in 2:nsteps){ # runs 3PG for one time step for every member  (would use forecasting code)
  
  #i is the same as mo - make them the same, climate matrices should be built before mo 
  # need to update times in addition to states

  site[3] = InitialYear
  site[4] = InitialMonth
  site[5] = start_age
  #Create array to hold DAPPER predictions for each ensemble
  x_star = array(NA, dim = c(nmembers,nstates))
  for(m in 1:nmembers){
    
    ### updates the states in site array each timestep 
    site[6] = x[i-1, m, 1] #"WFi" # can take LAI or WFi
    site[7] = x[i-1, m, 4] # WFR - is this correct?
    site[8] = x[i-1, m, 2] #"WSi"
    site[9] = x[i-1, m, 5] #"StemNoi"
    site[10] = x[i-1, m, 6] #"ASWi"
    site[20] = x[i-1, m, 3] # WCRi
    site[26] = x[i-1, m, 1] # same as WFi above
    
    DAPPER_states <- update_states(mo, plotnum, output_dim, pars, site, nopars, nosite, thin_event, met)
    
    #4) Fill x_star with temperatures from GLM
    x_star[m,] = DAPPER_states # all of the updated states for the current member iteration
  }
  
  InitialMonth = InitialMonth + 1
  if(InitialMonth > 12){
    InitialMonth = 1
    InitialYear = InitialYear + 1
  }
  start_age = start_age + 1/12
  
  #Corruption [nmembers x nstates] 
  NQt = rmvnorm(n=nmembers, sigma=as.matrix(Qt)) # matrix of diagaonal - error term (noise step in forecasting code)
  
  #Matrix Corrupted state estimate [nmembers x nstates]
  x_corr = x_star + NQt
  ## builing ^ starting point  
  
  
  #Obs for time step
  z_index = which(!is.na(z[i, ]))
  
  #if no observations at a time step then just propogate model uncertainity
  if(length(z_index) == 0 | i > start_forecast_step){
    x[i,,] = x_corr
    if(NO_UNCERT){
      x[i,,] = x_star
    }
    
  }else{
    print(i)
    #if observation then calucate Kalman adjustment
    zt = z[i,z_index] # will be time series of LAI from landsat
    z_states_t = z_states[i,z_index]
    yit = array(NA,dim=c(nmembers,length(zt)))
    
    #Assign which states have obs in the time step
    H <- array(0,dim=c(length(zt),nstates))
    for(j in 1:length(z_index)){
      H[j,z_states_t[j]] = 1
    }
    
    #Extract the data uncertainity for the data types present during the time-step
    if(length(z_index)>1){
      psi_t = diag(psi[z_index])
    }else{
      #Special case where there is only one data type during the time-step
      psi_t = psi[z_index]
    }
    
    #Ensemble mean
    ens_mean = apply(x_corr, 2, mean)
    #ens_mean_star = apply(x_star, 2, mean) #Adaptive noise estimation
    
    #Loop through ensemble members
    for(m in 1:nmembers){  
      
      #Ensemble specific deviation
      dit[m,] = x_corr[m,]-ens_mean
      #dit_star[m,] x_star[m,] - ens_mean_star #Adaptive noise estimation
      
      #Observational uncertainity
      N_psi = rmvnorm(n=1,sigma=as.matrix(psi_t))
      
      #Ensemble specific innovations
      yit[m,] = t(zt - crossprod(t(H),(x_corr[m,])) + t(N_psi))
      
      #Ensemble specific estimate and innovation covariance
      if(m == 1){
        Pit = tcrossprod(dit[m,])
        #Pit_star = tcrossprod(dit_star[m,])  #Adaptive noise estimation
        Sit = tcrossprod(yit[m,])
      }else{
        Pit = tcrossprod(dit[m,]) +  Pit
        #Pit_star = tcrossprod(dit_star[m,]) + Pit_star #Adaptive noise estimation
        Sit = tcrossprod(yit[m,]) +  Sit
      }
    }
    
    #estimate covariance
    Pt = Pit/nmembers
    #Pt_star = Pit_star/nmembers  #Adaptive noise estimation
    #Innovations covariance
    St = Sit/nmembers
    
    #Kalman gain
    Kt <- crossprod(t(crossprod(Pt,t(H))), solve(St, tol=1e-30))
    
    #Adaptive noise estimation
    beta = 0.55
    #Gammat = (1 - beta)*crossprod(t(crossprod(Pt,t(H))), solve(H, tol=1e-30))
    #Qt_hat = Gammat(St - H*Pt_star*H - t(N_psi))*Gammat
    #Q = alpha*Q+ (1-alpha)*Qt_hat
    #Ensemble specific updated state
    for(m in 1:nmembers){
      x[i,m,] = t((x_corr[m,]) + crossprod(t(Kt), yit[m,]))
      if(NO_UNCERT){
        x[i,m,] = x_star[m,]
      }
    }
    if(length(which(is.na(x[i,,]))) > 0){dies = i}
  }
}

#par(mfrow=c(2,3),mar = c(4,4,2,2),oma = c(3,3,2,2))
#hist(x[nsteps,,1],main='state 1',xlab='value')
#hist(x[nsteps,,2],main='state 2',xlab='value')
#hist(x[nsteps,,3],main='state 3',xlab='value')

par(mfrow=c(4,3))

if(!USE_OBS_CONTRAINT){
  z = z_obs
}

#c(1,4,7,10,13,16,19,22,25,28,29)
for(i in c(1,4,7,10,13,16,19,22,25,28,29)){
  #for(i in 1:nlayers_init){
  model = i
  if(length(which(z_states[1,] == i) > 0)){
    obs = which(z_states[1,] == i)
  }else{
    obs = NA
  }
  #ylim = range(c(x[,,model],c(z[,obs])),na.rm = TRUE)
  ylim = range(c(x[,,],c(z[,])),na.rm = TRUE)
  plot(x[,1,model],type='l',ylab='surface temperature (celsius)',xlab='time step (day)',main = the_depths_init[i],ylim=ylim)
  for(m in 2:nmembers){
    points(x[,m,model],type='l')
  }
  
  if(!is.na(obs)){
    tmp = z[,obs]
    tmp[is.na(tmp)] = -999
    points(tmp,col='red',pch=19,cex=2.0)
  }
}




