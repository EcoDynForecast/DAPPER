if (!"mvtnorm" %in% installed.packages()) install.packages("mvtnorm")
if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4")
if (!"glmtools" %in% installed.packages()) install.packages('glmtools', repos=c('http://cran.rstudio.com', 'http://owi.usgs.gov/R'))
library(mvtnorm)
library(glmtools)
library(ncdf4)


#### source required functions
source(prepare_update_states.R)
source(update_states.R)

#---------------------------------------------------------
###### MAIN EnKF CODE ######

#See table 1 in Rastetter et al. 2010
nmembers = 50
NO_UNCERT = FALSE
ADD_NOISE_TO_OBS = FALSE
USE_SYNTHETIC_DATA = FALSE
USE_OBS_COORD = FALSE
USE_OBS_CONTRAINT = TRUE

nsteps = 1 # nmonths
start_forecast_step = 1

#calibFolder = '/Users/quinn/Dropbox (VTFRS)/Research/SSC_forecasting/FCR-GLM/SCC/Calibration' # fix
pathDAPPER = '/Users/laurapuckett/Documents/Research/GIT/DAPPER/' 
EnKF_directory = '/Users/laurapuckett/Documents/Spring_2018/Research_Spring_2018/'
setwd(EnKF_directory)

# GET Prepare_update_states info
prep_update_states <- unlist(prepare_update_states())

# LAI, WS, WCR, WFR, Stem_Count, ASW
gamma_LAI = median_pars[52]
gamma_WS = median_pars[53]
gamma_WCR = median_pars[54]
gamma_WFR = median_pars[55] - median_pars[54] ### check this equation
gamma_Stem_Count = 
gamma_ASW = 

rho_LAI = 0 #median_pars[63]
rho_WS = 0  #median_pars[64]
rho_WCR =0 # median_pars[65]
rho_WFR = 0 #median_pars[66]
rho_Stem_Count
rho_ASW




#USE FIRST OBSERVATION AS THE INITIAL CONDITIONS
if(!USE_SYNTHETIC_DATA){
  #realdata = read.csv(file.path(calibFolder, 'Files', "FCR_CTD_wide_Dec14_Dec16.csv", fsep = .Platform$file.sep), header = TRUE)
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
for(i in 1:length(observedLAI)){
  obs_index[i] = which.min(abs(the_depths_init - observedLAI[i])) # how does this work?
}

#A matrix for knowing which state the observation corresponds to
z_states <- t(matrix(obs_index, nrow = length(obs_index), ncol = nsteps))

#Process error 
# fix - this needs a lot of work. Supposed to be sigma = gamma + rho*value

sigma = array(data = NA, dim = 15)
sigma[1] = # what should sigma Age be
sigma[2] = gamma_LAI
sigma[3] = gamma_WS # gamma_WS+ rho_WS * last stem measurement
sigma[4] = gamma_Stem_Density
sigma[5] = gamma_WCR
sigma[6] = gamma_WR
sigma[7] = gamma_foliage_production
sigma[8] = 0 # what should sigma total be?
sigma[9] = 0 # what should fSW be?
sigma[10] = gamma_ET
sigma[11] = gamma_Ctrans # is this total Ctrans?
sigma[12] = 0 # is this gamma_GEP - gamma_ET?
sigma[13] = 0 # runoff
sigma[14] = 0 # WUE_ctrans
sigma[15] = 0 # WUE_ET
Qt = diag(sigma) # fix this part


#Measurement error 
psi = rep(0.0001,length(obs_index)) # update this later

#Initial conditions
x <- array(NA,dim=c(nsteps,nmembers,nstates))
x[1,,] <- rmvnorm(n=nmembers, mean=init_states, sigma=as.matrix(Qt)) # create init_states at beginning


#Matrix to store ensemble specific deviations and innovations
dit = array(NA,dim=c(nmembers,nstates))
#dit_star = array(NA,dim=c(nmembers,nstates)) #Adaptive noise estimation

Lat
SoilClass
MaxASW
MinASW
FR

#loop through time steps
for(i in 2:nsteps){ # runs 3PG for one time step for every member  (would use forecasting code)
  
  #i is the same as mo - make them the same, climate matrices should be built before mo 
  InitialYear
  InitialMonth
  start_age
  #Create array to hold DAPPER predictions for each ensemble
  x_star = array(NA, dim = c(nmembers,nstates))
  for(m in 1:nmembers){
    
    # update site vector with previous states for each ensemble (site_in is site vector)
    # this needs to be moved out of the update_states code
   
    
    ### need to update this each timestep 
    site_in = c(PlantedYear, #PlantedYear # where it left off in last time-step
                PlantedMonth, #"PlantedMonth"
                InitialYear, #"InitialYear"
                InitialMonth, #"InitialMonth"
                start_age,
                WFi = , #"WFi" # can take LAI or WFi
                WRi = x[i-1,m,1], #"WRi"
                WSi = x[i-1,m,2], #"WSi"
                StemNum = , #"StemNoi"
                ASWi = , #"ASWi"
                Lat, #"Lat"
                FR = , #"FR"
                SoilClass, #"SoilClass"
                MaxASW, #"MaxASW"
                MinASW, #"MinASW"
                TotalMonths = 1,
                WFi_H = 0.001,
                WSi_H = 0.001,
                WRi_H = 0.001,
                WCRi = , # 
                IrrigRate = 0.0,
                Throughfall = 1.0,
                tmp_site_index,  
                WCRi_H = 0.0,
                Wbud_H = 0.0,
                LAI = , #
                LAI_h = 0.01
    )
    
    DAPPERstates <- unlist(update_states(mo,site_in))
    
    #4) Fill x_star with temperatures from GLM
    x_star[m,] = DAPPERstates # all of the updated states for the current member iteration
  }
  
  #Corruption [nmembers x nstates] 
  NQt = rmvnorm(n=nmembers, sigma=as.matrix(Qt)) # matrix of diagaonal - error term (noise step in forecasting code)
  
  #Matrix Corrupted state estimate [nmembers x nstates]
  x_corr = x_star + NQt
  ## builing ^ starting point  
  
  
  #Obs for time step
  z_index = which(!is.na(z[i,]))
  
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




