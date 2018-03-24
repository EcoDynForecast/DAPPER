EnKF_DAPPER <- function(psi, init_Fr, observed, Fr_uncert, LAI_uncert, run_name){ 
  #rm(list = ls())
  if (!"mvtnorm" %in% installed.packages()) install.packages("mvtnorm")
  if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4")
  if (!"glmtools" %in% installed.packages()) install.packages('glmtools', repos=c('http://cran.rstudio.com', 'http://owi.usgs.gov/R'))
  library(mvtnorm)
  library(glmtools)
  library(ncdf4)
  
  #### set working directory
  #pathDAPPER = '/Users/laurapuckett/Documents/Research/Current/DAPPER/' 
  #EnKF_directory = '/Users/laurapuckett/Documents/Research/Current/DAPPER/analysis_tools/'
  setwd(EnKF_directory)
  
  #### source the required functions
  source('prepare_update_states.R')
  source('update_states.R')
  
  #---------------------------------------------------------
  ###### MAIN EnKF CODE ######
  
  
  years = observed$Year
  startyear = min(years)
  startmonth = min(observed$Month[which(observed$Year == startyear)])
  endyear = max(years)
  start_age = observed$Age[which(observed$Year == startyear && observed$Month == startmonth)]
  nmonths = (endyear - startyear +1)*12
  observed_months = (observed$Year-startyear)*12 + observed$Month
  
  # create array to hold observation or NA for each month
  z <- array(data = NA, dim = c(nmonths,1))
  for (i in 1:nmonths){
    for (j in 1:length(observed_months))
    {
      if(i == observed_months[j] && j %% 1 == 0)
      {
        z[i, 1] <- observed$LAI[j]
        
      }
    }
  }
  
  # Retreieve variables from prepare_update_states() needed for update_states()
  #mo = 1
  prep_update_states <- prepare_update_states(plotnum, startyear, nmonths, start_age, endyear)
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
  
  #Initial states
  nstates <- 7 #nlayers_init
  init_LAI = 1.1 #realdata$LAI[1] # 
  init_WS = 49.56 # WSi
  init_WCR = 13.82 # WCRi
  init_WR = 2.68 # WRi was 0.16
  init_Stem_Count = 1511 #"StemNoi"
  init_ASW = ASW_max 
  #init_Fr = site[12] - function input for now
  init_states <- array(data = c(init_LAI, init_WS, init_WCR, init_WR, init_Stem_Count, init_ASW, init_Fr), dim = c(nstates,1))
  
  if(!USE_OBS_CONTRAINT){
    z_obs = z
    z[,] = NA
  }
  
  # Defining variables for Qt matrix
  # need to incorporate rho into sigma later
  # the values a re squared, because the variance, not std. dev., is what is used in rmvnorm()
  sigma = array(data = NA, dim = nstates)
  sigma[1] = (LAI_uncert)^2 #0.01*median_pars[52] # gamma_LAI
  sigma[2] = (median_pars[53] +  median_pars[64]*init_states[2])^2 # WS, gamma_WS+ rho_WS * last stem measurement
  sigma[3] = (median_pars[54])^2 # WCR
  sigma[4] = (median_pars[55])^2 ### check this one ## WFR
  sigma[5] = (.1)^2 # gamma_Stem_Count - look into this, should have a gamma
  sigma[6] = (.1)^2 # gamma_ASW
  sigma[7] = (Fr_uncert)^2 ##0.001
  
  Qt = diag(sigma) 
  #rho_LAI = 0 #median_pars[63]
  #rho_WS = 0  #median_pars[64]
  #rho_WCR =0 # median_pars[65]
  #rho_WFR = 0 #median_pars[66]
  #rho_Stem_Count = 0
  #rho_ASW = 0
  
  #Measurement error 
  #psi = 0.0001 # function input for now
  
  Qt_init = Qt
  x <- array(NA,dim=c(nmonths,nmembers,nstates))
  x[1,,] <- rmvnorm(n=nmembers, mean=init_states, sigma=as.matrix(Qt_init)) 
  x[1,which(x[1,,7] < 0),7] = 0
  x[1,which(x[1,,7] > 1 ),7] = 1
  
  InitialYear = startyear
  InitialMonth = startmonth
  
  #Matrix to store ensemble specific deviations and innovations
  dit = array(NA,dim=c(nmembers,nstates))
  #dit_star = array(NA,dim=c(nmembers,nstates)) #Adaptive noise estimation
  
  #loop through time steps
  for(i in 2:nmonths){ # runs 3PG for one time step for every member
    site[3] = InitialYear
    site[4] = InitialMonth
    site[5] = start_age
    #Create array to hold DAPPER predictions for each ensemble
    x_star = array(NA, dim = c(nmembers,nstates))
    x_corr = array(NA, dim = c(nmembers,nstates))
    for(m in 1:nmembers){
      
      ### updates the states in site array each timestep 
      site[6] = -99 #"WFi" # can take LAI or WFi
      site[7] = x[i-1, m, 4] # WR - is this correct?
      if(site[7] < 0 | is.na(site[7]) ){ # cheat for now, need to fix the problem for real later
        site[7] = 0
      }
      site[8] = x[i-1, m, 2] #"WSi"
      site[9] = x[i-1, m, 5] #"StemNoi"
      site[10] = x[i-1, m, 6] #"ASWi"
      site[20] = x[i-1, m, 3] # WCRi
      site[26] = x[i-1, m, 1] # same as WFi above
      site[12] = x[i-1,m,7] #FR
      
      DAPPER_states <- update_states(mo = i, plotnum, output_dim, pars, site, nopars, nosite, thin_event, met, exclude_hardwoods, median_pars, ASW_max)
      
      #Fill x_star with states from DAPPER
      x_star[m,] = c(DAPPER_states,site[12]) # all of the updated states for the current member iteration
      #NQt = rmvnorm(n=1, sigma=as.matrix(Qt))
      NQt = array(NA, dim = c(7))
      NQt[1] = rnorm(n=1, 0, sigma[1])
      NQt[2] = rnorm(n=1, 0, sigma[2])
      NQt[3] = rnorm(n=1, 0, sigma[3])
      NQt[4] = rnorm(n=1, 0, sigma[4])
      NQt[5] = rnorm(n=1, 0, sigma[5])
      NQt[6] = rnorm(n=1, 0, sigma[6])
      NQt[7] = rnorm(n=1, 0, sigma[7])
      x_corr[m,] = x_star[m,] + NQt
      
      if(x_corr[m, 7] >1)
      {
        x_corr[m,7] = 1
      }
      for(sta in 1:nstates){
        if(x_corr[m, sta] < 0){
          x_corr[m, sta] = 0
        }
      }
    }
    
    InitialMonth = InitialMonth + 1
    if(InitialMonth > 12){
      InitialMonth = 1
      InitialYear = InitialYear + 1
    }
    start_age = start_age + 1/12
    
    #Obs for time step
    z_index = which(!is.na(z[i, ]))
    
    #if no observations at a time step then just propogate model uncertainity
    if(length(which(!is.na(z[i, ]))) == 0){
      x[i,,] = x_corr
      if(NO_UNCERT){
        x[i,,] = x_star
      }
    }else{
      #   print(i)
      #if observation then calucate Kalman adjustment
      zt = z[i,z_index] # will be time series of LAI from landsat
      #z_states_t = z_states[i,z_index]
      yit = array(NA,dim=c(nmembers,length(zt)))
      
      H <- array(0,dim=c(length(zt),nstates))
      for(j in 1:length(z_index)){
        H[j,1] = 1
      }
      
      #Extract the data uncertainity for LAI at time-step
      psi_t = psi[z_index]
      
      
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
          Sit = tcrossprod(yit[m,])
        }else{
          Pit = tcrossprod(dit[m,]) +  Pit
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
      
      #Ensemble specific updated state
      for(m in 1:nmembers){
        x[i,m,] = t((x_corr[m,]) + crossprod(t(Kt), yit[m,]))
        if(NO_UNCERT){
          x[i,m,] = x_star[m,]
        }
        for(n in 1:nstates){
          if( x[i,m,n] < 0 ){
            x[i,m,n] = 0
          }
        }
      }
      if(length(which(is.na(x[i,,]))) > 0){dies = i}
    }
  }
  
  ###### save x as .rdata file ####
  saveRDS(x, file=paste(EnKF_directory, "/x.rds", sep = ""))
  ###### PLOTS #############
  
  par(mfrow = c(3,3), oma=c(0,0,2,0))
  # LAI
  plot(x[ , 1, 1],type='l',ylim=c(0,5), main = 'LAI')
  for (n in 2:nmembers)
  {
    points(x[ , n, 1],type='l') 
  }
  points(z, col = 'red')
  
  title(run_name, outer = TRUE, cex = 1.5)
  
  # WS
  ylim = range(c(x[ ,,2]))
  plot(x[ , 1, 2],type='l',ylim=ylim, main = 'WS')
  for (n in 2:nmembers)
  {
    points(x[ , n, 2],type='l') 
  }
  #Fr
  ylim = range(c(x[ ,,7]))
  plot(x[ , 1, 7],type='l',ylim=ylim, main = 'Fr')
  for (n in 2:nmembers)
  {
    points(x[ , n, 7],type='l') 
  }
  
  quant_WS = array(NA,dim=c(nmonths, 3))
  quant_LAI = array(NA, dim = c(nmonths, 3))
  quant_Fr = array(NA, dim = c(nmonths, 3))
  for (n in 1:nmonths){
    quant_WS[n, ] = quantile(x[n, , 2], c(0.10,0.5,0.90))
    quant_LAI[n, ] = quantile(x[n, , 1], c(0.10, 0.5, 0.90))
    quant_Fr[n, ] = quantile(x[n, , 7], c(0.10, 0.5, 0.90))
  }
  plot(quant_LAI[ , 1], type = 'l', ylim = range(quant_LAI, na.rm = TRUE), xlab = 'months', ylab = 'LAI')
  lines(quant_LAI[ ,2])
  lines(quant_LAI[ ,3])
  
  plot(quant_WS[ , 1],type='l',ylim=range(quant_WS, na.rm = TRUE),xlab = 'months',ylab = 'WS')
  lines(quant_WS[ ,2])
  lines(quant_WS[ ,3])
  
  plot(quant_Fr[ , 1],type='l',ylim=range(quant_Fr, na.rm = TRUE),xlab = 'months',ylab = 'Fr')
  lines(quant_Fr[ ,2])
  lines(quant_Fr[ ,3])
  
  plot(0,0) # place-holder
  
  hist(x[nmonths, , 3], main = 'WS at last step')
  hist(x[nmonths, , 7], main = 'Fr at last step')

  #### need to save the data as .rdata file
}



