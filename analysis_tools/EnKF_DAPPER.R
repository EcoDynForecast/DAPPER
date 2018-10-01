EnKF_DAPPER <- function(psi, init_Fr, observed, sigma_Fr, LAI_uncert, run_name){ 
  #rm(list = ls())
  if (!"mvtnorm" %in% installed.packages()) install.packages("mvtnorm")
  if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4")
  library(mvtnorm)
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
  plotnum = plotlist - 40000
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
  #ASW_max = prep_update_states$ASW_max
  
  #Initial states
  nstates <- 7 #nlayers_init
  init_LAI = observed$LAI[1] # 
  init_WS = site[8] # WSi
  init_WCR = site[20] # WCRi
  init_WR = site[7] # WRi 
  init_Stem_Count = site[9]
  init_ASW = site[10] # initASW = MaxASW
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
  sigma[1] = (LAI_uncert) #0.01*median_pars[52] # gamma_LAI
  sigma[2] = (median_pars[53] +  median_pars[64]*init_states[2]) # WS, gamma_WS+ rho_WS * last stem measurement
  sigma[3] = (median_pars[54]) # WCR
  sigma[4] = (median_pars[55]) ### check this one ## WFR
  sigma[5] = (.1) # gamma_Stem_Count - look into this, should have a gamma
  sigma[6] = (.1) # gamma_ASW
  sigma[7] = (sigma_Fr) #
  
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
  x[1,,] <- rmvnorm(n=nmembers, mean=init_states, sigma=as.matrix(Qt_init^2)) 
  x[1,which(x[1,,7] < 0),7] = 0
  x[1,which(x[1,,7] > 1 ),7] = 1
  
  InitialYear = startyear
  InitialMonth = startmonth
  
  #Matrix to store ensemble specific deviations and innovations
  dit = array(NA,dim=c(nmembers,nstates))
  #dit_star = array(NA,dim=c(nmembers,nstates)) #Adaptive noise estimation
  
  #site[11] = plot_lat # lat
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
      site[26] = x[i-1, m, 1] # LAI
      site[12] = x[i-1,m,7] #FR
      
      DAPPER_states <- update_states(mo = i, pixelnum = 1, output_dim, pars, site, nopars, nosite, thin_event, met, exclude_hardwoods, median_pars)
      
      #Fill x_star with states from DAPPER
      x_star[m,] = c(DAPPER_states,site[12]) # all of the updated states for the current member iteration
      NQt = array(NA, dim = c(7))
      NQt = rmvnorm(n=1, sigma=as.matrix(Qt^2))
      # rmvnorm wants variance, rnorm wants standard deviation
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
      ens_mean = apply(x_corr,2, mean)
      #ens_mean_star = apply(x_star, 2, mean) #Adaptive noise estimation
      
      #Loop through ensemble members
      for(m in 1:nmembers){  
        
        #Ensemble specific deviation
        dit[m,] = x_corr[m,]-ens_mean
        #dit_star[m,] x_star[m,] - ens_mean_star #Adaptive noise estimation
        
        #Ensemble specific estimate and innovation covariance
        if(m == 1){
          Pit = dit[m,] %*% t(dit[m,]) 
        }else{
          Pit = dit[m,] %*% t(dit[m,])  +  Pit
        }
      }
      
      #estimate covariance
      Pt <- Pit/nmembers
      
      #Kalman gain
      Kt <- Pt %*% t(H) %*% solve(H%*%Pt%*%t(H)+psi_t)
    
      #Update states array (transposes are necessary to convert between the dims here and the dims in the EnKF formulations)
      x[i,,] <- t(t(x_corr) + Kt%*%(D_mat - H%*%t(x_corr)))
      
      if(length(which(is.na(x[i,,]))) > 0){dies = i}
    }
  }
  
  ###### save x as .rdata file ####
  saveRDS(x, file=paste(EnKF_directory, "/x_", run_name, ".RDS", sep = ""))
  ###### PLOTS #############
  
  #LAI
  plot(startyear + (1/12)*1:nrow(x),x[ , 1, 1],type='l',ylim=c(0,5), main = 'LAI', ylab = 'LAI')
  for (n in 2:nmembers)
  {
    points(startyear + (1/12)*1:nrow(x),x[ , n, 1],type='l')
  }
  points(startyear + (1/12)*1:nrow(z), z, col = 'red')
  legend(x = "topright",inset = 0,
         legend = c("ensemble members", "observations"),     col=c("black","red"), lwd=5, cex= 0.75, horiz = FALSE, lty = c(1,0), pch = c(26,8))
  
  # NOT REAL DATA - making example diagram for paper
  # axis1 = c(1,2,3,4,5);
  # plot(axis1,x[8:12, 1, 1],type='l',ylim=c(1,4.5), main = 'LAI', xlab = 'month', ylab = 'LAI', xaxt='n')
  # axis(1,at=seq_along(axis1), labels=c("aug","sep","oct","nov","dec"), las=2)
  # for (n in 40:70)
  # {
  #   points(axis1,x[8:12, n, 1],type='l') 
  # }
  # points(c(3,4,5), c(3.2, 2.5, 2.2), col = 'red', pch = 8)
  # 
  # legend(x = "topright",inset = 0,
  #        legend = c("ensemble members", "observations"),     col=c("black","red"), lwd=5, cex= 0.75, horiz = TRUE, lty = c(1,0), pch = c(26,8))
  # 
  
  title(run_name, outer = TRUE, cex = 1.5)
  quant_Fr = array(NA, dim = c(nmonths, 3))
  for (n in 1:nmonths){
    quant_Fr[n, ] = quantile(x[n, , 7], c(0.10, 0.5, 0.90))
  }
  
  # LAI model quantiles
  plot(startyear + (1/12)*1:nrow(x) -(1/12), quant_Fr[ , 1],type='l',ylim=range(quant_Fr, na.rm = TRUE),xlab = 'months',ylab = 'Fr')
  lines(startyear + (1/12)*1:nrow(x) -(1/12), quant_Fr[ ,2])
  lines(startyear + (1/12)*1:nrow(x)-(1/12), quant_Fr[ ,3])
  
  #Fr as histogram
  hist(x[nmonths, , 7], main = 'Fr at last step')

  #### need to save the data as .rdata file
}



