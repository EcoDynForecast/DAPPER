### 3-PG Regional Runs for Quinn, Incorporating Parameter Uncertainty

#if (!"spatial.tools" %in% installed.packages()) install.packages("spatial.tools")
#if (!"readr" %in% installed.packages()) install.packages("readr")
if (!"compiler" %in% installed.packages()) install.packages("compiler")
#if (!"dplyr" %in% installed.packages()) install.packages("dplyr")
if (!"doParallel" %in% installed.packages()) install.packages("doParallel")
#if(!"tidyr" %in% installed.packages()) install.packages("tidyr")

#library(spatial.tools)
#library(readr)
#library(dplyr)
library(compiler)
library(doParallel)

enableJIT(1)

GCMs <- c("bcc-csm1-1-m","bcc-csm1-1","BNU-ESM","CanESM2","CCSM4","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-ESM2G","GFDL-ESM2M","HadGEM2-CC365","HadGEM2-ES365","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC-ESM-CHEM","MIROC-ESM","MIROC5","MRI-CGCM3","NorESM1-M", "metdata") # metdata is the Idaho reference dataset, indexed by rcpCase [1]

rcpCases <- c('baseline', 'rcp45','rcp85')

## Critical Parameters -- USER INPUT REQUIRED ####-----------------------------
# Working Directory 
setwd('/home/rqthomas/DAPER_regional_analysis')

# GCM
myGCM_list <- c(GCMs[3],GCMs[1:5]) # CCSM4, for example
myGCM_list <- c(GCMs[1],GCMs[1],GCMs[2],GCMs[2],GCMs[4],GCMs[4])
#myGCM_list <- c(GCMs[21],GCMs[1:20],GCMs[1:20],GCMs[1:20],GCMs[1:20])
# RCP
myRCP_list <- c(rcpCases[1],rep(rcpCases[3],5),rep(rcpCases[3],5)) # RCP8.5, for example
myRCP_list <- c(rcpCases[3],rcpCases[3],rcpCases[3],rcpCases[3],rcpCases[3],rcpCases[3])
#myRCP_list <- c(rcpCases[1],rep(rcpCases[3],80))

## Tuning Parameters
# Note that the metdata available years run from 1971-2011
yearStart_list <- c(1985,rep(1985,5),rep(2030,5)) # Change this as needed!
yearStart_list <- c(1985,2030,1985,2030,1985,2030)
#yearStart_list <- c(1985,rep(1985,20),rep(2020,20),rep(2045,20),rep(2070,20))

# Variability type
variation_type <- 'parameter_and_process' #parameter_and_process' # Choose from c('parameter_only', 'process_only', 'parameter_and_process')

Mort_uncert = FALSE
FR_uncert = FALSE
HOLD_CO2 = FALSE
precipModifier <- 1.0 # No toggle to precip
MaxFR <- FALSE # No fertilization
#variation_type <- 'parameter_only'
# Number of samples from parameter chain (will be ignored in the process_only case)
par_sample_size <- 500

# Number of samples from process error distribution (will be ignored in the parameter_only case)
proc_sample_size <- 10


# Number of processors to be used
nprocessors <- 24 # Make sure to change this to the max processors!

## Tuning Parameters
rotationAge <- 25

load("/home/rqthomas/hokieone/DAPER/chains/revision_base_chain1.1.2017-04-28.15.56.15.final.Rdata")
load('/home/rqthomas/DAPER_regional_analysis/FR_par_chain.Rdata')
code_library = "/home/rqthomas/DAPER_regional_analysis/r3pg_interface.so"
CO2 <- read.csv('CO2_Concentrations_from_CMIP5_1950-2099.csv')
Soils <- read.csv('Soil_Inputs_LPNR_Clipped_and_Imputed_v4.csv')


## 3-PG Wrapper ####-----------------------------------------------------------
run.3PG.quick <- function(myHUC, yearStart, yearEnd, startAge = 2, precipModifier = 1, MaxFR = FALSE, curr_pars,base_year = 1950, curr_FR_pars, FR_pars,HOLD_CO2 = FALSE, ...){ 

  pars <- curr_pars[1:npars_used_by_fortran]
  # Reference Parameters
  nomonths <- (yearEnd - yearStart + 1) * 12
  HUCIndex <- which(frostDays[, 1] == myHUC)[1] # For the Climate Variables
  HUCIndexSoils <- which(Soils[, 1] == myHUC)[1] # For the Soils variables
  
  ## Extracting Climate
  yearStartPosition <- (yearStart - base_year) * 12 + 1 # Selecting the first column to draw from
  
  # CO2
   if(HOLD_CO2 == TRUE){
    yearStart = yearStart_list[r]
    yearEnd <- yearStart + rotationAge - 1
    tmpCO2 <- CO2$CO2_Concentration_RCP85[is.element(CO2$Year, 1985:(1985 + rotationAge - 1))] # RCP8.5 should match historical readings in your timeframe
    monthlyCO2 <- c(matrix(rep(tmpCO2, 12), nrow = 12, byrow = T))
  }else{
  tmpCO2 <- CO2$CO2_Concentration_RCP85[is.element(CO2$Year, yearStart:yearEnd)] # RCP8.5 should match historical readings in your timeframe
  monthlyCO2 <- c(matrix(rep(tmpCO2, 12), nrow = 12, byrow = T))
  } 
  # All Climate
  col_index <- c(1 + yearStartPosition:(yearStartPosition + 12 * rotationAge - 1))
  met <- cbind(
    unlist(tMin[HUCIndex,col_index]),
    unlist(tMax[HUCIndex,col_index]),
    unlist(rain[HUCIndex,col_index])*precipModifier,
    unlist(solarRad[HUCIndex,col_index]/1000000*86400), # Converting from W/m^2_s to W/ha_day
    unlist(frostDays[HUCIndex,col_index]),
    monthlyCO2
  )
  met[is.na(met[,5]),5] <- 0
  if(sum(is.na(met)) > 0){return(rep(NA, 6))}
  
  
  ## Main inputs -----------------------------------------------------
  tmpSoils <- Soils[HUCIndexSoils,][1,]
  
  if(is.na(tmpSoils$aws0150)) {return(rep(NA, 6))}
  
  # FR -----------------------------------------------------
  if(MaxFR){
    tmpFR <- 1.0
  } else {
      tmpFR = 1/(1+exp((FR_pars[1]*tmpSoils$MAT-FR_pars[2] * tmpSoils$SIm)))
      #tmpFR <- 1/(1 + exp((FR_pars[1] - FR_pars[2] * tmpSoils$SIm)))
      #tmpFR <- FR_pars[1] + FR_pars[2] * tmpSoils$SIm
      #if(FR_uncert == TRUE){
      #  tmpFR = rnorm(1,tmpFR,FR_pars[3])
      #}
  }
  
  if(Mort_uncert == TRUE){
  pars[20] = rnorm(1,curr_pars[20],curr_pars[79])
  pars[40] = rnorm(1,curr_pars[40],curr_pars[80])
  }
   
  #if(is.nan(tmpFR)){print(c(tmpFR,FR_pars[1],FR_pars[2],tmpSoils$SIm))}
  #if(is.na(tmpFR)){print(c(tmpFR,FR_pars[1],FR_pars[2],tmpSoils$SIm))}
  #if(is.infinite(tmpFR)){print(c(tmpFR,FR_pars[1],FR_pars[2],tmpSoils$SIm))}

  # Initialization Inputs (+ Soil) -----------------------------------------------------
  site_in <- c(yearStart, # PlantedYear
               1, # PlantedMonth
               yearStart, # InitialYear
               0, # InitialMonth
               startAge, # 'EndAge', ACTUALLY THE STARTING AGE!!!
               .95, # WFi
               0.028, # WRi
               2, # WSi
               1235, # StemNoi
               tmpSoils$aws0150, # ASWi
               tmpSoils$LAT, # Lat
               tmpFR, # FR
               1, # SoilClass
               as.numeric(tmpSoils$aws0150), # MaxASW
               0, # MinASW
               nomonths-startAge*12, # TotalMonths
               0.001, # WFi_H
               0.001, # WSi_H
               0.001, # WRi_H
               2*0.2,
               IrrigRate = 0.0,
               Throughfall = 1.0,
               tmp_site_index = 0.0) # WCRi
  site <- array(site_in)
  site[is.na(site)] <- 0
  
  thin <- array(NA, dim = c(nomonths,2))
  thin[is.na(thin)] <- 0
  
  #--------------------------------------------------------------------------------
  tmp <- .Fortran( "r3pg_interface",
                   output_dim = as.integer(output_dim),
                   met = as.double(t(met)), 
                   pars = as.double(pars),
                   site = as.double(site),
                   thin = as.double(thin),
                   out_var = as.double(array(0,dim=c(nomonths,output_dim))),
                   nopars = as.integer(nopars),
                   nomet = as.integer(dim(met)[2]),
                   nosite = as.integer(length(site)),
                   nooutputs = as.integer(output_dim),
                   nomonths=as.integer(nomonths),
                   nometmonths = as.integer(dim(met)[1]),
                   nothin = as.integer(dim(thin)[1]),
                   exclude_hardwoods = as.integer(exclude_hardwoods))
  #-------------------------------------------------------------------------------
  
  # READ IN FORTRAN OUTPUT FILE PROCESS FORTRAN OUTPUT SO THAT IT CAN BE COMPARED TO THE DATA
  output0 <- data.frame(array(tmp$out_var, dim = c(nomonths, output_dim)))
  
  names(output0) <- c('Year','Month','StandAge','WF','WR','WS','StemNo','LAI','delLitter','avStemMass','Vob','BasArea','avDBH','ASW','GPPdm','NPP','ANPP','pF','pR','pS','conductance','Transpiration','InterceptedRain','EvapTransp','MinASW','MaxASW','runoff','Rain','monthlyIrrig','fNutr','fT','fFrost','genetics','fCalpha','fVPD','fSW','fAge','fComp','fCg','VPD','wSmax','WF_h','WR_h','WS_h','LAI_h', 'GPPdm_h','NPP_h','ANPP_h','pF_h','pR_h','pS_h','delLitter_h','Transpiration_h','CoarseRoots','TotalBiomass','SLA','NPPplusRABG')
  
  # Selecting the outputs to write  
  output <- output0[,c('StandAge','WS', 'WF', 'WR', 'CoarseRoots','LAI','SLA')]
  
  age_25 <- which(output$StandAge >= 24.9 & output$StandAge < 26.0)
  output <- colMeans(output[age_25,])
  
  # Final Asembly and Output
  finalOutput <- c(myHUC, output)
  
  finalOutput
}

## Model Runs ####-------------------------------------------------------------



if(variation_type == 'none'){
  nsamples <- 1
  njitters <- 1
  # Setting up the pass to Fortran
  niter <- dim(accepted_pars_thinned_burned)[1]
  npars <- 80
  npars_used_by_fortran <- 48
  set.seed(2016)
  pars_sample_index <- sample(seq(1,niter,1),nsamples,replace = TRUE)
  simulations_parm_proc <- array(NA,dim = c(nsamples,7))
  simulations_proc <- array(NA,dim = c(nsamples,5))
}

# Parameters
if(variation_type == 'parameter_only'){
  nsamples <- par_sample_size
  njitters <- 1
  # Setting up the pass to Fortran
  niter <- dim(accepted_pars_thinned_burned)[1]
  npars <- 80
  npars_used_by_fortran <- 48
  set.seed(2016)
  pars_sample_index <- sample(seq(1,niter,1),nsamples,replace = TRUE)
  simulations_parm_proc <- array(NA,dim = c(nsamples,7))
  simulations_proc <- array(NA,dim = c(nsamples,5))
}

if(variation_type == 'process_only'){
  nsamples <- 1
  njitters <- proc_sample_size
  # Setting up the pass to Fortran
  niter <- dim(accepted_pars_thinned_burned)[1]
  npars <- 80
  npars_used_by_fortran <- 48
  set.seed(2016)
  pars_sample_index <- sample(seq(1,niter,1),nsamples,replace = TRUE)
  simulations_parm_proc <- array(NA,dim = c(nsamples,7))
  simulations_proc <- array(NA,dim = c(nsamples,5))
}

if(variation_type == 'parameter_and_process'){
  nsamples <- par_sample_size
  njitters <- proc_sample_size
  # Setting up the pass to Fortran
  niter <- dim(accepted_pars_thinned_burned)[1]
  npars <- 80
  npars_used_by_fortran <- 48
  set.seed(2016)
  pars_sample_index <- sample(seq(1,niter,1),nsamples,replace = TRUE)
  simulations_parm_proc <- array(NA,dim = c(nsamples,7))
  simulations_proc <- array(NA,dim = c(nsamples,5))
}

for(r in 1:length(myRCP_list)){
  
  myRCP = myRCP_list[r]
  myGCM = myGCM_list[r]
  yearStart = yearStart_list[r]
  yearEnd <- yearStart + rotationAge - 1
  
  
  print(paste('Run ',r,': starting RCP ',myRCP,' for model ',myGCM,' for start year ',yearStart,sep="")) 
  
  npars <- 80
  curr_pars <- rep(NA,npars)
  for(i in 1:npars){
    curr_pars[i] <- median(accepted_pars_thinned_burned[,i]) # Taking median here, modify for bootstrapping
  }
  
  npars_used_by_fortran <- 48
  pars <- curr_pars[1:npars_used_by_fortran]
  
  noutput_variables <- 63
  output_dim <- noutput_variables  # NUMBER OF OUTPUT VARIABLES
  nomet <- 6  # NUMBER OF VARIABLES IN METEROLOGY (met)
  nopars <- length(pars)
  nomonths <- (yearEnd-yearStart+1)*12
  exclude_hardwoods <- 1
  npars <- 80
  
  
  RCP_index <- which(rcpCases == myRCP)
  base_year <- 1950
  if(RCP_index == 1){base_year <- 1979}

# Climate Data Loading
  frostDays <- read.csv(paste('/home/rqthomas/DAPER_regional_analysis/met_files/frostDays_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
  rain <- read.csv(paste('/home/rqthomas/DAPER_regional_analysis/met_files/rain_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
  solarRad <- read.csv(paste('/home/rqthomas/DAPER_regional_analysis/met_files/solarRad_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
  tMax <- read.csv(paste('/home/rqthomas/DAPER_regional_analysis/met_files/tMax_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
  tMin <- read.csv(paste('/home/rqthomas/DAPER_regional_analysis/met_files/tMin_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
  
myHUCsFinal <- intersect(Soils$HUC12, frostDays$HUC12)

#myHUCsFinal <- myHUCsFinal[which(myHUCsFinal == 31101030505 | myHUCsFinal == 30601050202 | myHUCsFinal == 111401070404 | myHUCsFinal == 20802031301)]

cl <- makeCluster(nprocessors)  
registerDoParallel(cl)
print('starting sample loop')
print(nsamples)

# for(ns in 1:nsamples){
 final_output_GCM_RCP <- foreach(ns = 1:nsamples, .export = ls()) %dopar% { # Will need to add objects,
  #final_output_GCM_RCP <- foreach(ns = 1:nsamples) %do% { # Will need to add objects,packages as needed
  # Load Fortran code
  dyn.load(code_library)

  par_ns_index <- pars_sample_index[ns]
  curr_pars <- accepted_pars_thinned_burned[par_ns_index,]
  #curr_FR_pars <- FR_par_chain[par_ns_index,]
  if(variation_type == 'process_only' | ns == 1){
    curr_pars <- rep(NA,npars)
    for(i in 1:npars){
      curr_pars[i] <- median(accepted_pars_thinned_burned[,i]) # Taking median here, modify for bootstrapping
    }
  }
  #FR_pars = c(FR_par_chain[par_ns_index,1],FR_par_chain[par_ns_index,2],FR_par_chain[par_ns_index,3]) 
  #if(variation_type == 'process_only' | ns == 1){ 
  #   FR_pars = c(median(FR_par_chain[,1]),median(FR_par_chain[,2]),median(FR_par_chain[,3])) 
  #}

  FR_pars = c(curr_pars[50],curr_pars[51])

  if(njitters == 1){
  bucket <- array(NA, dim = c(length(myHUCsFinal), 3, 1)) # columns for HUC12, WF, WS, WR, CoarseRoots, and TotalBiomass
  tmp_bucket <- array(NA, dim = c(length(myHUCsFinal), 6, 1)) # columns for HUC12, WF, WS, WR, CoarseRoots, and TotalBiomass
  }else{
  bucket <- array(NA, dim = c(length(myHUCsFinal), 3, (njitters+1))) # columns for HUC12, WF, WS, WR, CoarseRoots, and TotalBiomass
  tmp_bucket <- array(NA, dim = c(length(myHUCsFinal), 6, (njitters+1))) # columns for HUC12, WF, WS, WR, CoarseRoots, and TotalBiomass
  }
  tmp <- matrix(NA, nrow = length(myHUCsFinal), ncol = 8)
  for(HUCindex in 1:length(myHUCsFinal)){
    # bucket[HUCindex,,ns] 
    tmp[HUCindex,] <- run.3PG.quick(myHUC = myHUCsFinal[HUCindex], yearStart = yearStart, yearEnd = yearEnd, startAge = 2, precipModifier = 1, MaxFR = MaxFR, base_year = base_year, curr_pars = curr_pars,FR_pars = FR_pars, HOLD_CO2 = HOLD_CO2, Mort_uncert = Mort_uncert, FR_uncert = FR_uncert) 
  }
  tmp <- tmp[, -2]
  # c('HUC12','StandAge','WS', 'WF', 'WR', 'CoarseRoots')
  # c('HUC12', 'WS', 'WF', 'WR', 'CoarseRoots')
  
  # if(variation_type == 'parameter_only') # Possibly add a special case here?
  # HUC 12
if(njitters == 1){


  tmp_bucket[,1,] <- tmp[,1] 

  # WS
  tmp_bucket[,2,] <- tmp[,2]

  #WF
  #BECAUSE THE PROCESS ERROR WAS CALCULATED ON LAI RATHER THAN WF, THE PROCESS ERROR NEEDS TO BE APPLIED TO LAI AND THEN CONVERTED TO WF USING SLA
  #tmp_bucket[,3,] <- tmp[,6]
  #SLA = 0.1*tmp[,7]
  
  SLA = curr_pars[7] + (curr_pars[7] - curr_pars[8]) * exp(-0.693147181 * (25 / curr_pars[9])**2)
  tmp_bucket[,3,] <- tmp[,6]/(0.1*SLA)

  # WR
  tmp_bucket[,4,] <- tmp[,4]

  # CoarseRoots
  tmp_bucket[,5,]  <- tmp[,5]
  #Total Sums 
  tmp_bucket[,6,] <- rowSums(tmp_bucket[, -c(1,6),], c(3))
}else{
 
 tmp_bucket[,1,] <- tmp[,1]
 #tmp_bucket[,1,2:(njitters+1)] <- matrix(tmp[,1], nrow = length(myHUCsFinal), ncol = njitters)

  # WS
  tmp_bucket[,2,1] <- tmp[,2]
  tmp_bucket[,2,2:(njitters+1)] <- matrix(rnorm(njitters * length(myHUCsFinal), mean = tmp[,2], sd = curr_pars[53] + curr_pars[66] * tmp[,2]), ncol = njitters)

  #WF
  #BECAUSE THE PROCESS ERROR WAS CALCULATED ON LAI RATHER THAN WF, THE PROCESS ERROR NEEDS TO BE APPLIED TO LAI AND THEN CONVERTED TO WF USING SLA
  SLA = curr_pars[7] + (curr_pars[7] - curr_pars[8]) * exp(-0.693147181 * (25 / curr_pars[9])**2)

  #THIS WORKS BECAUSE SLA DEPENDES ON AGE AND ALL HUCS HAVE THE SAME AGE
  tmp_bucket[,3,1] = tmp[,6]/(0.1*SLA)
  tmp_bucket[,3,2:(njitters+1)] <- matrix(rnorm(njitters * length(myHUCsFinal), mean = tmp[,6], sd = curr_pars[57]), ncol = njitters)/(0.1*SLA)

  # WR
  tmp_bucket[,4,1] <- tmp[,4]
  tmp_bucket[,4,2:(njitters+1)] <- matrix(rnorm(njitters * length(myHUCsFinal), mean = tmp[,4], sd = curr_pars[58]), ncol = njitters)

  # CoarseRoots
  tmp_bucket[,5,1] <- tmp[,5]
  tmp_bucket[,5,2:(njitters+1)] <- matrix(rnorm(njitters * length(myHUCsFinal), mean = tmp[,5], sd = curr_pars[61]), ncol = njitters)


  # TotalBiomass
      tmp_bucket[,6,] <- apply(tmp_bucket[, -c(1,6),], c(3), rowSums, na.rm = T)
}
  bucket[,1,] = tmp_bucket[,1,]
  bucket[,2,] = tmp_bucket[,2,]
  bucket[,3,] = tmp_bucket[,6,]
  
  bucket
}
stopCluster(cl)

big_bucket <- array(unlist(final_output_GCM_RCP), dim = c(dim(final_output_GCM_RCP[[1]]),length(final_output_GCM_RCP)))

rm(frostDays, rain, solarRad, tMax, tMin)

save.image(file = paste('/work/newriver/rqthomas/DAPER_regional_analysis/3PG_runs_for_', variation_type, '_using_', nsamples, '_paramter_draws_and_', njitters, '_process_error_draws_based_on_', myGCM, '_and_', myRCP,'_and_startyear_of_',yearStart,'_FRset_',maxFR,'_HOLDCO2_',HOLD_CO2,'.Rdata', sep = ''))      

#library(tidyr)
#library(dplyr)
#library(readr)
#library(foreach)

#load(paste('3PG_runs_for_', variation_type, '_using_', nsamples, '_paramter_draws_and_', njitters, '_process_error_draws_based_on_', myGCM, '_and_', myRCP, '.Rdata'))
#load('3PG_runs_for_parameter_and_process_using_10_paramter_draws_and_10_process_error_draws_based_on_metdata_and_baseline.Rdata')

#summarise.by.HUC <- function(HUC_index, var_index = 2, method = 'all', summary_function = function(x){mean(x, na.rm = T)}, my_data = big_bucket, ...){
#  # my_HUC <- my_data[HUC_index,2,1,1]
#
#   tmp <- my_data[HUC_index, var_index, , ]
#
#  output <- NA
#
#  if(method == 'all'){
#    parameter_stats <- apply(tmp, 1, FUN = summary_function)
#    proc_stats <- apply(tmp, 2, FUN = summary_function)
#    par_and_proc_stats <- summary_function(c(tmp))
#
#    output <- c(parameter_stats, proc_stats, par_and_proc_stats)
#  }
#
#  if(method == 'par'){
#    parameter_stats <- apply(tmp, 1, FUN = summary_function)
#    output <- c(parameter_stats)
#  }
#
#  if(method == 'process'){
#    proc_stats <- apply(tmp, 2, FUN = summary_function)
#    output <- c(proc_stats)
#  }
#
#  if(method == 'par_and_proc'){
#    par_and_proc_stats <- summary_function(c(tmp))
#
#    output <- c(par_and_proc_stats)
#  }
#
#  output
#}



#summarise.by.HUC(sample(1:nrow(big_bucket), 1), method = 'par_and_proc')

#summarise.by.HUC(sample(1:nrow(big_bucket), 1), method = 'par_and_proc', summary_function = sd)

#test_output <- foreach(HUC_index = 1:nrow(big_bucket), .combine = 'rbind') %do% {i
#  summarise.by.HUC(HUC_index, var_index = 3,method = 'par_and_proc', summary_function = mean)
#}

#test_output_df <- tbl_df(cbind(big_bucket[,1,1,1], test_output))
#names(test_output_df) <- c('HUC12', 'par_mean')

#write.csv(test_output_df, paste('MEAN_3PG_runs_for_', variation_type, '_using_', nsamples, '_paramter_draws_and_', njitters, '_process_error_draws_based_on_', myGCM, '_and_', myRCP,'_and_startyear_of_',yearStart,'.csv',sep=''))

#test_output <- foreach(HUC_index = 1:nrow(big_bucket), .combine = 'rbind') %do% {
#  summarise.by.HUC(HUC_index, var_index = 3, method = 'par_and_proc', summary_function = sd)
#}

#test_output_df <- tbl_df(cbind(big_bucket[,1,1,1], test_output))
#names(test_output_df) <- c('HUC12', 'par_sd')

#write.csv(test_output_df, paste('SD_3PG_runs_for_', variation_type, '_using_', nsamples, '_paramter_draws_and_', njitters, '_process_error_draws_based_on_', myGCM, '_and_', myRCP,'_and_startyear_of_',yearStart,'.csv',sep=''))
}





