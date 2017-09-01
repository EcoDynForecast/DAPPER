### 3-PG Regional Runs for Quinn, Incorporating Parameter Uncertainty
rm(list = ls())
if (!"compiler" %in% installed.packages()) install.packages("compiler")
if (!"doParallel" %in% installed.packages()) install.packages("doParallel")

library(compiler)
library(doParallel)

enableJIT(1)

output_location = "/Users/quinn/Downloads/"
load("/Users/quinn/Dropbox (VTFRS)/Research/DAPPER/chains/SS_val6.1.2017-08-31.08.47.08.Rdata")
code_library = "/Users/quinn/Dropbox (VTFRS)/Research/DAPPER/source_code/r3pg_interface.so"
CO2 <- read.csv('/Users/quinn/Dropbox (VTFRS)/Research/DAPPER_inputdata/CO2/CO2_Concentrations_from_CMIP5_1950-2095.csv')
Soils <- read.csv('/Users/quinn/Documents/PINEMAP_big_files/parameter_run/Soil_Inputs_LPNR_Clipped_and_Imputed_v4.csv')
GCMs <- c("bcc-csm1-1-m","bcc-csm1-1","BNU-ESM","CanESM2","CCSM4","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-ESM2G","GFDL-ESM2M","HadGEM2-CC365","HadGEM2-ES365","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC-ESM-CHEM","MIROC-ESM","MIROC5","MRI-CGCM3","NorESM1-M", "metdata") # metdata is the Idaho reference dataset, indexed by rcpCase [1]
rcpCases <- c('baseline', 'rcp45','rcp85')

## Critical Parameters -- USER INPUT REQUIRED ####-----------------------------
useParallel = TRUE
HOLD_CO2 = FALSE
test4HUCS = FALSE
use_readr = TRUE
precipModifier <- 1.0 # No toggle to precip
MaxFR <- FALSE # No fertilization
# Number of samples from parameter chain (will be ignored in the process_only case)
sample_size <- 4
# Number of processors to be used
nprocessors <- 2 # Make sure to change this to the max processors!
rotationAge <- 25


PAR_UNCERT = TRUE
PROCESS_UNCERT = FALSE
CLIMATE_UNCERT = FALSE
variation_type <- 'median_climate'
yearStart_list <- c(1985)
myRCP_list <- rcpCases[3]



if(CLIMATE_UNCERT){
myGCM_list <- c(GCMs[1:20],GCMs[1:20])
}else{
myGCM_list <- c(GCMs[3],GCMs[3])
}

if(PAR_UNCERT | PROCESS_UNCERT){
  sample_size <- sample_size
}else{
  sample_size <- 1
}

if(use_readr){
  if (!"readr" %in% installed.packages()) install.packages("readr")
  library(readr)
}

## 3-PG Wrapper ####-----------------------------------------------------------
run.3PG.quick <- function(myHUC, yearStart, yearEnd, startAge = 2, precipModifier = 1, MaxFR = FALSE, curr_pars,base_year = 1950, curr_FR_pars, FR_pars,HOLD_CO2 = FALSE,PROCESS_UNCERT = TRUE, ...){ 
  
  pars <- curr_pars[1:npars_used_by_fortran]
  # Reference Parameters
  HUCIndex <- which(frostDays[, 1] == myHUC)[1] # For the Climate Variables
  HUCIndexSoils <- which(Soils[, 1] == myHUC)[1] # For the Soils variables
  
  ## Extracting Climate
  yearStartPosition <- (yearStart+startAge - base_year) * 12 + 1 # Selecting the first column to draw from
  nomonths <- (yearEnd - (yearStart+startAge) + 1) * 12
  # CO2
  if(HOLD_CO2 == TRUE){
    yearStart = yearStart_list[r]
    yearEnd <- yearStart+startAge + rotationAge - 1
    tmpCO2 <- CO2$CO2_Concentration_RCP85[is.element(CO2$Year, 1985:(1985 + rotationAge - 1))] # RCP8.5 should match historical readings in your timeframe
    monthlyCO2 <- c(matrix(rep(tmpCO2, 12), nrow = 12, byrow = T))
  }else{
    tmpCO2 <- CO2$CO2_Concentration_RCP85[is.element(CO2$Year, (yearStart+startAge):yearEnd)] # RCP8.5 should match historical readings in your timeframe
    monthlyCO2 <- c(matrix(rep(tmpCO2, 12), nrow = 12, byrow = T))
  } 
  # All Climate
  col_index <- c(1 + yearStartPosition:(yearStartPosition +nomonths-1))
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
    #tmpFR = 0.8
  }
  
  # Initialization Inputs (+ Soil) -----------------------------------------------------
  site_in <- c(yearStart, # PlantedYear
               1, # PlantedMonth
               yearStart, # InitialYear
               1, # InitialMonth
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
               tmp_site_index = 0.0,
               WCRi_H = 0.0,
               Wbud_H = 0.0,
               LAI = 0.5,
               LAI_h = 0.01) # WCRi
  site <- array(site_in)
  site[is.na(site)] <- 0
  nosite = length(site)
  
  nmonths = nomonths
  
  thin_event = array(0,dim=c(nmonths))
  
  
  output0 = array(NA,dim=c(nmonths,3))

  for(mo in 1:nmonths){

    tmp=.Fortran( "r3pg_interface",
                  output_dim=as.integer(output_dim),
                  met=as.double(t(met[mo,])),
                  pars=as.double(pars),
                  site = as.double(site),
                  thin_event = as.double(thin_event[mo]),
                  out_var=as.double(array(0,dim=c(1,output_dim))),
                  nopars=as.integer(nopars),
                  nomet=as.integer(dim(met)[2]),
                  nosite = as.integer(nosite),
                  nooutputs=as.integer(output_dim),
                  nomonths_plot=as.integer(1),
                  nothin = 1,
                  exclude_hardwoods = as.integer(1),
                  mo_start_end = as.integer(c(1,1)),
                  nmonths = 1
    )
    #--------------------------------------------------------------------------------
    
    output=array(tmp$out_var, dim=c(1,output_dim))
    
    if(length(which(is.nan(output))) == 0){
      
      if(output[2] == 12){
        site[3] = output[1]+1 #InitialYear
        site[4] = 1  #InitialMonth
      }else{
        site[3] = output[1] #InitialYear
        site[4] = output[2]+1  #InitialMonth	
      }
      site[5] = output[3] + (1.0/12.) #StartAge
      
      if(is.nan(output[26]) | is.na(output[26]) | output[26] < 0 ) {
        site[26] = 0.5 * SLA * 0.1
        site[8] = 1
        site[7] = 0.5
        site[9] = 1500
        site[20] = 0.30
        site[6] = 0.5
        site[27] = max(rnorm(1,output[9],curr_pars[52]),0.0)  #Hardwood LAI
        site[25] = output[26] #Hardwood Bud
        site[18] = rnorm(1,output[10],curr_pars[57]) #WS_H 
        site[24] = output[11]  #WCR_h
        site[19] = rnorm(1,output[12],curr_pars[55]) #WR_H
        site[10] = output[14] # ASW
        site[17] = output[23] #WF_H	
      }else{ 
        if(PROCESS_UNCERT){
          site[26] = max(rnorm(1,output[4],curr_pars[52]),0.1) #LAI
          site[8] = rnorm(1,output[5],(1.3+new_pars[53] +output[5]*new_pars[64]))  #WS
          site[20] = rnorm(1,output[6],curr_pars[54])   #WCR
          site[7] = rnorm(1,output[7],curr_pars[55])  #WRi
          site[9] = rnorm(1,output[8],curr_pars[56]) #StemNo
          site[6] = output[22] #WFi
          site[27] = max(rnorm(1,output[9],curr_pars[52]),0.0)  #Hardwood LAI
          site[25] = output[26] #Hardwood Bud
          site[18] = rnorm(1,output[10],curr_pars[57]) #WS_H 
          site[24] = output[11]  #WCR_h
          site[19] = rnorm(1,output[12],curr_pars[55]) #WR_H
          site[10] = output[14] # ASW
          site[17] = output[23] #WF_H	
        }else{
          site[26] = output[4]
          site[8] = output[5]
          site[20] = output[6]
          site[7] = output[7]
          site[9] = output[8]
          site[6] = output[22] #WFi    
          site[27] = output[9] #  #Hardwood LAI
          site[25] = output[26] #Hardwood Bud
          site[18] = output[10] #WS_H 
          site[24] = output[11]  #WCR_h
          site[19] = output[12] #WR_H
          site[10] = output[14] # ASW
          site[17] = output[23] #WF_H	
        }
      }
      
      total =  site[6] + site[8] + site[20] + site[7]
      output0[mo,1] =output[3]
      output0[mo,2] = site[8]
      output0[mo,3] =  total
      #print(c(mo,output[3],site[8],total))
    }else{
      site_in <- c(yearStart, # PlantedYear
                   1, # PlantedMonth
                   yearStart, # InitialYear
                   1, # InitialMonth
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
                   tmp_site_index = 0.0,
                   WCRi_H = 0.0,
                   Wbud_H = 0.0,
                   LAI = 0.5,
                   LAI_h = 0.01) # WCRi
      
      site = array(site_in)
      
      mo = 1
    }
    #-------------------------------------------------------------------------------
}
  age_25 <- which(output0[,1] >= 24.9 & output0[,1] < 26.0)
  output <- colMeans(output0[age_25,])
  
  # Final Asembly and Output
  finalOutput <- c(myHUC, colMeans(output0[age_25,]))
  
  finalOutput
}

## Model Runs ####-------------------------------------------------------------

nsamples <- sample_size
niter <- dim(accepted_pars_thinned_burned)[1]
npars <- 80
npars_used_by_fortran <- 48
set.seed(2016)
pars_sample_index <- sample(seq(1,niter,1),nsamples,replace = TRUE)


for(r in 1:length(myRCP_list)){
  
  myRCP = myRCP_list[r]
  myGCM = myGCM_list[r]
  yearStart = yearStart_list[r]
  yearEnd <- yearStart + rotationAge
  
  print(paste('Run ',r,': starting RCP ',myRCP,' for model ',myGCM,' for start year ',yearStart,sep="")) 
  
  npars <- 80
  curr_pars <- rep(NA,npars)
  for(i in 1:npars){
    curr_pars[i] <- median(accepted_pars_thinned_burned[,i]) # Taking median here, modify for bootstrapping
  }
  
  npars_used_by_fortran <- 48
  pars <- curr_pars[1:npars_used_by_fortran]
  
  noutput_variables <- 67
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
  if(use_readr){
  frostDays <- read_csv(paste('/Users/quinn/Documents/PINEMAP_big_files/parameter_run/Consolidated_Files_Imput/frostDays_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
  rain <- read_csv(paste('/Users/quinn/Documents/PINEMAP_big_files/parameter_run/Consolidated_Files_Imput/rain_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
  solarRad <- read_csv(paste('/Users/quinn/Documents/PINEMAP_big_files/parameter_run/Consolidated_Files_Imput/solarRad_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
  tMax <- read_csv(paste('/Users/quinn/Documents/PINEMAP_big_files/parameter_run/Consolidated_Files_Imput/tMax_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
  tMin <- read_csv(paste('/Users/quinn/Documents/PINEMAP_big_files/parameter_run/Consolidated_Files_Imput/tMin_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
  }else{
    frostDays <- read.csv(paste('/Users/quinn/Documents/PINEMAP_big_files/parameter_run/Consolidated_Files_Imput/frostDays_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
    rain <- read.csv(paste('/Users/quinn/Documents/PINEMAP_big_files/parameter_run/Consolidated_Files_Imput/rain_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
    solarRad <- read.csv(paste('/Users/quinn/Documents/PINEMAP_big_files/parameter_run/Consolidated_Files_Imput/solarRad_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
    tMax <- read.csv(paste('/Users/quinn/Documents/PINEMAP_big_files/parameter_run/Consolidated_Files_Imput/tMax_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
    tMin <- read.csv(paste('/Users/quinn/Documents/PINEMAP_big_files/parameter_run/Consolidated_Files_Imput/tMin_', RCP_index,'_', myGCM,'_imput.csv', sep = ''))
  }
  
  myHUCsFinal <- intersect(Soils$HUC12, frostDays$HUC12)
  if(test4HUCS){
  myHUCsFinal <- myHUCsFinal[which(myHUCsFinal == 31101030505 | myHUCsFinal == 30601050202 | myHUCsFinal == 111401070404 | myHUCsFinal == 20802031301)]
  }
  
  cl <- makeCluster(nprocessors)  
  registerDoParallel(cl)
  print('starting sample loop')
  print(nsamples)
  
  # for(ns in 1:nsamples){
  final_output_GCM_RCP <- foreach(ns = 1:nsamples, .export = ls()) %dopar% { # Will need to add objects,
    #final_output_GCM_RCP <- foreach(ns = 1:nsamples) %do% { # Will need to add objects,packages as needed
    # Load Fortran code
    par_ns_index <- pars_sample_index[ns]
    curr_pars <- accepted_pars_thinned_burned[par_ns_index,]
    if(PAR_UNCERT == FALSE){
      curr_pars <- rep(NA,npars)
      for(i in 1:npars){
        curr_pars[i] <- median(accepted_pars_thinned_burned[,i]) # Taking median here, modify for bootstrapping
      }
    }
    
    FR_pars = c(curr_pars[50],curr_pars[51])
    
    bucket <- array(NA, dim = c(length(myHUCsFinal), 3)) # columns for HUC12, WF, WS, WR, CoarseRoots, and TotalBiomass
    tmp_bucket <- array(NA, dim = c(length(myHUCsFinal), 3)) # columns for HUC12, WF, WS, WR, CoarseRoots, and TotalBiomass
    tmp <- matrix(NA, nrow = length(myHUCsFinal), ncol = 4)
    
    #if(!useParallel){
      dyn.load(code_library)
    #}
    for(HUCindex in 1:length(myHUCsFinal)){
      tmp[HUCindex,] <- run.3PG.quick(myHUC = myHUCsFinal[HUCindex], yearStart = yearStart, yearEnd = yearEnd, startAge = 2, precipModifier = 1, MaxFR = MaxFR, base_year = base_year, curr_pars = curr_pars,FR_pars = FR_pars, HOLD_CO2 = HOLD_CO2, PROCESS_UNCERT = PROCESS_UNCERT) 
    }
    
    tmp_bucket[,1] <- tmp[,1] 
    tmp_bucket[,2] <- tmp[,3]
    tmp_bucket[,3] <- tmp[,4]
    
    bucket = tmp_bucket
  }
  stopCluster(cl)
  
  big_bucket <- array(unlist(final_output_GCM_RCP), dim = c(dim(final_output_GCM_RCP[[1]]),length(final_output_GCM_RCP)))
  
  rm(frostDays, rain, solarRad, tMax, tMin)
  
  save.image(file = paste(output_location,'/3PG_runs_for_', variation_type, '_using_', nsamples, '_based_on_', myGCM, '_and_', myRCP,'_and_startyear_of_',yearStart,'_FRset_',MaxFR,'_HOLDCO2_',HOLD_CO2,'.Rdata', sep = ''))      
}




