prepare_update_states <- function(plotnum, startyear, nmonths, start_age, endyear){ 
  #rm(list = ls())
  #run_name = 'EnKF_trial_1'
  #restart_chain <- c('Duke_for_EnKF_1.1.2018-03-30.15.11.15.Rdata')
  #load(paste(working_directory,'/chains/',restart_chain,sep=''))
  #priors_file = 'default_priors.csv' ## fix (need to look at this)
  plots = FALSE
  single_GCM = TRUE
  fr_model = 1  # 1 = estimate FR for each plot, 2 = empirical FR model
  FR_fert_assumption = 0 #0 = assume fertilization plots have FR = 1, 1 = do not assume fertilization plots have FR = 1
  use_fol = TRUE  #TRUE= use allometric estimates of foliage biomass in fitting
  use_dk_pars = 1  #0 = do not use 3 specific parameters for the Duke site, 1 = use the 3 specific parameters
  nstreams = 19
  state_space = 1
  plotFR = NA
  windows_machine = FALSE
  
  #--OTHER INFO (WON'T CHANGE UNLESS MODIFY MODEL)-----
  npars_used_by_fortran = 48
  noutput_variables = 67
  process_model_pars = 51
  # npars calculated a few lines below, length(priors_in$parnames)

  
  #---- ENTER THE FORTRAN LIBRARY NAMES HERE ----------
  if(windows_machine){
    code_library_plot = paste(working_directory,'/source_code/r3pg_interface.dll',sep='')
  }else{
    code_library_plot = paste(working_directory,'/source_code/r3pg_interface.so',sep='')
  }
  
  setwd(paste(working_directory,'/scripts/',sep=''))
  source('prepare_obs.R')
  #source('prepare_state_space_obs.R')
  #source('assign_control_plots.R')
  source('prepare_met.R')
  #source('initialize_pars.R')
  #source('init_state_space.R')
  
  PARAMETER_UNCERT=FALSE
  PROCESS_UNCERT=FALSE
  HOLD_CO2_PARAMETERS=FALSE
  HOLD_NONCO2_PARAMETERS=FALSE
  rcp = 85

  adjust_rain = 1
  adjust_fert = 0
  adjust_CO2 = 0
  climate_data_index = 1
  
  setwd(working_directory) 
  
  #----------------------------------------------------
  
  
  priors_in = read.csv(paste(working_directory,'/priors/',priors_file,sep=''))
  npars = length(priors_in$parnames)
  ntotalpars = 125 ############# put this in is as fix for parameter_uncert = false, which was causing code to break
  # priormatrix = matrix(NA,npars,6)
  # priormatrix[,1] = priors_in$initial_value
  # priormatrix[,2] = priors_in$dist_par1
  # priormatrix[,3] = priors_in$dist_par2
  # priormatrix[,4] = priors_in$dist_type
  # priormatrix[,5] = priors_in$fit_par
  # priormatrix[,6] = priors_in$par_group
  
  parnames = priors_in$parnames
  
  #---  PREPARE OBSERVATIONS ---------------------------
  obs_list = prepare_obs(obs_set,FR_fert_assumption,use_fol)
  plotlist = obs_list$plotlist
  StudyName = obs_list$StudyName
  Treatment = obs_list$Treatment
  nplots= obs_list$nplots
  observations= obs_list$observations
  initdata= obs_list$initdata
  met_in = obs_list$met_in
  co2_in = obs_list$co2_in
  use_fol_state = obs_list$use_fol_state
  #-------------------------------------------------
  
  #-----INDEXING FOR THE PARAMETER VECTOR--------------------------
  # this helps speed up the analysis -----------------------------
  co2_in = read.csv(paste(input_directory,'/CO2/CO2_Concentrations_from_CMIP5_1950-2095.csv', sep = ""))
  
  #-----TURN OFF HARDWOOD SIMULATION (0 = HARDWOODS ARE SIMULATED)-------
  exclude_hardwoods = array(1,dim=nplots)
  exclude_hardwoods[which(initdata$PlotID >= 40000 & initdata$PlotID < 42000)]=0
  
  
  dyn.load(code_library_plot)
  
  # -------------------------------------------- the function -----------------------
  if(rcp == 85){
    load(paste(output_directory,'/Duke_all_climate_models_RCP85.Rdata',sep=''))
  }else if(rcp == 45){
    load(paste(output_directory,'/Duke_all_climate_models_RCP45.Rdata',sep=''))
  }
  
  final_pdf = paste(output_directory,'/',run_name,'.pdf',sep='')
  HOLD_CO2 = FALSE
  
  outfile = paste(output_directory,'/',run_name,'.Rdata',sep='')
  if(single_GCM){
    clim_model = c(3)
  }else{
    clim_model = seq(1,20,1)
  }
  
  

  
  print(all_studies)
  #----------------------------------------------------
  Age = -99
  lai  = -99
  stem  = -99
  stem_density = -99
  coarse_root= -99
  fine_root = -99
  fol= -99
  total= -99
  fSW= -99
  ET= -99
  Total_Ctrans= -99
  GPP= -99
  runoff= 0
  WUE_ctrans= -99
  WUE_ET= -99
  
  LAI_quant = array(NA,dim=c(3))
  stem_quant = array(NA,dim=c(3))
  stem_density_quant = array(NA,dim=c(3))
  coarse_root_quant = array(NA,dim=c(3))
  fine_root_quant = array(NA,dim=c(3))
  fol_quant = array(NA,dim=c(3))
  total_quant = array(NA,dim=c(3))
  fSW_quant = array(NA,dim=c(3))
  ET_quant = array(NA,dim=c(3))
  Total_Ctrans_quant = array(NA,dim=c(3))
  GPP_quant = array(NA,dim=c(3))
  runoff_quant = array(NA,dim=c(3))
  WUE_ctrans_quant = array(NA,dim=c(3))
  WUE_ET_quant = array(NA,dim=c(3))
  
  
  
  load(paste(working_directory,'/chains/',restart_chain,sep=''))
  median_pars = rep(NA,ntotalpars)
  for(p in 1:ntotalpars){
    median_pars[p] = median(accepted_pars_thinned_burned[,p])
  }
  new_pars = median_pars
  
  pars = new_pars[1:npars_used_by_fortran]
  FR = new_pars[npars+ plotnum] # fixed from hard-coded "77" to npars
  # delete this loop (m in nmodels)
  
  clim_model = clim_model
  
  #--- CREATE CLIMATE INPUT ARRAYS --------------------------------
  if(climate_data_index == 1){
    curr_met_in_frost = met_in_frost[which(met_in_model == clim_model),]
    curr_met_in_pr =  met_in_pr[which(met_in_model == clim_model),]
    curr_met_in_tasmax =  met_in_tasmax[which(met_in_model == clim_model),]
    curr_met_in_tasmin = met_in_tasmin[which(met_in_model == clim_model),]
    curr_met_in_rsds =  met_in_rsds[which(met_in_model == clim_model),]
    
    
    met = array(NA,c(6))
    yrnum = (startyear+start_age) - 1950 # years since 1950 + start_age
    endno = (((yrnum*12)+1)+nmonths)
    tmax = curr_met_in_tasmax[((yrnum*12)+2):endno]
    tmin = curr_met_in_tasmin[((yrnum*12)+2):endno]
    rain = curr_met_in_pr[((yrnum*12)+2):endno]
    solar = (curr_met_in_rsds[((yrnum*12)+2):endno])/1000000*86400
    frost = curr_met_in_frost[((yrnum*12)+2):endno]
    
    c02start = startyear # used to be startyear + start_age # WHY??
    co2plot = co2_in[which(co2_in$Year >= c02start & co2_in$Year <= endyear), ]
    if (rcp == 45){
      co2 = rep(co2plot[,2],each=12 )
    }
    if (rcp == 85){
      co2 = rep(co2plot[,3],each=12 )
    }
    
    met = array(NA,dim=c(nplots,6,length(tmax)))
    met[,1,] = t(tmin)  #Tmin
    met[,2,] = t(tmax)  # Tmax
    met[,3,] = t(rain) # Rain
    met[,3,] = met[,3,]*adjust_rain
    met[,4,] = t(solar) # SolarRad
    met[,5,] = t(frost) # FrostDays
    met[,6,] = t(co2)
    met[is.na(met)] <- 0
  }else{
    met_tmp = prepare_met(met_in,initdata,mo_start_end,co2_in,nplots,nmonths,months,years)
    met = array(NA,dim=c(nplots,6,length(met_tmp$tmin[1,])))
    met[,1,] = met_tmp$tmin
    met[,2,] = met_tmp$tmax
    met[,3,] = met_tmp$precip
    met[,4,] = met_tmp$ra
    met[,5,] = met_tmp$frost
    met[,6,] = met_tmp$co2
  }
  
  tmp_initdata = initdata[plotnum,]
  
  PlotID = initdata[plotnum,1]
  LAT_WGS84=initdata[plotnum,3]
  ASW_min = initdata[plotnum,8]
  ASW_max=initdata[plotnum,9]
  SoilClass = initdata[plotnum,10]
  SI = initdata[plotnum,11]
  Mean_temp = initdata[plotnum,26]
  
  
  PlantedYear = 0
  PlantedMonth = 0
  
  InitialYear = startyear # delete initial year
  InitialMonth = 1
  
  WFi=initdata[plotnum,13]
  WSi=initdata[plotnum,14]
  WRi=initdata[plotnum,15]
  WCRi=initdata[plotnum,32]
  
  WFi_H = initdata[plotnum,22]
  WSi_H = initdata[plotnum,23]
  WRi_H = initdata[plotnum,24]
  StemNum = initdata[plotnum,6] # PlantDensityHa
  
  #READ IN SITE DATA FROM FILE BASED ON plotlist
  Lat = LAT_WGS84
  ASWi = ASW_max
  MaxASW = ASW_max
  MinASW = ASW_min
  SoilClass=SoilClass
  
  #if(!is.na(plotFR)){
  #  FR = plotFR
  #}else{
  #  FR = 1/(1+exp((new_pars[49] + new_pars[50]*Mean_temp-new_pars[51]*SI)))
  #}
  
  if(initdata[plotnum,12] == 1) { FR = 1}
  
  
  tmp_site_index = 0
  if(PlotID > 40000 & PlotID < 41000 & use_dk_pars == 1){
    tmp_site_index = 1
  }
  
  SLA = 3.5754 + (5.4287 - 3.5754) * exp(-log(2) * (start_age / 5.9705)^2)
  SLA_h = 16.2
  #PlantedYear # where it left off in last time-step
  site_in = c(PlantedYear, #PlantedYear # where it left off in last time-step
              PlantedMonth, #"PlantedMonth"
              InitialYear, #"InitialYear"
              InitialMonth, #"InitialMonth"
              start_age,
              WFi , #"WFi" # can take LAI or WFi
              WRi, #"WRi"
              WSi, #"WSi"
              StemNum, #"StemNoi"
              ASWi = ASW_max, #"ASWi"
              Lat, #"Lat"
              FR, #"FR"
              SoilClass, #"SoilClass"
              MaxASW, #"MaxASW"
              MinASW, #"MinASW"
              TotalMonths = 1,
              WFi_H = 0.001,
              WSi_H = 0.001,
              WRi_H = 0.001,
              WCRi, # 
              IrrigRate = 0.0,
              Throughfall = 1.0,
              tmp_site_index,  
              WCRi_H = 0.0,
              Wbud_H = 0.0,
              LAI = -99, #
              LAI_h = 0.01
  )
  
  site = array(site_in)
  nosite = length(site_in)  # LENGTH OF SITE ARRAY
  
  
  #THIS DEALS WITH THINNING BASED ON PROPORTION OF STEMS REMOVED
  thin_event <- array(data = rep(0),dim = c(plotnum, nmonths))
  
  output_dim = noutput_variables  # NUMBER OF OUTPUT VARIABLES
  nomet = 6  # NUMBER OF VARIABLES IN METEROLOGY (met)
  nopars = length(pars)
  
  if(PlotID > 40000 & PlotID < 41000){
    pars[19] = new_pars[48]
  }
  return(list(output_dim = output_dim, pars = pars, site = site, nopars = nopars, nosite = nosite, median_pars = median_pars, met = met, thin_event = thin_event, exclude_hardwoods = exclude_hardwoods, ASW_max = ASW_max))
}    
