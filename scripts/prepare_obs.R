prepare_obs <- function(obs_set,FR_fert_assumption,use_fol){
  #-----READ IN OBSERVATED DATA------------------------

  initdata = NULL
  for(s in 1:length(all_studies)){
    d = read.csv(paste(input_directory,all_studies[s],'_plotlist.csv',sep=''))
    d$Treatment = as.factor(d$Treatment)
    d$StudyName = as.factor(d$StudyName)
    if(all_studies[s] != '/Duke/TIER4_Duke' & all_studies[s] != '/NC2/TIER4_NC2'  & all_studies[s] != '/Waycross/TIER4_Waycross' & all_studies[s] != '/SETRES/TIER4_SETRES'){
      d$Initial_LAI = -99
      d$Initial_LAI_code = -99  
      d$Initial_WR_code = -99
    }
    
    initdata = rbind(initdata,d)
  }
  
  initdata$Initial_WR[which(initdata$PlotID > 42000 & initdata$PlotID < 43000)] = initdata$Initial_WR[which(initdata$PlotID > 42000 & initdata$PlotID < 43000)]/0.54
  initdata = initdata[which(initdata$PlotID != 10007  &
                              initdata$PlotID != 10052  &
                              initdata$PlotID != 10157  &
                              initdata$PlotID != 10269  &
                              initdata$PlotID != 10456  &
                              initdata$PlotID != 10536  &
                              initdata$PlotID != 10776  &
                              initdata$PlotID != 10806  &
                              initdata$PlotID != 10813  &
                              initdata$PlotID != 10876  &
                              initdata$PlotID != 11302  &
                              initdata$PlotID != 11311  &
                              initdata$PlotID != 11392), ]
  
  initdata = initdata[which(initdata$MeasureLength > 1 | initdata$MeasureLength == -99), ]
  
  #initdata = initdata[which(initdata$ThinTreatment == 1 | initdata$ThinTreatment == 0), ]

  initdata = data.frame(PlotID = initdata$PlotID,SiteID = initdata$SiteID,LAT_WGS84=initdata$LAT_WGS84,
                        Planting_year = initdata$Planting_year,PlantMonth = initdata$PlantMonth,
                        PlantDensityHa = initdata$PlantDensityHa,Initial_ASW = initdata$Initial_ASW, ASW_min = initdata$ASW_min,
                        ASW_max=initdata$ASW_max,SoilClass = initdata$SoilClass,SI = initdata$SI,FR = initdata$FR,
                        Initial_WF = initdata$Initial_WF,Initial_WS = initdata$Initial_WS,Initial_WR = initdata$Initial_WR,
                        DroughtLevel = initdata$DroughtLevel,DroughtStart = initdata$DroughtStart, FertFlag=initdata$FertFlag,CO2flag = initdata$CO2flag,
                        CO2elev = initdata$CO2elev, ControlPlotID = initdata$ControlPlotID,
                        Initial_WF_H = initdata$Initial_WF_H,Initial_WS_H = initdata$Initial_WS_H,Initial_WR_H = initdata$Initial_WR_H, LONG_WGS84 = initdata$LONG_WGS84,
                        Mean_temp = initdata$mean_annual_temp,mean_precip = initdata$mean_annual_precip, MatchedFRPlotID = initdata$matched_FR_plotid, 
                        InitYear = initdata$Initial_year,InitMonth = initdata$Initial_month, StartAge = initdata$Initial_age,Initial_WCR = initdata$Initial_WCR,
                        IrrFlag = initdata$IrrFlag,IrrLevel = initdata$IrrLevel,Initial_LAI = initdata$Initial_LAI,Initial_LAI_code = initdata$Initial_LAI_code,
                        Initial_WR_code = initdata$Initial_WR_code, StudyName = initdata$StudyName, Treatment = initdata$Treatment, ThinTreatment = initdata$ThinTreatment)
  

  initdata$SoilClass = 1.0 

  
  if(FR_fert_assumption == 1){
    initdata$FR[which(initdata$PlotID > 20000 & initdata$PlotID < 30000)] = -99
    initdata$FR[which(initdata$PlotID > 30000 & initdata$PlotID < 40000)] = -99
    initdata$FR[which(initdata$PlotID > 40000 & initdata$PlotID < 41000)] = -99
    initdata$FR[which(initdata$PlotID > 41000 & initdata$PlotID < 42000)] = -99
    initdata$FR[which(initdata$PlotID > 42000 & initdata$PlotID < 43000)] = -99
    initdata$FR[which(initdata$PlotID > 43000 & initdata$PlotID < 44000)] = -99
  }else{
    initdata$FR[which(initdata$PlotID > 30000 & initdata$PlotID < 40000)] = -99
    #initdata$FR[which(initdata$PlotID > 40000 & initdata$PlotID < 41000)] = -99
    initdata$FR[which(initdata$PlotID > 41000 & initdata$PlotID < 42000)] = -99
  }
  
  observations = NULL
  for(s in 1:length(all_studies)){
    d = read.csv(paste(input_directory,all_studies[s],'_obs.csv',sep=''))
    d$Treatment = as.factor(d$Treatment)
    d$StudyName = as.factor(d$StudyName)
    d$PlotSizeHa = as.numeric(d$PlotSizeHa)
    d$ind_removed = as.numeric(d$ind_removed)
    d$ind_removed_prop = as.numeric(d$ind_removed_prop)
    if(all_studies[s] == '/FIA/VA_FIA'){  
      observations = rbind(observations,d[,1:55])
    }else{
      observations = rbind(observations,d) 
    }
  }
  
  observations$ROOT_FINE[which(observations$PlotID > 42000 & observations$PlotID < 43000 & observations$ROOT_FINE>-99)]  = observations$ROOT_FINE[which(observations$PlotID > 42000 & observations$PlotID < 43000 & observations$ROOT_FINE>-99)]/0.54
  observations = observations[which(observations$PlotID != 10007  &
                                      observations$PlotID != 10052  &
                                      observations$PlotID != 10157  &
                                      observations$PlotID != 10269  &
                                      observations$PlotID != 10456  &
                                      observations$PlotID != 10536  &
                                      observations$PlotID != 10776  &
                                      observations$PlotID != 10806  &
                                      observations$PlotID != 10813  &
                                      observations$PlotID != 10876  &
                                      observations$PlotID != 11302  &
                                      observations$PlotID != 11311  &
                                      observations$PlotID != 11392), ]
  
  observations = data.frame(PlotID = observations$PlotID,MonthMeas = observations$MonthMeas,YearMeas = observations$YearMeas,
                            AgeMeas=observations$AgeMeas,FOL=observations$FOL,WOODY=observations$WOODY,ROOT_FINE=observations$ROOT_FINE,
                            Nha=observations$Nha,ind_removed=observations$ind_removed,WFest_sd=observations$WFest_sd,
                            WSest_sd=observations$WSest_sd,GEP = observations$GEP,ET = observations$ET, Ctrans = observations$Ctrans,Ctrans_H = observations$Ctrans_H,
                            WOODY_H = observations$WOODY_H,FOL_H = observations$FOL_H,
                            LAI = observations$LAI,FOL_PROD = observations$FOL_PROD,FOL_PROD_H = observations$FOL_PROD_H,
                            FOL_PROD_TOTAL = observations$FOL_PROD_TOTAL, ROOT_COARSE = observations$ROOT_COARSE, LAI_H = observations$LAI_H,LAI_TOTAL = observations$LAI_TOTAL,
                            GEP_sd = observations$GEP_sd, ET_sd = observations$ET_sd,WRest_sd = observations$WRest_sd,Ctrans_sd = observations$Ctrans_sd, ROOT_PROD_TOTAL = observations$ROOT_PROD_TOTAL)
  
  observations$ROOT_FINE[which(observations$PlotID < 40000)] = -99
  observations$FOL[which(observations$PlotID > 42000 & observations$PlotID < 43000)] = -99
  if(use_fol == FALSE){
    observations$FOL[which(observations$PlotID < 42000 | observations$PlotID >= 43000)] = -99
  }
  
  
  observations$FOL[which(observations$PlotID > 100000)] = -99
  observations$ROOT_COARSE[which(observations$PlotID > 100000)] = -99
  observations$ROOT_FINE[which(observations$PlotID > 100000)] = -99
  observations$WOODY_H[which(observations$PlotID > 100000)] = -99
  observations$FOL_H[which(observations$PlotID > 100000)] = -99 
  observations = observations[which(observations$AgeMeas <= 30),]
  
  observations$ind_removed = -99
  
  #observations$LAI_TOTAL = -99

  
  met_in = NULL
  for(s in 1:length(all_studies)){
    d = read.csv(paste(input_directory,all_studies[s],'_met.csv',sep=''))
    met_in = rbind(met_in,d)
  }
  
  co2_in = read.csv(paste(input_directory,'/CO2/CO2_Concentrations_from_CMIP5_1950-2095.csv',sep=''))
  
  initdata = initdata[which(!is.na(initdata$ASW_max) & initdata$ASW_max > 0.0),]
  
  observations = observations[observations$PlotID %in% initdata$PlotID, ]
  

  #-----------------------------------------------------------------
  #  Remove foliage data for treatment plots because based on model
  for(i in 1:length(observations$FOL)){
    fert = initdata$FertFlag[which(initdata$PlotID == observations$PlotID[i])]
    irr = initdata$IrrFlag[which(initdata$PlotID == observations$PlotID[i])]
    drought = initdata$DroughtLevel[which(initdata$PlotID == observations$PlotID[i])]
    if(fert == 1 | irr == 1 | drought < 1){
      observations$FOL[i] = -99        
    }
  } 
  


  #-----ORGANIZE CONTROL PLOT ID (FOR USE WITH EXPERIMENTAL DATA)
  control_plot_index = which(initdata$FertFlag == 0 & (initdata$CO2flag == -99 | initdata$CO2flag == 0)  & initdata$DroughtLevel == 1.0 & initdata$IrrFlag == 0.0)
  control_plot_fert_index = which((initdata$CO2flag == -99 | initdata$CO2flag == 0) & initdata$DroughtLevel == 1.0 & initdata$IrrFlag == 0.0) 
  control_plot_drought_index = which(initdata$FertFlag == 0 & (initdata$CO2flag == -99 | initdata$CO2flag == 0))
  control_plot_co2_index = which(initdata$FertFlag == 0 & initdata$DroughtLevel == 1.0 & initdata$IrrFlag == 0.0)
  control_plot_co2_fert_index = which(initdata$DroughtLevel == 1.0 & initdata$IrrFlag == 0.0)
  control_plot_fert_drought_index = which((initdata$CO2flag == -99 | initdata$CO2flag == 0))
  control_plot_co2_drought_index = which((initdata$FertFlag == 0))
  flux_sites = which(initdata$PlotID == 40001 | initdata$PlotID == 41001)
  duke_nc2 = which(initdata$PlotID > 40000 & initdata$PlotID < 42000)
  duke_nc2_setres_wcross_pinemap = which(initdata$PlotID > 30000 & initdata$PlotID < 44000)
  not_duke_nc2 = which(initdata$PlotID < 40000 | initdata$PlotID >= 42000)
  setres = which(initdata$PlotID > 42000 & initdata$PlotID < 43000)
  pinemap_tier3 = which(initdata$PlotID > 30000 & initdata$PlotID < 40000)
  duke =  which(initdata$PlotID > 40000 & initdata$PlotID < 41000) 
  duke_wcross =  which(initdata$PlotID > 40000 & initdata$PlotID < 41000 | (initdata$PlotID > 43000)) 
  wcross =  which(initdata$PlotID > 43000 & initdata$PlotID < 44000) 
  pinemap_setres = which(initdata$PlotID > 30000 & initdata$PlotID < 40000 | (initdata$PlotID > 42000 & initdata$PlotID < 43000)) 
  pinemap_setres_fluxes = which(initdata$PlotID > 30000 & initdata$PlotID < 40000 | (initdata$PlotID > 42000 & initdata$PlotID < 43000) | (initdata$PlotID == 40001 | initdata$PlotID == 41001))  
  no_drought = which(initdata$DroughtLevel == 1.0 & initdata$IrrFlag == 0.0)
  duke_wcross_setress_nc2_pinemap =  which((initdata$PlotID > 30000 & initdata$PlotID < 40000) | (initdata$PlotID > 40000 & initdata$PlotID < 41000) | (initdata$PlotID >= 43000 & initdata$PlotID < 44000) | (initdata$PlotID > 42000 & initdata$PlotID < 43000) | (initdata$PlotID == 41001)) 
  not_new = which(initdata$PlotID < 50000)
  nc2 = which(initdata$PlotID == 41001)
  IMP = which(initdata$PlotID > 50000 | initdata$PlotID < 20000)
  test = which(initdata$PlotID >= 100000)
  
  all_plots =  seq(1,length(initdata$PlotID),1)
  
  #-----------------------------------------------------------------
  if(is.na(focal_plotID)){
    if(obs_set == 0){
      index = 1
    }else if(obs_set == 1){
      #index = control_plot_index
      #uses set_fitted_plots.R to remove the selected plots
      index = all_plots
    }else if(obs_set == 2){
      index=control_plot_fert_index 
    }else if(obs_set == 3){
      index = control_plot_drought_index 
    }else if(obs_set == 4){
      index = control_plot_co2_index
    }else if(obs_set == 5){
      #index = control_plot_co2_fert_index 
      #uses set_fitted_plots.R to remove the selected plots
      index =  all_plots
    }else if(obs_set == 6){
      #index = control_plot_fert_drought_index
      #uses set_fitted_plots.R to remove the selected plots
      index =  all_plots
    }else if(obs_set == 7){
      #index = control_plot_co2_drought_index
      #uses set_fitted_plots.R to remove the selected plots
      index =  all_plots
    }else if(obs_set == 8){
      index = flux_sites
    }else if(obs_set == 9){
      index =  duke_nc2
    }else if(obs_set == 10){
      index = not_duke_nc2
    }else if(obs_set == 11){
      index =  all_plots
    }else if(obs_set == 12){
      index = setres
    }else if(obs_set == 13){
      index = pinemap_tier3 
    }else if(obs_set == 14){
      index = duke
    }else if(obs_set == 15){
      index = pinemap_setres
    }else if(obs_set == 16){
      index = pinemap_setres_fluxes
    }else if (obs_set == 17){
      index =wcross
    }else if(obs_set == 18){
      index = no_drought
    }else if(obs_set == 19){
      index = all_plots
    }else if(obs_set == 20){
      index = nc2
    }else if(obs_set == 21){
      index = not_new
    }else if(obs_set == 22){
      index = duke_nc2_setres_wcross_pinemap
    }else if(obs_set == 23){
      index = test
    }
  }else{
    
    index = which(initdata$PlotID == focal_plotID )
  }
  
  #-------SELECT WHICH PLOTS TO USE---------------------------------
  plotlist = initdata$PlotID[index]  #[1:64] #initdata$PlotID[1:276] #initdata$PlotID[101:112]
  nplots = length(plotlist)
  observations =  observations[observations$PlotID %in% plotlist, ]
  initdata = initdata[initdata$PlotID %in% plotlist, ]
  
  use_fol_state = array(0,dim=nplots)
  for(plotnum in 1:nplots){
    tmp = observations[which(observations$PlotID == initdata$PlotID[plotnum]),]
    if(length(which(tmp$FOL != -99)) != 0 & length(which(tmp$LAI != -99))==0){
      use_fol_state[plotnum] = 1
    }
  }
  
  return(list(plotlist = plotlist,nplots= nplots,observations= observations,initdata= initdata,met_in=met_in,co2_in=co2_in,
              use_fol_state = use_fol_state))
  
}
