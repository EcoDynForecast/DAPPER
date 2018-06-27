obs_set = 1
#all_studies = c(
#  '/Duke/TIER4_Duke',
#  '/SETRES/TIER4_SETRES'
#)

all_studies = c(
  '/Duke/TIER4_Duke'
)

co2 = read.csv('/Users/quinn/Dropbox (VTFRS)/Research/DAPPER_inputdata/CO2/CO2_Concentrations_from_CMIP5_1950-2095.csv')
input_directory = "/Users/quinn/Dropbox/Research/DAPPER_inputdata/"


initdata = NULL
for(s in 1:length(all_studies)){
  d = read.csv(paste(input_directory,all_studies[s],'_plotlist.csv',sep=''))
  d$Treatment = as.character(d$Treatment)
  d$StudyName = as.character(d$StudyName)
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


initdata$ThinTreatment[which(initdata$PlotID >= 52001 & initdata$PlotID <= 52467)]  = initdata$Plot[which(initdata$PlotID >= 52001 & initdata$PlotID <= 52467)] 
initdata$ThinTreatment[which(initdata$PlotID >= 72001 & initdata$PlotID <= 72076)]  = initdata$Plot[which(initdata$PlotID >= 72001 & initdata$PlotID <= 72076)] 

initdata$Treatment[which(initdata$PlotID >= 52001 & initdata$PlotID <= 52467)]  = as.character(initdata$Plot[which(initdata$PlotID >= 52001 & initdata$PlotID <= 52467)]) 
initdata$Treatment[which(initdata$PlotID >= 72001 & initdata$PlotID <= 72076)]  = as.character(initdata$Plot[which(initdata$PlotID >= 72001 & initdata$PlotID <= 72076)]) 
initdata$Treatment[which(initdata$PlotID >= 10001 & initdata$PlotID < 20000)]  = as.character(initdata$ThinTreatment[which(initdata$PlotID >= 10001 & initdata$PlotID < 20000)]) 
initdata = initdata[which(initdata$ThinTreatment == 1 | initdata$ThinTreatment == 0), ]
initdata$Treatment = as.factor(initdata$Treatment)

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
                      Initial_WR_code = initdata$Initial_WR_code, StudyName = initdata$StudyName, Treatment = initdata$Treatment, ThinTreatment = initdata$ThinTreatment,Plot = initdata$Plot)


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

obs_uncertainity_proportion = 0.01 
#observations$WFest_sd[which(observations$WFest_sd == -99)] = observations$WF[which(observations$WFest_sd == -99)]*obs_uncertainity_proportion
observations$WSest_sd = observations$WOODY*0.025
observations$WRest_sd = observations$ROOT_FINE*0.1
observations$ROOT_COARSE_sd = observations$ROOT_COARSE*obs_uncertainity_proportion
observations$ROOT_FINE_sd= observations$ROOT_FINE*obs_uncertainity_proportion
observations$Nha_sd = observations$Nha*obs_uncertainity_proportion
observations$LAI_sd= observations$LAI*obs_uncertainity_proportion
observations$LAI_H_sd= observations$LAI_H*obs_uncertainity_proportion
observations$WOODY_H_sd= observations$WOODY_H*obs_uncertainity_proportion
#observations$GEP_sd = observations$GEP*0.1
#observations$ET_sd= observations$ET*obs_uncertainity_proportion
#observations$Ctrans_sd= observations$Ctrans*obs_uncertainity_proportion
#observations$Ctrans_H_sd = observations$Ctrans_H*obs_uncertainity_proportion
observations$ROOT_PROD_TOTAL_sd= observations$ROOT_PROD_TOTAL*obs_uncertainity_proportion
observations$FOL_PROD_TOTAL_sd  = observations$FOL_PROD_TOTAL*obs_uncertainity_proportion
observations$LAI_TOTAL_sd   = observations$LAI_TOTAL*obs_uncertainity_proportion


#observations = data.frame(PlotID = observations$PlotID,MonthMeas = observations$MonthMeas,YearMeas = observations$YearMeas,
#                          AgeMeas=observations$AgeMeas,FOL=observations$FOL,WOODY=observations$WOODY,ROOT_FINE=observations$ROOT_FINE,
#                          Nha=observations$Nha,ind_removed=observations$ind_removed,WFest_sd=observations$WFest_sd,
#                          WSest_sd=observations$WSest_sd,GEP = observations$GEP,ET = observations$ET, Ctrans = observations$Ctrans,Ctrans_H = observations$Ctrans_H,
#                          WOODY_H = observations$WOODY_H,FOL_H = observations$FOL_H,
#                          LAI = observations$LAI,FOL_PROD = observations$FOL_PROD,FOL_PROD_H = observations$FOL_PROD_H,
#                          FOL_PROD_TOTAL = observations$FOL_PROD_TOTAL, ROOT_COARSE = observations$ROOT_COARSE, LAI_H = observations$LAI_H,LAI_TOTAL = observations$LAI_TOTAL,
#                          GEP_sd = observations$GEP_sd, ET_sd = observations$ET_sd,WRest_sd = observations$WRest_sd,Ctrans_sd = observations$Ctrans_sd, 
#                          ROOT_PROD_TOTAL = observations$ROOT_PROD_TOTAL,Ctrans_H_sd = observations$Ctrans_H_sd)



observations = data.frame(PlotID = observations$PlotID,
                          MonthMeas = observations$MonthMeas,
                          YearMeas = observations$YearMeas,
                          AgeMeas=observations$AgeMeas,
                          FOL=observations$FOL,
                          WOODY=observations$WOODY,
                          ROOT_FINE=observations$ROOT_FINE,
                          Nha=observations$Nha,
                          ind_removed=observations$ind_removed,
                          WFest_sd=observations$WFest_sd,
                          WSest_sd=observations$WSest_sd,
                          GEP = observations$GEP,
                          ET = observations$ET, 
                          Ctrans = observations$Ctrans,
                          Ctrans_H = observations$Ctrans_H,
                          WOODY_H = observations$WOODY_H,
                          FOL_H = observations$FOL_H,
                          LAI = observations$LAI,
                          FOL_PROD = observations$FOL_PROD,
                          FOL_PROD_H = observations$FOL_PROD_H,
                          FOL_PROD_TOTAL = observations$FOL_PROD_TOTAL, 
                          ROOT_COARSE = observations$ROOT_COARSE, 
                          LAI_H = observations$LAI_H,
                          LAI_TOTAL = observations$LAI_TOTAL,
                          GEP_sd = observations$GEP_sd, 
                          ET_sd = observations$ET_sd,
                          WRest_sd = observations$WRest_sd,
                          Ctrans_sd = observations$Ctrans_sd, 
                          ROOT_PROD_TOTAL = observations$ROOT_PROD_TOTAL,
                          Ctrans_H_sd = observations$Ctrans_H_sd,
                          ROOT_COARSE_sd = observations$ROOT_COARSE_sd,
                          ROOT_FINE_sd= observations$ROOT_FINE_sd,
                          Nha_sd= observations$Nha_sd,
                          LAI_sd = observations$LAI_sd,
                          LAI_H_sd= observations$LAI_H_sd,
                          WOODY_H_sd= observations$WOODY_H_sd,
                          ROOT_PROD_TOTAL_sd= observations$ROOT_PROD_TOTAL_sd,
                          FOL_PROD_TOTAL_sd  = observations$FOL_PROD_TOTAL_sd,
                          LAI_TOTAL_sd   = observations$LAI_TOTAL_sd
                          )


#if((meas_month != 8) & tmp$PlotID[1] > 40000 & tmp$PlotID[1] < 42000){obs[19,plotnum,index]=-99}
#if((meas_month < 6 | meas_month > 9) & tmp$PlotID[1] > 40000 & tmp$PlotID[1] < 42000){obs[16,plotnum,index]=-99}



if(use_ctrans == 0){
  observations$Ctrans = -99
  observations$Ctrans_H = -99
}

if(use_et == 0){
  observations$ET = -99
}

if(use_gep == 0){
  observations$GEP = -99
}

if(use_ctrans_uncert == 0){
  observations$Ctrans_sd = 0.001
  observations$Ctrans_H_sd = 0.001
}

if(use_et_uncert == 0){
  observations$ET_sd = 0.001
}

if(use_gep_uncert == 0){
  observations$GEP_sd = 0.001
}

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
test = which(initdata$PlotID < 50000 | (initdata$PlotID >= 52001 & initdata$PlotID <=52467) | (initdata$PlotID >= 72001 & initdata$PlotID <= 72076))
duke_1plot = which(initdata$PlotID == 40001)
duke_control = which(initdata$PlotID > 40000 & initdata$PlotID < 40001 & initdata$FertFlag == 0 & (initdata$CO2flag == -99 | initdata$CO2flag == 0)) 
duke_control_co2 = which(initdata$PlotID > 40000 & initdata$PlotID < 40001 & initdata$FertFlag == 0) 
duke_control_nfert = which(initdata$PlotID > 40000 & initdata$PlotID < 40001 & (initdata$CO2flag == -99 | initdata$CO2flag == 0)) 
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
    index =  duke_1plot
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
  }else if(obs_set == 24){
    index = duke_control
  }else if(obs_set == 25){
    index = duke_control_nfert
  }else if(obs_set == 26){
    index = duke_control_co2
  }
}else{
  
  index = which(initdata$PlotID == focal_plotID )
}

#-------SELECT WHICH PLOTS TO USE---------------------------------
plotlist = initdata$PlotID[index]  #[1:64] #initdata$PlotID[1:276] #initdata$PlotID[101:112]
nplots = length(plotlist)
observations =  observations[observations$PlotID %in% plotlist, ]
initdata = initdata[initdata$PlotID %in% plotlist, ]

write.csv(observations,file='/Users/quinn/Dropbox (VTFRS)/Research/DAPPER_branch/DAPPER/example_run/example_observations.csv',row.names = FALSE)

write.csv(initdata,file='/Users/quinn/Dropbox (VTFRS)/Research/DAPPER_branch/DAPPER/example_run/example_plots.csv',row.names = FALSE)


####

met_in_list = NULL
for(s in 1:length(all_studies)){
  d = read.csv(paste(input_directory,all_studies[s],'_met.csv',sep=''))
  met_in_list = rbind(met_in_list,d)
}
years = unique(met_in_list$YEAR)
SiteID = unique(initdata$SiteID)

met_in = array(NA, dim=c(length(years)*12*length(SiteID),9))
for(sitenum in 1:length(SiteID)){
for(i in 1:length(years)){
  index = (sitenum-1)*length(years)*12 +  (i-1)*12 + 1
  curr_m = met_in_list[which(met_in_list$SiteID == SiteID[sitenum]),]
  met_in[index:(index+11),1] = SiteID[sitenum]
  met_in[index:(index+11),2] = years[i]
  met_in[index:(index+11),3] = seq(1,12,1)
  met_in[index:(index+11),4] = unlist(curr_m[i,15:26])
  met_in[index:(index+11),5] = unlist(curr_m[i,3:14])
  met_in[index:(index+11),6] = unlist(curr_m[i,27:38])
  met_in[index:(index+11),7] = unlist(curr_m[i,39:50])
  met_in[index:(index+11),8] = unlist(curr_m[i,51:62])
  met_in[index:(index+11),9] = rep(co2$CO2_Concentration_RCP45[which(co2$Year == years[i])],12)
}
}

colnames(met_in) = c('SiteID','Year','Month','Tmin','Tmax','Precip','RA','Frost','CO2')

write.csv(met_in,file='/Users/quinn/Dropbox (VTFRS)/Research/DAPPER_branch/DAPPER/example_run/example_met.csv',row.names = FALSE)








