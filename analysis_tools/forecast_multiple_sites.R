rm(list = ls())
#---CONTROL INFORMATION----------------------------
working_directory = '/Users/quinn/Dropbox (VTFRS)/Research/DAPPER'
input_directory = '/Users/quinn/Dropbox (VTFRS)/Research/DAPPER_inputdata/'
run_name = 'test'
#restart_chain = 'duke_state_space_without_trans_2.1.2017-07-21.13.19.13.Rdata'
restart_chain =  'SS_val6.1.2017-08-31.17.33.17.Rdata'
priors_file = 'default_priors.csv'
obs_set = 21 #14 #Select which plots are used in analysis.  See prepare_obs.R for number guide 
focal_plotID = NA #14 #Select which plots are used in analysis.  See prepare_obs.R for number guide 
val_set = 6
fr_model = 1  # 1 = estimate FR for each plot, 2 = empirical FR model
FR_fert_assumption = 0 #0 = assume fertilization plots have FR = 1, 1 = do not assume fertilization plots have FR = 1
use_fol = TRUE  #TRUE= use allometric estimates of foliage biomass in fitting
use_dk_pars = 1  #0 = do not use 3 specific parameters for the Duke site, 1 = use the 3 specific parameters
nstreams = 19
state_space = 1
plotFR = NA
windows_machine = FALSE

PARAMETER_UNCERT = TRUE
nsamples = 500

load(paste(working_directory,'/chains/',restart_chain,sep=''))

all_studies = c(
  '/SETRES/TIER4_SETRES',
  '/PINEMAP/TIER3_PINEMAP',
  '/NC2/TIER4_NC2',
  '/Duke/TIER4_Duke',
  '/Waycross/TIER4_Waycross',
  '/FMC_Thinning/TIER1_FMC_Thinning',
  #'/FBRC_AMERIFLU/TIER2_AMERIFLU',
  #'/FBRC_IMPAC/TIER1_IMPAC',
  #'/FBRC_IMPAC2/TIER2_IMPAC2',
  #'/FBRC_PPINES/TIER2_PPINES',
  #'/FBRC_VAR1/TIER2_VAR1',
  #'/FBRC_WPPINES/TIER2_WPPINES',
  #'/FMC_IMP_TIER1/TIER1_IMP',
  #'/FMC_IMP_TIER2/TIER2_IMP',
  #'/FPC_RS1/TIER1_RS1',
  #'/FPC_RS2/TIER1_RS2',
  #'/FPC_RS3/TIER1_RS3',
  #'/FPC_RS5/TIER1_RS5',
  #'/FPC_RS6/TIER1_RS6',
  #'/FPC_RS7/TIER1_RS7',
  #'/FPC_RS8/TIER1_RS8',
  '/FPC_RW18/TIER2_RW18'
  #'/FPC_RW19/TIER2_RW19',
  #'/FPC_RW20/TIER2_RW20',
  #'/PMRC_CPCD96_TIER1/TIER1_CPCD96',
  #'/PMRC_CPCD96_TIER2/TIER2_CPCD96',
  #'/PMRC_HGLOB87/TIER1_HGLOB87',
  #'/PMRC_SAGCD96_TIER1/TIER1_SAGCD96',
  #'/PMRC_SAGCD96_TIER2/TIER2_SAGCD96',
  #'/PMRC_SAGSP85_TIER1/TIER1_SAGSP85',
  #'/PMRC_SAGSP85_TIER2/TIER2_SAGSP85',
  #'/PMRC_WGCD01_TIER1/TIER1_WGCD01',
  #'/PMRC_WGCD01_TIER2/TIER2_WGCD01',
  #'/TAMU_GSSS/TIER1_GSSS'
  #'/FIA/VA_FIA'
)
#----------------------------------------------------

#---SELECT COMPONENTS THAT ARE ALLOWED TO HAVE UNCERTAINITY--
#plot_WSx1000 = FALSE  #include plot specific WSx1000 parameter
#plot_thinpower = FALSE #include plot specific thinpower parameter
#plot_mort_rate = FALSE #include plot specific mortality rate parameter

#----------------------------------------------------

#----------------------------------------------------

#--OTHER INFO (WON'T CHANGE UNLESS MODIFY MODEL)-----
npars_used_by_fortran = 48
noutput_variables = 67
process_model_pars = 51
npars =80

#---- ENTER THE FORTRAN LIBRARY NAMES HERE ----------
if(windows_machine){
  code_library_plot = paste(working_directory,'/source_code/r3pg_interface.dll',sep='')
}else{
  code_library_plot = paste(working_directory,'/source_code/r3pg_interface.so',sep='')
}

final_pdf = paste(working_directory,'/figures/',run_name,'.pdf',sep='')

setwd(paste(working_directory,'/scripts/',sep=''))
source('prepare_obs.R')
source('prepare_state_space_obs.R')
source('assign_control_plots.R')
source('prepare_met.R')
source('set_fitted_plots.R')
source('init_state_space.R')

setwd(working_directory)
#----------------------------------------------------

priors_in = read.csv(paste(working_directory,'/priors/',priors_file,sep=''))
npars = length(priors_in$parnames)
priormatrix = matrix(NA,npars,6)
priormatrix[,1] = priors_in$initial_value
priormatrix[,2] = priors_in$dist_par1
priormatrix[,3] = priors_in$dist_par2
priormatrix[,4] = priors_in$dist_type
priormatrix[,5] = priors_in$fit_par
priormatrix[,6] = priors_in$par_group

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

state_space_obs = prepare_state_space_obs()
mo_start_end=state_space_obs$mo_start_end
years = state_space_obs$years
months = state_space_obs$months
nmonths = state_space_obs$nmonths
obs = state_space_obs$obs
thin_event = state_space_obs$thin_event
obs_uncert = state_space_obs$obs_uncert
init_obs = state_space_obs$init_obs
init_uncert = state_space_obs$init_uncert


#----SET CONTROL PLOT INDEX---------------------------------------
# this assigns the control plot to match with the treatment plot

control_list = assign_control_plots(nplots,initdata,plotlist)
control_plot_index =  control_list$control_plot_index
matched_FR_plot_index = control_list$matched_FR_plot_index

#--- CREATE CLIMATE INPUT ARRAYS --------------------------------
met_tmp = prepare_met(met_in,initdata,mo_start_end,co2_in,nplots,nmonths,months,years)
met = array(NA,dim=c(nplots,6,length(met_tmp$tmin[1,])))
met[,1,] = met_tmp$tmin
met[,2,] = met_tmp$tmax
met[,3,] = met_tmp$precip
met[,4,] = met_tmp$ra
met[,5,] = met_tmp$frost
met[,6,] = met_tmp$co2

#-----INDEXING FOR THE PARAMETER VECTOR--------------------------
# this helps speed up the analysis -----------------------------

#index_guide = create_index_guide(npars,nplots)

#-----TURN OFF HARDWOOD SIMULATION (0 = HARDWOODS ARE SIMULATED)-------
exclude_hardwoods = array(1,dim=nplots)
exclude_hardwoods[which(initdata$PlotID >= 40000 & initdata$PlotID < 42000)]=0

init_pars = accepted_pars_thinned_burned[1,]
new_pars = init_pars


latent = obs
init_states = init_state_space()
latent[1,,] = init_states$init_state_space_LAI
latent[2,,]=init_states$init_state_space_WS
latent[3,,]=init_states$init_state_space_WCR
latent[4,,]=init_states$init_state_space_WR
latent[5,,]=init_states$init_state_space_stem_density
latent[6,,] = init_states$init_state_space_LAI_H
latent[7,,]=init_states$init_state_space_WS_H

age = init_states$age
obs_gap = init_states$obs_gap
obs_gap_next = init_states$obs_gap_next

#####
thinning_study_second_thin = read.csv(paste(input_directory,'/FMC_Thinning/TIER1_FMC_list_of_second_thin_plots.csv',sep=''))
thin_event = array(0,dim=c(nplots,nmonths))
for(plotnum in 1:nplots){
  tmp_initdata = initdata[which(initdata$PlotID == plotlist[plotnum]),]
  prev_nha = init_obs[5,plotnum]
  thin_occured = 0
  double_thin = 0
  if(length(which(thinning_study_second_thin$PlotID ==  tmp_initdata$Plot)) == 1 & plotlist[plotnum] > 10000 & plotlist[plotnum] < 20000){
    double_thin = 1
  }
  for(mo in (mo_start_end[plotnum,1]+1):mo_start_end[plotnum,2])
    if(obs[5,plotnum,mo] != -99){
      thin_event[plotnum,mo-1] =  prev_nha - obs[5,plotnum,mo] 
      if(plotlist[plotnum] > 20000 & plotlist[plotnum] < 41000){
        thin_event[plotnum,mo-1] = 0.0
      }
      if(plotlist[plotnum] > 42000 & plotlist[plotnum] < 50000){
        thin_event[plotnum,mo-1] = 0.0
      }
      if(plotlist[plotnum] > 10000 & plotlist[plotnum] < 20000 & tmp_initdata$ThinTreatment == 1){
        thin_event[plotnum,mo-1] = 0.0
      }
      if(plotlist[plotnum] >= 52001 & plotlist[plotnum] <= 52467 & tmp_initdata$ThinTreatment == 1){
        thin_event[plotnum,mo-1] = 0.0
      }
      if(plotlist[plotnum] >= 72001 & plotlist[plotnum] <= 72076 & tmp_initdata$ThinTreatment == 1){
        thin_event[plotnum,mo-1] = 0.0
      }
      if(plotlist[plotnum] > 10000 & plotlist[plotnum] < 20000 & tmp_initdata$ThinTreatment > 1 & thin_occured == 1 & double_thin == 0){
        thin_event[plotnum,mo-1] = 0.0
      }
      if(plotlist[plotnum] >= 52001 & plotlist[plotnum] <= 52467 & tmp_initdata$ThinTreatment > 1 & (thin_event[plotnum,mo-1] < 400 | thin_occured == 1)){
        thin_event[plotnum,mo-1] = 0.0
      }
      if(plotlist[plotnum] >= 72001 & plotlist[plotnum] <= 72076 & tmp_initdata$ThinTreatment > 1 & (thin_event[plotnum,mo-1] < 400 | thin_occured == 1)){
        thin_event[plotnum,mo-1] = 0.0
      }
      
      if(thin_event[plotnum,mo-1] < 200){
        thin_event[plotnum,mo-1]  = 0
      }else if(thin_occured == 0 & thin_event[plotnum,mo-1] >= 200){
        thin_occured = 1
      }
      
      prev_nha = obs[5,plotnum,mo]
    }
}

#--- ASSIGN VECTOR FLAGGING PLOTS AS USED IN FITTING OR NOT
fit_plot = set_fitted_plots(nplots,val_set,plotlist,initdata)
num_in95 = 0
num_out95 = 0
num_tot = 0


pdf(paste(working_directory,'/figures/',run_name,'.pdf',sep=''),width = 11,height = 11)
par(mfrow=c(4,4),mar = c(4,4,2,2),oma = c(3,3,2,2))

dyn.load(code_library_plot)



for(plotnum in 1:nplots){
  if((val_set > 0 & fit_plot[plotnum]==0) | val_set == 0){
    
    median_pars = rep(NA,npars)
    for(p in 1:npars){
      median_pars[p] = median(accepted_pars_thinned_burned[,p])
    }
    
    plot_nmonths = length(mo_start_end[plotnum,1]:mo_start_end[plotnum,2])
    age_model = array(-99,dim=c(nsamples,plot_nmonths))
    lai  = array(-99,dim=c(nsamples,plot_nmonths))
    stem  = array(-99,dim=c(nsamples,plot_nmonths))
    stem_density = array(-99,dim=c(nsamples,plot_nmonths))
    
    for(s in 1:nsamples){
      
      if(!PARAMETER_UNCERT){
        new_pars = median_pars
      }else{
        curr_sample = sample(seq(1,length(accepted_pars_thinned_burned[,1])),1)   
        new_pars = accepted_pars_thinned_burned[curr_sample,]
      }
      pars = new_pars[1:npars_used_by_fortran]
      
      new_FR = plotFR
      
      tmp_initdata = initdata[plotnum,]
      
      PlotID = initdata[plotnum,1]
      SiteID = initdata[plotnum,2]
      LAT_WGS84=initdata[plotnum,3]
      Planting_year = initdata[plotnum,4]
      PlantMonth = initdata[plotnum,5]
      PlantDensityHa = initdata[plotnum,6]
      Initial_ASW = initdata[plotnum,7]
      ASW_min = initdata[plotnum,8]
      ASW_max=initdata[plotnum,9]
      SoilClass = initdata[plotnum,10]
      SI = initdata[plotnum,11]
      FR = initdata[plotnum,12]
      Initial_WF = initdata[plotnum,13]
      Initial_WS = initdata[plotnum,14]
      Initial_WR = initdata[plotnum,15]
      DroughtLevel = initdata[plotnum,16]
      DroughtStart = initdata[plotnum,17]
      FertFlag=initdata[plotnum,18]
      CO2flag = initdata[plotnum,19]
      CO2elev = initdata[plotnum,20]
      ControlPlotID = initdata[plotnum,21]
      Initial_WF_H = initdata[plotnum,22]
      Initial_WS_H =initdata[plotnum,23]
      Initial_WR_H =initdata[plotnum,24]
      InitialYear = initdata[plotnum,29]
      InitialMonth = initdata[plotnum,30]
      StartAge = initdata[plotnum,31] 
      IrrFlag = initdata[plotnum,33] 
      Mean_temp = initdata[plotnum,26]
      
      
      PlantedYear = 0
      PlantedMonth = 0
      
      InitialYear = years[mo_start_end[plotnum]]
      InitialMonth = months[mo_start_end[plotnum]]
      
      WFi=initdata[plotnum,13]
      WSi=initdata[plotnum,14]
      WRi=initdata[plotnum,15]
      WCRi=initdata[plotnum,32]
      
      WFi_H = initdata[plotnum,22]
      WSi_H = initdata[plotnum,23]
      WRi_H = initdata[plotnum,24]
      
      StemNum = PlantDensityHa
      nomonths_plot = mo_start_end[plotnum,2] - mo_start_end[plotnum,1]+1
      
      #READ IN SITE DATA FROM FILE BASED ON PLOTNUM
      Lat = LAT_WGS84
      ASWi = ASW_max
      MaxASW = ASW_max
      MinASW = ASW_min
      SoilClass=SoilClass
      
      if(!is.na(plotFR)){
        FR = new_FR
      }else{
        FR = 1/(1+exp((new_pars[49] + new_pars[50]*Mean_temp-new_pars[51]*SI)))
      }
      
      if(initdata[plotnum,12] == 1) { FR = 1}
      
      IrrigRate = 0.0
      if(IrrFlag == 1){
        IrrigRate = (658/9)
      }
      
      tmp_site_index = 0
      if(PlotID > 40000 & PlotID < 41000 & use_dk_pars == 1){
        tmp_site_index = 1
      }
      
      SLA = 3.5754 + (5.4287 - 3.5754) * exp(-log(2) * (StartAge / 5.9705)^2)
      SLA_h = 16.2
      
      site_in = c(PlantedYear, #PlantedYear
                  PlantedMonth, #"PlantedMonth"
                  InitialYear, #"InitialYear"
                  InitialMonth, #"InitialMonth"
                  StartAge, #"EndAge"
                  WFi, #"WFi"
                  WRi, #"WRi"
                  WSi, #"WSi"
                  StemNum, #"StemNoi"
                  ASWi, #"ASWi"
                  Lat, #"Lat"
                  FR, #"FR"
                  SoilClass, #"SoilClass"
                  MaxASW, #"MaxASW"
                  MinASW, #"MinASW"
                  TotalMonths = 1,
                  WFi_H = Initial_WF_H,
                  WSi_H = Initial_WS_H,
                  WRi_H = Initial_WR_H,
                  WCRi,
                  IrrigRate = IrrigRate,
                  Throughfall = DroughtLevel,
                  tmp_site_index,  
                  WCRi_H = WSi_H*0.30,
                  Wbud_H = 0.0,
                  LAI = tmp_initdata$Initial_LAI,
                  LAI_h = WFi_H * SLA_h *0.1
      )
      
      site = array(site_in)
      
      #THIS DEALS WITH THINNING BASED ON PROPORTION OF STEMS REMOVED
      #thin_event = array(0,dim=c(nplots,nmonths))
      
      
      output_dim = noutput_variables  # NUMBER OF OUTPUT VARIABLES
      nosite = length(site_in)  # LENGTH OF SITE ARRAY
      nomet = 6  # NUMBER OF VARIABLES IN METEROLOGY (met)
      nopars = length(pars)
      
      #Wsx1000 (plot level variability in parameter)
      # pars[19] = new_pars[index_guide[5]+plotnum - 1]
      #thinpower (plot level variability in parameter)
      #  pars[20] = new_pars[index_guide[7]+plotnum - 1]
      #  pars[40] = new_pars[index_guide[9]+plotnum - 1]    
      #Read in Fortran code
      if(PlotID > 40000 & PlotID < 41000){
        pars[19] = new_pars[48]
      }
      
      mo_index = 0
      for(mo in mo_start_end[plotnum,1]:mo_start_end[plotnum,2]){
        mo_index = mo_index + 1
        
        
        tmp=.Fortran( "r3pg_interface",
                      output_dim=as.integer(output_dim),
                      met=as.double(met[plotnum,,mo]),
                      pars=as.double(pars),
                      site = as.double(site),
                      thin_event = as.double(thin_event[plotnum,mo]),
                      out_var=as.double(array(0,dim=c(1,output_dim))),
                      nopars=as.integer(nopars),
                      nomet=as.integer(dim(met)[2]),
                      nosite = as.integer(nosite),
                      nooutputs=as.integer(output_dim),
                      nomonths_plot=as.integer(1),
                      nothin = 1,
                      exclude_hardwoods = as.integer(exclude_hardwoods[plotnum]),
                      mo_start_end = as.integer(c(1,1)),
                      nmonths = 1
        )
        
        output=array(tmp$out_var, dim=c(1,output_dim))
        
        if(output[2] == 12){
          site[3] = output[1]+1 #InitialYear
          site[4] = 1  #InitialMonth
        }else{
          site[3] = output[1] #InitialYear
          site[4] = output[2]+1  #InitialMonth	
        }
        site[5] = output[3] + (1.0/12.) #StartAge
        site[26] = rnorm(1,output[4],new_pars[52]) #LAI
        if(is.na(site[26])) {site[26]=0.1}
        if(site[26] < 0.0) {site[26]=0.1}

        site[8] = rnorm(1,output[5],(1.3+new_pars[53] +output[5]*new_pars[64]))  #WS
        #site[8] = rnorm(1,output[5],0.15+new_pars[53] +output[5]*(0.015))  #WS
        if(is.na(site[8])) {site[8]=0.1}
        if(site[8]< 0.0) {site[8]=0.1}
        site[20] = rnorm(1,output[6],new_pars[54])   #WCR
        site[7] = rnorm(1,output[7],new_pars[55])  #WRi
        site[9] = rnorm(1,output[8],new_pars[56]) #StemNo
        
        site[27] = max(rnorm(1,output[9],new_pars[52]),0.0)  #Hardwood LAI
        site[25] = output[26] #Hardwood Bud
        site[18] = rnorm(1,output[10],new_pars[57]) #WS_H 
        site[24] = output[11]  #WCR_h
        site[19] = rnorm(1,output[12],new_pars[55]) #WR_H
        
        site[10] = output[14] # ASW
        
        site[6] = output[22] #WFi
        site[17] = output[23] #WF_H	
        
        age_model[s,mo_index] = output[3]
        lai[s,mo_index]  = site[26]
        stem[s,mo_index]  = site[8]
        stem_density[s,mo_index] = site[9]
      }
    }
    
    LAI_quant = array(NA,dim=c(length(age_model[1,]),3))
    stem_quant = array(NA,dim=c(length(age_model[1,]),3))
    stem_density_quant = array(NA,dim=c(length(age_model[1,]),3))
    
    
    modeled_age = age_model[1,]
    for(i in 1:length(modeled_age)){
      LAI_quant[i,] = quantile(lai[,i],c(0.025,0.5,0.975))
      stem_quant[i,] = quantile(stem[,i],c(0.025,0.5,0.975))
      stem_density_quant[i,] = quantile(stem_density[,i],c(0.025,0.5,0.975))
    }
    
    data_stream = 2
    obs[data_stream,plotnum,which(obs[data_stream,plotnum,]!=-99)]
    age[plotnum,which(obs[data_stream,plotnum,]!=-99)]
    
    #Calculate total  
    tmp_age_obs = age[plotnum,which(obs[data_stream,plotnum,]!=-99)]
    tmp_stem_obs = obs[data_stream,plotnum,which(obs[data_stream,plotnum,]!=-99)]
    tmp_stem_obs_uncert = obs_uncert[data_stream,plotnum,which(obs_uncert[data_stream,plotnum,]!=-99)]
    tmp_stem_obs_uncert = tmp_stem_obs*0.025
    xlim_range = c(0,30) #c(min(output[,3])-1,max(output[,3])+1)
    ylim_range = c(0,max(c(tmp_stem_obs,stem_quant),na.rm=TRUE)) #range(c(observed,modeled[data_stream,]))
    plot(modeled_age,stem_quant[,2],type='l',ylim=ylim_range, xlab = 'Stand Age',ylab = 'Stem Biomass (Mg/ha)',main=plotlist[plotnum])
    polygon(c(modeled_age,rev(modeled_age)),c(stem_quant[,1],rev(stem_quant[,3])),col="lightblue",border=NA)
    points(modeled_age,stem_quant[,2],type='l',col="blue",lwd=1)
    for(i in 1:length(tmp_age_obs)){
      index = which.min(abs(modeled_age-tmp_age_obs[i]))
      tmp_stem_obs_range = quantile(rnorm(1000,tmp_stem_obs[i],tmp_stem_obs_uncert[i]),c(0.025,0.975))
      points(tmp_age_obs[i],tmp_stem_obs[i],col='black',pch=20)
      segments(tmp_age_obs[i],tmp_stem_obs_range[1],tmp_age_obs[i],tmp_stem_obs_range[2])
      if((tmp_stem_obs_range[1] >= stem_quant[index,1] & tmp_stem_obs_range[1] <= stem_quant[index,3]) |
         (tmp_stem_obs_range[2] >= stem_quant[index,1] & tmp_stem_obs_range[2] <= stem_quant[index,3])){
        num_in95 = num_in95 + 1
        num_tot = num_tot + 1
      }else{
        num_out95 = num_out95 + 1
        num_tot = num_tot + 1      
      }
      
    }
    
  }
}

dev.off()
