
#---CONTROL INFORMATION----------------------------
working_directory =  getwd()
run_name = 'with_trans_importance'
restart_from_chain = FALSE
#restart_chain = 'duke_state_space_without_trans_2.1.2017-07-21.13.19.13.Rdata'
restart_chain = 'duke_state_space.1.2017-07-24.06.59.06.Rdata'
priors_file = 'default_priors.csv'
obs_set = 14 #14 #Select which plots are used in analysis.  See prepare_obs.R for number guide 
focal_plotID = 40001 #14 #Select which plots are used in analysis.  See prepare_obs.R for number guide 
fr_model = 1  # 1 = estimate FR for each plot, 2 = empirical FR model
FR_fert_assumption = 0 #0 = assume fertilization plots have FR = 1, 1 = do not assume fertilization plots have FR = 1
use_fol = TRUE  #TRUE= use allometric estimates of foliage biomass in fitting
use_dk_pars = 1  #0 = do not use 3 specific parameters for the Duke site, 1 = use the 3 specific parameters
nstreams = 18
state_space = 1
plotFR = 1

load(paste(working_directory,'/chains/',restart_chain,sep=''))


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
code_library_plot = paste(working_directory,'/source_code/r3pg_interface.so',sep='')

final_pdf = paste(working_directory,'/figures/',run_name,'.pdf',sep='')

setwd(paste(working_directory,'/scripts/',sep=''))
source('prepare_obs.R')
source('prepare_state_space_obs.R')
source('assign_control_plots.R')
source('prepare_met.R')
#source('initialize_pars.R')
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

init_pars = priormatrix[,1]
latent = obs


dyn.load(code_library_plot)

plotnum = 1

age = array(-99,dim=c(npars,nplots,nmonths))
lai  = array(-99,dim=c(npars,nplots,nmonths))
stem  = array(-99,dim=c(npars,nplots,nmonths))
stem_density = array(-99,dim=c(npars,nplots,nmonths))
coarse_root= array(-99,dim=c(npars,nplots,nmonths))
fine_root = array(-99,dim=c(npars,nplots,nmonths))
fol= array(-99,dim=c(npars,nplots,nmonths))
total= array(-99,dim=c(npars,nplots,nmonths))
fSW= array(-99,dim=c(npars,nplots,nmonths))

p_change = array(-99,dim=npars)

median_pars = rep(NA,npars)
for(p in 1:npars){
  median_pars[p] = median(accepted_pars_thinned_burned[,p])
}

for(p in 1:(npars_used_by_fortran+1)){
  print(p)
  
    new_pars = median_pars
    if(p <= npars_used_by_fortran & priormatrix[p,5] != 1){
      if(p == 28){
    new_pars[p] = new_pars[p]*0.97
      }else{
        new_pars[p] = new_pars[p]*1.03     
      }
    
    p_change[p] = new_pars[p] - median_pars[p]
    print(c(p,p_change[p]))
    }

    pars = new_pars[1:npars_used_by_fortran]
  
    new_FR = plotFR
  
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
  if(FertFlag == 1 | fr_model == 1){
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
              LAI = WFi* SLA * 0.1,
              LAI_h = WFi_H * SLA_h *0.1
  )
  
  site = array(site_in)
  
  #THIS DEALS WITH THINNING BASED ON PROPORTION OF STEMS REMOVED
  thin_event = array(0,dim=c(nplots,nmonths))
  
  
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
  for(mo in mo_start_end[plotnum,1]:mo_start_end[2]){
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
    
    output=array(tmp$out_var, dim=c(nomonths_plot,output_dim))
    
    if(output[2] == 12){
      site[3] = output[1]+1 #InitialYear
      site[4] = 1  #InitialMonth
    }else{
      site[3] = output[1] #InitialYear
      site[4] = output[2]+1  #InitialMonth	
    }
    site[5] = output[3] + (1.0/12.) #StartAge
    site[26] = output[4] #LAI
    if(site[26] < 0.0) {site[26]=0.1}
    site[8] = output[5] #WS
    site[20] = output[6]   #WCR
    site[7] = output[7]  #WRi
    site[9] = output[8] #StemNo
    
    site[27] = output[9]  #Hardwood LAI
    site[25] = output[26] #Hardwood Bud
    site[18] = output[10] #WS_H 
    site[24] = output[11]  #WCR_h
    site[19] = output[12] #WR_H
    
    site[10] = output[14] # ASW
    
    site[6] = output[22] #WFi
    site[17] = output[23] #WF_H	
    
    age[p,plotnum,mo] = output[3]
    lai[p,plotnum,mo]  = output[4]
    stem[p,plotnum,mo]  = output[5]
    stem_density[p,plotnum,mo] = output[8]
    coarse_root[p,plotnum,mo] = output[6]
    fine_root[p,plotnum,mo] = output[7]
    fol[p,plotnum,mo] = output[22]
    total[p,plotnum,mo] = output[22] + output[5] + output[6] + output[7]
    fSW[p,plotnum,mo] = output[49]
  }
}

delta_tot=rep(NA,npars_used_by_fortran)
sens_tot=rep(NA,npars_used_by_fortran)
var_p =rep(NA,npars_used_by_fortran)
par_importance =rep(NA,npars_used_by_fortran)
prior_var =rep(NA,npars_used_by_fortran)
par_importance_prior =rep(NA,npars_used_by_fortran)
for(p in 1:(npars_used_by_fortran)){
  delta_tot[p] = total[p,plotnum,mo_start_end[2]] - total[npars_used_by_fortran+1,plotnum,mo_start_end[2]]
  sens_tot[p] = (delta_tot[p]/p_change[p])^2
  var_p[p] = var(accepted_pars_thinned_burned[,p])
  par_importance[p] = sens_tot[p]*var_p[p]
  if(priormatrix[p,4]==1){
  prior_var[p] = (1/12)*(priormatrix[p,3] - priormatrix[p,2])^2
  }else if(priormatrix[p,4]==2){
    prior_var[p] = priormatrix[p,3]     
  }
  par_importance_prior[p] = sens_tot[p]*prior_var[p]
}

t = data.frame(parname = parnames[1:npars_used_by_fortran],imp = par_importance[1:npars_used_by_fortran],dy = delta_tot[1:npars_used_by_fortran],dydx = sens_tot[1:npars_used_by_fortran],var = var_p[1:npars_used_by_fortran],imp_prior = par_importance_prior)
t = t[which(priormatrix[1:npars_used_by_fortran,6]==0),]
t_sort = t[order(t$imp),]
sd_per_matrix = t_sort$imp
sd_per_prior_matrix = t_sort$imp_prior

sd_per_matrix_index = seq(1,length(t_sort$imp),1)

y = c(0,(length(sd_per_matrix_index)+1))
x = c(0,max(sd_per_matrix))
xlim=c(-10,(max(c(sd_per_matrix,sd_per_prior_matrix),na.rm=TRUE)+10))
xlim=c(-10,(max(c(sd_per_matrix),na.rm=TRUE)+1))
#xlim=c(-10,20000)
pdf(paste(working_directory,'/chains/',run_name,'.pdf',sep=''),width = 5,height = 8)
par(mfrow=c(1,1),mar = c(4,4,0,1.01),oma = c(0,10,0,0))
plot(x,y,col='white',bty='n',xaxs="i", yaxs="i",xaxt='n',yaxt='n',xlab='Posterior Importance',ylab='',xlim=xlim)
axis(1)
axis(2,at=sd_per_matrix_index,
     labels=t_sort$parname,
     las=2,tick = FALSE)
segments(0,0,0,max(sd_per_matrix_index),lwd=2)
points(sd_per_matrix,sd_per_matrix_index,pch=20)
for(i in 1:length(sd_per_matrix)){
  segments(0,sd_per_matrix_index[i],sd_per_matrix[i],sd_per_matrix_index[i])
}
points(sd_per_prior_matrix,sd_per_matrix_index,pch=20,col='gray')
for(i in 1:length(sd_per_matrix)){
  segments(0,sd_per_matrix_index[i],sd_per_prior_matrix[i],sd_per_matrix_index[i])
}

mtext('Parameter',side = 2, line = 13,at = 25)
dev.off()

p = 14
plot(total[p,1,],type='l',ylim=c(0,max(c(total[,1,]))))
points(total[(npars_used_by_fortran+1),1,],type='l',col='red')