print('--- INITIALIZING DAPPER---')

#VARIABLES THAT CODE LOOKS FOR BUT AREN'T USED ANYMORE
start_adapt = 1 #NOT USED ANYMORE NEED TO CUT
cost_type = 2 #NOTE USED # 1 = specified undertainity 2 = fitted uncertainity
high_freq_obs = 1 #NOT USED #Include the monthly flux observations in fitting
initial_conditions_npar_group = 1  #Parameter group for the initial condition (NOT TESTED)
obs_uncertainity = 1 #NOT USED #0 = no observational uncertainity; 1 = include observational uncertainity


#----MCMC TUNING INFORMATION-------------------------
jump_size_init = 0.1
lower_accept_bound=0.23
upper_accept_bound=0.44
nadapt = 1000
adaptfac = 1.5
std_fac = 2
#----------------------------------------------------
setwd(paste(DAPPER_directory,'/scripts/',sep=''))
source('prepare_obs.R')
source('prepare_state_space_obs.R')
source('set_fitted_plots.R')
source('assign_control_plots.R')
source('prepare_met.R')
source('create_index_guide.R')
source('initialize_pars.R')
source('init_state_space.R')
source('plot_output.R')

setwd(working_directory)
#--OTHER INFO (WON'T CHANGE UNLESS MODIFY MODEL)-----
npars_used_by_fortran = 48
noutput_variables = 68
process_model_pars = 51
#----------------------------------------------------

priors_in = read.csv(paste0(working_directory,'/',priors_file))

npars = length(priors_in$parnames)
priormatrix = matrix(NA,npars,6)
priormatrix[,1] = priors_in$initial_value
priormatrix[,2] = priors_in$dist_par1
priormatrix[,3] = priors_in$dist_par2
priormatrix[which(priors_in$dist_type == 'uniform'),4] = 1
priormatrix[which(priors_in$dist_type == 'normal'),4] = 2
priormatrix[which(priors_in$fit_par == 'TRUE'),5] = 0
priormatrix[which(priors_in$fit_par == 'FALSE'),5] = 1
priormatrix[,6] = priors_in$par_group
parnames = priors_in$parnames

index = 1
for(p in 1:npars){
  if(priormatrix[p,5] == 0){
    priormatrix[p,6] = index
    index = index + 1
  }else{
    priormatrix[p,6] = 10000
  }
}



#---- ENTER THE FORTRAN LIBRARY NAMES HERE ----------
if(windows_machine){
  code_library_iter = paste(DAPPER_directory,'/source_code/DAPPER_MCMC.dll',sep='')
  code_library_plot = paste(DAPPER_directory,'/source_code/r3pg_interface.dll',sep='')
}else{
  code_library_iter = paste(DAPPER_directory,'/source_code/DAPPER_MCMC.so',sep='')
  code_library_plot = paste(DAPPER_directory,'/source_code/r3pg_interface.so',sep='')
}

final_pdf = paste(working_directory,run_name,'.pdf',sep='')
#----------------------------------------------------

#---  PREPARE OBSERVATIONS ---------------------------
obs_list = prepare_obs(obs_set,FR_fert_assumption,use_fol)
plotlist = obs_list$plotlist
nplots= obs_list$nplots
observations= obs_list$observations
initdata= obs_list$initdata
met_in = obs_list$met_in
co2_in = obs_list$co2_in
use_fol_state = obs_list$use_fol_state
#-------------------------------------------------
data_uncertainity_npar_group = index
index = index + 1
if(plot_WSx1000){
  plot_WSx1000_pargroup = index
  index = plot_WSx1000_pargroup + 8
}else{
  plot_WSx1000_pargroup = 10000
}

FR_npar_group = index
if(fr_model == 1 & FR_separate_npar_groups == 1){ #when fitting FR there are more par groups
  npar_groups = FR_npar_group + 8
}else if(fr_model == 2 & FR_separate_npar_groups == 1){
  npar_groups = FR_npar_group + 4
}else if(fr_model == 1 & FR_separate_npar_groups == 2){
  npar_groups = FR_npar_group + nplots - 1
}else{
  npar_groups = FR_npar_group
}


plot_thinpower_pargroup = 10000 
plot_mort_rate_pargroup = npar_groups + 10000 
#--------------------------------


state_space_obs = prepare_state_space_obs()
mo_start_end=state_space_obs$mo_start_end
nobs = state_space_obs$nobs
years = state_space_obs$years
months = state_space_obs$months
nmonths = state_space_obs$nmonths
obs = state_space_obs$obs
thin_event = state_space_obs$thin_event
obs_uncert = state_space_obs$obs_uncert
init_obs = state_space_obs$init_obs
init_uncert = state_space_obs$init_uncert

if(length(which(all_studies == '/FMC_Thinning/TIER1_FMC'))){
  thinning_study_second_thin = read.csv(paste(input_directory,'/FMC_Thinning/TIER1_FMC_list_of_second_thin_plots.csv',sep=''))
}

thin_event = array(0,dim=c(nplots,nmonths))
for(plotnum in 1:nplots){
  tmp_initdata = initdata[which(initdata$PlotID == plotlist[plotnum]),]
  prev_nha = init_obs[5,plotnum]
  thin_occured = 0
  double_thin = 0
  if(length(which(all_studies == '/FMC_Thinning/TIER1_FMC'))){
    if(length(which(thinning_study_second_thin$PlotID ==  tmp_initdata$Plot)) == 1 & plotlist[plotnum] > 10000 & plotlist[plotnum] < 20000){
      double_thin = 1
    }
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

#if(length(which(plotlist != 41001))>0){
#  thin_event[which(plotlist != 41001),] =0.0
#}


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

index_guide = create_index_guide(npars,nplots)

#--INITIALIZE PARAMETER AND JUMP VECTOR--------------------------

init_list = initialize_pars(index_guide,priormatrix,jump_size_init,observations,initdata,fr_model,chain_number,process_model_pars)
init_pars = init_list$init_pars
new_pars = init_list$new_pars
jump_pars = init_list$jump_pars
prior_parameter1 = init_list$prior_parameter1
prior_parameter2 = init_list$prior_parameter2
prior_dist = init_list$prior_dist
fix_par = init_list$fix_par
par_group = init_list$par_group

#-----TURN OFF HARDWOOD SIMULATION (0 = HARDWOODS ARE SIMULATED)-------
exclude_hardwoods = array(1,dim=nplots)
exclude_hardwoods[which(initdata$PlotID >= 40000 & initdata$PlotID < 42000)]=0

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
#init_state_space_ASW=init_states$init_state_space_ASW

#---ASSIGN CHAIN NAME --------------------------------------------
tmp =  strsplit(as.character(Sys.time()), " ")
date = tmp[[1]][1]
time = strsplit(as.character(tmp[[1]][2]), ":")
#chain_file_name = paste(run_name,chain_number,date,time[[1]][1],time[[1]][2],time[[1]][1],'Rdata',sep=".")
chain_file_name = paste(run_name,'Rdata',sep=".")
final_chain_file_name = paste(working_directory,chain_file_name,sep="")

#-----------------------------------------------------------------

#----INITIAL VARIABLES FOR MCMC---------------------------------
pnow = -1e30
global_accept = 0
local_accept = 0
current_like = -1e30

if(restart_from_chain){ #start from a previously run chain
  #if the chain number is > 1 then start from a random location along previous chain
  #if the chain number is 1 then start from the end of the previous chain
  if(chain_number == 1){
    rm(jump_pars)
    load(paste(working_directory,restart_chain,sep=''))
    max_iter = dim(accepted_pars_thinned_burned)[1]
    init_pars=accepted_pars_thinned_burned[max_iter,]
  }else{
    rm(jump_pars)
    load(paste(working_directory,restart_chain,sep=''))
    max_iter = dim(accepted_pars_thinned_burned)[1]
    s = sample(1,seq(1,max_iter,1))
    init_pars=accepted_pars_thinned_burned[s,]
  }
}

#--- ASSIGN VECTOR FLAGGING PLOTS AS USED IN FITTING OR NOT
fit_plot = set_fitted_plots(nplots,val_set,plotlist,initdata)

#---ONLY RETURN BACK TO R A SAMPLED CHAIN (SAVES MEMORY DEMANDS)
sample_index = seq(burn,niter,thin_interval)
accepted_pars_thinned_burned = array(-999,dim=c(length(sample_index),length(init_pars)))
like_chain = array(-999,dim=c(length(sample_index)))
tracked_plot = array(-99, dim=c(length(sample_index),nstreams,nmonths))
#----------------------------------------------------------------


#------RUN FORTRAN CODE---------------------------------
if(only_create_plot == FALSE){
  #----LOAD FORTRAN CODE------------------------------------
  dyn.load(code_library_iter)
  ptm <- proc.time()
  print('--- RUNNING MCMC---')
  fortran_output=.Fortran("DAPPER_MCMC"
                          ,nopars = as.integer(length(init_pars))
                          ,nplots	= as.integer(nplots)
                          ,initdata_dim =  as.integer(dim(initdata[,1:34])[2])
                          ,control_plot_index = as.integer(control_plot_index)
                          ,index_guide = as.integer(index_guide)
                          ,nosamples = as.integer(length(sample_index))
                          ,sample_index = as.integer(sample_index)
                          ,exclude_hardwoods = as.integer(exclude_hardwoods)
                          ,fit_plot = as.integer(fit_plot)
                          ,matched_FR_plot_index = as.integer(matched_FR_plot_index)
                          ,par_group =  as.integer(par_group)
                          ,control_pars = c(as.integer(niter),
                                            as.integer(start_adapt),
                                            as.integer(cost_type),
                                            as.integer(fr_model),
                                            as.integer(high_freq_obs),
                                            as.integer(obs_uncertainity),
                                            as.integer(use_dk_pars),
                                            as.integer(use_age_edc),
                                            as.integer(use_fr_edc),
                                            as.integer(use_sm_edc),
                                            as.integer(state_space))
                          ,npar_groups = as.integer(npar_groups)
                          ,data_uncertainity_npar_group  = as.integer(data_uncertainity_npar_group)
                          ,nstreams = as.integer(nstreams)
                          ,nmonths = as.integer(nmonths)
                          ,years = as.integer(years)
                          ,months =as.integer(months)
                          ,mo_start_end = as.integer(mo_start_end)
                          ,met = as.double(met)
                          ,initdata =data.matrix(initdata[,1:34])
                          ,obs = as.double(obs)
                          ,thin_event =as.double(thin_event)
                          ,init_pars = as.double(init_pars)
                          ,prior_parameter1 = as.double(prior_parameter1)
                          ,prior_parameter2 = as.double(prior_parameter2)
                          ,prior_dist = as.double(prior_dist)
                          ,fix_par = as.double(fix_par)
                          ,obs_uncert = as.double(obs_uncert)
                          ,latent = as.double(latent)
                          ,jump_pars = as.double(jump_pars)
                          ,pnow = as.double(pnow)
                          ,accepted_pars_thinned_burned = as.double(accepted_pars_thinned_burned)
                          ,like_chain = as.double(like_chain)
                          ,current_like = as.double(current_like)
                          ,init_obs = as.double(init_obs)
                          ,init_uncert = as.double(init_uncert)
                          ,tracked_plot = as.double(tracked_plot)
                          ,tracked_plotnum = as.integer(tracked_plotnum)
                          ,obs_gap = as.integer(obs_gap)
                          ,obs_gap_next = as.integer(obs_gap_next)
                          ,use_fol_state = as.integer(use_fol_state))
  
  print(proc.time() - ptm)
  
  
  #-----------------------------------------------------------------
  
  #-----PROCESS FORTRAN OUTPUT--------------------------------------
  
  accepted_pars_thinned_burned=array(fortran_output$accepted_pars_thinned_burned, dim=c(length(sample_index),length(init_pars)))
  jump_pars = array(fortran_output$jump_pars,dim=c(as.integer(length(init_pars))))
  pnow = fortran_output$pnow
  like_chain = array(fortran_output$like_chain,dim=c(length(sample_index)))
  current_like = fortran_output$current_like  
  latent = array(fortran_output$latent,dim=c(nstreams,nplots,nmonths))
  tracked_plot = array(fortran_output$tracked_plot, dim=c(length(sample_index),nstreams,nmonths))
  
  #-----------------------------------------------------------------
  
  #---SAVE CHAIN AS AN R BINARY-------------------------------------
  save(accepted_pars_thinned_burned,jump_pars,pnow,index_guide,obs,months,years,
       control_plot_index,current_like,like_chain,fit_plot,initdata,
       latent,age,tracked_plot,tracked_plotnum,file = final_chain_file_name)
  #-----------------------------------------------------------------
}else{
  load(paste(working_directory,restart_chain,sep=''))
}
#------PLOT FINAL RESULTS FROM CHAIN------------------------------
if(create_plot){
  print('--- CREATING PLOTS---')
  plot_output()
}
#-----------------------------------------------------------------
