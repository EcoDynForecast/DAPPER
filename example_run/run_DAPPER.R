rm(list = ls())
#---CONTROL INFORMATION----------------------------
working_directory = "/Users/quinn/Dropbox (VTFRS)/Research/DAPPER_branch/DAPPER/example_run/"
DAPPER_directory =  "/Users/quinn/Dropbox (VTFRS)/Research/DAPPER_branch/DAPPER/"
input_directory = "/Users/quinn/Dropbox/Research/DAPPER_inputdata_public/"
niter = 1
chain_number = 1
burn =  1
thin_interval = 1
print_debug = 1
run_name = 'Your_duke_assimilation'
restart_from_chain = FALSE
restart_chain =  NA
priors_file = 'example_priors.csv'
validation_set_file = NA #file name of .csv that defines plot numbers for plots that are not fit but are compared to the obs.
met_file = 'example_met.csv'
plot_file = 'example_plots.csv'
observations_file = 'example_observations.csv'
create_plot = TRUE
only_create_plot = FALSE
focal_plotID = NA #30001 #Setting a value here causes only a single plot to be simulated and fit
fr_model = 1  # 1 = estimate FR for each plot, 2 = empirical FR model
FR_fert_assumption = 0 #0 = assume fertilization plots have FR = 1, 1 = do not assume fertilization plots have FR = 1
FR_separate_npar_groups = 2  #Assigns a different parameter group to groups of FR values: 0 = all one group, 1 = separate groups, 2 = all plots separate
use_fol = TRUE  #TRUE= use allometric estimates of foliage biomass in fitting
use_dk_pars = 1  #0 = do not use 3 specific parameters for the Duke site, 1 = use the 3 specific parameters
use_age_edc = 0  #0 = do not use an ecological constraint on the age function (see code); 1 = use the constraint
use_sm_edc = 0  #0 = do not use an ecological constraint on the soil moisture function (see code); 1 = use the constraint
use_fr_edc = 0   #0 = do note use an ecological constraint on the SI - FR function (see code); 1 = use the constraint
nstreams = 19
use_gep = 1
use_et = 1  
use_ctrans =1 
use_gep_uncert = 1
use_et_uncert = 1
use_ctrans_uncert = 1
state_space = 1
tracked_plotnum = 1
windows_machine = FALSE
#----------------------------------------------------
#obs_set = 14 #Select which plots are used in analysis.  See prepare_obs.R for number guide 
#all_studies = c(
#  '/Duke/TIER4_Duke'
#)

#---SELECT COMPONENTS THAT ARE ALLOWED TO HAVE UNCERTAINITY--
plot_WSx1000 = FALSE  #include plot specific WSx1000 parameter
plot_thinpower = FALSE #include plot specific thinpower parameter
plot_mort_rate = FALSE #include plot specific mortality rate parameter

setwd(paste(DAPPER_directory,'/scripts/',sep=''))
source('prepare_DAPPER_MCMC.R')