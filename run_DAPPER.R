rm(list = ls())
#---CONTROL INFORMATION----------------------------
working_directory =  '/Users/quinn/Dropbox (VTFRS)/Research/DAPPER/'
input_directory = '/Users/quinn/Dropbox (VTFRS)/Research/DAPPER_inputdata/'
niter = 100 #150000 #10000 #500000#15000000 #5000000 # 5000000
chain_number = 1
burn = 1 #75000 #25000 #15000000 # 100000 #2500000 #10000000 #500000 #2500000
thin_interval = 1 # 10 #10 #10 #10 # 50 #500 #100 #100 #100
run_name = 'test_old3'
restart_from_chain = FALSE
restart_chain =  'state_space.1.2017-06-25.16.52.16.Rdata'
priors_file = 'dukectrans_priors.csv'
create_plot = TRUE
only_create_plot = FALSE
obs_set = 14 #Select which plots are used in analysis.  See prepare_obs.R for number guide 
focal_plotID = NA #30001 #Setting a value here causes only a single plot to be simulated and fit
val_set = 0  #Select which plots are withheld from fitting.  0 includes all plot
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
use_gep_uncert = 0
use_et_uncert = 0
use_ctrans_uncert = 0
state_space = 1
tracked_plotnum = 1
windows_machine = FALSE
#----------------------------------------------------
all_studies = c(
  #'/SETRES/TIER4_SETRES',
  #'/PINEMAP/TIER3_PINEMAP',
  #'/NC2/TIER4_NC2',
  '/Duke/TIER4_Duke'
  #'/Waycross/TIER4_Waycross',
  #'/FMC_Thinning/TIER1_FMC_Thinning',
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
  #'/FPC_RW18/TIER2_RW18'
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

#---SELECT COMPONENTS THAT ARE ALLOWED TO HAVE UNCERTAINITY--
plot_WSx1000 = FALSE  #include plot specific WSx1000 parameter
plot_thinpower = FALSE #include plot specific thinpower parameter
plot_mort_rate = FALSE #include plot specific mortality rate parameter

setwd(paste(working_directory,'/scripts/',sep=''))
source('prepare_DAPPER_MCMC.R')
