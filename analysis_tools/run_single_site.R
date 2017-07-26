rm(list = ls())
#---CONTROL INFORMATION----------------------------
working_directory = '/Users/quinn/Dropbox (VTFRS)/Research/DAPPER'
input_directory = '/Users/quinn/Dropbox (VTFRS)/Research/DAPPER_inputdata/'
run_name = 'test2'
restart_from_chain = FALSE
restart_chain = 'state_space.1.2017-06-25.16.52.16.Rdata'
priors_file = 'default_priors.csv'
obs_set = 14 #14 #Select which plots are used in analysis.  See prepare_obs.R for number guide 
focal_plotID = 40001 #14 #Select which plots are used in analysis.  See prepare_obs.R for number guide 
fr_model = 1  # 1 = estimate FR for each plot, 2 = empirical FR model
FR_fert_assumption = 0 #0 = assume fertilization plots have FR = 1, 1 = do not assume fertilization plots have FR = 1
use_fol = TRUE  #TRUE= use allometric estimates of foliage biomass in fitting
use_dk_pars = 1  #0 = do not use 3 specific parameters for the Duke site, 1 = use the 3 specific parameters
nstreams = 18
state_space = 1
plotFR = 0.5
#----------------------------------------------------
all_studies = c(
  #'/SETRES/TIER4_SETRES',
  #'/PINEMAP/TIER3_PINEMAP',
  #'/NC2/TIER4_NC2',
  '/Duke/TIER4_Duke',
  #'/FMC_Thinning/TIER1_FMC_Thinning',
  '/FBRC_AMERIFLU/TIER2_AMERIFLU',
  '/FBRC_IMPAC/TIER1_IMPAC',
  #'/FBRC_IMPAC2/TIER2_IMPAC2',
  '/FBRC_PPINES/TIER2_PPINES',
  '/FBRC_VAR1/TIER2_VAR1',
  '/FBRC_WPPINES/TIER2_WPPINES',
  '/FMC_IMP/TIER2_IMP',
  '/FPC_RS1/TIER1_RS1',
  '/FPC_RS2/TIER1_RS2',
  '/FPC_RS3/TIER1_RS3',
  '/FPC_RS5/TIER1_RS5',
  '/FPC_RS6/TIER1_RS6',
  '/FPC_RS7/TIER1_RS7',
  '/FPC_RS8/TIER1_RS8',
  '/FPC_RW18/TIER2_RW18',
  '/FPC_RW19/TIER2_RW19',
  #'/FPC_RW20/TIER2_RW20'
  '/PMRC_CPCD96_TIER1/TIER1_CPCD96',
  '/PMRC_CPCD96_TIER2/TIER2_CPCD96',
  '/PMRC_HGLOB87/TIER1_HGLOB87',
  '/PMRC_SAGCD96_TIER1/TIER1_SAGCD96',
  '/PMRC_SAGCD96_TIER2/TIER2_SAGCD96',
  '/PMRC_SAGSP85_TIER1/TIER1_SAGSP85',
  '/PMRC_SAGSP85_TIER2/TIER2_SAGSP85',
  '/PMRC_WGCD01_TIER1/TIER1_WGCD01',
  '/PMRC_WGCD01_TIER2/TIER2_WGCD01',
  '/TAMU_GSSS/TIER1_GSSS'
)



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

new_pars = init_pars

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
            LAI = tmp_initdata$Initial_LAI,
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



tmp=.Fortran( "r3pg_interface",
              output_dim=as.integer(output_dim),
              met=as.double(met[plotnum,,]),
              pars=as.double(pars),
              site = as.double(site),
              thin_event = as.double(thin_event[plotnum,]),
              out_var=as.double(array(0,dim=c(nomonths_plot,output_dim))),
              nopars=as.integer(nopars),
              nomet=as.integer(dim(met)[2]),
              nosite = as.integer(nosite),
              nooutputs=as.integer(output_dim),
              nomonths_plot=as.integer(nomonths_plot),
              nothin = 1,
              exclude_hardwoods = as.integer(exclude_hardwoods[plotnum]),
              mo_start_end = as.integer(mo_start_end[plotnum,]),
              nmonths = nmonths
)

output=array(tmp$out_var, dim=c(nomonths_plot,output_dim))

age = array(-99,dim=c(nplots,nmonths))
age[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]] = output[,3]


final_pdf = paste(working_directory,'/figures/',run_name,'_single_site_',focal_plotID,'_.pdf',sep='')
pdf(final_pdf,height = 3.4252*3.0,width = 3.4252*3.0)
plotnum = 1
#par(mfrow=c(4,4))
par(mfrow=c(4,5),mar = c(4,4,2,2),oma = c(3,3,2,2))
xlim_range = c(0,30) #c(min(output[,3])-1,max(output[,3])+1)
modeled = output[,6]
modeled_y = output[,3]
observed = obs[2,plotnum,which(obs[2,plotnum,]!=-99)]
observed_y =age[plotnum,which(obs[2,plotnum,]!=-99)]
ylim_range = range(c(observed,modeled))
plot(modeled_y,modeled,pch=20,xlim=xlim_range,ylim=ylim_range,col='white',xlab='',ylab='',main='',xaxt='n', yaxt ='n',ann=FALSE)
legend('topleft',legend=c('data','model'),lty=c(1,1),col=c('gray','black'),cex=1,bty='n')
legend('bottomleft',legend=c('fNutr','fT','fFrost','fCalpha','fVPD','fSW','fAge'),lty=c(1,1,1,1,1,1,1),col=c('black','gray','brown','blue','red','green','purple'),cex=1.0,bty='n')

ylab=rep(NA,nstreams)
ylab[1] = 'LAI'
ylab[2] = 'Stem'
ylab[3] = 'Coarse root'
ylab[4] = 'Fine root'
ylab[5] = 'Stem density'
ylab[6] = 'LAI Hard'
ylab[7] = 'Stem Hard'
ylab[8] = 'Coarse root Hard'
ylab[9] = 'Fine root Hard'
ylab[10] = 'Stem density Hard'
ylab[11] = 'ASW'    
ylab[12] = 'GEP'
ylab[13] = 'NEE'
ylab[14] = 'ET'
ylab[15] ='Ctrans Pine'
ylab[16] ='Ctrans Hard'

modeled = array(NA,dim=c(nstreams,length(output[,4])))

modeled[1,] = output[,4] #!Pine LAI
modeled[2,] = output[,5] #WS
modeled[3,] = output[,6] #WCR
modeled[4,] = output[,7] + output[,12] #WR
modeled[5,] = output[,8] #Stem Density
modeled[6,] = output[,9] #LAI_H        
modeled[7,] = output[,10] #Hardwood Stem   
modeled[8,] = output[,11] #Hardwood Coarse roots
modeled[9,] = output[,12] #Hardwood Fine roots 
modeled[10,] = output[,13] #Hardwood Stem density 
modeled[11,] = output[,14] # ASW        
modeled[12,] = output[,15] # GEP        
modeled[13,] = output[,16] # NEE
modeled[14,] = output[,17] # ET  
modeled[15,] = output[,18] # Ctrans Pine           
modeled[16,] = output[,19] # Ctrans Hardwood     


modeled_age = output[,3]

for(data_stream in 1:nstreams){
  
  if(data_stream != 8 & data_stream != 9 & data_stream != 10 & data_stream != 11 & data_stream != 13 & data_stream != 17 & data_stream != 18){
    observed = obs[data_stream,plotnum,which(obs[data_stream,plotnum,]!=-99)]
    observed_y =age[plotnum,which(obs[data_stream,plotnum,]!=-99)]
    xlim_range = c(0,30) #c(min(output[,3])-1,max(output[,3])+1)
    ylim_range = c(0,max(c(observed,modeled[data_stream,]))) #range(c(observed,modeled[data_stream,]))
    plot(modeled_age,modeled[data_stream,],type='l',xlim=xlim_range,ylim=ylim_range,col='black',xlab='plot age (yr)',ylab=ylab[data_stream])
    points(observed_y,observed,col='gray',pch=20)
    if(data_stream == 1) {title(paste(plotlist[plotnum],StudyName[plotnum],sep=' '))}
    if(data_stream == 2) {title(paste('Treat: ',Treatment[plotnum],sep=''))}
  }
}

plot((output[,1]+(output[,2]-1)/12),output[,14]/output[,39],ylim=c(0,1),type='l',ylab='Avail Soil Water / Max Avail Soil Water',xlab='Year')
plot((output[,1]+(output[,2]-1)/12),output[,43],type='l',ylab='Env. modifier',xlab='Year',ylim=c(0,1.5))
points((output[,1]+(output[,2]-1)/12),output[,44],type='l',col='gray')
points((output[,1]+(output[,2]-1)/12),output[,45],type='l',col='brown')
points((output[,1]+(output[,2]-1)/12),output[,47],type='l',col='blue')
points((output[,1]+(output[,2]-1)/12),output[,48],type='l',col='red')
points((output[,1]+(output[,2]-1)/12),output[,49],type='l',col='green')
points((output[,1]+(output[,2]-1)/12),output[,50],type='l',col='purple')
#legend('topright',legend=c('fNutr','fT','fFrost','fCalpha','fVPD','fSW','fAge'),lty=c(1,1,1,1,1,1,1),col=c('black','gray','brown','blue','red','green','purple'),cex=0.5)

dev.off()

