
#---CONTROL INFORMATION----------------------------
working_directory = '/Users/quinn/Dropbox (VTFRS)/Research/DAPPER'
input_directory = '/Users/quinn/Dropbox (VTFRS)/Research/DAPPER_inputdata/'
run_name = 'process'
#restart_chain = 'duke_state_space_without_trans_2.1.2017-07-21.13.19.13.Rdata'
restart_chain = 'Duke_without_Ctrans.1.2017-08-09.07.32.07.Rdata'
priors_file = 'default_priors.csv'
fr_model = 1  # 1 = estimate FR for each plot, 2 = empirical FR model
FR_fert_assumption = 0 #0 = assume fertilization plots have FR = 1, 1 = do not assume fertilization plots have FR = 1
use_fol = TRUE  #TRUE= use allometric estimates of foliage biomass in fitting
use_dk_pars = 1  #0 = do not use 3 specific parameters for the Duke site, 1 = use the 3 specific parameters
nstreams = 19
state_space = 1
plotFR = 1
windows_machine = FALSE

load(paste(working_directory,'/chains/',restart_chain,sep=''))

outfile = paste(working_directory,'/figures/',run_name,'.Rdata',sep='')

plotlist = c(30001,30018,30049,30041)
sitelist = c(30001,30002,30003,30004)
statelist = c('FL','GA','VA','OK')

rcp = 85
plotlist = c(30018) #c(30018)
plotnum = 3
focal_plotID = plotlist
plotSI = c(15.49341)
plotMaxASW = c(141)
climate_data_index = c(1,1)
clim_model = 1 #seq(1,20,1)
startyear = c(1985,2030)
adjust_rain = c(1,1)
adjust_fert =c(0,0)
adjust_CO2 = c(0,0) 
PROCESS_UNCERT = TRUE
PARAMETER_UNCERT = TRUE
HOLD_CO2 = FALSE
nsamples = 100


all_studies = c(
  '/PINEMAP/TIER3_PINEMAP'
)
#----------------------------------------------------

#---SELECT COMPONENTS THAT ARE ALLOWED TO HAVE UNCERTAINITY--
if(rcp == 85){
  load('/Users/quinn/Dropbox (VTFRS)/Research/PINEMAP/DAPER_development/DAPER_master/DAPER_simple/analysis/Tier3_all_climate_models_RCP85.Rdata')
}else if(rcp == 45){
  load(paste(working_directory,'analysis/Tier3_all_climate_models_RCP45.Rdata',sep=''))
}
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

#-----INDEXING FOR THE PARAMETER VECTOR--------------------------
# this helps speed up the analysis -----------------------------
co2_in = read.csv('/Users/quinn/Dropbox/Research/PINEMAP/DAPER_development/DAPER_master/DAPER_simple/climate_data/CO2_Concentrations_from_CMIP5_1950-2095.csv')
#index_guide = create_index_guide(npars,nplots)

#-----TURN OFF HARDWOOD SIMULATION (0 = HARDWOODS ARE SIMULATED)-------
exclude_hardwoods = array(1,dim=nplots)
exclude_hardwoods[which(initdata$PlotID >= 40000 & initdata$PlotID < 42000)]=0


dyn.load(code_library_plot)



nmodels = length(clim_model)
nperiods = length(startyear)
plotnum = 1


age = array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
lai  = array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
stem  = array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
stem_density = array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
coarse_root= array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
fine_root = array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
fol= array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
total= array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
fSW= array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
ET= array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
Total_Ctrans= array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
GPP= array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
runoff= array(0,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
WUE_ctrans= array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))
WUE_ET= array(-99,dim=c(nsamples,nperiods,nmodels,nplots,nmonths))

median_pars = rep(NA,npars)
for(p in 1:npars){
  median_pars[p] = median(accepted_pars_thinned_burned[,p])
}

start_age = 2
WFi = 1 #"WFi"
WSi = 2 #"WSi"
WCRi = 0.6
WRi = 0.5 #"WRi"
StemNum = 1500 #"StemNoi"

for(s in 1:nsamples){
  if(PARAMETER_UNCERT == FALSE){
    new_pars = median_pars
  }else{
    curr_sample = sample(seq(1,length(accepted_pars_thinned_burned[,1])),1)   
    new_pars = accepted_pars_thinned_burned[curr_sample,]
  }
  pars = new_pars[1:npars_used_by_fortran]

  for(m in 1:nmodels){
    
    curr_model = clim_model[m]
    
    for(yearnum in 1:nperiods){
    
      #--- CREATE CLIMATE INPUT ARRAYS --------------------------------
      if(climate_data_index[yearnum] == 1){
        
      
        curr_met_in_frost = met_in_frost[which(met_in_model == curr_model & met_in_plotnum == plotnum),]
        curr_met_in_pr =  met_in_pr[which(met_in_model == curr_model & met_in_plotnum == plotnum),]
        curr_met_in_tasmax =  met_in_tasmax[which(met_in_model == curr_model & met_in_plotnum == plotnum),]
        curr_met_in_tasmin = met_in_tasmin[which(met_in_model == curr_model & met_in_plotnum == plotnum),]
        curr_met_in_rsds =  met_in_rsds[which(met_in_model == curr_model & met_in_plotnum == plotnum),]
        
        nyears = 26 - start_age
        nmonths = nyears*12
        met = array(NA,c(nmonths,6))
        yrnum = (startyear[yearnum]+start_age) - 1950
        endno = (((yrnum*12)+1)+nmonths)
        tmax = curr_met_in_tasmax[((yrnum*12)+2):endno]
        tmin = curr_met_in_tasmin[((yrnum*12)+2):endno]
        rain = curr_met_in_pr[((yrnum*12)+2):endno]
        solar = (curr_met_in_rsds[((yrnum*12)+2):endno])/1000000*86400
        frost = curr_met_in_frost[((yrnum*12)+2):endno]
        
        c02start = startyear[yearnum]+start_age
        endyear = c02start + nyears-1
        co2plot = co2_in[which(co2_in >= c02start & co2_in <= endyear),]
        if (rcp == 45){
          co2 = rep(co2plot[,2],each=12 )
        }
        if (rcp == 85){
          co2 = rep(co2plot[,3],each=12 )
        }
        
        if(adjust_CO2[yearnum]>0){
          co2[]=adjust_CO2[yearnum]
        }
        
        if(HOLD_CO2 == TRUE){
          c02start = 1985+start_age
          endyear = c02start + nyears-1
          co2plot = co2_in[which(co2_in >= c02start & co2_in <= endyear),]
          co2 = rep(co2plot[,2],each=12)
          year = rep(year_month,each=12)
          month = rep(c(1,2,3,4,5,6,7,8,9,10,11,12),length(year_month))
          if(adjust_CO2[yearnum]>0){
            co2[]=adjust_CO2[yearnum]
          }
        }
        met = array(NA,dim=c(nplots,6,length(tmax)))
        met[,1,] = t(tmin)  #Tmin
        met[,2,] = t(tmax)  # Tmax
        met[,3,] = t(rain) # Rain
        met[,3,] = met[,3,]*adjust_rain[yearnum]
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
      
      InitialYear = startyear[yearnum]
      InitialMonth = 1
      
      WFi=initdata[plotnum,13]
      WSi=initdata[plotnum,14]
      WRi=initdata[plotnum,15]
      WCRi=initdata[plotnum,32]
      
      WFi_H = initdata[plotnum,22]
      WSi_H = initdata[plotnum,23]
      WRi_H = initdata[plotnum,24]
      
      #READ IN SITE DATA FROM FILE BASED ON PLOTNUM
      Lat = LAT_WGS84
      ASWi = ASW_max
      MaxASW = ASW_max
      MinASW = ASW_min
      SoilClass=SoilClass
      
      #if(FertFlag == 1 | fr_model == 1){
      #  FR = new_FR
      #}else{
      #  FR = 1/(1+exp((new_pars[49] + new_pars[50]*Mean_temp-new_pars[51]*SI)))
      #}
      
      FR = plotFR
      
      if(initdata[plotnum,12] == 1) { FR = 1}
      
      
      tmp_site_index = 0
      if(PlotID > 40000 & PlotID < 41000 & use_dk_pars == 1){
        tmp_site_index = 1
      }
      
      SLA = 3.5754 + (5.4287 - 3.5754) * exp(-log(2) * (start_age / 5.9705)^2)
      SLA_h = 16.2
      
      site_in = c(PlantedYear, #PlantedYear
                  PlantedMonth, #"PlantedMonth"
                  InitialYear, #"InitialYear"
                  InitialMonth, #"InitialMonth"
                  start_age,
                  WFi , #"WFi"
                  WRi, #"WRi"
                  WSi, #"WSi"
                  StemNum, #"StemNoi"
                  ASWi = MaxASW, #"ASWi"
                  Lat, #"Lat"
                  FR, #"FR"
                  SoilClass, #"SoilClass"
                  MaxASW, #"MaxASW"
                  MinASW, #"MinASW"
                  TotalMonths = 1,
                  WFi_H = 0.001,
                  WSi_H = 0.001,
                  WRi_H = 0.001,
                  WCRi,
                  IrrigRate = 0.0,
                  Throughfall = 1.0,
                  tmp_site_index,  
                  WCRi_H = 0.0,
                  Wbud_H = 0.0,
                  LAI = -99,
                  LAI_h = 0.01
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
      
      for(mo in 1:nmonths){
 
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
            site[27] = output[9]  #Hardwood LAI
            site[25] = output[26] #Hardwood Bud
            site[18] = output[10] #WS_H 
            site[24] = output[11]  #WCR_h
            site[19] = output[12] #WR_H
            site[10] = output[14] # ASW
            site[17] = output[23] #WF_H	
          }else{ 
            if(PROCESS_UNCERT){
              site[26] = max(rnorm(1,output[4],new_pars[52]),0.1) #LAI
              site[8] = rnorm(1,output[5],new_pars[53]) #WS
              site[20] = rnorm(1,output[6],new_pars[54])   #WCR
              site[7] = rnorm(1,output[7],new_pars[55])  #WRi
              site[9] = rnorm(1,output[8],new_pars[56]) #StemNo
              site[6] = output[22] #WFi
              site[27] = max(rnorm(1,output[9],new_pars[52]),0.0)  #Hardwood LAI
              site[25] = output[26] #Hardwood Bud
              site[18] = rnorm(1,output[10],new_pars[57]) #WS_H 
              site[24] = output[11]  #WCR_h
              site[19] = rnorm(1,output[12],new_pars[55]) #WR_H
              site[10] = output[14] # ASW
              site[17] = output[23] #WF_H	
            }else{
              site[26] = output[4]
              site[8] = output[5]
              site[20] = output[6]
              site[7] = output[7]
              site[9] = output[8]
              site[6] = output[22] #WFi    
              site[27] = output[9]  #Hardwood LAI
              site[25] = output[26] #Hardwood Bud
              site[18] = output[10] #WS_H 
              site[24] = output[11]  #WCR_h
              site[19] = output[12] #WR_H
              site[10] = output[14] # ASW
              site[17] = output[23] #WF_H	
            }
          }
          
          
          
          age[s,yearnum,m,plotnum,mo] = output[3]
          lai[s,yearnum,m,plotnum,mo]  = site[26]
          stem[s,yearnum,m,plotnum,mo]  = site[8]
          stem_density[s,yearnum,m,plotnum,mo] = site[9]
          coarse_root[s,yearnum,m,plotnum,mo] = site[20]
          fine_root[s,yearnum,m,plotnum,mo] = site[7]
          fol[s,yearnum,m,plotnum,mo] = output[22]
          total[s,yearnum,m,plotnum,mo] = fol[s,yearnum,m,plotnum,mo] +  stem[s,yearnum,m,plotnum,mo] +  fine_root[s,yearnum,m,plotnum,mo] + coarse_root[s,yearnum,m,plotnum,mo]
          fSW[s,yearnum,m,plotnum,mo] = output[49]
          ET[s,yearnum,m,plotnum,mo]= output[17]
          Total_Ctrans[s,yearnum,m,plotnum,mo]= output[18] + output[19]
          GPP[s,yearnum,m,plotnum,mo]= output[22]*0.5
          runoff[s,yearnum,m,plotnum,mo]= output[40] 
          #if(m > mo_start_end[plotnum,1]){
          #runoff[s,plotnum,mo]= runoff[s,plotnum,mo-1] + output[40] 
          #}else{
          #  runoff[s,plotnum,mo]= output[40]    
          #}
          WUE_ctrans[s,yearnum,m,plotnum,mo]= GPP[s,yearnum,m,plotnum,mo]/ET[s,yearnum,m,plotnum,mo]
          WUE_ET[s,yearnum,m,plotnum,mo]= GPP[s,yearnum,m,plotnum,mo]/Total_Ctrans[s,yearnum,m,plotnum,mo]
        }else{
          mo = 1
          site_in = c(PlantedYear, #PlantedYear
                      PlantedMonth, #"PlantedMonth"
                      InitialYear, #"InitialYear"
                      InitialMonth, #"InitialMonth"
                      start_age, 
                      WFi, #"WFi"
                      WRi, #"WRi"
                      WSi, #"WSi"
                      StemNum, #"StemNoi"
                      ASWi = MaxASW, #"ASWi"
                      Lat, #"Lat"
                      FR, #"FR"
                      SoilClass, #"SoilClass"
                      MaxASW, #"MaxASW"
                      MinASW, #"MinASW"
                      TotalMonths = 1,
                      WFi_H = 0.001,
                      WSi_H = 0.001,
                      WRi_H = 0.001,
                      WCRi,
                      IrrigRate = IrrigRate,
                      Throughfall = 1.0,
                      tmp_site_index,  
                      WCRi_H = 0.0,
                      Wbud_H = 0.0,
                      LAI = -99,
                      LAI_h = 0.01
          )
          
          site = array(site_in)
        }
      }
    }
  }
}



nmonths = length(age[1,1,1,1,])

LAI_quant = array(NA,dim=c(nperiods,nmonths,3))
stem_quant = array(NA,dim=c(nperiods,nmonths,3))
stem_density_quant = array(NA,dim=c(nperiods,nmonths,3))
coarse_root_quant = array(NA,dim=c(nperiods,nmonths,3))
fine_root_quant = array(NA,dim=c(nperiods,nmonths,3))
fol_quant = array(NA,dim=c(nperiods,nmonths,3))
total_quant = array(NA,dim=c(nperiods,nmonths,3))
fSW_quant = array(NA,dim=c(nperiods,nmonths,3))
ET_quant = array(NA,dim=c(nperiods,nmonths,3))
Total_Ctrans_quant = array(NA,dim=c(nperiods,nmonths,3))
GPP_quant = array(NA,dim=c(nperiods,nmonths,3))
runoff_quant = array(NA,dim=c(nperiods,nmonths,3))
WUE_ctrans_quant = array(NA,dim=c(nperiods,nmonths,3))
WUE_ET_quant = array(NA,dim=c(nperiods,nmonths,3))

modeled_age = age[1,1,1,1,]
for(yearnum in 1:nperiods){
for(i in 1:length(modeled_age)){
  LAI_quant[yearnum,i,] = quantile(lai[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  stem_quant[yearnum,i,] = quantile(stem[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  stem_density_quant[yearnum,i,] = quantile(stem_density[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  coarse_root_quant[yearnum,i,] = quantile(coarse_root[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  fine_root_quant[yearnum,i,] = quantile(fine_root[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  fol_quant[yearnum,i,] = quantile(fol[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  total_quant[yearnum,i,] = quantile(total[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  fSW_quant[yearnum,i,] = quantile(fSW[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  ET_quant[yearnum,i,] = quantile(ET[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  Total_Ctrans_quant[yearnum,i,] = quantile(Total_Ctrans[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  GPP_quant[yearnum,i,] = quantile(GPP[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  runoff_quant[yearnum,i,] = quantile(runoff[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  WUE_ctrans_quant[yearnum,i,] = quantile(WUE_ctrans[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
  WUE_ET_quant[yearnum,i,] = quantile(WUE_ET[,yearnum,,1,i],c(0.025,0.5,0.975),na.rm=TRUE)
}
}

save(lai,stem,stem_density,coarse_root,fine_root,fol,total,fSW,ET,Total_Ctrans,GPP,runoff,WUE_ctrans,WUE_ET,age,nperiods,file = outfile)

pdf(paste(working_directory,'/figures/',run_name,'.pdf',sep=''),width = 11,height = 11)

par(mfrow=c(3,4),mar = c(4,4,2,2),oma = c(3,3,2,2))

plot(modeled_age,LAI_quant[1,,2],type='l',ylim=range(LAI_quant),xlab = 'Stand Age',ylab = 'Leaf Area Index')
polygon(c(modeled_age,rev(modeled_age)),c(LAI_quant[1,,1],rev(LAI_quant[1,,3])),col="lightblue",border=NA)
polygon(c(modeled_age,rev(modeled_age)),c(LAI_quant[2,,1],rev(LAI_quant[2,,3])),col=adjustcolor("salmon",alpha=0.4),border=NA)
points(modeled_age,LAI_quant[1,,2],type='l',col="blue",lwd=1)
points(modeled_age,LAI_quant[2,,2],type='l',col="red",lwd=1)


plot(modeled_age,stem_quant[1,,2],type='l',ylim=range(stem_quant), xlab = 'Stand Age',ylab = 'Stem Biomass (Mg/ha)')
polygon(c(modeled_age,rev(modeled_age)),c(stem_quant[1,,1],rev(stem_quant[1,,3])),col="lightblue",border=NA)
polygon(c(modeled_age,rev(modeled_age)),c(stem_quant[2,,1],rev(stem_quant[2,,3])),col=adjustcolor("salmon",alpha=0.4),border=NA)
points(modeled_age,stem_quant[1,,2],type='l',col="blue",lwd=1)
points(modeled_age,stem_quant[2,,2],type='l',col="red",lwd=1)


plot(modeled_age,total_quant[1,,2],type='l',ylim=range(total_quant), xlab = 'Stand Age',ylab = 'Total Biomass (Mg/ha)')
polygon(c(modeled_age,rev(modeled_age)),c(total_quant[1,,1],rev(total_quant[1,,3])),col="lightblue",border=NA)
polygon(c(modeled_age,rev(modeled_age)),c(total_quant[2,,1],rev(total_quant[2,,3])),col=adjustcolor("salmon",alpha=0.4),border=NA)
points(modeled_age,total_quant[1,,2],type='l',col="blue",lwd=1)
points(modeled_age,total_quant[2,,2],type='l',col="red",lwd=1)

plot(modeled_age,stem_density_quant[1,,2],type='l',ylim=range(stem_density_quant), xlab = 'Stand Age',ylab = 'Stem Density (ind/ha)')
polygon(c(modeled_age,rev(modeled_age)),c(stem_density_quant[1,,1],rev(stem_density_quant[1,,3])),col="lightblue",border=NA)
polygon(c(modeled_age,rev(modeled_age)),c(stem_density_quant[2,,1],rev(stem_density_quant[2,,3])),col=adjustcolor("salmon",alpha=0.4),border=NA)
points(modeled_age,stem_density_quant[1,,2],type='l',col="blue",lwd=1)
points(modeled_age,stem_density_quant[2,,2],type='l',col="red",lwd=1)

plot(modeled_age,fSW_quant[1,,2],type='l',ylim=range(fSW_quant),xlab = 'Stand Age',ylab = 'fSW')
polygon(c(modeled_age,rev(modeled_age)),c(fSW_quant[1,,1],rev(fSW_quant[1,,3])),col="lightblue",border=NA)
polygon(c(modeled_age,rev(modeled_age)),c(fSW_quant[2,,1],rev(fSW_quant[2,,3])),col=adjustcolor("salmon",alpha=0.4),border=NA)
points(modeled_age,fSW_quant[1,,2],type='l',col="blue",lwd=1)
points(modeled_age,fSW_quant[2,,2],type='l',col="red",lwd=1)

plot(modeled_age,ET_quant[1,,2],type='l',ylim=range(ET_quant),xlab = 'Stand Age',ylab = 'ET')
polygon(c(modeled_age,rev(modeled_age)),c(ET_quant[1,,1],rev(ET_quant[1,,3])),col="lightblue",border=NA)
polygon(c(modeled_age,rev(modeled_age)),c(ET_quant[2,,1],rev(ET_quant[2,,3])),col=adjustcolor("salmon",alpha=0.4),border=NA)
points(modeled_age,ET_quant[1,,2],type='l',col="blue",lwd=1)
points(modeled_age,ET_quant[2,,2],type='l',col="red",lwd=1)

plot(modeled_age,Total_Ctrans_quant[1,,2],type='l',ylim=range(Total_Ctrans_quant),xlab = 'Stand Age',ylab = 'Transpiration')
polygon(c(modeled_age,rev(modeled_age)),c(Total_Ctrans_quant[1,,1],rev(Total_Ctrans_quant[1,,3])),col="lightblue",border=NA)
polygon(c(modeled_age,rev(modeled_age)),c(Total_Ctrans_quant[2,,1],rev(Total_Ctrans_quant[2,,3])),col=adjustcolor("salmon",alpha=0.4),border=NA)
points(modeled_age,Total_Ctrans_quant[1,,2],type='l',col="blue",lwd=1)
points(modeled_age,Total_Ctrans_quant[2,,2],type='l',col="red",lwd=1)

plot(modeled_age,GPP_quant[1,,2],type='l',ylim=range(GPP_quant),xlab = 'Stand Age',ylab = 'GPP')
polygon(c(modeled_age,rev(modeled_age)),c(GPP_quant[1,,1],rev(GPP_quant[1,,3])),col="lightblue",border=NA)
polygon(c(modeled_age,rev(modeled_age)),c(GPP_quant[2,,1],rev(GPP_quant[2,,3])),col="lightgreen",col=adjustcolor("salmon",alpha=0.4),border=NA)
points(modeled_age,GPP_quant[1,,2],type='l',col="blue",lwd=1)
points(modeled_age,GPP_quant[2,,2],type='l',col="red",lwd=1)

plot(modeled_age,runoff_quant[1,,2],type='l',ylim=range(runoff_quant),xlab = 'Stand Age',ylab = 'Runoff')
polygon(c(modeled_age,rev(modeled_age)),c(runoff_quant[1,,1],rev(runoff_quant[1,,3])),col="lightblue",border=NA)
polygon(c(modeled_age,rev(modeled_age)),c(runoff_quant[2,,1],rev(runoff_quant[2,,3])),col="lightgreen",col=adjustcolor("salmon",alpha=0.4),border=NA)
points(modeled_age,runoff_quant[1,,2],type='l',col="blue",lwd=1)
points(modeled_age,runoff_quant[2,,2],type='l',col="red",lwd=1)

plot(modeled_age,WUE_ctrans_quant[1,,2],type='l',ylim=range(WUE_ctrans_quant),xlab = 'Stand Age',ylab = 'WUE (Transpiration)')
polygon(c(modeled_age,rev(modeled_age)),c(WUE_ctrans_quant[1,,1],rev(WUE_ctrans_quant[1,,3])),col="lightblue",border=NA)
polygon(c(modeled_age,rev(modeled_age)),c(WUE_ctrans_quant[2,,1],rev(WUE_ctrans_quant[2,,3])),col=adjustcolor("salmon",alpha=0.4),border=NA)
points(modeled_age,WUE_ctrans_quant[1,,2],type='l',col="blue",lwd=1)
points(modeled_age,WUE_ctrans_quant[2,,2],type='l',col="red",lwd=1)

plot(modeled_age,WUE_ET_quant[1,,2],type='l',ylim=range(WUE_ET_quant),xlab = 'Stand Age',ylab = 'WUE (ET)')
polygon(c(modeled_age,rev(modeled_age)),c(WUE_ET_quant[1,,1],rev(WUE_ET_quant[1,,3])),col="lightblue",border=NA)
polygon(c(modeled_age,rev(modeled_age)),c(WUE_ET_quant[2,,1],rev(WUE_ET_quant[2,,3])),col=adjustcolor("salmon",alpha=0.4),border=NA)
points(modeled_age,WUE_ET_quant[1,,2],type='l',col="blue",lwd=1)
points(modeled_age,WUE_ET_quant[2,,2],type='l',col="red",lwd=1)

#par(mfrow=c(1,1),mar = c(4,4,2,2),oma = c(3,3,2,2))
#plot(density(stem[,1,length(modeled_age)]),xlim=c(0,300),ylim=c(0,1),col='red')
#points(density(fol[,1,length(modeled_age)]),type='l',col='green')
#points(density(fine_root[,1,length(modeled_age)]),type='l',col='blue')
#points(density(coarse_root[,1,length(modeled_age)]),type='l',col='orange')
#points(density(total[,1,length(modeled_age)]),type='l',col='black')
dev.off()