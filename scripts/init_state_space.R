init_state_space <- function(){
  init_state_space_LAI = array(-99,dim=c(nplots,nmonths))
  init_state_space_WS = array(-99,dim=c(nplots,nmonths))
  init_state_space_WR = array(-99,dim=c(nplots,nmonths))
  init_state_space_WCR = array(-99,dim=c(nplots,nmonths))
  init_state_space_stem_density = array(-99,dim=c(nplots,nmonths))
  init_state_space_LAI_H= array(-99,dim=c(nplots,nmonths))
  init_state_space_WS_H= array(-99,dim=c(nplots,nmonths))
  init_state_space_WCR_H= array(-99,dim=c(nplots,nmonths))
  init_state_space_WR_H= array(-99,dim=c(nplots,nmonths))
  init_state_space_stem_density_H= array(-99,dim=c(nplots,nmonths))
  
  
  
  gap = array(0,dim=c(nstreams))
  prev_mo= array(0,dim=c(nstreams))
  obs_gap = array(-99,dim=c(nstreams,nplots,nmonths))
  obs_gap_next = array(-99,dim=c(nstreams,nplots,nmonths))
  prev_obs= array(0,dim=c(5))
  
  
  init_state_space_ASW = array(-99,dim=c(nplots,nmonths))
  age = array(-99,dim=c(nplots,nmonths))
  
  dyn.load(code_library_plot)
  
  for(plotnum in 1:nplots){

    new_pars = init_pars
    
    pars = new_pars[1:npars_used_by_fortran]
    
    index1 = index_guide[3]+plotnum-1
    new_FR = new_pars[index1]
    
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
    
    if(use_fol_state[plotnum] == 0){
    SLA = 3.5754 + (5.4287 - 3.5754) * exp(-log(2) * (StartAge / 5.9705)^2)
    SLA_h = 16.2
    LAI = WFi* SLA * 0.1
    LAI_h = WFi_H * SLA_h *0.1
    }else{
      LAI = -99
      LAI_h = -99
    }

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
                LAI = LAI,
                LAI_h = LAI_h
    )
    
    site = array(site_in)
    
    #THIS DEALS WITH THINNING BASED ON PROPORTION OF STEMS REMOVED
    #thin_event = array(0,dim=c(nplots,nmonths))
    
    
    output_dim = noutput_variables  # NUMBER OF OUTPUT VARIABLES
    nosite = length(site_in)  # LENGTH OF SITE ARRAY
    nomet = 6  # NUMBER OF VARIABLES IN METEROLOGY (met)
    nopars = length(pars)
    
    #Wsx1000 (plot level variability in parameter)
    pars[19] = new_pars[index_guide[5]+plotnum - 1]
    #thinpower (plot level variability in parameter)
    pars[20] = new_pars[index_guide[7]+plotnum - 1]
    pars[40] = new_pars[index_guide[9]+plotnum - 1]    
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
    
    age[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]]  = output[,3]
    
    if(use_fol_state[plotnum]== 0){
    init_state_space_LAI[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]] = output[,4]
    }else{
      init_state_space_LAI[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]] = output[,22]
    }
    init_state_space_WS[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]] = output[,5]
    init_state_space_WCR[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]] = output[,6]
    init_state_space_WR[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]] = output[,7]
    init_state_space_stem_density[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]]= output[,8]
    
    
    init_state_space_LAI_H[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]] = output[,9]
    init_state_space_WS_H[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]] = output[,10]
    init_state_space_WCR_H[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]] = output[,11]
    init_state_space_WR_H[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]] = output[,12]
    init_state_space_stem_density_H[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]]= output[,13]
    
    init_state_space_ASW[plotnum,mo_start_end[plotnum,1]:mo_start_end[plotnum,2]] = output[,14]
    
    if(state_space == 0 | state_space == 1){
      for(mo in (mo_start_end[plotnum,1]+1):mo_start_end[plotnum,2]){
        if(months[mo] != 8 & PlotID < 40000 & PlotID >= 41000){
          init_state_space_LAI[plotnum,mo]=-99
          obs[1,plotnum,mo] = -99
        }
        if(obs[1,plotnum,mo] == -99){init_state_space_LAI[plotnum,mo]=-99}
        if(obs[2,plotnum,mo] == -99){init_state_space_WS[plotnum,mo]=-99}
        if(obs[3,plotnum,mo] == -99){init_state_space_WCR[plotnum,mo]=-99}
        if(obs[4,plotnum,mo] == -99){init_state_space_WR[plotnum,mo]=-99}
        if(obs[5,plotnum,mo] == -99){init_state_space_stem_density[plotnum,mo]=-99}
        if(obs[6,plotnum,mo] == -99){init_state_space_LAI_H[plotnum,mo]=-99}
        if(obs[7,plotnum,mo] == -99){init_state_space_WS_H[plotnum,mo]=-99}
        
        if(obs[1,plotnum,mo] != -99){init_state_space_LAI[plotnum,mo]=obs[1,plotnum,mo]}
        if(obs[2,plotnum,mo] != -99){init_state_space_WS[plotnum,mo]=obs[2,plotnum,mo]}
        if(obs[3,plotnum,mo] != -99){init_state_space_WCR[plotnum,mo]=obs[3,plotnum,mo]}
        if(obs[4,plotnum,mo] != -99){init_state_space_WR[plotnum,mo]=obs[4,plotnum,mo]}
        if(obs[5,plotnum,mo] != -99){init_state_space_stem_density[plotnum,mo]=obs[5,plotnum,mo]}
        if(obs[6,plotnum,mo] != -99){init_state_space_LAI_H[plotnum,mo]=obs[6,plotnum,mo]}
        if(obs[7,plotnum,mo] != -99){init_state_space_WS_H[plotnum,mo]=obs[7,plotnum,mo]}
      }
    }else if(state_space == 2){
      for(mo in (mo_start_end[plotnum,1]+1):mo_start_end[plotnum,2]){
        if(obs[1,plotnum,mo] != -99){init_state_space_LAI[plotnum,mo:mo_start_end[plotnum,2]]=obs[1,plotnum,mo]}
        if(obs[2,plotnum,mo] != -99){init_state_space_WS[plotnum,mo:mo_start_end[plotnum,2]]=obs[2,plotnum,mo]}
        if(obs[3,plotnum,mo] != -99){init_state_space_WCR[plotnum,mo:mo_start_end[plotnum,2]]=obs[3,plotnum,mo]}
        if(obs[4,plotnum,mo] != -99){init_state_space_WR[plotnum,mo:mo_start_end[plotnum,2]]=obs[4,plotnum,mo]}
        if(obs[5,plotnum,mo] != -99){init_state_space_stem_density[plotnum,mo:mo_start_end[plotnum,2]]=obs[5,plotnum,mo]}
        if(obs[6,plotnum,mo] != -99){init_state_space_LAI_H[plotnum,mo:mo_start_end[plotnum,2]]=obs[6,plotnum,mo]}
        if(obs[7,plotnum,mo] != -99){init_state_space_WS_H[plotnum,mo:mo_start_end[plotnum,2]]=obs[7,plotnum,mo]}
      }
    }
    
    
    for(data_stream in 1:7){
      prev_mo[data_stream] = (mo_start_end[plotnum,1]+1)
      gap[data_stream] = 0
      for(mo in (mo_start_end[plotnum,1]+1):mo_start_end[plotnum,2]){
        gap[data_stream] =  gap[data_stream] + 1
        if(obs[data_stream,plotnum,mo] != -99){
          if(state_space == 1){
            obs_gap[data_stream,plotnum,mo] = gap[data_stream]
            obs_gap_next[data_stream,plotnum,prev_mo[data_stream]] = gap[data_stream]
            prev_mo[data_stream] = mo
            gap[data_stream] = 0
          }else if(state_space == 0){
            obs_gap[data_stream,plotnum,mo] = 1 #gap[data_stream]
          }
        }
        if(state_space == 2){
          obs_gap[data_stream,plotnum,mo] = 1
          obs_gap_next[data_stream,plotnum,] = 1
        }
      }
    }
    
  }
  return(list(init_state_space_LAI=init_state_space_LAI,init_state_space_WS=init_state_space_WS,init_state_space_WR=init_state_space_WR,
              init_state_space_WCR = init_state_space_WCR,init_state_space_stem_density=init_state_space_stem_density,
              init_state_space_LAI_H=init_state_space_LAI_H,init_state_space_WS_H=init_state_space_WS_H,init_state_space_WR_H=init_state_space_WR_H,
              init_state_space_WCR_H = init_state_space_WCR_H,init_state_space_stem_density_H=init_state_space_stem_density_H,
              init_state_space_ASW=init_state_space_ASW,age=age,obs_gap=obs_gap,obs_gap_next = obs_gap_next))
}
