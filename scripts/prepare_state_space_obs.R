prepare_state_space_obs <- function(){
  
  obs_uncertainity_proportion = 0.01 
  
  
  earliestYear = min(observations$YearMeas)
  earliestMonth = min(observations$MonthMeas[which(observations$YearMeas==earliestYear)])
  latestYear = max(observations$YearMeas)
  latestMonth = max(observations$MonthMeas[which(observations$YearMeas==latestYear)])
  nmonths =(latestYear-earliestYear)*12+1 + (latestMonth-earliestMonth)
  
  years = rep(NA,nmonths)
  months = rep(NA,nmonths)
  curr_year = earliestYear
  curr_month = earliestMonth
  for(mo in 1:nmonths){
    years[mo] = curr_year
    months[mo] = curr_month
    curr_month = curr_month + 1
    if(curr_month == 13){
      curr_month = 1
      curr_year = curr_year + 1
    }
  }
  
  plots = unique(observations$PlotID)
  init_obs=array(-99,dim=c(nstreams,nplots))
  init_uncert=array(-99,dim=c(nstreams,nplots))
  nplots = length(plots)
  mo.start = rep(0,nplots)
  mo.end = rep(0,nplots)
  for(plotnum in 1:nplots){
    tmp = observations[which(observations$PlotID == plots[plotnum]),]
    tmp2 = initdata[which(initdata$PlotID == plots[plotnum]),]
    plot_earliestYear = min(tmp2$InitYear)
    plot_latestYear = max(tmp$YearMeas)
    plot_earliestMonth= min(tmp2$InitMonth)
    plot_latestMonth = max(tmp$MonthMeas[which(tmp$YearMeas==plot_latestYear)])
    plot_earliest_index = (plot_earliestYear - earliestYear)*12+1 + (plot_earliestMonth-earliestMonth)
    plot_latest_index = (plot_latestYear - earliestYear)*12+1 + (plot_latestMonth-earliestMonth)
    mo.start[plotnum] = plot_earliest_index
    mo.end[plotnum] = plot_latest_index
    priormatrix[7,]  = c(5.4287,5.529,0.44,2,0,6)#SLA0 
    priormatrix[8,]  = c(3.5754,3.875,0.11,2,0,7) #SLA1 
    priormatrix[9,]  = c(5.9705,5.9705,2.15,2,0,8) #tSLA
    SLA = 3.5754 + (5.4287 - 3.5754) * exp(-log(2) * (tmp2$StartAge / 5.9705)^2)
    SLA_H = 16.2
    if(use_fol_state[plotnum] == 0){
      init_obs[1,plotnum] = tmp2$Initial_LAI
      if(tmp2$Initial_LAI == -99){
        init_obs[1,plotnum] = tmp2$Initial_WF*SLA * 0.1
      }
      if(tmp2$Initial_LAI_code == 0){
        init_uncert[1,plotnum] = init_obs[1,plotnum]*0.25
      }else{
        init_uncert[1,plotnum] = init_obs[1,plotnum]*obs_uncertainity_proportion
      }      
    }else{
      init_obs[1,plotnum] = tmp2$Initial_WF
      init_uncert[1,plotnum] = init_obs[1,plotnum]*obs_uncertainity_proportion
    }
    init_obs[2,plotnum] = tmp2$Initial_WS
    init_obs[3,plotnum] = tmp2$Initial_WCR
    init_obs[4,plotnum] = tmp2$Initial_WR
    if(tmp2$Initial_WR_code == 1){
      init_uncert[4,plotnum] = init_obs[4,plotnum]*obs_uncertainity_proportion
    }else{
      init_uncert[4,plotnum] = init_obs[4,plotnum]*0.25
    }
    init_obs[5,plotnum] = tmp2$PlantDensityHa
    if(length(which(tmp$LAI_H != -99)) >0 |(length(which(tmp$FOL_H != -99)))){
      init_obs[6,plotnum] = tmp$LAI_H[which(tmp$LAI_H != -99 & tmp$MonthMeas == 8)][1]  # tmp2$Initial_WF_H * SLA_H * 0.1
      if(is.na(init_obs[6,plotnum])){
        init_obs[6,plotnum] = tmp2$Initial_WF_H * SLA_H * 0.1
      }
    }else{
      init_obs[6,plotnum] = 0.0
    }
    
    init_obs[7,plotnum] = tmp2$Initial_WS_H
    
    init_uncert[2,plotnum] = init_obs[2,plotnum]*0.01
    init_uncert[3,plotnum] = init_obs[3,plotnum]*obs_uncertainity_proportion
    init_uncert[5,plotnum] = init_obs[5,plotnum]*0.01
    init_uncert[6,plotnum] = init_obs[6,plotnum]*obs_uncertainity_proportion
    init_uncert[7,plotnum] = init_obs[7,plotnum]*obs_uncertainity_proportion
  }
  
  obs = array(-99,dim=c(nstreams,nplots,nmonths))
  obs_uncert = array(-99,dim=c(nstreams,nplots,nmonths))
  thin_event = array(0,dim=c(nplots,nmonths))
  
  for(plotnum in 1:nplots){
    
    tmp = observations[which(observations$PlotID == plots[plotnum]),]
    
    ##### STOCKS ##################
    
    # PINE LAI
    if(use_fol_state[plotnum] == 0){
      if(length(tmp$LAI[which(tmp$LAI != -99)]) >0){
        for(i in 2:(length(tmp$LAI[which(tmp$LAI != -99)]))){
          meas_year = tmp$YearMeas[which(tmp$LAI != -99)][i]
          meas_month = tmp$MonthMeas[which(tmp$LAI != -99)][i]
          index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
          obs[1,plotnum,index] = tmp$LAI[which(tmp$LAI != -99)][i]
          #ONLY USE AUGUST LAI VALUES AT DUKE
          #if(meas_month != 8 & tmp$PlotID[1] > 40000 & tmp$PlotID[1] < 41000){obs[1,plotnum,index]=-99}
          obs_uncert[1,plotnum,index] = tmp$LAI[which(tmp$LAI != -99)][i]*obs_uncertainity_proportion
        }
      }
    }else{
      #Pine foliage
      if(length(tmp$FOL[which(tmp$FOL != -99)]) >0){
        for(i in 2:(length(tmp$FOL[which(tmp$FOL != -99)]))){
          meas_year = tmp$YearMeas[which(tmp$FOL != -99)][i]
          meas_month = tmp$MonthMeas[which(tmp$FOL != -99)][i]
          index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
          obs[1,plotnum,index] = tmp$FOL[which(tmp$FOL != -99)][i]
          obs_uncert[1,plotnum,index] = tmp$WFest_sd[i]
          if(tmp$WFest_sd[i] == -99){
            obs_uncert[1,plotnum,index] = tmp$FOL[which(tmp$FOL != -99)][i]*obs_uncertainity_proportion
          }
        }
      }
    }
    
    #STEM
    if(length(tmp$WOODY[which(tmp$WOODY != -99)]) >0){
      for(i in 2:(length(tmp$WOODY[which(tmp$WOODY != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$WOODY != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$WOODY != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[2,plotnum,index] = tmp$WOODY[which(tmp$WOODY != -99)][i]
        obs_uncert[2,plotnum,index] = tmp$WSest_sd[i]
        if(tmp$WSest_sd[i] == -99){
          obs_uncert[2,plotnum,index] = tmp$WOODY[which(tmp$WOODY != -99)][i]*0.01
        }
      }
    }
    #Coarse Roots
    if(length(tmp$ROOT_COARSE[which(tmp$ROOT_COARSE != -99)]) >0){
      for(i in 2:(length(tmp$ROOT_COARSE[which(tmp$ROOT_COARSE != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$ROOT_COARSE != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$ROOT_COARSE != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[3,plotnum,index] = tmp$ROOT_COARSE[which(tmp$ROOT_COARSE != -99)][i]
        obs_uncert[3,plotnum,index] = tmp$ROOT_COARSE[which(tmp$ROOT_COARSE != -99)][i]*obs_uncertainity_proportion
      }
    }
    
    #Fine roots
    if(length(tmp$ROOT_FINE[which(tmp$ROOT_FINE != -99)]) >0){
      for(i in 2:(length(tmp$ROOT_FINE[which(tmp$ROOT_FINE != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$ROOT_FINE != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$ROOT_FINE != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[4,plotnum,index] = tmp$ROOT_FINE[which(tmp$ROOT_FINE != -99)][i]
        obs_uncert[4,plotnum,index] = tmp$ROOT_FINE[which(tmp$ROOT_FINE != -99)][i]*obs_uncertainity_proportion
      }
    }
    
    #Stem Density
    if(length(tmp$Nha[which(tmp$Nha != -99)]) >0){
      for(i in 2:(length(tmp$Nha[which(tmp$Nha != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$Nha != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$Nha != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[5,plotnum,index] = tmp$Nha[which(tmp$Nha != -99)][i]
        obs_uncert[5,plotnum,index] = tmp$Nha[which(tmp$Nha != -99)][i]*0.01
      }
    }
    #Hardwood LAI
    if(length(tmp$LAI_H[which(tmp$LAI_H != -99)]) >0){
      for(i in 1:(length(tmp$LAI_H[which(tmp$LAI_H != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$LAI_H != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$LAI_H != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[6,plotnum,index] = tmp$LAI_H[which(tmp$LAI_H != -99)][i]
        #if(obs[6,plotnum,index]  == 0.0){
        #  obs[6,plotnum,index] = -99
        #}
        if(meas_month != 8 & tmp$PlotID[1] > 40000 & tmp$PlotID[1] < 42000){
          obs[6,plotnum,index]=-99
        }
        
        obs_uncert[6,plotnum,index] = tmp$LAI_H[which(tmp$LAI_H != -99)][i]*obs_uncertainity_proportion
        if(obs_uncert[6,plotnum,index] == 0.0){obs_uncert[6,plotnum,index]= 0.001}
      }
    }
    #Hardwood stem
    if(length(tmp$WOODY_H[which(tmp$WOODY_H != -99)]) >0){
      for(i in 1:(length(tmp$WOODY_H[which(tmp$WOODY_H != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$WOODY_H != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$WOODY_H != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[7,plotnum,index] = tmp$WOODY_H[which(tmp$WOODY_H != -99)][i]
        obs_uncert[7,plotnum,index] = tmp$WOODY_H[which(tmp$WOODY_H != -99)][i]*obs_uncertainity_proportion
      }
    } 
    #Hardwood coarse roots
    
    obs[8,plotnum,] = -99
    obs_uncert[8,plotnum,index] = -99
    
    #Hardwood fine roots
    obs[9,plotnum,] = -99
    obs_uncert[9,plotnum,index] = -99
    
    
    #Hardwood stem density
    obs[10,plotnum,] = -99
    obs_uncert[10,plotnum,index] = -99
    
    #Avialable soil water
    obs[11,plotnum,] = -99
    obs_uncert[11,plotnum,index] = -99
    
    ########  FLUXES ##################
    #GEP
    if(length(tmp$GEP[which(tmp$GEP != -99)]) >0){
      for(i in 1:(length(tmp$GEP[which(tmp$GEP != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$GEP != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$GEP != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[12,plotnum,index] = tmp$GEP[which(tmp$GEP != -99)][i]
        obs_uncert[12,plotnum,index] = tmp$GEP[which(tmp$GEP != -99)][i]*obs_uncertainity_proportion
      }
    }
    #ET
    if(length(tmp$ET[which(tmp$ET != -99)]) >0){
      for(i in 1:(length(tmp$ET[which(tmp$ET != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$ET != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$ET != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[14,plotnum,index] = tmp$ET[which(tmp$ET != -99)][i]
        obs_uncert[14,plotnum,index] = tmp$ET[which(tmp$ET != -99)][i]*obs_uncertainity_proportion
      }
    } 
    #Ctrans
    if(length(tmp$Ctrans[which(tmp$Ctrans != -99)]) >0){
      for(i in 1:(length(tmp$Ctrans[which(tmp$Ctrans != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$Ctrans != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$Ctrans != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[15,plotnum,index] = tmp$Ctrans[which(tmp$Ctrans != -99)][i]
        obs_uncert[15,plotnum,index] = tmp$Ctrans_sd[which(tmp$Ctrans != -99)][i]
      }
    } 
    
    #Hardwood transpiration
    if(length(tmp$Ctrans_H[which(tmp$Ctrans_H != -99)]) >0){
      for(i in 1:(length(tmp$Ctrans_H[which(tmp$Ctrans_H != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$Ctrans_H != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$Ctrans_H != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[16,plotnum,index] = tmp$Ctrans_H[which(tmp$Ctrans_H != -99)][i]
        if((meas_month != 8) & tmp$PlotID[1] > 40000 & tmp$PlotID[1] < 42000){obs[16,plotnum,index]=-99}
        obs_uncert[16,plotnum,index] = tmp$Ctrans_sd[which(tmp$Ctrans_H != -99)][i]
      }
    }  
    
    #Total Root Production
    if(length(tmp$ROOT_PROD_TOTAL[which(tmp$ROOT_PROD_TOTAL != -99)]) >0){
      for(i in 1:(length(tmp$ROOT_PROD_TOTAL[which(tmp$ROOT_PROD_TOTAL != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$ROOT_PROD_TOTAL != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$ROOT_PROD_TOTAL != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[17,plotnum,index] = tmp$ROOT_PROD_TOTAL[which(tmp$ROOT_PROD_TOTAL != -99)][i]
        obs_uncert[17,plotnum,index] = tmp$ROOT_PROD_TOTAL[which(tmp$ROOT_PROD_TOTAL != -99)][i]*obs_uncertainity_proportion
      }
    }
    #Total Foliage Production
    if(length(tmp$FOL_PROD_TOTAL[which(tmp$FOL_PROD_TOTAL != -99)]) >0){
      for(i in 1:(length(tmp$FOL_PROD_TOTAL[which(tmp$FOL_PROD_TOTAL != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$FOL_PROD_TOTAL != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$FOL_PROD_TOTAL != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[18,plotnum,index] = tmp$FOL_PROD_TOTAL[which(tmp$FOL_PROD_TOTAL != -99)][i]
        obs_uncert[18,plotnum,index] = tmp$FOL_PROD_TOTAL[which(tmp$FOL_PROD_TOTAL != -99)][i]*obs_uncertainity_proportion
      }
    }
    
    #Total LAI
    if(length(tmp$LAI_TOTAL[which(tmp$LAI_TOTAL != -99)]) >0){
      for(i in 1:(length(tmp$LAI_TOTAL[which(tmp$LAI_TOTAL != -99)]))){
        meas_year = tmp$YearMeas[which(tmp$LAI_TOTAL != -99)][i]
        meas_month = tmp$MonthMeas[which(tmp$LAI_TOTAL != -99)][i]
        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
        obs[19,plotnum,index] = tmp$LAI_TOTAL[which(tmp$LAI_TOTAL != -99)][i]
        if((meas_month != 8) & tmp$PlotID[1] > 40000 & tmp$PlotID[1] < 42000){obs[19,plotnum,index]=-99}
        obs_uncert[19,plotnum,index] = tmp$LAI_TOTAL[which(tmp$LAI_TOTAL != -99)][i]*obs_uncertainity_proportion
      }
    }
    
    #### EXTRA STOCKS #######
    
    #   #Hardwood foliage
    #    if(length(tmp$FOL_H[which(tmp$FOL_H != -99)]) >0){
    #     for(i in 1:(length(tmp$FOL_H[which(tmp$FOL_H != -99)]))){
    #        meas_year = tmp$YearMeas[which(tmp$FOL_H != -99)][i]
    #       meas_month = tmp$MonthMeas[which(tmp$FOL_H != -99)][i]
    #        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
    #       obs[19,plotnum,index] = tmp$FOL_H[which(tmp$FOL_H != -99)][i]
    #        if(obs[19,plotnum,index]  == 0.0){
    #         obs[19,plotnum,index] = -99
    #        }
    #        obs_uncert[19,plotnum,index] = tmp$FOL_H[which(tmp$FOL_H != -99)][i]*obs_uncertainity_proportion
    #      }
    #    }
    #   #Pine foliage
    #    if(length(tmp$FOL[which(tmp$FOL != -99)]) >0){
    #      for(i in 1:(length(tmp$FOL[which(tmp$FOL != -99)]))){
    #        meas_year = tmp$YearMeas[which(tmp$FOL != -99)][i]
    #        meas_month = tmp$MonthMeas[which(tmp$FOL != -99)][i]
    #        index = (meas_year - earliestYear)*12+1 + (meas_month-earliestMonth)
    #        obs[20,plotnum,index] = tmp$FOL[which(tmp$FOL != -99)][i]
    #        obs_uncert[20,plotnum,index] = tmp$FOL[which(tmp$FOL != -99)][i]*obs_uncertainity_proportion
    #      }
    #    }
    
    
    
  }
  
  mo_start_end=cbind(mo.start,mo.end)
  return(list(mo_start_end = mo_start_end,init_obs=init_obs,init_uncert=init_uncert,years = years,months= months,nmonths= nmonths,obs= obs,thin_event=thin_event,obs_uncert=obs_uncert))
}
