prepare_met <- function(met_in,initdata,mo_start_end,nplots,nmonths,months,years){
  
  co2elev_year = 1996
  co2elev_month = 8
  
  plots = unique(initdata$PlotID)
  tmin = array(-99,dim=c(nplots,nmonths))
  tmax = array(-99,dim=c(nplots,nmonths))
  precip = array(-99,dim=c(nplots,nmonths))
  ra = array(-99,dim=c(nplots,nmonths))
  frost = array(-99,dim=c(nplots,nmonths))
  co2 = array(-99,dim=c(nplots,nmonths))
  for(plotnum in 1:nplots){
    tmp = met_in[which(met_in$SiteID == initdata$SiteID[plotnum]),]
    tmp2 = initdata[which(initdata$PlotID == initdata$PlotID[plotnum]),]
    for(mo in mo_start_end[plotnum,1]:mo_start_end[plotnum,2]){
      curr_year_index = which(tmp$Year == years[mo] & tmp$Month == months[mo])
      tmin[plotnum,mo] = met_in$Tmin[curr_year_index]
      tmax[plotnum,mo] = met_in$Tmax[curr_year_index]
      precip[plotnum,mo] = met_in$Precip[curr_year_index]
      ra[plotnum,mo] = met_in$RA[curr_year_index]
      frost[plotnum,mo] =met_in$Frost[curr_year_index]
      co2[plotnum,mo] = met_in$CO2[curr_year_index]
      
      if(years[mo] >= tmp2$DroughtStart){
        precip[plotnum,mo] = precip[plotnum,mo]#*DroughtLevel
      }
      if(tmp2$CO2flag == 1  & years[mo] == co2elev_year & months[mo] >= co2elev_month){
        co2[plotnum,mo]=tmp2$CO2elev
      }else if(tmp2$CO2flag == 1 & years[mo] > co2elev_year){
        co2[plotnum,mo]=tmp2$CO2elev
      }else if(tmp2$CO2flag == 1 & tmp2$PlotID > 40006){
        co2[plotnum,mo]=tmp2$CO2elev
      }
    }
  }
  
  return(list(tmin=tmin,tmax=tmax,precip=precip,ra=ra,frost=frost,co2=co2))
}




