prepare_met <- function(met_in,initdata,mo_start_end,co2_in,nplots,nmonths,months,years){
  
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
      curr_year_index = which(tmp$YEAR == years[mo])
      co2[plotnum,mo]= co2_in$CO2_Concentration_RCP85[which(co2_in$Year == years[mo])]
      if(months[mo] == 1){
        tmin[plotnum,mo] = tmp$TMIN1[curr_year_index]
        tmax[plotnum,mo] = tmp$TMAX1[curr_year_index]
        precip[plotnum,mo] = tmp$PRECIP1[curr_year_index]
        ra[plotnum,mo] = tmp$RA1[curr_year_index]
        frost[plotnum,mo] = tmp$FROST1[curr_year_index]
      }
      if(months[mo] == 2){
        tmin[plotnum,mo] = tmp$TMIN2[curr_year_index]
        tmax[plotnum,mo] = tmp$TMAX2[curr_year_index]
        precip[plotnum,mo] = tmp$PRECIP2[curr_year_index]
        ra[plotnum,mo] = tmp$RA2[curr_year_index]
        frost[plotnum,mo] = tmp$FROST2[curr_year_index]
      }
      if(months[mo] == 3){
        tmin[plotnum,mo] = tmp$TMIN3[curr_year_index]
        tmax[plotnum,mo] = tmp$TMAX3[curr_year_index]
        precip[plotnum,mo] = tmp$PRECIP3[curr_year_index]
        ra[plotnum,mo] = tmp$RA3[curr_year_index]
        frost[plotnum,mo] = tmp$FROST3[curr_year_index] 
      }
      if(months[mo] == 4){
        tmin[plotnum,mo] = tmp$TMIN4[curr_year_index]
        tmax[plotnum,mo] = tmp$TMAX4[curr_year_index]
        precip[plotnum,mo] = tmp$PRECIP4[curr_year_index]
        ra[plotnum,mo] = tmp$RA4[curr_year_index]
        frost[plotnum,mo] = tmp$FROST4[curr_year_index]
      }
      if(months[mo] == 5){
        tmin[plotnum,mo] = tmp$TMIN5[curr_year_index]
        tmax[plotnum,mo] = tmp$TMAX5[curr_year_index]
        precip[plotnum,mo] = tmp$PRECIP5[curr_year_index]
        ra[plotnum,mo] = tmp$RA5[curr_year_index]
        frost[plotnum,mo] = tmp$FROST5[curr_year_index]   
      }
      if(months[mo] == 6){
        tmin[plotnum,mo] = tmp$TMIN6[curr_year_index]
        tmax[plotnum,mo] = tmp$TMAX6[curr_year_index]
        precip[plotnum,mo] = tmp$PRECIP6[curr_year_index]
        ra[plotnum,mo] = tmp$RA6[curr_year_index]
        frost[plotnum,mo] = tmp$FROST6[curr_year_index]   
      }
      if(months[mo] == 7){
        tmin[plotnum,mo] = tmp$TMIN7[curr_year_index]
        tmax[plotnum,mo] = tmp$TMAX7[curr_year_index]
        precip[plotnum,mo] = tmp$PRECIP7[curr_year_index]
        ra[plotnum,mo] = tmp$RA7[curr_year_index]
        frost[plotnum,mo] = tmp$FROST7[curr_year_index]
      }
      if(months[mo] == 8){
        tmin[plotnum,mo] = tmp$TMIN8[curr_year_index]
        tmax[plotnum,mo] = tmp$TMAX8[curr_year_index]
        precip[plotnum,mo] = tmp$PRECIP8[curr_year_index]
        ra[plotnum,mo] = tmp$RA8[curr_year_index]
        frost[plotnum,mo] = tmp$FROST8[curr_year_index] 
      }
      if(months[mo] == 9){
        tmin[plotnum,mo] = tmp$TMIN9[curr_year_index]
        tmax[plotnum,mo] = tmp$TMAX9[curr_year_index]
        precip[plotnum,mo] = tmp$PRECIP9[curr_year_index]
        ra[plotnum,mo] = tmp$RA9[curr_year_index]
        frost[plotnum,mo] = tmp$FROST9[curr_year_index]    
      }
      if(months[mo] == 10){
        tmin[plotnum,mo] = tmp$TMIN10[curr_year_index]
        tmax[plotnum,mo] = tmp$TMAX10[curr_year_index]
        precip[plotnum,mo] = tmp$PRECIP10[curr_year_index]
        ra[plotnum,mo] = tmp$RA10[curr_year_index]
        frost[plotnum,mo] = tmp$FROST10[curr_year_index] 
      }
      if(months[mo] == 11){
        tmin[plotnum,mo] = tmp$TMIN11[curr_year_index]
        tmax[plotnum,mo] = tmp$TMAX11[curr_year_index]
        precip[plotnum,mo] = tmp$PRECIP11[curr_year_index]
        ra[plotnum,mo] = tmp$RA11[curr_year_index]
        frost[plotnum,mo] = tmp$FROST11[curr_year_index] 
      }
      if(months[mo] == 12){
        tmin[plotnum,mo] = tmp$TMIN12[curr_year_index]
        tmax[plotnum,mo] = tmp$TMAX12[curr_year_index]
        precip[plotnum,mo] = tmp$PRECIP12[curr_year_index]
        ra[plotnum,mo] = tmp$RA12[curr_year_index]
        frost[plotnum,mo] = tmp$FROST12[curr_year_index]   
      }
      
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




