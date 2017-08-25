plot_output <- function(){
  
  pdf(final_pdf,height = 3.4252*3.0,width = 3.4252*3.0)
  
  ylab = rep(NA,7)
  ylab[1] = 'LAI'
  ylab[2] = 'Stem'
  ylab[3] = 'Coarse root'
  ylab[4] = 'Fine root'
  ylab[5] = 'Stem density'
  ylab[6] = 'LAI Hard'
  ylab[7] = 'Stem Hard'
  
  for(data_stream in 1:7){
    tmp = tracked_plot[1,data_stream,]
    tmp[which(is.nan(tmp))] = -99
    tmp[which(is.na(tmp))] = -99
    tmp1 = tmp
    tmp2 = tmp
    index = 1
    for(mo in mo_start_end[tracked_plotnum,1]:mo_start_end[tracked_plotnum,2]){
      tmp[mo] = median(tracked_plot[,data_stream,mo],na.rm=TRUE)
      tmp1[mo] = quantile(tracked_plot[,data_stream,mo],0.025,na.rm=TRUE)
      tmp2[mo] = quantile(tracked_plot[,data_stream,mo],0.975,na.rm=TRUE)
    }
    
    #iter = 1
    #plot(age[tracked_plotnum,mo_start_end[tracked_plotnum,1]:mo_start_end[tracked_plotnum,2]],tracked_plot[iter,data_stream,],type='l')
    #points(age[tracked_plotnum,mo_start_end[tracked_plotnum,1]:mo_start_end[tracked_plotnum,2]],
    #       obs[data_stream,tracked_plotnum,mo_start_end[tracked_plotnum,1]:mo_start_end[tracked_plotnum,2]],col='red')
    #iter = iter + 5
    
    plot(age[tracked_plotnum,mo_start_end[tracked_plotnum,1]:mo_start_end[tracked_plotnum,2]],tmp[mo_start_end[tracked_plotnum,1]:mo_start_end[tracked_plotnum,2]],type='l',ylim=c(0,max(c(tmp2,1),na.rm=TRUE)),main=paste('data stream:',data_stream,sep=' '),col='white',xlab='age',ylab=ylab[data_stream])
    for(mo in mo_start_end[tracked_plotnum,1]:mo_start_end[tracked_plotnum,2]){
      if(latent[data_stream,tracked_plotnum,mo] != -99 & !is.nan(latent[data_stream,tracked_plotnum,mo]) ){
        points(age[tracked_plotnum,mo],tmp[mo],ylim=c(0,max(tmp2)),main=paste('data stream:',data_stream,sep=' '))
        segments(age[tracked_plotnum,mo],tmp1[mo],age[tracked_plotnum,mo],tmp2[mo])
      }
      if(obs[data_stream,tracked_plotnum,mo] != -99){
        points(age[tracked_plotnum,mo],obs[data_stream,tracked_plotnum,mo],ylim=c(0,max(tmp2)),col='red')
      }
    }
    
    points(age[tracked_plotnum,mo_start_end[tracked_plotnum,1]],
           tmp[mo_start_end[tracked_plotnum,1]],col='red')
    segments(age[tracked_plotnum,mo_start_end[tracked_plotnum,1]],tmp1[mo_start_end[tracked_plotnum,1]],age[tracked_plotnum,mo_start_end[tracked_plotnum,1]],tmp2[mo_start_end[tracked_plotnum,1]])
    
  }
  
  max_iter = dim(accepted_pars_thinned_burned)[1]
  
  new_pars = rep(NA,dim(accepted_pars_thinned_burned)[2])
  
  for(i in 1:dim(accepted_pars_thinned_burned)[2]){
    new_pars[i] = median(accepted_pars_thinned_burned[,i])
  }
  #new_pars = accepted_pars_thinned_burned[dim(accepted_pars_thinned_burned)[1],]
  
  FR1_new = new_pars[49]
  FR2_new = new_pars[50]
  FR3_new = new_pars[51]
  
  pars=array(NA,npars-9)
  pars = new_pars[1:npars_used_by_fortran]
  pars=as.array(pars)
  pars[is.na(pars)] <- 0
  
  fol_max = max(observations$FOL)+5
  wood_max = max(observations$WOODY)+10
  Nha_max = max(observations$Nha)+200
  fol_min = 0 #min(observations$FOL)
  wood_min = 0 #min(observations$WOODY)
  Nha_min = 0 #min(observations$Nha)
  age_max = max(observations$AgeMeas)
  age_min = min(observations$AgeMeas)
  
  val_RMSE = rep(NA,nplots)
  stem_predicted  = NULL
  stem_observed  = NULL
  stem_plotid  = NULL
  
  root_predicted  = NULL
  root_observed  = NULL
  root_plotid  = NULL
  
  coarse_root_predicted  = NULL
  coarse_root_observed = NULL
  coarse_root_plotid  = NULL
  
  fol_predicted = NULL
  fol_observed = NULL
  fol_plotid = NULL
  
  lai_predicted = NULL
  lai_observed  = NULL
  lai_plotid = NULL
  
  nha_predicted = NULL
  nha_observed = NULL
  nha_plotid = NULL
  
  gep_predicted = NULL
  gep_observed = NULL
  gep_plotid = NULL
  gep_year= NULL
  gep_month = NULL
  
  et_predicted = NULL
  et_observed = NULL
  et_plotid = NULL
  et_year= NULL
  et_month = NULL
  
  ctrans_pine_predicted = NULL
  ctrans_pine_observed = NULL
  ctrans_pine_plotid = NULL
  ctrans_pine_year= NULL
  ctrans_pine_month = NULL
  
  ctrans_hard_predicted = NULL
  ctrans_hard_observed = NULL
  ctrans_hard_plotid = NULL
  ctrans_hard_year= NULL
  ctrans_hard_month = NULL
  
  for(plotnum in 1:nplots){
    
    tmp_initdata = initdata[plotnum,]
    
    pars = new_pars[1:npars_used_by_fortran]
    
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
    
    index1 = index_guide[3]+plotnum-1
    new_FR = new_pars[index1]
    control_plotnum = control_plot_index[plotnum]
    index = index_guide[3] + control_plotnum - 1
    control_FR = new_pars[index]
    
    PlantedYear = Planting_year
    PlantedMonth = PlantMonth
    
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
    
    if(use_fol_state[plotnum] == 1){
      LAI = -99
    }else{
      LAI = WFi* SLA * 0.1
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
                LAI = tmp_initdata$Initial_LAI,
                LAI_H = init_obs[6,plotnum]
    )
    
    site = array(site_in)
    
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
    
    dyn.load(code_library_plot)
    
    output = NULL
    for(mo in (mo_start_end[plotnum,1]+1):mo_start_end[plotnum,2]){
      
      site[3] = years[mo]
      site[4]  = months[mo]
    
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
      
      output_mo=array(tmp$out_var, dim=c(1,output_dim))
      
      modeled = rep(NA,nstreams)
      #--- PROCESS MODEL OUTPUT      
      if(use_fol_state[plotnum] == 1){	
      modeled[1] = output_mo[22] #output(4) !Pine Foliage
      }else{
        modeled[1] = output_mo[4] #output(4) !Pine LAI
      }
      
      
      modeled[2] = output_mo[5] #WS
      modeled[3] = output_mo[6] #+ output_mo[43] !WCR
      modeled[4] = output_mo[7] #+ output_mo[58] !WR
      modeled[5] = output_mo[8] #Stem Density
      modeled[6] = output_mo[9] #LAI_H        
      modeled[7] = output_mo[10] #Hardwood Stem   
      modeled[8] = output_mo[11] #Hardwood Coarse roots
      modeled[9] = output_mo[12] #Hardwood Fine roots 
      modeled[10] = output_mo[13] #Hardwood Stem density 
      modeled[11] = output_mo[14] # ASW        
      modeled[12] = output_mo[15] # GEP        
      modeled[13] = output_mo[16] # NEE
      modeled[14] = output_mo[17] # ET  
      modeled[15] = output_mo[18] # Ctrans Pine           
      modeled[16] = output_mo[19] # Ctrans Hardwood   
      modeled[19] =  output_mo[4]  + output_mo[9]
      
      output = rbind(output,output_mo)
      
      site[5] =  output_mo[3]+ (1.0/12.) #Age
      
      if(state_space > 0 & fit_plot[plotnum] == 1){
        if(obs[1,plotnum,mo] != -99.0){
          if(use_fol_state[plotnum] == 1){
            site[6] = obs[1,plotnum,mo]
            site[26] = -99
          }else{
            site[6] = output_mo[22]
            site[26] = obs[1,plotnum,mo]
          }
        }else{
          if(use_fol_state[plotnum] == 1){
            site[6] = modeled[1]
            site[26] = -99
          }else{
            site[6] = output_mo[22]
            site[26] = modeled[1]
          }
        }
        if(obs[2,plotnum,mo] != -99.0){
          site[8] = obs[2,plotnum,mo] #WSi
        }else{
          site[8] = modeled[2]
        }   
        
        if(obs[3,plotnum,mo] != -99.0){
          site[20] = obs[3,plotnum,mo] #WCR
        }else{
          site[20] = modeled[3]
        }  
        
        if(obs[4,plotnum,mo] != -99.0){
          site[7] = obs[4,plotnum,mo]* (modeled[4]/(modeled[9]+modeled[4])) 
          site[19] =obs[4,plotnum,mo]* (modeled[9]/(modeled[9]+modeled[4]))
        }else{
          site[7] = modeled[4]
          site[19] = modeled[9]
        }        
        
        if(obs[5,plotnum,mo] != -99.0){
          site[9] = obs[5,plotnum,mo] #StemNo
        }else{
          site[9] = modeled[5]
        }
        
        if(obs[6,plotnum,mo] != -99.0){
          site[27] = obs[6,plotnum,mo] # Hardwood LAI
          site[25] =  output_mo[26]
        }else{
          site[27] = modeled[6]
          site[25] =  output_mo[26]
        }
        
        if(obs[7,plotnum,mo] != -99.0){
          site[18] = obs[7,plotnum,mo] # Hardwood WS
        }else{
          site[18] = modeled[7]
        }
        
        
      }else{
        if(use_fol_state[plotnum] == 1){
          site[6] = modeled[1]
          site[26] = -99
        }else{
          site[6] = output_mo[22]
          site[26] = modeled[1]
        }
        site[8] =  modeled[2] #WS
        site[20] = modeled[3] #WCR
        site[7] = modeled[4] #WRi
        site[9] = modeled[5] #StemNo
        site[27] = modeled[6]
        site[18] = modeled[7]
        site[19] = modeled[9] #WR_H
        site[25] =  output_mo[26]
      }
      
      #Other state variables
      
      site[10] = modeled[11] #ASW
      site[17] = output_mo[23] #WF_H
      site[24] =  output_mo[11]#WCR_h
      
    }
    
    if(fit_plot[plotnum] == 0){
      for(mo in 1:nomonths_plot){
        if(obs[2,plotnum,mo] != -99){
          index = mo_start_end[plotnum,1]+mo
          val_RMSE[plotnum] = sqrt((obs[2,plotnum,mo] - output[mo,6])^2)
        }
      }
    }
    
    if(nplots > 1){
      # 1 - year, 2-month, 3-stand age, 4 - WF, 5 - WR, 6-WS, 7-StemNo
      for(mo in 2:(nomonths_plot-1)){
        index = mo_start_end[plotnum,1]+mo
        if(obs[1,plotnum,index] != -99 & use_fol_state[plotnum] ==  0){
          lai_predicted = c(lai_predicted,output[mo,4])
          lai_observed = c(lai_observed,obs[1,plotnum,index])
          lai_plotid = c(lai_plotid,PlotID)
        }
        if(obs[2,plotnum,index] != -99){
          stem_predicted = c(stem_predicted,output[mo,5])
          stem_observed = c(stem_observed,obs[2,plotnum,index])
          stem_plotid = c(stem_plotid,PlotID)
        }
        
        if(obs[3,plotnum,index] != -99){
          coarse_root_predicted = c(coarse_root_predicted,output[mo,6])
          coarse_root_observed = c(coarse_root_observed,obs[3,plotnum,index])
          coarse_root_plotid = c(coarse_root_plotid,PlotID)
        }
        if(obs[4,plotnum,index] != -99){
          root_predicted = c(root_predicted,output[mo,7] + output[mo,12])
          root_observed = c(root_observed,obs[4,plotnum,index])
          root_plotid = c(root_plotid,PlotID)
        }
        if(obs[5,plotnum,index] != -99){
          nha_predicted = c(nha_predicted,output[mo,8])
          nha_observed = c(nha_observed,obs[5,plotnum,index])
          nha_plotid = c(nha_plotid,PlotID)
        }
        
        if(obs[12,plotnum,index] != -99){
          gep_predicted = c(gep_predicted,output[mo,15])
          gep_observed = c(gep_observed,obs[12,plotnum,index])
          gep_plotid = c(gep_plotid,PlotID)
          gep_year = c(gep_year,years[index])
          gep_month = c(gep_month,months[index])
        }
        if(obs[14,plotnum,index] != -99){
          et_predicted = c(et_predicted,output[mo,17])
          et_observed = c(et_observed,obs[14,plotnum,index])
          et_plotid = c(et_plotid,PlotID)
          et_year = c(et_year,years[index])
          et_month = c(et_month,months[index])
        }
        if(obs[15,plotnum,index] != -99){
          ctrans_pine_predicted = c(ctrans_pine_predicted,output[mo,18])
          ctrans_pine_observed = c(ctrans_pine_observed,obs[15,plotnum,index])
          ctrans_pine_plotid = c(ctrans_pine_plotid,PlotID)
          ctrans_pine_year = c(ctrans_pine_year,years[index])
          ctrans_pine_month = c(ctrans_pine_month,months[index])
        }
        if(obs[16,plotnum,index] != -99){
          ctrans_hard_predicted = c(ctrans_hard_predicted,output[mo,19])
          ctrans_hard_observed = c(ctrans_hard_observed,obs[16,plotnum,index])
          ctrans_hard_plotid = c(ctrans_hard_plotid,PlotID)
          ctrans_hard_year = c(ctrans_hard_year,years[index])
          ctrans_hard_month = c(ctrans_hard_month,months[index])
        }
      }
    }
    
    #par(mfrow=c(4,4))
    par(mfrow=c(4,5),mar = c(4,4,2,2),oma = c(3,3,2,2))
    xlim_range = c(0,30) #c(min(output[,3])-1,max(output[,3])+1)
    modeled = output[which(!is.nan(output[,6])),6]
    modeled_y = output[which(!is.nan(output[,6])),3]
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
    
    if(use_fol_state[plotnum] ==  0){
      modeled[1,] = output[,4] #Pine LAI
      ylab[1] = 'LAI'
    }else{
      modeled[1,] = output[,22] #Pine Foliage
      ylab[1] = 'Foliage'
    }
    
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
    
    for(data_stream in 1:(nstreams-1)){
      if(data_stream != 8 & data_stream != 9 & data_stream != 10 & data_stream != 11 & data_stream != 13 & data_stream != 17 & data_stream != 18){
        observed = obs[data_stream,plotnum,which(obs[data_stream,plotnum,]!=-99)]
        observed_y =age[plotnum,which(obs[data_stream,plotnum,]!=-99)]
        xlim_range = c(0,30) #c(min(output[,3])-1,max(output[,3])+1)
        ylim_range = c(0,max(c(observed,modeled[data_stream,]),na.rm=TRUE)) #range(c(observed,modeled[data_stream,]))
        plot(modeled_age[which(!is.nan(modeled[data_stream,]))],modeled[data_stream,which(!is.nan(modeled[data_stream,]))],type='l',xlim=xlim_range,ylim=ylim_range,col='black',xlab='plot age (yr)',ylab=ylab[data_stream])
        points(observed_y,observed,col='gray',pch=20)
        if(data_stream == 1) {title(paste(plotlist[plotnum],tmp_initdata$StudyName,sep=' '))}
        if(data_stream == 2) {title(paste('Treat: ',tmp_initdata$Treatment,sep=''))}
        if(data_stream == 3) {title(paste('plotnum: ',plotnum,sep=''))} 
        if(data_stream == 4) {title(paste('Fit site?: ',fit_plot[plotnum],sep=''))} 
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
    
    index1 = index_guide[3]+plotnum-1
    plot(accepted_pars_thinned_burned[,index1],xlab='iterations',ylab='FR',type='l',col='red',ylim=c(0,1))
    index1 = index_guide[3]+control_plot_index[plotnum]-1
    points(accepted_pars_thinned_burned[,index1],xlab='iterations',ylab='FR',type='l',col='black')
    legend('bottomright',legend=c('FR','FR (control plot)'),lty=c(1,1),col=c('red','black'),cex=1.0,bty='n')
    #legend('topright',legend=c('FR:current','FR: control','pFS2','pFS20'),lty=c(1,1,2,2),col=c('red','black','green','blue'),cex=0.5)
    
    index1 = index_guide[5]+plotnum-1
    plot(accepted_pars_thinned_burned[,index1],xlab='iterations',ylab='wsx1000',type='l',lty=1,col='black')
    index1 = index_guide[7]+plotnum-1
    plot(accepted_pars_thinned_burned[,index1],xlab='iterations',ylab='thinpower',type='l',lty =1,col='black')
    index1 = index_guide[9]+plotnum-1
    plot(accepted_pars_thinned_burned[,index1],xlab='iterations',ylab='mort_rate',type='l',lty =1,col='black')
  }
  
  
  par(mfrow=c(4,4),mar=c(5, 4, 4, 2), oma=c(0,0,1,0))
  plot(like_chain,type = 'l',ylab = 'Likelihood',xlab='iteration')
  
  for(i in 1:npars){
    plot(accepted_pars_thinned_burned[1:max_iter,i],type='l',ylab = parnames[i],xlab = 'iteration',main=paste(parnames[i]))
  }
  par(mfrow=c(4,4),mar=c(5, 4, 4, 2), oma=c(0,0,1,0))
  hist(like_chain,ylab='frequency')
  for(i in 1:npars){
    hist(accepted_pars_thinned_burned[1:max_iter,i],xlab = parnames[i],ylab = 'frequency',main = parnames[i])
  }
  
  par(mfrow=c(3,3),mar=c(5, 4, 4, 3), oma=c(0,0,1,0))
  #plot(stem_estimated[,4],stem_estimated_predicted,xlab='observed',ylab='predicted',ylim=c(0,400),xlim=c(0,400))
  #abline(0,1)
  plot(stem_observed,stem_predicted,xlab='observed',ylab='predicted',ylim=c(0,400),xlim=c(0,400),main=c('Stem Biomass'))
  abline(0,1)
  
  
  # if(length(which(initdata[,16] < 1.0)) > 0){
  #   points(stem_max_observed[which(initdata[,16] < 1.0)],stem_max_predicted[which(initdata[,16] < 1.0)],col='blue')
  # }
  # if(length(which(fit_plot == 0))>0){
  #   points(stem_observed[which(fit_plot == 0)],stem_max_predicted[which(fit_plot == 0)],col='red')
  #   RMSE = sqrt(mean((stem_max_predicted[which(fit_plot == 0)]- stem_max_observed[which(fit_plot == 0)])^2))
  #   SSres = sum((stem_max_predicted[which(fit_plot == 0)]- stem_max_observed[which(fit_plot == 0)])^2)
  #   SStot = sum((stem_max_observed[which(fit_plot == 0)] - mean(stem_max_observed[which(fit_plot == 0)]))^2)
  #   nonfit_R2 = 1-(SSres/SStot)
  #   SSres = sum((stem_max_predicted[which(fit_plot == 1)]- stem_max_observed[which(fit_plot == 1)])^2)
  #   SStot = sum((stem_max_observed[which(fit_plot == 1)] - mean(stem_max_observed[which(fit_plot == 1)]))^2)
  #   fit_R2 = 1-(SSres/SStot)
  #   legend('topleft',c(paste('Val. R2 = ',round(nonfit_R2,2),sep=''),paste('Fit R2 = ',round(fit_R2,2),sep=''),paste('Val. RMSE = ',(round(RMSE,2)),sep='')),bty='n')
  # }else{
  #   SSres = sum((stem_max_predicted[which(fit_plot == 1)]- stem_max_observed[which(fit_plot == 1)])^2)
  #   SStot = sum((stem_max_observed[which(fit_plot == 1)] - mean(stem_max_observed[which(fit_plot == 1)]))^2)
  #   fit_R2 = 1-(SSres/SStot)
  #   legend('topleft',c(paste('fit R2 = ',fit_R2,sep='')),bty='n')
  # }
  
  abline(0,1)
  plot(lai_observed,lai_predicted,xlab='observed',ylab='predicted',ylim=c(0,6),xlim=c(0,6),main=c('LAI'))
  abline(0,1)
  plot(coarse_root_observed,coarse_root_predicted,xlab='observed',ylab='predicted',ylim=c(0,75),xlim=c(0,75),main=c('Coarse Root Biomass'))
  abline(0,1)
  plot(root_observed,root_predicted,xlab='observed',ylab='predicted',ylim=c(0,10),xlim=c(0,10),main=c('Fine Root Biomass'))
  abline(0,1)
  plot(nha_observed,nha_predicted,xlab='observed',ylab='predicted',ylim=c(0,2000),xlim=c(0,2000),main=c('Stem Density'))
  abline(0,1)
  plot(gep_observed,gep_predicted,xlab='observed',ylab='predicted',ylim=c(0,12),xlim=c(0,12),main=c('GEP'))
  abline(0,1)
  plot(et_observed,et_predicted,xlab='observed',ylab='predicted',ylim=c(0,200),xlim=c(0,200),main=c('ET'))
  abline(0,1)
  plot(ctrans_pine_observed,ctrans_pine_predicted,xlab='observed',ylab='predicted',ylim=c(0,100),xlim=c(0,100),main=c('Ctrans Pine'))
  abline(0,1)
  plot(ctrans_hard_observed,ctrans_hard_predicted,xlab='observed',ylab='predicted',ylim=c(0,100),xlim=c(0,100),main=c('Ctrans Hard'))
  abline(0,1)
  dev.off()
  
  fname = paste(working_directory,'/figures/',sep='')
  write.table(data.frame(PlotID = initdata[,1],fit_plot,FertFlag=initdata[,18],IrrFlag=initdata[,33],DroughtLevel=initdata[,16],
                         CO2flag=initdata[,19],InitFR = initdata[,12],Mean_temp = initdata[,26],mean_precip = initdata[,27]),
              paste(fname,run_name,'_predicted_vs_observed.csv',sep=""),sep=",",col.names = TRUE,row.names = FALSE)
  write.table(data.frame(stem_plotid,stem_observed,stem_predicted),
              paste(fname,run_name,'_stem_predicted_vs_observed.csv',sep=""),sep=",",col.names = TRUE,row.names = FALSE)
  write.table(data.frame(lai_plotid,lai_observed,lai_predicted),
              paste(fname,run_name,'_LAI_predicted_vs_observed.csv',sep=""),sep=",",col.names = TRUE,row.names = FALSE)
  write.table(data.frame(coarse_root_plotid,coarse_root_observed,coarse_root_predicted),
              paste(fname,run_name,'_coarse_root_predicted_vs_observed.csv',sep=""),sep=",",col.names = TRUE,row.names = FALSE)
  write.table(data.frame(root_plotid,root_observed,root_predicted),
              paste(fname,run_name,'_fine_root_predicted_vs_observed.csv',sep=""),sep=",",col.names = TRUE,row.names = FALSE)
  write.table(data.frame(nha_plotid,nha_observed,nha_predicted),
              paste(fname,run_name,'_nha_predicted_vs_observed.csv',sep=""),sep=",",col.names = TRUE,row.names = FALSE)
  write.table(data.frame(gep_plotid,gep_year,gep_month,gep_observed,gep_predicted),
              paste(fname,run_name,'_GEP_predicted_vs_observed.csv',sep=""),sep=",",col.names = TRUE,row.names = FALSE)
  write.table(data.frame(et_plotid,et_year,et_month,et_observed,et_predicted),
              paste(fname,run_name,'_ET_predicted_vs_observed.csv',sep=""),sep=",",col.names = TRUE,row.names = FALSE)
  write.table(data.frame(ctrans_pine_plotid,ctrans_pine_year,ctrans_pine_month,ctrans_pine_observed,ctrans_hard_predicted,ctrans_hard_plotid,ctrans_hard_year,ctrans_hard_month,ctrans_hard_observed,ctrans_hard_predicted),
              paste(fname,run_name,'_Ctrans_predicted_vs_observed.csv',sep=""),sep=",",col.names = TRUE,row.names = FALSE)
}
