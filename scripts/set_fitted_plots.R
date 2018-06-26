set_fitted_plots <- function(nplots,val_set,plotlist,initdata){
  fit_plot = rep(1,nplots)
  if(!is.na(validation_set_file)){
    fit_set_all = read.csv(paste0(working_directory,'/',validation_set_file))  
    
    if(val_set > 0){
      not_fit_set = fit_set_all$PlotID[which(fit_set_all[,1] == val_set)]
      not_fit_set_index = rep(NA,length(not_fit_set))
      for(i in 1:length(not_fit_set)){
        not_fit_set_index[i] = which(plotlist == not_fit_set[i])
      }
      fit_plot[not_fit_set_index] = 0
    }
    
    if(obs_set == 1){
      not_fit_set = which((initdata$FertFlag == 1 | initdata$DroughtLevel < 1.0 | initdata$IrrFlag == 1 | initdata$CO2flag == 1))
      fit_plot[not_fit_set] = 0
    } 
    
    if(obs_set == 7){
      not_fit_set = which((initdata$FertFlag == 1))
      fit_plot[not_fit_set] = 0
    }
    
    if(obs_set == 5){
      not_fit_set = which(initdata$DroughtLevel < 1.0 | initdata$IrrFlag == 1)
      fit_plot[not_fit_set] = 0
    }
    
    if(obs_set == 6){
      not_fit_set = which((initdata$CO2flag == 1)) 
      fit_plot[not_fit_set] = 0
    }
    
    if(obs_set == 19){
      not_fit_set = which((initdata$CO2flag == 1 & initdata$FertFlag == 1) | (initdata$DroughtLevel < 1.0 & initdata$FertFlag == 1) | (initdata$IrrFlag == 1.0 & initdata$FertFlag == 1))
      fit_plot[not_fit_set] = 0
    }
  }
  return(fit_plot)
}
