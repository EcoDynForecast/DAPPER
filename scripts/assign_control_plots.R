assign_control_plots <- function(nplots,initdata,plotlist){
  
  control_plot_index = array(NA,dim=c(nplots))
  for(plotnum in 1:nplots){
    plot_init <- initdata[plotnum,]
    control_plot_index[plotnum] = which(plotlist==plot_init$ControlPlotID)
  }
  matched_FR_plot_index = array(99,dim=c(nplots))
  for(plotnum in 1:nplots){
    plot_init <- initdata[plotnum,]
    if(length(which(plotlist==plot_init$MatchedFRPlotID)) > 1){
      matched_FR_plot_index[plotnum] = which(plotlist==plot_init$MatchedFRPlotID)
    }
  }
  return(list(control_plot_index = control_plot_index,matched_FR_plot_index= matched_FR_plot_index))
}