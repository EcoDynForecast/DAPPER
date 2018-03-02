initialize_pars <- function(index_guide,priormatrix,jump_size_init,observations,initdata,fr_model,chain_number,process_model_pars){
  
  plots = unique(observations$PlotID)
  
  max_pars = max(index_guide)
  init_pars = array(NA,dim=c(max_pars))
  new_pars = array(NA,dim=c(max_pars))
  jump_pars = array(NA,dim=c(max_pars))
  prior_parameter1= array(NA,dim=c(max_pars))
  prior_parameter2= array(NA,dim=c(max_pars)) 
  prior_dist= array(NA,dim=c(max_pars))
  fix_par= array(NA,dim=c(max_pars))
  par_group= array(NA,dim=c(max_pars))
  
  plot_index = array(NA,dim=c(max_pars))
  meas_index = array(NA,dim=c(max_pars))
  
  init_pars[index_guide[1]:index_guide[2]] = priormatrix[,1]
  prior_parameter1[index_guide[1]:index_guide[2]]= priormatrix[,2] 
  prior_parameter2[index_guide[1]:index_guide[2]]= priormatrix[,3] 
  prior_dist[index_guide[1]:index_guide[2]] =priormatrix[,4] 
  fix_par[index_guide[1]:index_guide[2]] = priormatrix[,5]
  jump_pars[index_guide[1]:index_guide[2]]=priormatrix[,1]*jump_size_init
  par_group[index_guide[1]:index_guide[2]]= priormatrix[,6]
  
  #FR
  index = index_guide[3]
  for(plotnum in 1:nplots){
    plot_obs <- observations[which(observations$PlotID == plots[plotnum]),]
    plot_init <- initdata[plotnum,]
    plot_index[index]=plotnum
    #USE DIFFERENT INITIALIZATIONS FOR DIFFERENT CHAINS
    if(chain_number == 1 | chain_number > 4){
      init_pars[index]=0.80
      if(plot_init$FertFlag == 1){
        init_pars[index]=0.99
      }
    }else if(chain_number == 2){
      init_pars[index]=0.70
      if(plot_init$FertFlag == 1){
        init_pars[index]=0.85
      }
    }else if(chain_number == 3){
      init_pars[index]=0.75
      if(plot_init$FertFlag == 1){
        init_pars[index]=0.90
      }
    }else if(chain_number == 4){
      init_pars[index]=0.85
      if(plot_init$FertFlag == 1){
        init_pars[index]=0.97
      }
      
    }
    prior_parameter1[index]= 0.1
    prior_parameter2[index]= 1
    prior_dist[index] =1 
    fix_par[index] = 0
    if(fr_model == 2 & plot_init$FertFlag == 0.0 & plot_init$FR == -99){
      fix_par[index] = 1
    }
    if(plot_init$FR == 1.0){
      fix_par[index] = 1
    }
    
    jump_pars[index] = 0.01
    if(FR_separate_npar_groups == 0){
      par_group[index] = FR_npar_group
    }else if(fr_model == 1 & FR_separate_npar_groups == 1){
      if(plot_init$PlotID > 10000 &  plot_init$PlotID < 20000){  #Thinning
        par_group[index] = FR_npar_group
      }else if(plot_init$PlotID > 20000 &  plot_init$PlotID < 30000){ #RW18
        par_group[index] = FR_npar_group +1
      }else if(plot_init$PlotID > 30000 &  plot_init$PlotID < 40000){ #TEIR3
        par_group[index] = FR_npar_group + 2
      }else if(plot_init$PlotID > 40000 &  plot_init$PlotID < 41000){  #DUKE
        par_group[index] = FR_npar_group + 3
      }else if(plot_init$PlotID > 41000 &  plot_init$PlotID < 42000){ #NC2
        par_group[index] = FR_npar_group + 4
      }else if(plot_init$PlotID > 42000 &  plot_init$PlotID < 43000){ #SETRES
        par_group[index] = FR_npar_group + 5
      }else if(plot_init$PlotID > 43000 &  plot_init$PlotID < 44000){ #WAYCROSS
        par_group[index] = FR_npar_group + 6
      }else if(plot_init$PlotID > 72001 & plot_init$PlotID < 72076 ){ 
        par_group[index] = FR_npar_group + 7
      }else if(plot_init$PlotID >= 52001 & plot_init$PlotID < 52467){ 
        par_group[index] = FR_npar_group + 8
      }
    }else if(fr_model == 2 & FR_separate_npar_groups == 1){
      if(plot_init$PlotID > 10000 &  plot_init$PlotID < 20000){  #Thinning
        par_group[index] = FR_npar_group
      }else if(plot_init$PlotID > 20000 &  plot_init$PlotID < 30000){ #RW18
        par_group[index] = FR_npar_group
      }else if(plot_init$PlotID > 30000 &  plot_init$PlotID < 40000){ #TEIR3
        par_group[index] = FR_npar_group + 1
      }else if(plot_init$PlotID > 40000 &  plot_init$PlotID < 41000){  #DUKE
        par_group[index] = FR_npar_group + 2
      }else if(plot_init$PlotID > 41000 &  plot_init$PlotID < 42000){ #NC2
        par_group[index] = FR_npar_group
      }else if(plot_init$PlotID > 42000 &  plot_init$PlotID < 43000){ #SETRES
        par_group[index] = FR_npar_group + 3
      }else if(plot_init$PlotID > 43000 &  plot_init$PlotID < 44000){ #WAYCROSS
        par_group[index] = FR_npar_group + 4
      }else if(plot_init$PlotID > 50000){ 
        par_group[index] = FR_npar_group
      }
    }else if(fr_model == 1 & FR_separate_npar_groups == 2){
      par_group[index] = FR_npar_group + plotnum - 1
    }
    index = index + 1
  }
  
  #WSx1000 INIT
  index = index_guide[5]
  for(plotnum in 1:nplots){
    plot_init <- initdata[plotnum,]
    plot_index[index] = plotnum
    init_pars[index] = rnorm(1,init_pars[19],25)
    prior_parameter1[index]= 0.1 #prior_parameter1[19]
    prior_parameter2[index]= 400 #prior_parameter2[19]
    prior_dist[index] = 1 #prior_dist[19]
    fix_par[index] = 0
    if((plot_init$PlotID > 10000 &  plot_init$PlotID < 20000) | (plot_WSx1000 == FALSE)){  #Thinning
      par_group[index] = plot_WSx1000_pargroup
    }else if(plot_init$PlotID > 20000 &  plot_init$PlotID < 30000){ #RW18
      par_group[index] = plot_WSx1000_pargroup +1
    }else if(plot_init$PlotID > 30000 &  plot_init$PlotID < 40000){ #TEIR3
      par_group[index] = plot_WSx1000_pargroup + 2
    }else if(plot_init$PlotID > 40000 &  plot_init$PlotID < 41000){  #DUKE
      par_group[index] = plot_WSx1000_pargroup + 3
    }else if(plot_init$PlotID > 41000 &  plot_init$PlotID < 42000){ #NC2
      par_group[index] = plot_WSx1000_pargroup + 4
    }else if(plot_init$PlotID > 42000 &  plot_init$PlotID < 43000){ #SETRES
      par_group[index] = plot_WSx1000_pargroup + 5
    }else if(plot_init$PlotID > 43000 &  plot_init$PlotID < 44000){ #WAYCROSS
      par_group[index] = plot_WSx1000_pargroup + 6
    }else if(plot_init$PlotID > 7200 & plot_init$PlotID <= 72076 ){ 
      par_group[index] = plot_WSx1000_pargroup + 7
    }else if(plot_init$PlotID > 5200 & plot_init$PlotID <= 52467){ 
      par_group[index] = plot_WSx1000_pargroup + 8
    }
    if(!plot_WSx1000){
      fix_par[index] = 1
    }
    jump_pars[index] = init_pars[19]*jump_size_init
    index = index + 1
  }
  
  #Thinpower INIT
  index = index_guide[7]
  for(plotnum in 1:nplots){
    plot_index[index] = plotnum
    init_pars[index] = init_pars[20]
    prior_parameter1[index]= prior_parameter1[20]
    prior_parameter2[index]= prior_parameter2[20]
    prior_dist[index] =prior_dist[20]
    fix_par[index] = 0
    par_group[index] = plot_thinpower_pargroup
    if(!plot_thinpower){
      fix_par[index] = 1
    }
    jump_pars[index] = init_pars[20]*jump_size_init
    index = index + 1
  }
  #Mort_rate INIT
  index = index_guide[9]
  for(plotnum in 1:nplots){
    plot_index[index] = plotnum
    init_pars[index] = init_pars[40]
    prior_parameter1[index] = 0.0 
    prior_parameter2[index]= prior_parameter2[40]
    prior_dist[index] = prior_dist[40]
    fix_par[index] = 0
    par_group[index] = plot_mort_rate_pargroup
    if(!plot_mort_rate){
      fix_par[index] = 1
    }
    jump_pars[index] = init_pars[40]*jump_size_init
    index = index + 1
  }
  
  
  
  
  return(list(init_pars=init_pars,new_pars=new_pars,jump_pars=jump_pars,prior_parameter1=prior_parameter1,prior_parameter2=prior_parameter2,prior_dist=prior_dist,fix_par=fix_par,par_group=par_group))
}
