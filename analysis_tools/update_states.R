update_states <- function(mo,plotnum,output_dim,pars,site,nopars,nosite, thin_event, met){ # fix - need to alter to accept initial states as input 
  #---CONTROL INFORMATION----------------------------
  #rm(list = ls())
  
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

  # below matchings have been confirmed
  LAI  = output[4]
  WS  = output[5]
  WCR = output[6]
  WR = output[7]
  Stem_Count = output[8]
  ASW = output[14]
  
  DAPPERstates <- c(LAI = LAI, WS = WS, WCR = WCR, WR = WR, Stem_Count = Stem_Count, ASW =  ASW) # LAI, WS, WCR, WFR, Stem_Count, ASW, FR
  return(DAPPERstates)
  
}  

  