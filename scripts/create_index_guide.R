create_index_guide <- function(npars,nplots){
  #-----INDEXING FOR THE PARAMETER VECTOR--------------------------
  pars_start = 1
  pars_end = npars
  FR_start = pars_end + 1
  FR_end = FR_start + nplots -1
  wSx1000_start = FR_end + 1 
  wSx1000_end = wSx1000_start + nplots -1 
  thinPower_start =  wSx1000_end + 1 
  thinPower_end =  thinPower_start + nplots -1 
  mort_rate_start =  thinPower_end + 1 
  mort_rate_end =  mort_rate_start + nplots -1 
    
  index_guide = c(pars_start,pars_end,FR_start,FR_end,wSx1000_start,wSx1000_end,thinPower_start,thinPower_end,
                  mort_rate_start,mort_rate_end)
  return(index_guide)
}