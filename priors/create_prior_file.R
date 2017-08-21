
parnames = c(
  'pFS2',
  'pFS20',
  'StemConst',
  'StemPower',
  'pRx',
  'pRn',
  'SLA0',
  'SLA1',
  'tSLA',
  'k',
  'fullCanAge',
  'MaxIntcptn',
  'LAImaxIntcptn',
  'alpha',
  'MaxCond',
  'LAIgcx',
  'CoeffCond',
   'BLcond (tmp DK thinPower)',
  'wSx100',
  'thinPower',
  'mF',
  'mR',
  'mS',
  'Rttover',
  'm0 (TMP DK CR)',
  'fN0',
  'TMin',
  'Topt',
  'Tmax',
  'kF',
  'MaxAge',
   'nAge',
  'rAge',
  'y',
  'FCalpha700_h',
  'FCalpha700',
  'fCg700',
  'SWconst1',
  'SWpower1',
  'mort_rate',
  'Alpha_h', 
  'pFS_h',
  'pR_h',
  'mort_rate_h',
  'SLA_h',
  'pCRS',
  'FCpFS700',
  'Leaftover (Temp DK wSx100)',
  'FR Par 1',
  'FR Par 2',
  'FR Par 3',
  'LAI_SD',
  'WS_SD',
  'WCR_SD',
  'WR_SD',
  'Stem Density_SD',
  'Hardwood_stem_SD',
  'GEP_SD',
  'ET_SD',
  'Ctrans_SD',
  'Root_production_SD',
  'Foliage_production_SD',
  'LAI_SD_2',
  'WS_SD_2',
  'WCR_SD_2',
  'WR_SD_2',
  'Stem Density_SD_2',
  'Hardwood_stem_SD_2',
  'GEP_SD_2',
  'ET_SD_2',
  'Ctrans_SD_2',
  'Root_production_SD_2',
  'Foliage_production_SD_2',
  'WsX1000_SD',
  'ThinPower_SD',
  'Mort_rate_SD',
  'Foliage_SD')


priormatrix = matrix(NA,77,6)


priormatrix[1,]  = c(.40,0.08,2.0,1,0,1) #pFS2 
priormatrix[2,]  = c(.30,0.1,1.0,1,0,2) # pFS20 
priormatrix[3,]  = c(.021754,.021754,0.005,2,0,3) #StemConst 
priormatrix[4,]  = c(2.77,2.77,0.2,2,0,4)  #StemPower 
priormatrix[5,]  = c(0.5,0.0,2.0,1,0,5) #pRx 
priormatrix[6,]  = c(0.07,0.01,1.0,1,1,100) #pRn 
priormatrix[7,]  = c(5.4287,5.529,0.44,2,0,6)#SLA0 
priormatrix[8,]  = c(3.5754,3.875,0.11,2,0,7) #SLA1 
priormatrix[9,]  = c(5.9705,5.9705,2.15,2,0,8) #tSLA 
priormatrix[10,] = c(.56,0.55,0.70,1,1,100) #k 
priormatrix[11,] = c(3,2,5,1,1,100) #fullCanAge  
priormatrix[12,] = c(.2,0,1,1,0,100) #MaxIntcptn 
priormatrix[13,] = c(2.1,2,6,1,1,100) #LAImaxIntcptn 
priormatrix[14,] = c(0.04,0.02,.08,1,0,9) #alpha 
priormatrix[15,] = c(0.02,0.005,0.03,1,0,10) #MaxCond 
priormatrix[16,] = c(4.5,2.0,5.0,1,0,11)  #LAIgcx
priormatrix[17,] = c(0.037,0.0408,0.0028,2,0,12)  #CoeffCond 
priormatrix[18,] = c(1.0,0.0,10,1,1,13)  # HardPineCondRatio
priormatrix[19,] = c(250,10,400,1,0,14) #wSx1000 c(170,100,450,1,0) 3PGPAPERS
priormatrix[20,] = c(1.5,1.5,0.1,2,0,15)  #thinPower
priormatrix[21,] = c(0.1,0.05,0.3,1,0,100) #mF
priormatrix[22,] = c(0.0,0.0,1.0,1,1,100)  #AVIALABLE mR
priormatrix[23,] = c(0.5,0.2,1,1,0,16) #mS
priormatrix[24,] = c(.04,0.0167,0.2,1,0,17)  #Rttover  c(.04,0.0167,0.0417,1,0,17) 
priormatrix[25,] = c(0.20,0.15,0.35,1,0,18) #pCR for DUKE CAN BE m0  VAGUE 3PG PAPERS 
priormatrix[26,] = c(1.5,1.5,0.1,2,0,15) # AVIALABLE mS for DUKE CAN BE fN0 VAGUE 3PG PAPERS
priormatrix[27,] = c(-2.44,4,2,2,0,19)  #Tmin 
priormatrix[28,] = c(23.72,25,2,2,0,20) #Topt 
priormatrix[29,] = c(40.51,38,2,2,0,21) #Tmax 
priormatrix[30,] = c(0.17,0.1780282,0.0162,2,0,22)  #kF 
priormatrix[31,] = c(199,16,200,1,0,23)  #MaxAge  
priormatrix[32,] = c(3.54,1,4,1,0,24) #nAge
priormatrix[33,] = c(0.5,0.01,3,1,0,25) #rAge 
priormatrix[34,] = c(.48,.3,0.65,1,0,26) #y 
priormatrix[35,] = c(1.50,1.00,2.5,1,0,27)#fCalpha700_h  
priormatrix[36,] = c(1.3,1.00,1.8,1,0,28)  #fCalpha  
priormatrix[37,] = c(1,0.2,1.8,1,1,100) #fCg700 
priormatrix[38,] = c(0.65,0.01,1.8,1,0,29) #SWconst 
priormatrix[39,] = c(8,1.0,13,1,0,30)##SWpower
priormatrix[40,] = c(9.8e-4,2e-04,0.001666667,1,0,31) #mort_rate c(9.8e-4,2e-04,0.001666667,1,0,31)
priormatrix[41,] = c(.02,0.005,0.08,1,0,32)  #alpha_h 
priormatrix[42,] = c(1.87,0.2,10.0,1,0,33) #pFS_h 
priormatrix[43,] = c(0.16,0.0,2,1,1,34)  #pR_h 
priormatrix[44,] = c(9.8e-4,0.0000,0.0025,1,1,100)#mort_rate_h 
priormatrix[45,] = c(16.2,16,3.8,2,0,35)  #SLA_h 
priormatrix[46,] = c(0.3,0.15,0.35,1,0,36) #pCRS 
priormatrix[47,] = c(0.84,0.50,1.5,1,0,37) #fCpFS700 
priormatrix[48,] = c(230,230,25,2,0,38) #Duke wSx1000 CAN ALSO BE Leaftover
priormatrix[49,] = c(0.0,-1.0,1.2,1,1,100) #FR Par 1 
priormatrix[50,] = c(0.02,0.0,1.0,1,0,39) #FR Par 2 
priormatrix[51,] = c(0.144,0.0,1.0,1,0,40) #FR Par 3 
priormatrix[52,] = c(10,0.001,100,1,0,41) #LAI_SD 
priormatrix[53,] = c(20,0.001,100,1,0,42) #WS_SD
priormatrix[54,] = c(2,0.001,100,1,0,43) #WCR_SD
priormatrix[55,] = c(5,0.001,100,1,0,44) #WR_SD
priormatrix[56,] = c(100,0.001,1000,1,0,45) #Stem Density_SD
priormatrix[57,] = c(1,0.001,100,1,0,52) #Hardwood stem SD
priormatrix[58,] = c(6,0.001,100,1,0,45) #GEP_SD
priormatrix[59,] = c(10,0.001,100,1,0,47) #ET_SD 
priormatrix[60,] = c(5,0.001,100,1,1,100) #Ctrans_SD
priormatrix[61,] = c(1,0.001,100,1,0,50) #Root production SD
priormatrix[62,] = c(1,0.001,100,1,0,49)#Foliage production SD
priormatrix[63,] = c(0.0,0.0,1000,1,1,100) #LAI_SD_2
priormatrix[64,] = c(0.0,0.0,1000,1,1,100)  #WS_SD_2
priormatrix[65,] = c(0.0,0.0,1000,1,1,100)#WCR_SD_2
priormatrix[66,] = c(0.0,0.0,1000,1,1,100)  #WR_SD_2
priormatrix[67,] = c(0.0,0.0,1000,1,1,100)  #Stem Density_SD_2
priormatrix[68,] = c(0.0,0.0,1000,1,1,100)  #Hardwood stem SD_2
priormatrix[69,] = c(0.1,0.0,1000,1,0,54)   #GEP_SD_2 
priormatrix[70,] = c(0.1,0.0,1000,1,0,55)  #ET_SD_2
priormatrix[71,] = c(0.1,0.0,1000,1,1,100)  #Ctrans_SD_2
priormatrix[72,] = c(0.0,0.0,1000,1,1,100)  #Root production SD_2
priormatrix[73,] = c(0.0,0.0,1000,1,1,100)  #Foliage production SD_2
priormatrix[74,] = c(1000,0.0,1000,1,1,100) #plot level uncertainity in WSx1000
priormatrix[75,] = c(0.1,0.0,2,1,1,100) #plot level uncertainity in ThinPower
priormatrix[76,] = c(0.01,0.0,2,1,1,100) #plot level uncertainity in mort_rate
priormatrix[77,] = c(0.2,0.0,100,1,1,100)  #Foliage SD


priors = cbind.data.frame(parnames,initial_value = priormatrix[,1],dist_par1 = priormatrix[,2], dist_par2 = priormatrix[,3],dist_type=priormatrix[,4],fit_par =priormatrix[,5],par_group =priormatrix[,5] ,row.names=NULL)
write.csv(priors,'/Users/quinn/Dropbox/Research/DAPPER/priors/default_priors.csv')


####

datafilename='/Users/quinn/Dropbox/Research/DAPER/chains/duke_state_space.1.2017-07-24.06.59.06.Rdata'

priormatrix = matrix(NA,77,6)
load(datafilename)

parnames = c(
  'pFS2',
  'pFS20',
  'StemConst',
  'StemPower',
  'pRx',
  'pRn',
  'SLA0',
  'SLA1',
  'tSLA',
  'k',
  'fullCanAge',
  'MaxIntcptn',
  'LAImaxIntcptn',
  'alpha',
  'MaxCond',
  'LAIgcx',
  'CoeffCond',
  'BLcond (tmp DK thinPower)',
  'wSx100',
  'thinPower',
  'mF',
  'mR',
  'mS',
  'Rttover',
  'm0 (TMP DK CR)',
  'fN0',
  'TMin',
  'Topt',
  'Tmax',
  'kF',
  'MaxAge',
  'nAge',
  'rAge',
  'y',
  'FCalpha700_h',
  'FCalpha700',
  'fCg700',
  'SWconst1',
  'SWpower1',
  'mort_rate',
  'Alpha_h', 
  'pFS_h',
  'pR_h',
  'mort_rate_h',
  'SLA_h',
  'pCRS',
  'FCpFS700',
  'Leaftover (Temp DK wSx100)',
  'FR Par 1',
  'FR Par 2',
  'FR Par 3',
  'LAI_SD',
  'WS_SD',
  'WCR_SD',
  'WR_SD',
  'Stem Density_SD',
  'Hardwood_stem_SD',
  'GEP_SD',
  'ET_SD',
  'Ctrans_SD',
  'Root_production_SD',
  'Foliage_production_SD',
  'LAI_SD_2',
  'WS_SD_2',
  'WCR_SD_2',
  'WR_SD_2',
  'Stem Density_SD_2',
  'Hardwood_stem_SD_2',
  'GEP_SD_2',
  'ET_SD_2',
  'Ctrans_SD_2',
  'Root_production_SD_2',
  'Foliage_production_SD_2',
  'WsX1000_SD',
  'ThinPower_SD',
  'Mort_rate_SD',
  'Foliage_SD')


for(p in 1:77){
  priormatrix[p,1] = mean(accepted_pars_thinned_burned[,p])
  priormatrix[p,2] = mean(accepted_pars_thinned_burned[,p])
  priormatrix[p,3] = sd(accepted_pars_thinned_burned[,p])
  priormatrix[p,4] = 2
  priormatrix[p,5] = 0
  if(priormatrix[p,3] == 0.0){
    priormatrix[p,5] = 1
  }
  priormatrix[p,6] = p
}

priors = cbind.data.frame(parnames,initial_value = priormatrix[,1],dist_par1 = priormatrix[,2], dist_par2 = priormatrix[,3],dist_type=priormatrix[,4],fit_par =priormatrix[,5],par_group =priormatrix[,5] ,row.names=NULL)

write.csv(priors,'/Users/quinn/Dropbox/Research/DAPER/priors/post_duke_posteriors.csv')

#########
datafilename='/Users/quinn/Dropbox/Research/DAPER/chains/revision_base_chain1.1.2017-04-27.11.25.11.final.Rdata'

priormatrix = matrix(NA,80,6)
load(datafilename)

parnames = c('pFS2','pFS20','StemConst','StemPower','pRx','pRn','SLA0','SLA1','tSLA','k',
             'fullCanAge','MaxIntcptn','LAImaxIntcptn','alpha','MaxCond','LAIgcx','CoeffCond',
             'BLcond (tmp DK thinPower)','wSx100','thinPower','mF','mR','mS','Rttover','m0 (TMP DK CR)','fN0','TMin','Topt','Tmax','kF','MaxAge',
             'nAge','rAge','y','FCalpha700_h','FCalpha700', 'fCg700','SWconst1','SWpower1','mort_rate','Alpha_h', 
             'pFS_h','pR_h','mort_rate_h','SLA_h','pCRS','FCpFS700','Leaftover (Temp DK wSx100)','FR Par 1','FR Par 2','FR Par 3',
             'WF_SD','WS_SD','WR_SD','WCR_SD','Stem Density_SD','GEP_SD','ET_SD','Ctrans_SD','LAI_SD','Foliage_production_SD',
             'Root_production_SD','Hardwood_foliage_SD','Hardwood_stem_SD',
             'WF_SD_2','WS_SD_2','WR_SD_2','WCR_SD_2','Stem Density_SD_2','GEP_SD_2','ET_SD_2',
             'Ctrans_SD_2','LAI_SD_2','Foliage_production_SD_2','Root_production_SD_2','Hardwood_foliage_SD_2','Hardwood_stem_SD_2',
             'WsX1000_SD','ThinPower_SD','Mort_rate_SD')

for(p in 1:80){
  priormatrix[p,1] = mean(accepted_pars_thinned_burned[,p])
  priormatrix[p,2] = mean(accepted_pars_thinned_burned[,p])
  priormatrix[p,3] = sd(accepted_pars_thinned_burned[,p])
  priormatrix[p,4] = 1
  priormatrix[p,5] = 0
  if(priormatrix[p,3] == 0.0){
    priormatrix[p,5] = 1
  }
  priormatrix[p,6] = p
}

priors = cbind.data.frame(parnames,initial_value = priormatrix[,1],dist_par1 = priormatrix[,2], dist_par2 = priormatrix[,3],dist_type=priormatrix[,4],fit_par =priormatrix[,5],par_group =priormatrix[,5] ,row.names=NULL)

write.csv(priors,'BG_posteriors.csv')
