module DAPPER_plot_mod

contains


subroutine likelihood(plotnum &
   						,nopars         & 
                        ,nplots			&
                        ,initdata_dim &
                       	,control_plot_index &
                       	,index_guide  &
                       	,nosamples &
                       	,sample_index &
                    	,exclude_hardwoods &
                    	,fit_plot &
                    	,matched_FR_plot_index &
                    	,par_group &
                    	,control_pars &
                    	,npar_groups & 
                    	,nstreams &
                    	,nmonths &
                    	,years  &
                    	,months &
                    	,mo_start_end &
						,new_pars &
 						,met &
                       	,initdata &
                       	,obs &
                       	,latent &
                       	,thin_event &
                       	,pred_new &
                       	,prob_new &
                       	,obs_uncert &
                       	,obs_gap &
                       	,fix_par &
                       	,use_fol_state)

	use R3PG_MODEL_MOD, only: R3PG_MODEL
 	use prob_functions_mod

 	implicit none

   	integer, intent(in) :: plotnum &
   						,nopars         &
                        ,nplots			&
                        ,initdata_dim &
                       	,control_plot_index(nplots) &
                       	,index_guide(10)  &
                       	,nosamples &
                       	,sample_index(nosamples) &
                    	,exclude_hardwoods(nplots) &
                    	,fit_plot(nplots) &
                    	,matched_FR_plot_index(nplots) &
                    	,par_group(nopars) &
                    	,control_pars(11) &
                    	,npar_groups & 
                    	,nstreams &
                    	,nmonths &
                    	,years(nmonths) &
                    	,months(nmonths) &
                    	,mo_start_end(nplots,2) &
                    	,obs_gap(nstreams,nplots,nmonths) &
                    	,use_fol_state(nplots)
                		

                        
 	double precision, intent(in)::met(nplots,6,nmonths) &
                       			,initdata(nplots,initdata_dim) &
                       			,obs(nstreams,nplots,nmonths) &
                       			,latent(nstreams,nplots,nmonths) &
                       			,thin_event(nplots,nmonths) &
                       			,obs_uncert(nstreams,nplots,nmonths) &
                       			,fix_par(nopars)

   			
   	double precision, intent(inout)::prob_new(nplots),pred_new(nstreams,nplots,nmonths), new_pars(nopars)	
 	integer :: nomet = 6
 	integer :: nosite = 25
 	integer :: output_dim = 64
 	integer :: npars_used_by_fortran = 48
    double precision :: output(64), modeled_old(nstreams)
     	
    double precision ::pars(48),LL(nstreams),modeled(nstreams),SD1(nstreams),SD2(nstreams)       			

   	double precision :: PlotID,SiteID,LAT_WGS84,Planting_year,PlantMonth, &
   		PlantDensityHa,Initial_ASW,ASW_min,ASW_max,SoilClass,SI,FR,Initial_WF,Initial_WS, &
   		Initial_WR,DroughtLevel, DroughtStart,FertFlag,CO2flag,CO2elev,ControlPlotID, &
    	FR1_new,FR2_new,FR3_new,Initial_WF_H,Initial_WS_H,Initial_WR_H,Initial_WCR,IrrFlag
    	
   	double precision :: new_FR,control_FR, prob_FR, prob_pFS,prob_plot_params   
   	double precision :: site(27), annual_fol_prod, annual_root_prod

   	integer::control_plotnum, index, index2,i, num,y_index, start_index,first_full_year
   	
   	integer :: PlantedYear, PlantedMonth,InitialYear, InitialMonth,EndMonth, StartAge, EndYear
	double precision :: WFi, WRi,WSi,WCRi, StemNum,measurement_year,measurement_month,WFi_H 
	double precision :: LAI_i, WRi_H,WSi_H,IrrigRate,ThroughFall
    double precision ::  new_pFS2,	new_pFS20, tmp, MeanTemp, MeanPrecip,Foliage_SD
	integer :: mo,fr_model,high_freq_obs,use_dk_pars,data_stream,state_space, time_since_last_measure(5)


	!INITIALIZE ARRAYS
	
	fr_model = control_pars(4)
 	high_freq_obs = control_pars(5)
 	use_dk_pars = control_pars(7)
 	state_space = control_pars(11)
 	
   	pars = new_pars(1:npars_used_by_fortran)
 
	PlotID = initdata(plotnum,1)
   	SiteID = initdata(plotnum,2)
   	LAT_WGS84=initdata(plotnum,3)
   	PlantDensityHa = initdata(plotnum,6)
   	Initial_ASW = initdata(plotnum,7) 
   	ASW_min = initdata(plotnum,8)
   	ASW_max=initdata(plotnum,9)
  	SoilClass = initdata(plotnum,10)
  	SI = initdata(plotnum,11)
  	FR = initdata(plotnum,12)
  	Initial_WF = initdata(plotnum,13)
  	Initial_WS = initdata(plotnum,14)
  	Initial_WR = initdata(plotnum,15)
  	DroughtLevel = initdata(plotnum,16)  
  	DroughtStart = initdata(plotnum,17) 
  	FertFlag=initdata(plotnum,18)
  	CO2flag = initdata(plotnum,19)
  	CO2elev = initdata(plotnum,20)
  	ControlPlotID = initdata(plotnum,21)
  	Initial_WF_H = initdata(plotnum,22)
  	Initial_WS_H = initdata(plotnum,23)
  	Initial_WR_H = initdata(plotnum,24)
  	MeanTemp = initdata(plotnum,26)
    MeanPrecip = initdata(plotnum,27)
  	InitialYear = initdata(plotnum,29)
  	InitialMonth = initdata(plotnum,30)
  	StartAge = initdata(plotnum,31) 
  	Initial_WCR = initdata(plotnum,32)
  	IrrFlag = initdata(plotnum,33)

	  	  	  		
  	FR1_new = new_pars(49)
  	FR2_new = new_pars(50)
  	FR3_new = new_pars(51)
  	
	SD1(1) = new_pars(52)
	SD1(2) = new_pars(53)
	SD1(3) = new_pars(54)
	SD1(4) = new_pars(55)
	SD1(5) = new_pars(56)
	SD1(6) = SD1(1)
	SD1(7) = new_pars(57)
	SD1(8) = -99
	SD1(9) = -99
	SD1(10) = -99
	SD1(11)  = -99
 	SD1(12) = new_pars(58)
	SD1(13) = -99
	SD1(14) = new_pars(59)
	SD1(15) = new_pars(60)
	SD1(16) = SD1(15)
	SD1(17) = new_pars(61)
	SD1(18) = new_pars(62)	
	Foliage_SD = new_pars(77)
		
	SD2(1) = new_pars(63)
	SD2(2) = new_pars(64)
	SD2(3) = new_pars(65)
	SD2(4) = new_pars(66)
	SD2(5) = new_pars(67)
	SD2(6) = SD2(1)
	SD2(7) = new_pars(68)
	SD2(8) = -99
	SD2(9) = -99
	SD2(10) = -99
	SD2(11)  = -99
 	SD2(12) = new_pars(69)
	SD2(13) = -99
	SD2(14) = new_pars(70)
	SD2(15) = new_pars(71)
	SD2(16) = SD2(15)
	SD2(17) = new_pars(72)
	SD2(18) = new_pars(73)	
	
	if(use_fol_state(plotnum) == 1) then
    	SD1(1) = Foliage_SD
	else
    	SD1(1) = new_pars(52)
    endif
    
	!-------------------------------------------------------------------------------------  
    !SET FR FOR EXPERIMENTAL TREATMENTS BASED ON THE CONTROL PLOT FR
    index = index_guide(3)+plotnum-1
    if(DroughtLevel<1 .AND. FertFlag==0) then
    	index2 = index_guide(3) + control_plot_index(plotnum) -1
       	new_pars(index) = new_pars(index2)
    end if
    if(PlotID >= 40000 .AND. FertFlag == 0) then  !Duke FACE plots
    	index2 = index_guide(3) + control_plot_index(plotnum) -1
       	new_pars(index) = new_pars(index2)
    end if
    if(IrrFlag == 1 .AND. FertFlag == 0) then
    	index2 = index_guide(3) + control_plot_index(plotnum) -1
        new_pars(index) = new_pars(index2)
    endif
    if(IrrFlag == 1 .AND. FertFlag == 1) then
    	index2 = index_guide(3) + control_plot_index(plotnum)-1
        new_pars(index) = new_pars(index2)
    endif
    if(DroughtLevel<1 .AND. FertFlag==1) then
    	index2 = index_guide(3) +control_plot_index(plotnum) -1
        new_pars(index) = new_pars(index2)
	end if
    if(FR == 1) then
    	new_pars(index) = 1.0
    end if
    !match FR for plots that are validation plot
    if(fit_plot(plotnum) == 0) then 
    	index2 = index_guide(3) + matched_FR_plot_index(plotnum) -1
     	new_pars(index) = new_pars(index2)
    endif
    
	!NOTE: INDEX GUIDE IS THE VECTOR THAT IDS WHERE DIFFERENT PARAMETERS ARE IN THE LONG PARAMETER VECTOR
  	new_FR = new_pars(index_guide(3)+plotnum-1)
  	control_plotnum = control_plot_index(plotnum)
  	control_FR = new_pars(index_guide(3) + control_plotnum - 1)
  	  	
  	! THIS IS A TEMPORARY SOLUTION SO THAT IRRIGATION IS TO THE SOIL RATHER THAN THE TOP OF THE CANOPY
  	IrrigRate = 0.0
  	if(IrrFlag == 1) then
		IrrigRate = 650.0/9.0
    endif
	!-------------------------------------------------------------------------------------  

	!-------------------------------------------------------------------------------------  
	!---PROBABILITY OF FR 
  	!---THIS IS WHERE THE N FERTILIZATION TREATMENT SIMULATED BY FORCING THE ESTIMATED FR TO BE GREATER THAN THE CONTROL
    prob_FR = 0.0
    if(new_FR < 0.0 .OR. new_FR > 1.0) then
        prob_FR = -1e30
    endif
    if(FertFlag == 1) then
    	if(new_FR<control_FR) then
    	      prob_FR=-1e30
    	endif
    endif	
    if(fit_plot(control_plotnum) == 0) then 
  		prob_FR = log(uniform_pdf(new_FR,0D+00,1D+00))  
  	endif
   	!-------------------------------------------------------------------------------------  
  	
  	
  	!----- PLOT LEVEL PARAMETERS THAT ARE FROM A REGIONAL DISTRIBUTION -------------------    
    prob_plot_params = 0.0
    !WsX1000 Variation
    index = index_guide(5) + plotnum -1
    if(fix_par(index) == 1) then
    	new_pars(index) = new_pars(19)
    else
    	prob_plot_params = prob_plot_params + &
        	log(normal_pdf(new_pars(index), new_pars(19), new_pars(74)))
	endif
    !ThinPower Variation
    index = index_guide(7) + plotnum -1
    if(fix_par(index) == 1) then
     	new_pars(index) = new_pars(20)
    else
     	prob_plot_params = prob_plot_params + &
        	log(normal_pdf(new_pars(index), new_pars(20), new_pars(75)))
    endif
    !MortRate Variation
    index = index_guide(9) + plotnum -1
    if(fix_par(index) == 1) then
    	new_pars(index) = new_pars(40)
    else
    	prob_plot_params = prob_plot_params + &
      		log(normal_pdf(new_pars(index), new_pars(40), new_pars(76)))
      	prob_plot_params = prob_plot_params + &
      		log(normal_pdf(new_pars(index), new_pars(40), new_pars(76)))
	endif 
	
	!Wsx1000 (plot level variability in parameter)
   	pars(19) = new_pars(index_guide(5)+plotnum - 1)
   	!thinpower (plot level variability in parameter)
   	pars(20) = new_pars(index_guide(7)+plotnum - 1)
   	!Mort_rate (plot level variability in parameter)
   	pars(40) = new_pars(index_guide(9)+plotnum - 1)
   	!-------------------------------------------------------------------------------------  
   	
  
  	!---LIKELIHOOD CALCULATIONS (PLOT LEVEL) -------------
    LL(:) = 0.0
	modeled(:) = 0.0
	
	
   	!---INITIAL CONDITIONS -------------
   	! --A BIT COMPILICATED
	site(1) =  PlantedYear !PlantedYear
    site(2) =  PlantedMonth !"PlantedMonth"
    site(3) =  months(mo_start_end(plotnum,1)) !"InitialYear"
    site(4) =  years(mo_start_end(plotnum,1))
    site(5) =  StartAge
    site(7) =  latent(4,plotnum,mo_start_end(plotnum,1))-Initial_WR_H !"WRi" 
    site(8) =  latent(2,plotnum,mo_start_end(plotnum,1))!"WSi" 
    site(9) =  latent(5,plotnum,mo_start_end(plotnum,1))!"StemNoi"
    site(10) = ASW_max
    site(11) = LAT_WGS84 !"Lat"
    site(12) = new_FR !"FR" 
    site(13) = SoilClass !"SoilClass"
    site(14) = ASW_max !"MaxASW"
    site(15) = ASW_min !"MinASW"
    site(16) = 1
    site(17) = Initial_WF_H
    site(18) = latent(7,plotnum,mo_start_end(plotnum,1))
    site(19) = Initial_WR_H
    site(20) = latent(3,plotnum,mo_start_end(plotnum,1)) !WCR
    site(21) = IrrigRate
    site(22) = DroughtLevel !ALLOWS FOR THROUGHFALL EXCLUSION
    site(23) = 0 !THIS IS A TEMP SOLUTION FOR SEPARATING THE SELF THINNINGAT DUKE  
    if(PlotID > 40000 .AND. PlotID < 41000 .AND. use_dk_pars == 1) then
    	site(23) = 1
    endif
	site(24) = latent(7,plotnum,mo_start_end(plotnum,1))*0.2
	site(25) = 0.0  !Bud C
	site(27) = latent(6,plotnum,mo_start_end(plotnum,1))
	
	if(use_fol_state(plotnum) == 1) then
		site(6) = latent(1,plotnum,mo_start_end(plotnum,1))
		site(26) = -99
	else
		site(6) = Initial_WF
		site(26) = latent(1,plotnum,mo_start_end(plotnum,1))
	endif
	
	
	!---START LOOPING THOUGHT MONTHS -------------
					
    do mo = (mo_start_end(plotnum,1)+1),mo_start_end(plotnum,2)
    	
    	site(3) = years(mo)
        site(4)  = months(mo)
        
       	call R3PG_MODEL(output_dim,&
                met(plotnum,:,mo), &
                pars, &
                site, &
                thin_event(plotnum,mo),&
                nopars,&
                6, &
				nosite,&
				output_dim,&
				1,&
				1,&
				1,&
				exclude_hardwoods(plotnum),&
				OUTPUT)
				
		 if(use_fol_state(plotnum) .EQ. 1) then	
			modeled(1) = output(22) !output(4) !Pine LAI
         else
         	modeled(1) = output(4) !output(4) !Pine LAI
         endif
         
         !--- PROCESS MODEL OUTPUT       
         modeled(2) = output(5) !WS
         modeled(3) = output(6) !+ output(43) !WCR
         modeled(4) = output(7) !+ output(58) !WR
         modeled(5) = output(8) !Stem Density
         modeled(6) = output(9) !LAI_H        
         modeled(7) = output(10) !Hardwood Stem   
         modeled(8) = output(11) !Hardwood Coarse roots
         modeled(9) = output(12) !Hardwood Fine roots 
         modeled(10) = output(13) !Hardwood Stem density 
         modeled(11) = output(14) ! ASW        
         modeled(12) = output(15) ! GEP        
         modeled(13) = output(16) ! NEE
         modeled(14) = output(17) ! ET  
         modeled(15) = output(18) ! Ctrans Pine           
         modeled(16) = output(19) ! Ctrans Hardwood     
         
        !ANNUAL PRODUCTION AND TURNOVER FLUXES
		if(output(2) > 1 .AND. mo == 1) then
    		first_full_year = 0
        elseif(output(2) == 12) then
        	first_full_year = 1
        endif
        if(output(2)==1 .OR. mo == 1) then
            modeled(17) = output(20) !root production
    		modeled(18) =  output(21)  !Fol production

    	else
    	    modeled(17) = modeled(17) + output(20) !root production
			modeled(18) = modeled(18) + output(21)  !Fol production

    	endif 
    	
    	! ASSIGN THE PREDATIONS THAT ARE USED IN THE LATENT STATE CALCULATIONS
                           
         pred_new(:,plotnum,mo) =  modeled(:)
         
         
        !-------  COMPARE MODEL PREDICTIONS TO LATENT STATES AND COMPARED LATENT STATES TO OBSERVATIONS
		!STATE-SPACE DATA STREAM
		do data_stream=1,7
			if(data_stream == 4) then 
				if(latent(data_stream,plotnum,mo) .NE. -99) then
					!Special case for fine roots because the pine and hardwood are added together
					pred_new(data_stream,plotnum,mo) = modeled(4) + modeled(9)
			 		LL(data_stream)  = LL(data_stream)  + &
						log(normal_pdf(latent(data_stream,plotnum,mo),& !added together
				  		pred_new(data_stream,plotnum,mo),(SD1(data_stream)+pred_new(data_stream,plotnum,mo)*SD2(data_stream)) &
                                                *obs_gap(data_stream,plotnum,mo)))
                    if(obs(data_stream,plotnum,mo) .NE. -99) then
						LL(data_stream)  = LL(data_stream)  + &
							log(normal_pdf(latent(data_stream,plotnum,mo),&  !added together
							obs(data_stream,plotnum,mo),obs_uncert(data_stream,plotnum,mo)))
					endif	
				endif				
			else
					if(latent(data_stream,plotnum,mo) .NE. -99) then
						LL(data_stream)  = LL(data_stream)  + log(normal_pdf(latent(data_stream,plotnum,mo),&
				  		pred_new(data_stream,plotnum,mo),(SD1(data_stream)+modeled(data_stream)*SD2(data_stream)) &
                                                *obs_gap(data_stream,plotnum,mo)))
                    	if(obs(data_stream,plotnum,mo) .NE. -99) then                           
						LL(data_stream)  = LL(data_stream)  + log(normal_pdf(latent(data_stream,plotnum,mo),&
							obs(data_stream,plotnum,mo),obs_uncert(data_stream,plotnum,mo)))
						endif		
					endif
	
    			endif
    	enddo
    	
    	!FLUXES
		do data_stream = 12,16			
			if(data_stream .NE. 13 .AND. obs(data_stream,plotnum,mo) .NE. -99) then
				 LL(data_stream)  = LL(data_stream)  + log(normal_pdf(latent(data_stream,plotnum,mo),&
				  	pred_new(data_stream,plotnum,mo),SD1(data_stream)+modeled(data_stream)*SD2(data_stream)))
				LL(data_stream)  = LL(data_stream)  + log(normal_pdf(latent(data_stream,plotnum,mo),&
					obs(data_stream,plotnum,mo),obs_uncert(data_stream,plotnum,mo)))	
    		endif
		end do
				    
		!ANNUAL PRODUCTION AND TURNOVER FLUXES
    	if(first_full_year == 1 .AND. months(mo)==12) then 
    		do data_stream = 17,18
        		if(obs(data_stream,plotnum,mo) .NE. -99) then
					LL(data_stream)  = LL(data_stream)  + log(normal_pdf(latent(data_stream,plotnum,mo),&
						modeled(data_stream),SD1(data_stream)+modeled(data_stream)*SD2(data_stream)))
					LL(data_stream)  = LL(data_stream)  + log(normal_pdf(latent(data_stream,plotnum,mo),&
						obs(data_stream,plotnum,mo),obs_uncert(data_stream,plotnum,mo)))
    			endif
        	end do
        endif

		!Updating inputs for the next month simulation
		!This is not internally updated when doing only a single month simulation
		! COMPLICATED BECAUSE THERE ARE 3 OPTIONS FOR SIMULATIONS
		
		site(5) =  output(3)+ (1.0/12.) !StartAge !StartAge

		if(state_space > 0) then
			if(latent(1,plotnum,mo) .NE. -99.0) then
				if(use_fol_state(plotnum) == 1) then
					site(6) = latent(1,plotnum,mo)
					site(26) = -99
				else
					site(6) = output(22)
					site(26) = latent(1,plotnum,mo)
				endif
        	else
        		if(use_fol_state(plotnum) == 1) then
					site(6) = modeled(1)
					site(26) = -99
				else
					site(6) = output(22)
					site(26) = modeled(1)
				endif
    
        	endif
        	if(latent(2,plotnum,mo) .NE. -99.0) then
        		site(8) = latent(2,plotnum,mo) !WSi
        	else
        		site(8) = modeled(2)
        	endif    
        	
         	if(latent(3,plotnum,mo) .NE. -99.0) then
        		site(20) = latent(3,plotnum,mo) !WCR
        	else
        		site(20) = modeled(3)
        	endif   
        	
        	if(latent(4,plotnum,mo) .NE. -99.0) then
        		site(7) = latent(4,plotnum,mo)* (modeled(4)/(modeled(9)+modeled(4))) 
        		site(19) =latent(4,plotnum,mo)* (modeled(9)/(modeled(9)+modeled(4)))
        	else
        		site(7) = modeled(4)
        		site(19) = modeled(9)
        	endif        
        	
        	if(latent(5,plotnum,mo) .NE. -99.0) then
        		site(9) = latent(5,plotnum,mo) !StemNo
        	else
        		site(9) = modeled(5)
        	endif
        	
        	if(latent(6,plotnum,mo) .NE. -99.0) then
        		site(27) = latent(6,plotnum,mo) ! Hardwood LAI
        		site(27) = modeled(6)
        		site(25) =  output(26)
        		!site(25) =  output(26)*(latent(6,plotnum,mo)/pred_new(6,plotnum,mo))
        	else
        		site(27) = modeled(6)
        		site(25) =  output(26)
        	endif
        	
        	if(latent(7,plotnum,mo) .NE. -99.0) then
        		site(18) = latent(7,plotnum,mo) ! Hardwood WS
        	else
        		site(18) = modeled(7)
        	endif

		 
        else
        	if(use_fol_state(plotnum) == 1) then
					site(6) = modeled(1)
					site(26) = -99
				else
					site(6) = output(22)
					site(26) = modeled(1)
			endif
        	site(26) = modeled(1)  !LAI Pine
        	site(8) =  modeled(2) !WS
        	site(20) = modeled(3) !WCR
			site(7) = modeled(4) !WRi
			site(9) = modeled(5) !StemNo
			site(27) = modeled(6)
			site(18) = modeled(7)
			site(19) = modeled(9) !WR_H
			site(25) =  output(26)
		endif

		!Other state variables

        site(10) = modeled(11) !ASW
        site(17) = output(23) !WF_H
  		site(24) =  output(11)!WCR_h
   		
   
	end do !END LOOPING THROUGH OUTPUT MONTHS AND SEARCHING FOR CORRESPONDING OBSERVATIONS
	
	! SUM UP OVER TIME AND STATE-STREAMS
	prob_new(plotnum)= LL(1) + LL(2) + LL(3) + LL(4) + LL(5) + LL(6) + LL(7) + &
					   LL(12) + LL(14) + LL(15) + LL(16) + &
					   LL(17)+ LL(18) + prob_FR + prob_plot_params

end subroutine likelihood

end module DAPPER_plot_mod




