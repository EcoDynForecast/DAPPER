subroutine DAPPER_MCMC( &
						nopars         &
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
                    	,data_uncertainity_npar_group &
                    	,nstreams &
                    	,nmonths &
                    	,years &
                    	,months &
                    	,mo_start_end &
                    	,met &
                       	,initdata &
                       	,obs &
                       	,thin_event &
                       	,init_pars &
                       	,prior_parameter1 &
         				,prior_parameter2 &
         				,prior_dist &
         				,fix_par &
         				,obs_uncert &
         				,latent &
						,jump_pars &
						,pnow &
						,accepted_pars_thinned_burned &
						,like_chain &
						,current_like &
						,init_obs &
						,init_uncert &
						,tracked_plot &
						,tracked_plotnum &
						,obs_gap &
						,obs_gap_next &
						,use_fol_state &
						)
                       	
	use DAPPER_plot_mod, only: likelihood
  	use prob_functions_mod, only: idum,randn,normal_sample,normal_pdf,uniform_pdf,std
  	use iso_fortran_env, only: int64
  	

  	! subroutine specificially deals with the calling of the fortran code model by
  	! R

  	implicit none
  	! declare input variables
  
    integer, intent(in) :: nopars         &
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
                    	,data_uncertainity_npar_group &
                    	,nstreams &
                    	,nmonths &
                    	,years(nmonths) &
                    	,months(nmonths) &
                    	,mo_start_end(nplots,2) &
                    	,tracked_plotnum &
                    	,obs_gap(nstreams,nplots,nmonths) &
                    	,obs_gap_next(nstreams,nplots,nmonths) &
                    	,use_fol_state(nplots)
                    	
                    	
                    	
	double precision, intent(in) ::met(nplots,nmonths) &
                       			,initdata(nplots,initdata_dim) &
                       			,obs(nstreams,nplots,nmonths) &
                       			,thin_event(nplots,nmonths) &
                       			,init_pars(nopars) &
                       			,prior_parameter1(nopars) &
         						,prior_parameter2(nopars) &
         						,prior_dist(nopars) &
         						,fix_par(nopars) &
         						,obs_uncert(nstreams,nplots,nmonths) &
         						,init_obs(nstreams,nplots) &
         						,init_uncert(nstreams,nplots)
         						
         						

	double precision, intent(inout) :: latent(nstreams,nplots,nmonths), &
									 	jump_pars(nopars),pnow,accepted_pars_thinned_burned(nosamples,nopars), &
									 	like_chain(nosamples), current_like,tracked_plot(nosamples,nstreams,nmonths)
	
	double precision :: prob_new(nplots), pred(nstreams,nplots,nmonths), &
						pred_new(nstreams,nplots,nmonths),pred_state_increment(nstreams,nplots,nmonths), &
						prob_plot_params(nplots),prob_plot_params_all, latent_old(nstreams,nplots,nmonths)
  	
   	integer :: plotnum, index2, num,mo,data_stream
   	integer :: max_fitted_pars = 77
   	integer :: nooutputs = 64
   	integer :: nadapt = 100
   
   	integer :: p,pass,index,iter,count1,niter,start_adapt,cost_type,fr_model,pgroup,use_dk_pars
   	integer :: use_age_edc,use_fr_edc,use_sm_edc,high_freq_obs
   	integer :: global_accept(npar_groups),local_accept(npar_groups)
   	double precision :: PlotID,DroughtLevel,FertFlag,FR,prior,pnew,new_pars(nopars),IrrFlag, &
   		local_accept_rate, r, z, like_new,FR1_new, FR2_new, FR3_new,SI,MeanTemp,MeanPrecip
   	double precision :: adaptfac = 1.25	
   	double precision :: upper_accept_bound = 0.44
    double precision :: lower_accept_bound = 0.23
    double precision :: std_fac = 2
    
   	integer :: dt(8), sample,time1, time2,time3,index1
   	double precision :: current_pars(nopars),t,SD1(nstreams),SD2(nstreams),tmp_sd
   	integer :: n, fp, adapt_index,i,nOut,saved_samples,thin_interval,thin_index,nPRINT,state_space
   	double precision  :: npar, rn,fadapt,adapt_adjust_rate,obs_uncertainity,rr(nplots),Foliage_SD
   	CHARACTER(LEN=30) :: Format 
	double precision  :: V,Vv, diff
	integer :: values(1:8), k
	integer, dimension(:), allocatable :: seed

	call date_and_time(values=values)

	call random_seed(size=k)
	allocate(seed(1:k))
	seed(:) = values(8)
	call random_seed(put=seed)
	call random_number(r)
   
    ! default values

 	!call random_seed(size = n)
  	!allocate(seed(n))
  	!call random_seed(get=seed)
  	!print *, seed

    !call system_clock(time1,time2,time3)
    !CALL RANDOM_SEED()
    !CALL RANDOM_NUMBER(r)
    !print *, r
 	! set seed value outside of the function, idum must be a negative number
 	idum=time1+time2+time3

 	niter = control_pars(1)
 	start_adapt = control_pars(2)
 	cost_type = control_pars(3)
 	fr_model = control_pars(4)
 	high_freq_obs = control_pars(5)
 	obs_uncertainity = control_pars(6)
 	use_dk_pars = control_pars(7)
	use_age_edc = control_pars(8)
    use_fr_edc = control_pars(9)
    use_sm_edc = control_pars(10)
    state_space = control_pars(11)

	print *, 'Start fitting'
	
  	local_accept(:) = 0.0
  	current_pars(:)=init_pars(:)
  	adapt_index = 1
  	nOut = nosamples
  	fadapt = 0.50
  	thin_index = 1
	nPRINT = 1000
 
   	sample = 1
   	global_accept(:) = 0

  	!---START LOOP THROUGH MCMC ITERATIONS
  	do iter = 1,(niter)
		do pgroup = 1,npar_groups
			!---PICK NEW PARAMETER VALUES
               do p = 1,nopars
                     pass = 0
                      count1 = 1
                      if(par_group(p)==pgroup .AND. fix_par(p) == 0 .AND. iter > 0) then
                       	do while(pass==0)
                             if(count1 > 1000000) then
                                 write(*,'(A,I5)') 'stuck on parameter ',p
                                 write(*,*) current_pars(p),jump_pars(p),prior_parameter1(p),prior_parameter2(p)
                                 stop
                              endif
                             if(count1 > 50000) then
                                 !write(*,'(A,I5)') 'stuck and adjusting parameter ',p
                                 !write(*,*) current_pars(p),jump_pars(p),prior_parameter1(p),prior_parameter2(p)
                                 jump_pars(p) = jump_pars(p)/adaptfac
                              endif

                         	new_pars(p)=normal_sample(current_pars(p),jump_pars(p)) 
                         	if(prior_dist(p)==1) then
                                 if(new_pars(p) >= prior_parameter1(p) .AND. new_pars(p) <= prior_parameter2(p)) then
                                 pass = 1
                                 end if
                          	else
                                 pass = 1
                         	end if
                                 count1 = count1 + 1
                      	end do
                    else
                         new_pars(p) = current_pars(p)
                    endif
            end do

			! FR set parameters to new values
        	FR1_new = new_pars(49)
        	FR2_new = new_pars(50)
        	FR3_new = new_pars(51)
        	! Deal with the different options for FR models (SI vs FR model is fr_model == 2)
          if(fr_model == 2) then 
                do plotnum = 1,nplots
                        MeanTemp = initdata(plotnum,26)
                        MeanPrecip = initdata(plotnum,27)
                        SI = initdata(plotnum,11) 
                         FertFlag=initdata(plotnum,18)
                        index = index_guide(3)+plotnum-1
                        if(FertFlag == 0) then
                                new_pars(index) = 1/(1+exp((FR1_new +FR2_new*MeanTemp-FR3_new*SI)))
                                if(new_pars(index) > 1.0) then
                                        new_pars(index) = 1.0
                                elseif(new_pars(index) < 0.0) then
                                        new_pars(index) = 0.0
                                endif
                        endif
                enddo
        endif	
        
        !1 = LAI Pine
		!2 = Stem Pine
		!3 = Coarse root pine
		!4 = Fine root pine
		!5 = Pine stem density
		!6 = LAI Hardwood
		!7 = Stem Hardwood
		!8 = Coarse root hardwood
		!9 = Fine root hardwood
		!10 = Hardwood stem density
		!11 = Available soil water
		!12 = GEP
		!13 = NEE  !NOT USED PLACE
		!14 = ET
		!15 = Ctrans Pine
		!16 = Ctrans Hardwood
		!17 = Fine root production
		!18 = Foliage production
		!19 = Pine foliage
		!20 = Hardwood foliage
		!21 = Total LAI
		!22 = Bud_h
        
        !assign the process error parameters to the correct data stream
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
		
		!SD2 is the process error term that scales with the size of the prediction
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
    			
			like_new = 0
   			prob_new(:) = 0.0
   			latent_old = latent
   			!$OMP PARALLEL	
			!$OMP DO PRIVATE(plotnum,V,Vv,tmp_sd,data_stream,mo)
    		do plotnum = 1,nplots
    		
    			!Deals with the fact that some plots using LAI as the state and other use Foliage biomass
    		    if(use_fol_state(plotnum) == 1) then
    				SD1(1) = Foliage_SD
    			else
    				SD1(1) = new_pars(52)
    			endif
    			
    			!Direct sampling of the latent states
    			if(pgroup == data_uncertainity_npar_group) then
        !STATE-SPACE WITHOUT GAP-FILLING
   			   		if(state_space == 1) then
    			   		do data_stream=1,7
    			   			!update initial value
    			   			!Equation 9.38a Clark 2007
  			   			
                            tmp_sd =(SD1(data_stream)+pred(data_stream,plotnum,mo_start_end(plotnum,1)) &
                                                        *SD2(data_stream))
    			   			V = 1/(1/init_uncert(data_stream,plotnum) +  &
    			   				 1/(tmp_sd*obs_gap_next(data_stream,plotnum,mo_start_end(plotnum,1))))
    			   			!Equation 9.38b Clark 2007 
    			   			Vv = V*((init_obs(data_stream,plotnum)/init_uncert(data_stream,plotnum) + &
    			   				latent(data_stream,plotnum,mo_start_end(plotnum,1))/ &
    			   				(tmp_sd*obs_gap_next(data_stream,plotnum,mo_start_end(plotnum,1)))))
    			   			latent(data_stream,plotnum,mo_start_end(plotnum,1)) = normal_sample(Vv,V)
    			   			!NEED TO FIX
    			   			!latent(data_stream,plotnum,mo_start_end(plotnum,1)) = init_obs(data_stream,plotnum)
    	  					if(latent(data_stream,plotnum,mo_start_end(plotnum,1)) < 0) then
    	  						latent(data_stream,plotnum,mo_start_end(plotnum,1)) = 0.00001
    	  					endif    	  					
    	  				end do
    					do mo = (mo_start_end(plotnum,1)+1),(mo_start_end(plotnum,2))
    						do data_stream=1,7
    	  						if(obs(data_stream,plotnum,mo) .NE. -99) then
    	  						!Equation 9.35 Clark 2007
    	  							if(obs_gap_next(data_stream,plotnum,mo) .NE. -99) then
    	  								!OBSERVATION WITH A PREVIOUS AND NEXT OBSERVATION
                                    	tmp_sd =(SD1(data_stream)+pred(data_stream,plotnum,mo)*SD2(data_stream))
    	  								V = 1/(1/(tmp_sd*obs_gap(data_stream,plotnum,mo)) &
    	  								 	+ 1/(tmp_sd*obs_gap_next(data_stream,plotnum,mo)) &
    	  								 	+ 1/obs_uncert(data_stream,plotnum,mo))
            							Vv = V*(pred(data_stream,plotnum,mo)/ &
            							 	(tmp_sd*obs_gap(data_stream,plotnum,mo)) &
            							 	+latent(data_stream,plotnum,mo)/ &
            							 	(tmp_sd*obs_gap_next(data_stream,plotnum,mo)) + &
            								obs(data_stream,plotnum,mo)/obs_uncert(data_stream,plotnum,mo))
            							latent(data_stream,plotnum,mo) = normal_sample(Vv,V) 
            						else
            						!LAST OBSERVATION
            					       tmp_sd = (SD1(data_stream)+ &
                                        	pred(data_stream,plotnum,mo)*SD2(data_stream))
    	  								V = 1/(1/(tmp_sd*obs_gap(data_stream,plotnum,mo)) &
    	  					 				+ 1/obs_uncert(data_stream, plotnum,mo))
     	  								Vv = V *(pred(data_stream,plotnum,mo) &
    	  									/(tmp_sd*obs_gap(data_stream,plotnum,mo)) &
     	  									+ obs(data_stream,plotnum,mo) &
     	  									/obs_uncert(data_stream,plotnum,mo))
    	  								latent(data_stream,plotnum,mo) = normal_sample(Vv,V)
    	  							endif	
            					endif
    	  					end do
    	  					do data_stream=12,18
    	  						if(obs(data_stream,plotnum,mo) .NE. -99) then
    	  							V = 1/(1/(SD1(data_stream)+pred(data_stream,plotnum,mo)*SD2(data_stream)) &
    	  							+ 1/obs_uncert(data_stream,plotnum,mo))
    	  							Vv = V*(latent(data_stream,plotnum,mo)/ &
    	  								(SD1(data_stream)+pred(data_stream,plotnum,mo)*SD2(data_stream)) + &
            							 obs(data_stream,plotnum,mo)/obs_uncert(data_stream,plotnum,mo))
									latent(data_stream,plotnum,mo) = normal_sample(Vv,V)  	 					 		
    	  						endif 
 								if(latent(data_stream,plotnum,mo) < 0) then
    	  							latent(data_stream,plotnum,mo) = 0.00001
    	  						endif   	  							
    	  					end do
        	  			end do
    				!STATE-SPACE WITH GAP-FILLING
        	  		elseif(state_space == 2) then
        	  		    	do data_stream=1,7
    			   			!update initial value
    			   			!Equation 9.38a Clark 2007
                            tmp_sd =(SD1(data_stream)+pred(data_stream,plotnum,mo_start_end(plotnum,1)) &
                                                        *SD2(data_stream))
    			   			V = 1/(1/init_uncert(data_stream,plotnum) + 1/tmp_sd)
    			   			!Equation 9.38b Clark 2007 
    			   			Vv = V*(init_obs(data_stream,plotnum)/init_uncert(data_stream,plotnum) &
    			   				+ latent(data_stream,plotnum,mo_start_end(plotnum,1))/tmp_sd)
    			   			latent(data_stream,plotnum,mo_start_end(plotnum,1)) = normal_sample(Vv,V)
    	  					if(latent(data_stream,plotnum,mo_start_end(plotnum,1)) < 0) then
    	  						latent(data_stream,plotnum,mo_start_end(plotnum,1)) = 0.00001
    	  					endif
    	  					
    	  					!LAST OBSERVATION
  
    	  					if(obs(data_stream,plotnum,mo_start_end(plotnum,2)) .NE. -99) then
                                tmp_sd = (SD1(data_stream)+ &
                                	pred(data_stream,plotnum,mo_start_end(plotnum,2))*SD2(data_stream))
    	  						V = 1/(1/tmp_sd &
    	  					 		+ 1/obs_uncert(data_stream, plotnum,mo_start_end(plotnum,2)))
     	  						Vv = V *(pred(data_stream,plotnum,mo_start_end(plotnum,2)) &
    	  							/tmp_sd &
     	  							+ obs(data_stream,plotnum,mo_start_end(plotnum,2)) &
     	  							/obs_uncert(data_stream,plotnum,mo_start_end(plotnum,2)))
    	  						latent(data_stream,plotnum,mo_start_end(plotnum,2)) = normal_sample(Vv,V)
    	  					else
    	  					
                                tmp_sd = (SD1(data_stream)+ &
                                	pred(data_stream,plotnum,mo_start_end(plotnum,2))*SD2(data_stream))
    	  						V = 1/(1/tmp_sd)
     	  						Vv = V *(pred(data_stream,plotnum,mo_start_end(plotnum,2)) &
    	  							/tmp_sd)
    	  						
    	  						latent(data_stream,plotnum,mo_start_end(plotnum,2)) = normal_sample(Vv,V)    	  					
    	  					endif
    	  				 	if(latent(data_stream,plotnum,mo_start_end(plotnum,2)) < 0) then
    	  							latent(data_stream,plotnum,mo_start_end(plotnum,2)) = 0.00001
    	  						endif   	
    	  					  	  				
    	  				end do
    					do mo = (mo_start_end(plotnum,1)+1),(mo_start_end(plotnum,2)-1)
    						do data_stream=1,7
    	  						if(obs(data_stream,plotnum,mo) .NE. -99) then
    	  						!Equation 9.35 Clark 2007
    	  								!OBSERVATION WITH A PREVIOUS AND NEXT OBSERVATION
    	  								diff = latent(data_stream,plotnum,mo+1) - &
    	  										 (pred(data_stream,plotnum,mo) - latent(data_stream,plotnum,mo))
                                    	tmp_sd =(SD1(data_stream)+pred(data_stream,plotnum,mo)*SD2(data_stream))
    	  								V = 1/(1/tmp_sd &
    	  								 	+ 1/tmp_sd &
    	  								 	+ 1/obs_uncert(data_stream,plotnum,mo))
            							Vv = V*(pred(data_stream,plotnum,mo)/ (tmp_sd) &
            							 	+latent(data_stream,plotnum,mo)/(tmp_sd) &
            								+ obs(data_stream,plotnum,mo)/obs_uncert(data_stream,plotnum,mo))
            							latent(data_stream,plotnum,mo) = normal_sample(Vv,V) 
            						else
            						diff = latent(data_stream,plotnum,mo+1) - &
    	  										 (pred(data_stream,plotnum,mo) - latent(data_stream,plotnum,mo))
                                    	tmp_sd =(SD1(data_stream)+pred(data_stream,plotnum,mo)*SD2(data_stream))
    	  								V = 1/(1/tmp_sd + 1/tmp_sd)
            							Vv = V*(pred(data_stream,plotnum,mo)/tmp_sd &
            							 	+latent(data_stream,plotnum,mo)/tmp_sd)
            							 latent(data_stream,plotnum,mo) = normal_sample(Vv,V) 
    	  							endif	
    	  						if(latent(data_stream,plotnum,mo) < 0) then
    	  							latent(data_stream,plotnum,mo) = 0.00001
    	  						endif
    	  					end do
    	  					do data_stream=12,18
    	  						if(obs(data_stream,plotnum,mo) .NE. -99) then
    	  							V = 1/(1/(SD1(data_stream)+pred(data_stream,plotnum,mo)*SD2(data_stream)) &
    	  							+ 1/obs_uncert(data_stream,plotnum,mo))
    	  							Vv = V*(latent(data_stream,plotnum,mo)/ &
    	  								(SD1(data_stream)+pred(data_stream,plotnum,mo)*SD2(data_stream)) + &
            							 obs(data_stream,plotnum,mo)/obs_uncert(data_stream,plotnum,mo))
									latent(data_stream,plotnum,mo) = normal_sample(Vv,V)  	 					 		
    	  						endif 
 								if(latent(data_stream,plotnum,mo) < 0) then
    	  							latent(data_stream,plotnum,mo) = 0.00001
    	  						endif   	  							
    	  					end do
        	  			end do
					!NON-STATE SPACE
        	  		else
  						do data_stream=1,7
  							!Initial condition
  						    tmp_sd =(SD1(data_stream)+pred(data_stream,plotnum,mo_start_end(plotnum,1)) &
                            	*SD2(data_stream))
    			   			V = 1/(1/init_uncert(data_stream,plotnum))
    			   			!Equation 9.38b Clark 2007 
    			   			Vv = V*((init_obs(data_stream,plotnum)/init_uncert(data_stream,plotnum)))
    			   			latent(data_stream,plotnum,mo_start_end(plotnum,1)) = normal_sample(Vv,V)
    			   		
    	  					if(latent(data_stream,plotnum,mo_start_end(plotnum,1)) < 0) then
    	  						latent(data_stream,plotnum,mo_start_end(plotnum,1)) = 0.00001
    	  					endif
  						
  							do mo = (mo_start_end(plotnum,1)+1),(mo_start_end(plotnum,2))
    	  						if(obs(data_stream,plotnum,mo) .NE. -99) then
                                                                tmp_sd =(SD1(data_stream)+ &
                                                                        pred(data_stream,plotnum,mo)*SD2(data_stream))
    	  							V = 1/(1/(tmp_sd*obs_gap(data_stream,plotnum,mo)) &
    	  								+ 1/obs_uncert(data_stream,plotnum,mo))
    	  							Vv = V*(latent(data_stream,plotnum,mo)/ &
    	  								(tmp_sd*obs_gap(data_stream,plotnum,mo)) + &
            							obs(data_stream,plotnum,mo)/obs_uncert(data_stream,plotnum,mo))
									latent(data_stream,plotnum,mo) = normal_sample(Vv,V) 	
									if(latent(data_stream,plotnum,mo) < 0) then
    	  								latent(data_stream,plotnum,mo) = 0.00001
    	  							endif 					 		
    	  						endif 					
    	  					end do 
    	  				end do    
    	  				do data_stream=12,18
  							do mo = (mo_start_end(plotnum,1)+1),(mo_start_end(plotnum,2))
    	  						if(obs(data_stream,plotnum,mo) .NE. -99) then
    	  						V = 1/(1/(SD1(data_stream)+pred(data_stream,plotnum,mo)*SD2(data_stream)) &
    	  							+ 1/obs_uncert(data_stream,plotnum,mo))
    	  							Vv = V*(latent(data_stream,plotnum,mo)/ &
    	  								(SD1(data_stream)+pred(data_stream,plotnum,mo)*SD2(data_stream)) + &
            							 obs(data_stream,plotnum,mo)/obs_uncert(data_stream,plotnum,mo))
									latent(data_stream,plotnum,mo) = normal_sample(Vv,V) 	
									if(latent(data_stream,plotnum,mo) < 0) then
    	  								latent(data_stream,plotnum,mo) = 0.00001
    	  							endif 					 		
    	  						endif 					
    	  					end do 
    	  				end do     	  			  		 	  			  			
					endif
    	  		endif
    		    		   	
   				!---LOOP THROUGH PLOTS TO CALCULATE LIKELIHOOD OF EACH PLOT 	
 				if(fit_plot(plotnum) == 1) then
  					call likelihood(plotnum &
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
        			else
        			prob_new(plotnum) = 0.0
        		endif
 			end do 			
 			!$OMP END DO
	 		!$OMP END PARALLEL	
	 			
			!---SUM UP LIKELIHOOD ACROSS PLOTS	
			like_new = sum(prob_new)
					
			!--CALCULATE PROBABILITY OF ALL PRIORS (MULTIPLY INDIVIDUAL PROABILITIES)
   			prior = 0
   
   			do p = 1,max_fitted_pars
   				if(fix_par(p) .NE. 1) then
     				if(prior_dist(p)==1) then !uniform distribution
       					prior = prior + log(uniform_pdf(new_pars(p), prior_parameter1(p), prior_parameter2(p)))
     				else if(prior_dist(p)==2 ) then !normal distribution
      					prior = prior + log(normal_pdf(new_pars(p), prior_parameter1(p), prior_parameter2(p)))
     				else if(prior_dist(p)==3) then !Holding for new distribution
       					!prior = prior + dlnorm(new_pars[p], meanlog=prior_parameter1[p], sdlog=prior_parameter2[p], log = TRUE)
     				else 
       					print *, 'PRIOR DISTRIBUTION NOT SPECIFIED!'
     				endif
     			endif
  			end do
  	

   			!---ECOLOGICAL CONSTRAINTS
   			!if(new_pars(5)<new_pars(6))then  !Root allocation decreases with FR
   			!		prior=-1e30
   			!endif

   			if(new_pars(7)<new_pars(8)) then  !SLA decreases as the forest ages
     			prior=-1e30
   			end if
   
   			if(new_pars(28) < new_pars(27) .OR. new_pars(29) < new_pars(28)) then  !T OPT has to be between MAX and MIN
     			prior=-1e30
   			end if
   
   			if(new_pars(1) < new_pars(2)) then  !Allocation to foliage declines as the forest ages
     			prior=-1e30
   			end if 

   			if(use_sm_edc == 1) then
   			! Ensure that the soil moisture modifer covers most of the 0 - 1 range
   				if((1.0D0 / (1.0D0 + ((1.0D0 - 0.0) / new_pars(38)) ** new_pars(39))) > 0.10) then
      				prior=-1e30
   				end if 
   				if((1.0D0 / (1.0D0 + ((1.0D0 - 1.0) / new_pars(38)) ** new_pars(39))) < 0.99) then
      				prior=-1e30        
   				end if
   			endif
   			! Ensure that a SI of 40m has add FR near 1 
   			if(use_fr_edc == 1) then
   				if(fr_model == 2) then
       				if((1/(1+exp((FR1_new + FR2_new*10-FR3_new*40)))) < 0.99) then
       					!if((FR1_new + FR2_new*10 +FR3_new*40) < 0.99) then
           				!prior=-1e30
       				endif
       				if((1/(1+exp((FR1_new + FR2_new*30-FR3_new*40)))) < 0.99) then
       					!if((FR1_new + FR2_new*30 +FR3_new*40) < 0.99) then
           				!prior=-1e30
       				endif
    			endif
    		endif

    		if(use_age_edc == 1) then
    			!Ensure that the age effect does not get too strong at young ages
    			if((1.0 / (1.0 + ((10/new_pars(31)) / new_pars(33)) ** new_pars(32))) <0.95) then
        		prior=-1e30
    			endif 
    		endif
   
			!---SUM UP THE LOG LIKELIHOOD FOR ALL THE PLOTS, THE PRIORS, AND THE PROB OF LATENT STATES
 
			!PNEW IS THE PROBABILITY FROM THE CURRENT ITERATION, SUMMED ACROSS ALL PLOTS

             pnew = like_new  + prior
             
			if(pnew < -1e30) then
      			pnew=-1e30
   			endif

			!---ACCEPTANCE CRITERION
			CALL RANDOM_NUMBER(z)
   			!z = randn(0)
   			r = exp(pnew-pnow)
   			
   			
   			! METROLOPIS HASTING STEP
  			!---IF THE CURRENT SET OF PARAMETERS IS MORE PROBABLE THAN THE PREVIOUS SET THEN ALWAYS ACCEPT
			!---IF THE CURRENT SET OF PARAMETERS IS LESS PROBABLE THAN THE PREVIOUS SET THEN ACCEPT IN PROPORTION TO HOW MUCH WORSE
			if(pgroup == data_uncertainity_npar_group .AND. pnew>-1e30) then
				!Always accept the latent states because they are directly sampled
				global_accept(pgroup) = global_accept(pgroup) + 1
        		local_accept(pgroup) = local_accept(pgroup) + 1
        		current_pars = new_pars
        		current_like = like_new
        		pred = pred_new
        		pnow = pnew
        	elseif(pgroup == data_uncertainity_npar_group .AND. pnew<=-1e30) then
        		!If the latent state cause the model to break then don't accept
        		latent = latent_old     	
        	else
    			if(z<r .AND. pnew>-1e30) then
    				global_accept(pgroup) = global_accept(pgroup) + 1
        			local_accept(pgroup) = local_accept(pgroup) + 1
        			current_pars = new_pars
        			current_like = like_new
        			pred = pred_new  		
        			pnow = pnew
    			else
        			current_pars = current_pars
        			current_like = current_like
        			pred = pred
        			pnow = pnow
    			end if		
    		end if

        	!---ADAPT THE JUMP SIZE USED IN RANDOM WALK TO ACHIEVE A GOAL ACCEPTANCE RATE
        	if(mod(iter,nadapt) == 0 .AND. iter .GE. start_adapt) then
				local_accept_rate = dble(local_accept(pgroup))/dble(nadapt)
        		if (local_accept_rate < lower_accept_bound) then
                	do p = 1,nopars
                		if(par_group(p) == pgroup .AND. fix_par(p) == 0) then
                    		jump_pars(p)=jump_pars(p)/adaptfac
                    	endif
                	end do
        		elseif (local_accept_rate > upper_accept_bound) then        
                	do p = 1, nopars
                		if(par_group(p) == pgroup  .AND. fix_par(p) == 0) then
                    		jump_pars(p)=jump_pars(p)*adaptfac
                    	endif
                	end do
        		end if !conditional if acceptance rate low or high
        		local_accept(pgroup) = 0
        		if(mod(iter,10000) == 0) then
        			write(*,*)"*********************"
           			write(*,*)'Pgroup = ',pgroup
           			Format = "(A,I8,A,I8)"
           			write(*,Format)"Total Accepted = ",global_accept(pgroup)
           			Format = "(A,F4.2)"
           			write(*,Format)"Local Acceptance rate ",local_accept_rate
           			Format = "(A,F10.2)"
           			write(*,Format)"Current Log Likelihood = ",pnow
           			Format = "(A,I8)"
          			write(*,Format)"Current iteration = ",iter
         		endif
     		endif
     	
		end do

    	if(iter == sample_index(sample)) then
    		tracked_plot(sample,:,:) = latent(:,tracked_plotnum,:) 
        	accepted_pars_thinned_burned(sample,:) = current_pars(:)
        	like_chain(sample) = pnow 
        	sample = sample + 1
    	endif

	end do
	
	
	


  ! return back to the subroutine then

end subroutine DAPPER_MCMC
