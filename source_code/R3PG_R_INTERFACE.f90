subroutine r3pg_interface(output_dim,met,pars,site,thin_event,out_var,nopars,nomet,nosite &
                      ,nooutputs,nomonths_plot,nothin,exclude_hardwoods,mo_start_end,nmonths)

  use R3PG_MODEL_MOD, only: R3PG_MODEL

  ! subroutine specificially deals with the calling of the fortran code model by
  ! R

  implicit none
  ! declare input variables
  
  
  integer, intent(in) :: nopars         & ! number of paremeters in vector
                        ,output_dim     & !
                        ,nosite			&
                        ,nomet          & ! number of meteorological fields
                        ,nooutputs        & ! number of model pools
                        ,nomonths_plot        &   ! number of days in simulation
                        ,nmonths     &
                        ,nothin 		&
                        ,exclude_hardwoods &
                        ,mo_start_end(2)

   double precision, intent(inout) :: met(nomet,nmonths)   & ! met drivers, note reverse of needed
                       ,pars(nopars)         & ! number of parameters
                       ,site(nosite)        ! number of site descriptors
 
    double precision, intent(in) :: thin_event(nmonths)

  ! output declaration
   double precision, intent(out), dimension(nomonths_plot,output_dim) :: out_var

  ! local variables

  double precision, dimension(nooutputs) :: OUTPUT
  
   integer :: mo,mo_index


  ! zero initial conditions
 	OUTPUT = 0.
 	out_var = 0.
    mo_index = 0
        
    do mo = mo_start_end(1),mo_start_end(2)
      	mo_index = mo_index + 1

     	! call the models
		call R3PG_MODEL(output_dim,met(:,mo),pars,site,thin_event(mo),nopars,nomet, &
					nosite,nooutputs,1,1,nothin,exclude_hardwoods,OUTPUT)
					
     	! now allocate the output the our 'output' variable
     	out_var(mo_index,:)  = OUTPUT
     	
		if(OUTPUT(2) == 12) then
			site(3) = OUTPUT(1)+1 !InitialYear
			site(4) = 1  !InitialMonth
		else
			site(3) = OUTPUT(1) !InitialYear
			site(4) = OUTPUT(2)+1  !InitialMonth	
		endif	
		site(5) = OUTPUT(3) + (1.0/12.) !StartAge
		site(26) = OUTPUT(4) !LAI
		site(8) = OUTPUT(5) !WS
		site(20) = OUTPUT(6)   !WCR
		site(7) = OUTPUT(7)  !WRi
		site(9) = OUTPUT(8) !StemNo
		
		site(27) = output(9)  !Hardwood LAI
		site(25) = OUTPUT(26) !Hardwood Bud
		site(18) = OUTPUT(10) !WS_H 
		site(24) = OUTPUT(11)  !WCR_h
		site(19) = OUTPUT(12) !WR_H

		site(10) = OUTPUT(14) ! ASW
		
		site(6) = OUTPUT(22) !WFi
		site(17) = OUTPUT(23) !WF_H	
		
		


    enddo
   
  ! return back to the subroutine then
  return

end subroutine r3pg_interface
