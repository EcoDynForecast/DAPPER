module R3PG_MODEL_MOD
       
contains

	!-----------------------------------------------------
	!-  main r3PG code
	!-----------------------------------------------------

subroutine R3PG_MODEL(output_dim,met,pars,site,thin_event,nopars,nomet, &
					nosite,nooutputs,nomonths,nometmonths,nothin,exclude_hardwoods,OUTPUT)

	implicit none

  	! declare input variables
  	integer, intent(in) :: nopars         & ! number of paremeters in vector
                        ,output_dim     & !
                        ,nosite			&
                        ,nomet          & ! number of meteorological fields
                        ,nooutputs      & ! number of model fluxes
                        ,nomonths         &  ! number of days in simulation
                        ,nometmonths    &
                        ,nothin			&
                        ,exclude_hardwoods

  	double precision, intent(in) :: met(nomet,nomonths)   & ! met drivers, note reverse of needed
                       ,pars(nopars)        & ! number of parameters
                       ,site(nosite)        &  ! number of site descriptors
                       ,thin_event(nomonths)

  	! output declaration
  	double precision,  dimension(nomonths,nooutputs), intent(inout) :: OUTPUT 
  	
  		INTEGER, DIMENSION(1:12), parameter :: daysInMonth = &
                (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    !DP = double precision
    INTEGER, parameter :: DP=SELECTED_REAL_KIND(8)
                
	double precision, parameter :: ln2 = 0.693147181D0  
 	double precision, parameter :: Pi = 3.141592654D0
	!----------------------------------------------------------------
	! 3PG PARAMETERS

	double precision :: pFS2,pFS20,stemConst,stemPower,pRx,pRn,  &              
        gammaFx,gammaF0,tgammaF,rttover,                  &
        Tmin,Topt,Tmax,kF,SWconst0,SWpower0,              &
        m0,fN0,MaxAge,nAge,rAge,                           &
        gammaN1,gammaN0,tgammaN,                  &
        wSx1000,thinPower,mF,mR,mS,SLA0,SLA1,tSLA,        &
        k,fullCanAge,MaxIntcptn,LAImaxIntcptn,            &
        alpha,Y,MinCond,MaxCond,LAIgcx,                      & 
        CoeffCond,BLcond,tBB,             &
        rho0,rho1,tRho,aH,nHB,nHN,aV,nVB,nVN,             &
        Qa,Qb,gDM_mol,molPAR_MJ,pfsConst,pfsPower,        & 
        Density,genetics,competition,formFactor,volRatio, &    
        tempAlpha, fCalpha700,fCalpha700_h,fCg700,SWconst1,SWconst2, &
        SWpower1,SWpower2,mort_rate,						 &
        alpha_h, pfs_h, pR_h, mort_rate_h,SLA_h,pCRS,maxavDBH, &
        fCpFS700,Leaftover,HardPineCondRatio

	!----------------------------------------------------------------
	! STAND VARIABLES
   
    INTEGER :: InitialYear,InitialMonth,     &
                PlantedYear,PlantedMonth,soilClass
    double precision  :: StartAge,StandAge,ASWi,ASW,StemNoi, & 
                StemNo,WFi,WF1,WF2,WRi,WR,WSi,WS, & 
                avDBH,TotalW,BasArea,LAI, &
                FR,Lat,MaxASW,MinASW,SLA,CanCover,  &
                WF_h,WS_h,WR_h,WCR_h,LAI_h,WBud_h, &
                WFi_H,WRi_H,WSi_H, WCR, WCRi,ThroughFall, IrrigRate, LAI_i,LAI_i_h	
                
	!---------------------------------------------------------------- 
	! WEATHER DATA VARIABLES

    double precision  :: DayLength,FrostDays,SolarRad,Tav,VPD,Rain,CO2
    
	!----------------------------------------------------------------
	! INTERMEDIATE RESULTS VARIABLES
    double precision  :: m,alphaC,L_epsilon,RAD,PAR,                    & 
                lightIntcptn,fAge,fT,fFrost,fVPD,fSW,fNutr, fcomp,         &
                PhysMod,CanCond,Transp,EvapTransp,                 &
                AvStemMass,wSmax,delWF,delWR,delWS,                  &
                GPP,NPP,pR,pS,pF,pFS,incrWF,incrWR,incrWS,        &
                lossWF,lossWR,gammaF,gammaN,SWconst,SWpower,  &
                APAR,APARu,GPPmolc,GPPdm,water_balance_error,  &
                Litter,TotalLitter,delStemNo,stemNoTemp,delStemN,        &
                canCovComp, fCalphax,fCalphax_h,fCg0, fCalpha,fCalpha_h, fCg, &
                conductance,LAIratio,runoff,delWCR,c_allometry,nlc,Transp_pine, 				&
                !Hardwoods
                lightIntcptn_h,APAR_h,alphaC_h,GPPmolc_h,GPPdm_h,NPP_h,Transp_hard, &
				pS_h,pF_h,delWF_h,delWR_h,delWCR_h,delWS_h,Littfall_h,delLitter_h, &
				delRoots_h,TotalW_h,TotalLitter_h,delWBud_h,leaved_out, &
				delWBud_to_WF_h,Tmin_met,fCpFS0,fCpFS,delStems,delLitter, &
				delRoots,Littfall,Vob,GPPc,GPPc_h
	!----------------------------------------------------------------     
	! 3PG Variables

    INTEGER :: month,year,metMonth,dayofyr,calyear,calmonth,TotalMonths
    double precision  :: RelAge,MoistRatio,Intcptn,age,LAImaxRatio

    INTEGER :: output_index
    LOGICAL :: clm	
    INTEGER :: i
    	
    !----------------------------	
	! SET PARAMETERS
	!----------------------------
	pFS2 = pars(1)
	pFS20 = pars(2)
	StemConst = pars(3)
	StemPower = pars(4)
	pRx = pars(5)
	pRn = pars(6)
	SLA0 = pars(7)
	SLA1 = pars(8)
	tSLA = pars(9)
	k = pars(10)
	fullCanAge = pars(11)
	MaxIntcptn = pars(12)
	LAImaxIntcptn = pars(13)
	alpha = pars(14)
	MaxCond = pars(15)
	LAIgcx = pars(16)
	CoeffCond = pars(17)
	BLcond = pars(18)
	wSx1000 = pars(19)
	thinPower = pars(20)
	mF = pars(21)
	mR = pars(22)
	mS = pars(23)
	Rttover = pars(24)
	m0 = pars(25)
	fN0 = pars(26)
	Tmin = pars(27)
	Topt = pars(28)
	Tmax = pars(29)
	kF = pars(30)
	MaxAge = pars(31)
	nAge = pars(32)
	rAge = pars(33)
	y = pars(34)
	fCalpha700_h = pars(35)
	fCalpha700= pars(36)
	fCg700= pars(37)
	SWconst1 = pars(38)
	SWpower1 = pars(39)
	mort_rate = pars(40)
	alpha_h = pars(41)
	pfs_h = pars(42)
	pR_h = pars(43)
	mort_rate_h = pars(44)
	SLA_h = pars(45)
	pCRS = pars(46)
	fCpFS700 = pars(47)
	Leaftover = pars(48)
	
    mR = 0.0
    fN0 = 0.0
    m0 = 0.0
    BLcond = 0.1 
	molPAR_MJ = 2.3
	Qa = -90
	Qb = 0.8
	gDM_mol = 24
	
	HardPineCondRatio = pars(18)

	!Unique parameters for the Duke Forest
    if(site(23) == 1) then
    	!mort_rate = pars(22)
        !thinPower = pars(18)
        !wSx1000 = pars(48)
        !mS = pars(26)
        pCRS = pars(25)
    endif
    
    !-----------------------------
    !-- Set site level variables
    !-----------------------------
	PlantedYear = site(1) 
	PlantedMonth = site(2)
	InitialYear = site(3)
	InitialMonth = site(4) 
	StandAge = site(5)
	WFi = site(6)
	WRi = site(7) 
	WSi = site(8) 
	StemNoi = site(9) 
	ASWi = site(10)
	Lat = site(11) 
	FR = site(12)
	SoilClass = site(13) 
	MaxASW = site(14)
	MinASW = site(15)
	TotalMonths = site(16)
	WF_H = site(17)
	WS_H = site(18) 
	WR_H = site(19) 
	WCRi = site(20)
	IrrigRate = site(21)
	ThroughFall = site(22)
	WCR_h = site(24)
	WBud_h = site(25)
	LAI_i = site(26)
	LAI_i_h = site(27)
	
	genetics = 1
    clm = .FALSE.
    
	!------------------------------------------------
	!-  INITIALIZE STAND STRUCTURE 
	!------------------------------------------------

	! INITILIZE SOIL WATER
    SWconst = SWconst1
    SWpower = SWpower1
    
	! Initial ASW must be between min and max ASW
    if(MinASW > MaxASW) then 
        MinASW = MaxASW
    endif
    if(ASWi <= MinASW) then 
        ASWi = MinASW
    elseif(ASWi >= MaxASW) then 
        ASWi = MaxASW
    endif
    ASW = ASWi

	! ----  INITILIZE STAND FACTORS
	pfsPower = Log(pFS20 / pFS2) / Log(20.0D0 / 2.0D0)
    pfsConst = pFS2 / 2.0D0 ** pfsPower
   
    WF1 = WFi
    WF2 = 0.0
    WS = WSi
    WCR =WCRi
    WR = WRi 
    StemNo = StemNoi
    
    if(LAI_i .NE. -99) then  !OVERWRITE INITIAL FOLIAGE IF INITIAL LAI IS SUPPLIED
    	LAI = LAI_i
    	SLA = SLA1 + (SLA0 - SLA1) * Exp(-ln2 * (StandAge / tSLA)**2)
    	WFi = LAI_i / (SLA * 0.1D0)
    else
    	if(clm) then
    		LAI  = ((SLA0*0.1D0)*(exp(tSLA*0.1D0*((WF1+WF2)))-1))/(tSLA*0.1D0)
        	SLA = LAI / (WF1+WF2)
    	else
    		SLA = SLA1 + (SLA0 - SLA1) * Exp(-ln2 * (StandAge / tSLA)**2)
        	LAI = (WF1+WF2) * SLA * 0.1D0
    	endif        
    endif
    

    !Hardwood LAI   
    if(LAI_i_h .NE. -99) then
    	LAI_h = LAI_i_h
    	WF_h  = LAI_h/ (SLA_h * 0.1D0) 
	endif
    	
    if((InitialMonth) .LE. 3 .AND. WBud_h == 0.0) then
		WBud_h = WF_H
		WF_h = 0.0
    endif

    LAI_h = WF_h * SLA_h * 0.1D0  

    AvStemMass = WS * 1000.0D0 / StemNo     !kg/tree 
    avDBH = (AvStemMass / stemConst) ** (1.0D0 / stemPower)
    BasArea = (((avDBH / 200.0D0) ** 2.0D0) * Pi) * StemNo 
    
 	TotalLitter = 0
 	leaved_out = 0
 	output_index = 1
 	
 	!--------------------------
 	!-- Initialize time
 	!--------------------------
 	calmonth = InitialMonth
 	metmonth = 1
 	calyear = InitialYear
 	
 	
 	do i = 1,(InitialMonth-1)
    	if (i == 1) then 
            dayofyr = -15
        end if
    	dayofyr = dayofyr + daysInMonth(i)
    end do
    
	
 	!--------------------------
 	!---  START SIMULATION
 	!--------------------------
  	DO month = 1,nomonths
  	  	
		!-----------------------------------
		! MONTHLY CLIMATE
		!-----------------------------------
		
		if (calmonth == 1) then 
            dayofyr = -15
        end if
        
        dayofyr = dayofyr + daysInMonth(calmonth)
        DayLength = 86400.0D0 * getDayLength(Lat,dayofyr) !Seconds
 		! MET: index 1 = Tmin, 2 = T max, 3 = Rain, 4 = SolarRad, 5 = FrostDays
 		SolarRad = met(4,metmonth)
		Tav = (met(1,metmonth) + met(2,metmonth)) / 2.0D0
		VPD = getVPD(met(2,metmonth),met(1,metmonth))
		FrostDays = met(5,metmonth)
		Rain = met(3,metmonth)
		CO2 = met(6,metmonth)
        Tmin_met = met(1,metmonth)
        
		!-----------------------------------
		! ASSIGN SILVICULTURE
		!-----------------------------------
       	!the value of fertility should be less than 1.0
       	if (FR > 1.0) then
           	FR = 1.0D0
       	end if
    
		!-----------------------------------
		! ENVIRONMENTAL MODIFIERS
		!-----------------------------------
		!calculate temperature response function
    	if ((Tav <= Tmin) .OR. (Tav >= Tmax)) then 
        	fT = 0
    	else 
        	fT = ((Tav - Tmin) / (Topt - Tmin)) * &
                  ((Tmax - Tav) / (Tmax - Topt)) ** &
                  ((Tmax - Topt) / (Topt - Tmin))
    	end if
		
    	!calculate VPD modifier
    	fVPD = Exp(-CoeffCond * VPD)
    	
    	!calculate soil water modifier
   	    MoistRatio = ASW / MaxASW
    	fSW = 1.0D0 / (1.0D0 + ((1.0D0 - MoistRatio) / SWconst) ** SWpower)
    	
    	!calculate soil nutrition modifier
    	fNutr = fN0 + (1.0D0 - fN0) * FR
        
    	!calculate frost modifier
        if(Tmin_met < 0) then
            fFrost = 1.0D0 - ((1.0D0-exp(kF*Tmin_met)) * (FrostDays/ 30.0D0))
        else
           fFrost = 1.0D0
        endif 

        ! calculate age modifier
        RelAge = StandAge/MaxAge
        if(nAge == 0) then
          fAge = 1
        else
    	  fAge = (1.0D0 / (1.0D0 + (RelAge / rAge) ** nAge))
    	end if

        ! Calculate atmosphere CO2 modifier
    	fCalphax = fCalpha700 / (2 - fCalpha700)
    	fCg0 = fCg700 / (2 * fCg700 - 1)
    	fCalpha = fCalphax * CO2 / (350D0 * (fCalphax - 1) + CO2)
        fCg = fCg0 / (1 + (fCg0 - 1) * CO2 / 350D0)
   
    	! competition equation
		fcomp =  1.0

        !-----------------------------------
        !--  DO WATER BALANCE
        !-----------------------------------
        LAIratio = (LAI+LAI_h) / LAIgcx
        conductance = MaxCond * MIN(1.0D0, LAIratio) * fCg * fVPD* fSW * fAge

        if (conductance == 0 .OR. conductance < 0.0001) then
        	CanCond = 0.0001D0
        else
            CanCond = conductance
        end if

        Transp = REAL(daysInMonth(calmonth),8) * &
        getTranspiration(SolarRad, VPD, DayLength, BLcond,CanCond,Qa,Qb)*fFrost  !rainfall interception
        Transp_pine = Transp*(HardPineCondRatio*LAI/(HardPineCondRatio*LAI + LAI_h))
        Transp_hard = Transp*(LAI_h/(HardPineCondRatio*LAI + LAI_h))
        if (LAImaxIntcptn > 0) then
            LAImaxRatio = (LAI+LAI_h) / LAImaxIntcptn
            Intcptn = MaxIntcptn * MIN(1.0D0, LAImaxRatio)
    	else
            Intcptn = MaxIntcptn
        end if

        EvapTransp = Transp + Intcptn * Rain
                
    	if(calmonth >= 3 .AND. calmonth <= 11) then !Irrigation occurs in growing season
        	ASW = ASW + Rain*ThroughFall + IrrigRate - EvapTransp
        else
        	ASW = ASW + Rain*ThroughFall - EvapTransp
        endif
        
        water_balance_error = 0

        if (ASW < MinASW) then
            water_balance_error = MinASW - ASW
            ASW = MinASW
        else if (ASW >= MaxASW) then
        	runoff = ASW - MaxASW
            ASW = MaxASW
        end if

		!-----------------------------------   
		! DETERMINE NPP
		!-----------------------------------
   		
		CanCover = 1
        if (fullCanAge > 0) then 
            CanCover = StandAge / fullCanAge
        end if                  
        if (CanCover > 1) then 
            CanCover = 1
        end if   
        
        lightIntcptn = (1.0D0 - (Exp(-k * LAI)))
        
    	! Calculate PAR, APAR, APARu and GPP
    	RAD = SolarRad * REAL(DaysInMonth(calmonth), 8)        !MJ/m^2
    	PAR = RAD * molPAR_MJ                                  !mol/m^2
    	APAR = PAR * lightIntcptn* CanCover
    	APARu = APAR  
    	! alphac has a new multiplier genetics
    	alphaC = alpha * fNutr * fT * fFrost * genetics * fCalpha * fVPD* fSW * fAge * fcomp
    	GPPmolc = APARu * alphaC             !mol/m**2    	
        GPPc = (GPPmolc * 12.0) / 100.0D0      !ton/ha	
    	GPPdm = (GPPc * 2.0)
    	NPP = GPPdm * y !assumes respiratory rate is constant

    
 
    	!---UNDERSTORY HARDWOOD
        fCalphax_h = fCalpha700_h / (2 - fCalpha700_h)
        fCalpha_h = fCalphax_h * CO2 / (350D0 * (fCalphax_h - 1) + CO2) 
    	lightIntcptn_h = (1.0D0 - (Exp(-k * LAI_h)))
    	APAR_h = (PAR - APAR) *lightIntcptn_h  !use the light left over by the understory pines\
    	alphaC_h = alpha_h * fNutr * fT * fFrost * fVPD * fSW * fCalpha_h
    	
    	GPPmolc_h = APAR_h * alphaC_h             !mol/m**2
    	GPPc_h = (GPPmolc_h * 12.0) / 100.0D0      !ton/ha	
    	GPPdm_h = (GPPc_h * 2.0)
    	NPP_h = GPPdm_h * y !assumes respiratory rate is constant
    	
		!-----------------------------------
		!- DETERMINE BIOMASS CHANGES
		!-----------------------------------

 		!calculate partitioning coefficients
    	m = m0 + (1.0D0 - m0) * FR
    	    	
    	fCpFS0 = fCpFS700 / (2 * fCpFS700 - 1)
    	fCpFS = fCpFS0 / (1 + (fCpFS0 - 1) * CO2 / 350D0)

        pFS = (pfsConst * avDBH ** pfsPower)*fCpFS
        
        if(clm) then
          pFS = pFS2 
        endif 
     	
    	pR = pRx
        !pR = pRn + (pRx - pRn)*(1-FR)
    	!pRx * pRn / (pRn + (pRx - pRn) * fVPD * fSW * m)   	

        c_allometry = 1 + pR + (1/pFS)*(1+pCRS)
    	nlc = NPP/c_allometry
        delWF = nlc
        delWR = nlc*pR
        delWS = nlc*(1/pFS)
        delWCR = nlc*(1/pFS)*pCRS
        

    	!calculate litterfall & root turnover
    	delRoots = Rttover * WR
    	
        WF1 = WF1 + delWF	
        
        delLitter = 0.0
        if(calmonth > 9 .OR. calmonth < 5) then
        delLitter = WF1*mF
        WF1 = WF1 -  delLitter
        endif
        
        WF2 = 0.0
        
  

        
        !if(calmonth == 10) then
    	!	WF2 = WF2+WF1
    	!	WF1 = 0.0
    	!endif
    	
    	!delLitter = WF2*mF
    	!WF2 = WF2 - delLitter
        
        !if(calmonth == 11) then
    	!	delLitter = WF2
    	!else
    	!	delLitter = 0.0
    	!endif
  
        !if(delLitter <= WF2) then     
        !    WF2 = WF2 - delLitter
        !else 
        !    WF2 = 0
        !    WF1 = WF1 - (delLitter - WF2)
        !endif

        WR = WR + delWR - delRoots
        WS = WS + delWS
        WCR = WCR + delWCR
    	TotalW = WF1 + WF2 + WR + WS
    	TotalLitter = TotalLitter + delLitter 

    	!----UNDERSTORY HARDWOOD
		pR_h = pR 
		!ASSUME SIMILAR LEAF TO FINE ROOT ALLOCATION BETWEEN PINES AND HARDWOOD BECAUSE FINE ROOT
		! DATA IS NOT SEPARATED
        c_allometry = 1+ pR_h + (1/pFS_h)*(1+0.2)
        nlc = NPP_h/c_allometry
		!---Phenology 
		delWBud_h = nlc    	
        delWR_h = nlc * pR_h  
        delWS_h =nlc * (1/pFS_h) 
        delWCR_h = nlc * (1/pFS_h)* 0.2 
             
		if(calmonth == 4) then
				delWBud_to_WF_h = WBud_h
		else
			    delWBud_h = nlc 
    	        delWBud_to_WF_h = 0.0
			    delLitter_h = 0.0
		endif
		if(calmonth == 10) then
			delLitter_h = WF_h
		else
			delLitter_h = 0.0
		endif
		
    	!calculate root turnover
    	delRoots_h = Rttover * WR_h

		WBud_h = WBud_h + delWBud_h - delWBud_to_WF_h
        WF_h = WF_h + delWBud_to_WF_h - delLitter_h
        WR_h = WR_h + delWR_h - delRoots_h
        WS_h = WS_h + delWS_h
        WCR_h = WCR_h + delWCR_h
        
 
    	TotalW_h = WF_h + WR_h + WS_h + WCR_h
    	TotalLitter_h = TotalLitter_h + delLitter_h 
    	
		!-----------------------------------
		!--  DO THINNING
		!-----------------------------------
	
		if(thin_event(month) > 0.0) then
        	delStemN = thin_event(month)
        	WF1 = WF1  - mS*delStemN * (WF1 / StemNo)
        	WF2 = WF2  - mS*delStemN * (WF2 / StemNo)
        	WR = WR  - mS*delStemN * (WR / StemNo)
        	WS = WS - mS*delStemN * (WS / StemNo)
        	WCR = WCR - mS*delStemN * (WCR / StemNo)
       		WS_h = WS_h - WS_h*(delStemN/StemNo)
     		WF_h = WF_h - WF_h*(delStemN/StemNo)
     		WR_h = WR_h - WR_h*(delStemN/StemNo)
     		WCR_h = WCR_h - WCR_h*(delStemN/StemNo) 	
        	StemNo = StemNo - delStemN
        endif

		!-----------------------------------
		!--  CALCULATE MORTALITY
		!-----------------------------------
		!-- Self-thinning mortality
  		wSmax = wSx1000 * (1000.0D0 / StemNo) ** thinPower
    	AvStemMass = WS * 1000.0D0 / StemNo
    	
    	delStems = 0
        if(clm) then
            WS = WS - mS*(WS*mort_rate)
            WF1 = WF1 - (WF1*mort_rate)
            WF2 = WF2 - (WF2*mort_rate)            
            WR = WR - (WR*mort_rate)
            WCR = WCR - mS*(WCR*mort_rate)
            StemNo =  StemNo - StemNo*mort_rate
            AvStemMass = WS * 1000.0D0 / StemNo
            wSmax = wSx1000 * (1000.0D0 / StemNo) ** thinPower
        else
    	    if (wSmax < AvStemMass) then
        		delStems = GetMortality(StemNo, WS,mS,wSx1000,thinPower)
                ! ASSUMES TREES THAT DIED THROUGH SELF-THNNING DIDN"T HAVE FOLIAGE
                ! ASSUMES TREES THAT DIED ARE SMALLER THAN THE AVERAGE SIZE TREE
                !WF1 = WF1 - mS * delStems*(WF1/StemNo)
            	!WF2 = WF2 - mS * delStems*(WF2/StemNo) 
            	!delLitter = delLitter + mS * delStems*(WF1/StemNo) + mS * delStems*(WF2/StemNo) 
        		WR = WR - mS * delStems * (WR / StemNo)
                WS = WS - mS * delStems * (WS / StemNo)
                WCR = WCR - mS * delStems * (WCR / StemNo)
        		StemNo = StemNo - delStems
        		wSmax = wSx1000 * (1000.0D0 / StemNo) ** thinPower       	
        		AvStemMass = WS * 1000.0D0 / StemNo
        	
     		end if
     	    !--- Non-density dependent mortality
     	    delStems = StemNo*mort_rate
  			WF1 = WF1 - delStems*(WF1/StemNo)
            WF2 = WF2 - delStems*(WF2/StemNo) 
            delLitter = delLitter + delStems*(WF1/StemNo) + delStems*(WF2/StemNo) 
        	WR = WR - delStems * (WR / StemNo)
            WS = WS - delStems * (WS / StemNo)
            WCR = WCR - delStems * (WCR / StemNo)
        	StemNo = StemNo - delStems
        	wSmax = wSx1000 * (1000.0D0 / StemNo) ** thinPower       	
        	AvStemMass = WS * 1000.0D0 / StemNo
     	endif
     	
     	  
     	
     	!---Understory hardwood mortality
     	!---ASSUME SAME AS PINES BECAUSE HARDWOOD STEM DENSITY DATA NOT AVIALABLE
     	WS_h = WS_h - (WS_h*mort_rate_h)
     	WF_h = WF_h - (WF_h*mort_rate_h)
     	WR_h = WR_h - (WR_h*mort_rate_h)
     	WCR_h = WCR_h - (WCR_h*mort_rate_h)
     	if(exclude_hardwoods == 1) then
     		WS_h = 0.00
     		WF_h = 0.00
     		WR_h = 0.00
     		WCR_h = 0.00
        endif
        
		!-----------------------------------
		!-- DO VOLUME CALCULATIONS
		!----------------------------------- 
	

        !avStemMass updated in calculateMortality (before this sub)
        avDBH = (AvStemMass / stemConst) ** (1.0D0 / stemPower)
        BasArea = (((avDBH / 200.) **2.) * Pi) * StemNo
        if(clm) then
             LAI  = ((SLA0*0.1D0)*(exp(tSLA*0.1D0*((WF1+WF2)))-1))/(tSLA*0.1D0)
             SLA = LAI / (WF1+WF2)
        else
             SLA = SLA1 + (SLA0 - SLA1) * Exp(-ln2 * (StandAge / tSLA)**2)
             LAI = (WF1+WF2) * SLA * 0.1D0
        endif
        !Hardwood LAI
        LAI_h = WF_h * SLA_h * 0.1D0 
        

		!-----------------------------------
		!-- UPDATE STATES AND FLUXES
		!----------------------------------- 

		!TIME VARIABLES
		OUTPUT(output_index,1) = calyear
		OUTPUT(output_index,2) = calmonth
		OUTPUT(output_index,3) = StandAge
		
		!STATE VARIABLES OR VARIABLES THAT HAVE OBSERVATIONS
		OUTPUT(output_index,4) = LAI
		OUTPUT(output_index,5) = WS
		OUTPUT(output_index,6) = WCR
		OUTPUT(output_index,7) = WR
		OUTPUT(output_index,8) = StemNo
		OUTPUT(output_index,9) = LAI_h
		OUTPUT(output_index,10) = WS_h
		OUTPUT(output_index,11) = WCR_h
		OUTPUT(output_index,12) = WR_h
		OUTPUT(output_index,13) = -99!StemNo_h
		OUTPUT(output_index,14) = ASW
		OUTPUT(output_index,15) = GPPdm + GPPdm_h
		OUTPUT(output_index,16) = -99 !NEE
		OUTPUT(output_index,17) = EvapTransp
		OUTPUT(output_index,18) = Transp_pine
		OUTPUT(output_index,19) = Transp_hard	
		OUTPUT(output_index,20) = delWR + delWR_h	
		OUTPUT(output_index,21) = delLitter + delLitter_h
		OUTPUT(output_index,22) = WF1 + WF2
		OUTPUT(output_index,23) = WF_h		
		OUTPUT(output_index,24) = LAI + LAI_h
		OUTPUT(output_index,25) = -99  !WBud
		OUTPUT(output_index,26) = WBud_h
		
		
		!OTHER VARIABLES
		OUTPUT(output_index,27) = delLitter
		OUTPUT(output_index,28) = delLitter_h		
		OUTPUT(output_index,29) = AvStemMass
		OUTPUT(output_index,30) = BasArea
		OUTPUT(output_index,31) = avDBH
		OUTPUT(output_index,32) = NPP
		OUTPUT(output_index,33) = NPP - NPP * pR 
		OUTPUT(output_index,34) = pF
		OUTPUT(output_index,35) = pS
		OUTPUT(output_index,36) = conductance
		OUTPUT(output_index,37) = Intcptn * Rain		
		OUTPUT(output_index,38) = MinASW
		OUTPUT(output_index,39) = MaxASW
		OUTPUT(output_index,40) = runoff		
		OUTPUT(output_index,41) = Rain	
		OUTPUT(output_index,42) = water_balance_error
		OUTPUT(output_index,43) = fNutr
		OUTPUT(output_index,44) = fT 
		OUTPUT(output_index,45) = fFrost
		OUTPUT(output_index,46) = genetics
		OUTPUT(output_index,47) = fCalpha
		OUTPUT(output_index,48) = fVPD
		OUTPUT(output_index,49) = fSW
		OUTPUT(output_index,50) = fAge
		OUTPUT(output_index,51) = fcomp
		OUTPUT(output_index,52) = fCg
		OUTPUT(output_index,53) = VPD
		OUTPUT(output_index,54) = wSmax	
		OUTPUT(output_index,55) = GPPc_h
		OUTPUT(output_index,56) = NPP_h
		OUTPUT(output_index,57) = NPP_h - NPP_h * pR_h 
		OUTPUT(output_index,58) = pF_h
		OUTPUT(output_index,59) = pS_h
		OUTPUT(output_index,60) = WF1 + WF2 + WS + WCR + WR
		OUTPUT(output_index,61) = SLA	
		OUTPUT(output_index,62) = NPP + GPPdm*(1-y)*(WR+WCR)/(WR+WCR+WF1+WF2+WS)
		OUTPUT(output_index,63) = delWR
		OUTPUT(output_index,64) = delWR_h
        OUTPUT(output_index,65) = APAR
        OUTPUT(output_index,66) = APARu
        
        if(output_index > 1) then
        	OUTPUT(output_index,67) = WS - OUTPUT(output_index-1,67)
        else
            OUTPUT(output_index,67) = WS - WSi
        endif


		!-----------------------------
		!--UPDATE TIME VARIABLES
		!------------------------------
        output_index = output_index + 1
        metMonth = metMonth + 1
        calmonth = calmonth + 1
                
        if (calmonth > 12) then
        	calmonth = 1
            calyear = calyear+1 
        end if
        
       	if (metMonth > nometmonths) then
       		metMonth = 1
        end if
                
		StandAge = StandAge + (1.0D0 / 12.0D0)

	END DO

END SUBROUTINE R3PG_MODEL

!-----------------------------------------------------
!  FUNCTIONS
!-----------------------------------------------------

!-----------------------------------------------------
!-  CALCULATE TRANSPIRATION
!-----------------------------------------------------
    
double precision FUNCTION getTranspiration(Q,VPD,h,gBL,gC,Qa,Qb)
	implicit none
    double precision, INTENT(IN) :: Q,VPD,h,gBL,gC,Qa,Qb
    !The following are constants in the PM formula (Landsberg & Gower, 1997)
    double precision, parameter :: e20 = 2.2D0            ! rate of change of saturated VP with T at 20C
    double precision, parameter :: rhoAir = 1.2D0         ! density of air, kg/m3
    double precision, parameter :: lambda = 2460000.0D0   ! latent heat of vapourisation of H2O (J/kg)
    double precision, parameter :: VPDconv = 0.000622D0   ! convert VPD to saturation deficit = 18/29/1000
    double precision :: netRad, defTerm, div, Etransp
    
    netRad = Qa + Qb * (Q * 10.0D0 ** 6.0D0 / h)        ! Q in MJ/m2/day --> W/m2
    defTerm = rhoAir * lambda * (VPDconv * VPD) * gBL
    div = (1.0D0 + e20 + gBL / gC)                        ! changed in 3PG+
    Etransp = (e20 * netRad + defTerm) / div               ! in J/m2/s; changed in 3PG+
    getTranspiration = Etransp / lambda * h             ! converted to kg/m2/day
    
End Function


!-----------------------------------------------------
!-  CALCULATE MORTALITY
!-----------------------------------------------------

double precision  FUNCTION getMortality(OldN,OldW,mS,wSx1000,thinPower)
	implicit none
    double precision , INTENT(IN) :: OldN,OldW,mS,wSx1000,thinPower
    double precision , PARAMETER :: accuracy = 1.0D0 / 1000.0D0
    INTEGER :: i
    double precision  :: fN,dfN,dN,n,x1,x2

    n = OldN / 1000.0D0
    x1 = 1000.0D0 * mS * OldW / OldN
    i = 0
    DO
        i = i + 1
        if (n < 0) EXIT !added in 3PG+
        x2 = wSx1000 * n ** (1.0D0 - thinPower)
        fN = x2 - x1 * n - (1.0D0 - mS) * OldW
        dfN = (1.0D0 - thinPower) * x2 / n - x1
        dN = -fN / dfN
        n = n + dN
        if (Abs(dN) <= accuracy .Or. i >= 5) EXIT
    END DO
    
    getMortality = OldN - 1000.0D0 * n
  
END FUNCTION getMortality

!-----------------------------------------------------
!-  CALCULATE VAPOR PRESSURE DEFICIT
!-----------------------------------------------------

double precision  FUNCTION getVPD(Tx, Tn)
	implicit none
	!Gets daily 'mean' VPD in mBar - based on daily max and min temperatures only
    double precision , INTENT(IN) :: Tx,Tn
    double precision :: VPDx, VPDn
    
	!THIS ASSUMES A CONSTANT RELATIVE HUMIDITY
	VPDx = 6.1078D0 * EXP(17.269D0 * Tx / (237.3D0 + Tx))
    VPDn = 6.1078D0 * EXP(17.269D0 * Tn / (237.3D0 + Tn))
    getVPD = (VPDx - VPDn) / 2.0D0

END FUNCTION

!-----------------------------------------------------
!-  CALCULATE ARC COS
!-----------------------------------------------------

double precision  FUNCTION Acos(x)
	implicit none
    double precision  :: x
    double precision, parameter :: Pi = 3.141592654D0
    
    Acos = ATAN(-x / SQRT(-x * x + 1.0D0)) + Pi / 2.0D0  !pjs note: arcos(x) = arctan((sqrt(1-x**2))/x)?

END FUNCTION Acos

!-----------------------------------------------------
!-  CALCULATE DAYLENGTH
!-----------------------------------------------------

!gets fraction of day when sun is "up"
double precision  function getDayLength(Lat, dayOfYear)
	implicit none
    double precision , INTENT(IN) :: Lat
    INTEGER, INTENT(IN) :: dayOfYear
    double precision  :: SLAt, cLat, sinDec, cosH0
    double precision, parameter :: Pi = 3.141592654D0

    SLAt = SIN(Pi * Lat / 180.0D0)
    cLat = COS(Pi * Lat / 180.0D0)
    sinDec = 0.4D0 * SIN(0.0172D0 * REAL(dayOfYear - 80,8))
    cosH0 = -sinDec * SLAt / (cLat * SQRT(1.0D0 - (sinDec) ** 2.0D0))
    if(cosH0 > 1) then
        getDayLength = 0.0D0
    elseif(cosH0 < -1) then
        getDayLength = 1.0D0
    else
        getDayLength = Acos(cosH0) / Pi
    end if
    

End Function      

end module R3PG_MODEL_MOD
