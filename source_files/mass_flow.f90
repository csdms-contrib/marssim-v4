      !     ####################################################################################################
    !     this is one of several source files for the marssim landform evolution model
    !     copyright (c) 2020 alan d. howard
    !     developer can be contacted by ah6p@virginia.edu or ahoward@psi.edu and department of environmental sciences, 
    !     p.o. box 400123, university of virginia, charlottesville, va 22904-4123
    !     this program is free software; you can redistribute it and/or modify it under the terms of the gnu general public license
    !       as published by the free software foundation; either version 3 of the
    !       license, or (at your option) any later version.
    !     this program is distributed in the hope that it will be useful, but without any warranty;
    !        without even the implied warranty of merchantability or fitness for a particular purpose. see the gnu
    !        general public license for more details.
    !      you should have received a copy of the gnu general public license along with this program; if not, write to
    !        the free software foundation, inc., 51 franklin street, fifth floor, boston, ma 02110-1301 usa.
    !        a web link:   http://www.gnu.org/licenses/gpl-3.0.txt
    !     ####################################################################################################
	!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
       ! glacier model additions 01-2014 to 04-2020 (Orkan Umurhan and Alan Howard)
	SUBROUTINE READ_MASS_FLOW_PARAMETERS()
		USE ERODE_GLOBALS
		USE MASS_FLOW_GLOBALS
        IMPLICIT NONE
		INTEGER :: WE_D0_MASS_FLOW,VARIABLE_MASS_FLOW_RATE_USE,PARAMS
        CHARACTER (80) :: TEXTLINE
! The mass flow subroutines route mass flows downstream. The routed material is considered to
!    be regolith, and all regolith is considered susceptible to flow unless there is a 
!    bingham_maximum_thickness greater than zero specified.  Two versions of mass flow are
!    considered, (1) Bingham-Law mass flow, which has a critical yield stress, parameterized in the
!    DELTA-45 parameter and (2) Glen's Law mass (glacial) flow, governed by the GLEN_LAW_ALPHA exponent.
!    The deformation rate for both versions is governed by the FLOW_DIFFUSIVITY parameter.  Flow rate
!    is also governed by the product of the sime of the topographic gradient times GRAVITY time the
!    FLOW_VOLUME_DENSITY (ice or wet regolith), as well as the depth of the regolith cover.
!    The means of calculating the flow depth(regolith thickness) depends upon the FLOW_SCHEME as detailed
!    below.  If the parameter VARIABLE_FLOW_VOLUME_RATE is set greater than one, then the time increment for
!    the simulation is determined to be governed by the maximum rate of surface elevation change due to the
!    flow according to the TARGET_MASS_FLOW_CHANGE_RATE parameter, subject to the simulation maximum and minimum
!    time increments.
! The mass flow may erode the "bedrock" base of the mass flow if the MASS_FLOW_EROSION_RATE is greater than
!    zero, subject to the MAXIMUM_FLOW_DEPTH_EROSION limit. ERoded bedrock becomes regolith. The rate of erosion is
!    specified in the Weathering module.  
! Glens flow law may also be governed by an upper layer of thickness CRITICAL_FLOW_THICKNESSS which does not
!    deform but is convected by the flow - that might be, for example, a debris cover or ice too cold to flow or 
!    is seasonal snow or low density ice.
        !READ  (PARAMS,22) TEXTLINE
        PARAMS=326
        OPEN(PARAMS,FILE='MASS_FLOW_PARAMETERS.PRM',ACTION='READ')
        read(params,22) textline
22      format(a)        
        READ (PARAMS,*) WE_D0_MASS_FLOW    !  (0 = none, 1=Bingham flow, any other integer > 0 is Glen Law)
        READ(PARAMS,*) GLEN_LAW_ALPHA  ! Exponent for Glen's Law glacier-like mass flow
		READ(PARAMS,*) FLOW_DIFFUSIVITY
		READ(PARAMS,*) DELTA_HB ! Units of meters - critical yield stress divided by (FLOW_VOLUME_density*gravity)
		READ(PARAMS,*) MASS_FLOW_EROSION_RATE
		READ(PARAMS,*) bingham_maximum_thickness 
		READ(PARAMS,*) FLOW_VOLUME_DENSITY   
        READ(PARAMS,*) FLUX_EXPONENT
		READ(PARAMS,*) DEPTH_EXPONENT
		READ(PARAMS,*) SINE_EXPONENT
		READ(PARAMS,*) MASS_FLOW_CRITICAL_VALUE
        READ(PARAMS,*) CRITICAL_FLOW_THICKNESSS
		READ(PARAMS,*) CRITICAL_SCOUR_THICKNESS
		READ(PARAMS,*) MAXIMUM_FLOW_DEPTH_EROSION
		READ(PARAMS,*) FLOW_SCHEME
        read(params,*) VARIABLE_MASS_FLOW_RATE_USE
		READ(PARAMS,*) TARGET_MASS_FLOW_CHANGE_RATE
        if (VARIABLE_MASS_FLOW_RATE_USE>0) then
            use_variable_FLOW_VOLUME_rate=.true.
        else
            use_variable_FLOW_VOLUME_rate=.false.
        endif
        MASS_FLOW_CRITICAL_VALUE=MASS_FLOW_CRITICAL_VALUE*FLOW_VOLUME_DENSITY*GRAVITY
        IF (WE_D0_MASS_FLOW> 0) THEN
            USE_FLOW_VOLUME = .TRUE.
            CUMULATIVE_FLOW_VOLUME_FLUX=0.0
            NUM_FLUX=0.0
            IF (WE_D0_MASS_FLOW == 1) THEN
                THIS_IS_BINGHAM_FLOW = .TRUE.
                ELSE
!                   DEFAULT SETTING IS TO USE GLEN LAW FORMULATION  !  added by (orkan) March 25
                THIS_IS_BINGHAM_FLOW  = .FALSE.
                ENDIF
            ELSE
                USE_FLOW_VOLUME = .FALSE.
            ENDIF
            SCHEME1 =.FALSE.
            SCHEME2 =.FALSE.
            SCHEME3 = .FALSE.
            IF (FLOW_SCHEME==1) SCHEME1=.TRUE.
            IF (FLOW_SCHEME==2) SCHEME2=.TRUE.
            IF (FLOW_SCHEME==3) SCHEME3=.TRUE.
        ! flow_scheme actions:
        !     if 0 or >3, use mixed scheme
        !     if 1, use local flow thickness
        !     if 2, use flow source thickness
        !     if 3, use average flow thickness
        ! end glacier model additions 02/26/2014 (orkan)
        !************* WRITE MASS FLOW PARAMETERS***************
        IF (USE_FLOW_VOLUME) THEN
            WRITE(OUTHIST,2201)
2201        FORMAT('**************** MASS FLOW PARAMETERS **************')
            WRITE (OUTHIST,2203) FLOW_DIFFUSIVITY,FLOW_VOLUME_DENSITY
2203        FORMAT('FLOW DIFFUSIVITY=',G12.5,' MASS FLOW DENSITY=',G12.5)
            WRITE(OUTHIST,2206) MASS_FLOW_EROSION_RATE, CRITICAL_SCOUR_THICKNESS, &
                MASS_FLOW_CRITICAL_VALUE
2206            FORMAT(' MASS FLOW EROSION RATE=',G12.5,' CRITICAL THICKNESS FOR SCOUR=' &
                ,G12.5,' MASS FLOW CRITICAL VALUE=',G12.5)
            IF (use_variable_FLOW_VOLUME_rate) THEN
                WRITE(OUTHIST,2208)TARGET_MASS_FLOW_CHANGE_RATE
2208            FORMAT(' TIME INTERVAL LIMITED TO MAXIMUM FLOW ELEVATION CHANGE=',G12.5)
            ENDIF
            IF (WE_D0_MASS_FLOW == 1) THEN
                WRITE(OUTHIST,2202)
2202            FORMAT('MODELING ASSUMES BINGHAM RHEOLOGY')
                WRITE (OUTHIST,2204) DELTA_HB,BINGHAM_MAXIMUM_THICKNESS 
2204            FORMAT(' CRITICAL VALUE FOR FLOW=',G12.5,' MAXIMUM FLOW THICKNESS=',G12.5)
                WRITE(OUTHIST,2207) fLOW_SCHEME
2207            FORMAT(' SCHEME FOR DETERMINING FLOW THICKNESS=',I6)
            ELSE
                WRITE(OUTHIST,2205)
2205            FORMAT('MODELING ASSUMES GLEN FLOW LAW')
                WRITE(OUTHIST,2209) GLEN_LAW_ALPHA,CRITICAL_FLOW_THICKNESSS, &
                    MAXIMUM_FLOW_DEPTH_EROSION
2209                FORMAT('GLEN LAW ALPHA=',G12.5,' CRITICAL THICKNESS TO ALLOW FLOW=',G12.5 &
                    ,' MAXIMUM DEPTH OF GLACIAL EROSION=',G12.5)
          ENDIF          
        ENDIF
        return
        end ! SUBROUTINE READ_MASS_FLOW_PARAMETERS
        !
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DO_BINGHAM_MASS_FLOW()
        USE ERODE_GLOBALS
		USE MASS_FLOW_GLOBALS
        IMPLICIT NONE
        !      ********************************************************************
        !        This subroutine models Bingham rheology mass flows and returns the net elevation
        !         change in ERODE_SLOPE(i,j)
        !   MODIFIES:  CFW, CFN, CFNE, CFNW, ERODE_SLOPE, IS_ROCK_SURFACE, NOMINAL_ERODE_SLOPE
        !              SLOPETOADD, SLGRAVTOADD
        !      ********************************************************************
        REAL (4) :: RAPID_CREEP
        EXTERNAL RAPID_CREEP
        INTEGER :: I,J,I_NEW,J_NEW,LOCAL_DIRECTION,I_OLD,J_OLD,IE,IW,IPATH,JS,JN
        integer :: ii_max, jj_max, ii_flux_max, jj_flux_max,i_regmax,j_regmax,i_regmin,j_regmin
        REAL (4) :: GRADIENT,SUMSLOPE
        REAL (4) :: REGEMIN,REGEMAX,REGEAVG,REGENUM
        REAL (4) :: ROCEMIN,ROCEMAX,ROCEAVG,ROCENUM
        real (4) :: xregmin,xregmax,xregaverage
        REAL (4) :: FLOW_THICKNESS,b_alpha,bsine,flow_param,xchange,shearnot,yeild_parameter
        real (4) :: smult,shearfac
		real (4) :: ti,ta,to,sheari,shearo,sheara,flowi,flowa,flowo
          !!!!!!!!write(outhist,4833)
4833      format('bingham called')
    
        !      ********************************************************************
        !       This is the main loop determining rates of Bingham mass flows
        !      *******************************************************************
	    MAX_FLOW_VOLUME_FLUX = -1.0
        MAX_FLOW_VOLUME_CHANGE = -1.0    !  added by (orkan) March 30
		FLOW_VOLUME_CHANGE=0.0
        FLOW_VOLUME_FLUX = 0.0
        AVG_FLOW_VOLUME_CHANGE = 0.0
        AVG_FLOW_VOLUME_FLUX = 0.0
        FRAC_BINGHAM_ACTIVE = 0.0
        FLOWTOADD=0.0
        FLGRAVTOADD=0.0
        XREGAVERAGE=0.0
        II_MAX=0
        JJ_MAX=0
        II_FLUX_MAX=0
        JJ_FLUX_MAX=0
        I_REGMIN=0
        I_REGMAX=0
        J_REGMIN=0
        J_REGMAX=0
        XREGMIN=1.0E25
        XREGMAX=-1.0E25
        ! b_alpha is the product of flow density and gravity
		b_alpha=FLOW_VOLUME_DENSITY*gravity
        ! yeild_parameter is the yield stress for the Bingham flow
		yeild_parameter=DELTA_HB*b_alpha
        ! The solution iterates through all cells and evaluates
        !  the mass flux to or from that cell to the adjacent
        !  cells to the north, northwest,northeast and west.
        !  Fluxes to the other four directions are calculated
        !  to and from the adjacent cells.
            DO  J = 1, MY
					JN=J-1
					IF (JN<1) THEN
						IF(IS_Y_PERIODIC) THEN
							JN=MY
						ELSE
							JN=2
						ENDIF
					ENDIF
					JS=J+1
					IF (JS>MY) THEN
						IF (IS_Y_PERIODIC) THEN
							JS=1
						ELSE
							JS=MY-1
						ENDIF
					ENDIF
                DO  I = 1, MX
					IW=I-1
					IF (IW<1) THEN
						IF(IS_X_PERIODIC) THEN
							IW=MX
						ELSE
							IW=2
						ENDIF
					ENDIF
					IE=I+1
					IF (IE>MX) THEN
						IF (IS_X_PERIODIC) THEN
							IE=1
						ELSE
							IE=MX-1
						ENDIF
					ENDIF
                    !      ********************************************************************
                    !       Determine driving topographic gradients.  Using sin(theta) rather than tan(theta)
                    !      ********************************************************************
                    GRADIENT=(ELEVATION(I,J)-ELEVATION(IW,J))/CELL_SIZE
					BSINE=SQRT(GRADIENT**2/(GRADIENT**2+1.0))
					CFW(I,J)=BSINE
					IF (GRADIENT<0.0) CFW(I,J)=-CFW(I,J)
	                GRADIENT=(ELEVATION(I,J)-ELEVATION(IE,JN))/(1.414*CELL_SIZE)    
					BSINE=SQRT(GRADIENT**2/(GRADIENT**2+1.0))
					CFNE(I,J)=BSINE
					IF (GRADIENT<0.0) CFNE(I,J)=-CFNE(I,J)
                    GRADIENT=(ELEVATION(I,J)-ELEVATION(I,JN))/CELL_SIZE
					BSINE=SQRT(GRADIENT**2/(GRADIENT**2+1.0))
                    CFN(I,J)=BSINE
					IF (GRADIENT<0.0) CFN(I,J)=-CFN(I,J)
                    GRADIENT=(ELEVATION(I,J)-ELEVATION(IW,JN))/(1.414*CELL_SIZE) 
					BSINE=SQRT(GRADIENT**2/(GRADIENT**2+1.0))
					CFNW(I,J)=BSINE
					IF (GRADIENT<0.0) CFNW(I,J)=-CFNW(I,J)
				ENDDO
			ENDDO
			DO  J = 1, MY
					JN=J-1
					IF (JN<1) THEN
						IF(IS_Y_PERIODIC) THEN
							JN=MY
						ELSE
							JN=2
						ENDIF
					ENDIF
					JS=J+1
					IF (JS>MY) THEN
						IF (IS_Y_PERIODIC) THEN
							JS=1
						ELSE
							JS=MY-1
						ENDIF
					ENDIF
                DO  I = 1, MX
					IW=I-1
					IF (IW<1) THEN
						IF(IS_X_PERIODIC) THEN
							IW=MX
						ELSE
							IW=2
						ENDIF
					ENDIF
					IE=I+1
					IF (IE>MX) THEN
						IF (IS_X_PERIODIC) THEN
							IE=1
						ELSE
							IE=MX-1
						ENDIF
                    ENDIF
                    ! to,ti,and ta are the local, adjacent, and average flow thicknesses
                    ! Which thickness is used depends upon the scheme used
                    !    if scheme=1 the local regolith thickness is used
                    !    if scheme=2 the regolith thickness is the location where the flow
                    !       is coming from which may be the local or adjacent cell
                    !    if scheme=3 the regolith is the average regolith thickness of the
                    !       local and adjacent cells
                    !    If scheme=0 or >3 then a mixed scheme is used.  
                    !       If the flow is coming from the adjacent cell
                    !         and the local cells exceeds the critical yield stress then
                    !         the average regolith thickness is used, otherwise the 
                    !         remote xell thickness is used
                    !       if the flow is going to the adjacent cell and the
                    !         local cell exceeds the critical yield stress then
                    !         the average thickness is used otherwise the local
                    !         thickness is used
					TO=MAX1(0.0,REGOLITH(I,J))
					TI=MAX1(0.0,REGOLITH(I,JN))
					TA=(TO+TI)/2.0
					IF (BINGHAM_MAXIMUM_THICKNESS>0.0) THEN
                        TO=MIN1(TO,BINGHAM_MAXIMUM_THICKNESS)
						TI=MIN1(TI,BINGHAM_MAXIMUM_THICKNESS)
						TA=MIN1(TA,BINGHAM_MAXIMUM_THICKNESS)
                    ENDIF
					
					SHEARFAC=ABS(CFN(I,J))*B_ALPHA
					SHEARI=SHEARFAC*TI
					SHEARO=SHEARFAC*TO
					SHEARA=SHEARFAC*TA
					FLOWI=SHEARI-yeild_parameter
					FLOWO=SHEARO-yeild_parameter
					FLOWA=SHEARA-yeild_parameter
					IF (CFN(I,J)>=0.0) THEN
					    IF (SCHEME1) THEN
					       FLOW_PARAM=FLOWO
					    ELSE
					       IF (SCHEME2) THEN
						       FLOW_PARAM=FLOWO
						   ELSE
							   IF (SCHEME3) THEN
							        FLOW_PARAM=FLOWA
							   ELSE
						           IF (FLOWO>0.0) THEN
							           FLOW_PARAM=FLOWA
						           ELSE
							           FLOW_PARAM=FLOWO
						           ENDIF
							   ENDIF
						   ENDIF
						ENDIF
					ELSE
                        IF (SCHEME1) THEN
					       FLOW_PARAM=FLOWO
					    ELSE
					       IF (SCHEME2) THEN
						       FLOW_PARAM=FLOWI
						   ELSE
							   IF (SCHEME3) THEN
							        FLOW_PARAM=FLOWA
							   ELSE
						           IF (FLOWO>0.0) THEN
							           FLOW_PARAM=FLOWA
						           ELSE
							           FLOW_PARAM=FLOWI
						           ENDIF
							   ENDIF
						   ENDIF
						ENDIF						
					ENDIF
					IF (FLOW_PARAM>0.0) THEN
                        FRAC_BINGHAM_ACTIVE=FRAC_BINGHAM_ACTIVE+0.5
					XCHANGE=CROSS_WEIGHTING*FLOW_DIFFUSIVITY*(FLOW_PARAM**3/3.0 + &
					   yeild_parameter*FLOW_PARAM**2/2.0)/SHEARFAC**2
                    SMULT=SIGN(1.0,CFN(I,J))
					FLOW_VOLUME_CHANGE(I,J)=FLOW_VOLUME_CHANGE(I,J)-SMULT*XCHANGE 
					FLOW_VOLUME_CHANGE(I,JN)=FLOW_VOLUME_CHANGE(I,JN)+SMULT*XCHANGE
					FLOW_VOLUME_FLUX(I,J)=FLOW_VOLUME_FLUX(I,J)+XCHANGE
					FLOW_VOLUME_FLUX(I,JN)=FLOW_VOLUME_FLUX(I,JN)+XCHANGE
                    FLGRAVTOADD=FLGRAVTOADD+ABS(XCHANGE*(ELEVATION(I,J)-ELEVATION(I,JN)))
                    FLOWTOADD=FLOWTOADD+ABS(XCHANGE)
					ENDIF
					TI=MAX1(0.0,REGOLITH(IW,J))
					TA=(TO+TI)/2.0
					IF (BINGHAM_MAXIMUM_THICKNESS>0.0) THEN
						TI=MIN1(TI,BINGHAM_MAXIMUM_THICKNESS)
						TA=MIN1(TA,BINGHAM_MAXIMUM_THICKNESS)
                    ENDIF					
					SHEARFAC=ABS(CFW(I,J))*B_ALPHA
					SHEARI=SHEARFAC*TI
					SHEARO=SHEARFAC*TO
					SHEARA=SHEARFAC*TA
					FLOWI=SHEARI-yeild_parameter
					FLOWO=SHEARO-yeild_parameter
					FLOWA=SHEARA-yeild_parameter
					IF (CFW(I,J)>=0.0) THEN
					    IF (SCHEME1) THEN
					       FLOW_PARAM=FLOWO
					    ELSE
					       IF (SCHEME2) THEN
						       FLOW_PARAM=FLOWO
						   ELSE
							   IF (SCHEME3) THEN
							        FLOW_PARAM=FLOWA
							   ELSE
						           IF (FLOWO>0.0) THEN
							           FLOW_PARAM=FLOWA
						           ELSE
							           FLOW_PARAM=FLOWO
						           ENDIF
							   ENDIF
						   ENDIF
						ENDIF
					ELSE
                        IF (SCHEME1) THEN
					       FLOW_PARAM=FLOWO
					    ELSE
					       IF (SCHEME2) THEN
						       FLOW_PARAM=FLOWI
						   ELSE
							   IF (SCHEME3) THEN
							        FLOW_PARAM=FLOWA
							   ELSE
						           IF (FLOWO>0.0) THEN
							           FLOW_PARAM=FLOWA
						           ELSE
							           FLOW_PARAM=FLOWI
						           ENDIF
							   ENDIF
						   ENDIF
						ENDIF						
                    ENDIF
					! Here the flow flux is calculated using the relationship for
                    !   depth-integrated Bingham flow velocity
					IF (FLOW_PARAM>0.0) THEN
                         FRAC_BINGHAM_ACTIVE=FRAC_BINGHAM_ACTIVE+0.5
					XCHANGE=CROSS_WEIGHTING*FLOW_DIFFUSIVITY*(FLOW_PARAM**3/3.0 + &
					   yeild_parameter*FLOW_PARAM**2/2.0)/SHEARFAC**2                         
                    SMULT=SIGN(1.0,CFW(I,J))
					FLOW_VOLUME_CHANGE(I,J)=FLOW_VOLUME_CHANGE(I,J)-SMULT*XCHANGE 
					FLOW_VOLUME_CHANGE(IW,J)=FLOW_VOLUME_CHANGE(IW,J)+SMULT*XCHANGE
					FLOW_VOLUME_FLUX(I,J)=FLOW_VOLUME_FLUX(I,J)+XCHANGE
					FLOW_VOLUME_FLUX(IW,J)=FLOW_VOLUME_FLUX(IW,J)+XCHANGE
                    FLGRAVTOADD=FLGRAVTOADD+ABS(XCHANGE*(ELEVATION(I,J)-ELEVATION(IW,J)))
                    FLOWTOADD=FLOWTOADD+ABS(XCHANGE)
					ENDIF
					TI=MAX1(0.0,REGOLITH(IE,JN))
					TA=(TO+TI)/2.0
					IF (BINGHAM_MAXIMUM_THICKNESS>0.0) THEN
						TI=MIN1(TI,BINGHAM_MAXIMUM_THICKNESS)
						TA=MIN1(TA,BINGHAM_MAXIMUM_THICKNESS)
                    ENDIF					
					SHEARFAC=ABS(CFNE(I,J))*B_ALPHA
					SHEARI=SHEARFAC*TI
					SHEARO=SHEARFAC*TO
					SHEARA=SHEARFAC*TA
					FLOWI=SHEARI-yeild_parameter
					FLOWO=SHEARO-yeild_parameter
					FLOWA=SHEARA-yeild_parameter
						IF (CFNE(I,J)>=0.0) THEN
					    IF (SCHEME1) THEN
					       FLOW_PARAM=FLOWO
					    ELSE
					       IF (SCHEME2) THEN
						       FLOW_PARAM=FLOWO
						   ELSE
							   IF (SCHEME3) THEN
							        FLOW_PARAM=FLOWA
							   ELSE
						           IF (FLOWO>0.0) THEN
							           FLOW_PARAM=FLOWA
						           ELSE
							           FLOW_PARAM=FLOWO
						           ENDIF
							   ENDIF
						   ENDIF
						ENDIF
					ELSE
                        IF (SCHEME1) THEN
					       FLOW_PARAM=FLOWO
					    ELSE
					       IF (SCHEME2) THEN
						       FLOW_PARAM=FLOWI
						   ELSE
							   IF (SCHEME3) THEN
							        FLOW_PARAM=FLOWA
							   ELSE
						           IF (FLOWO>0.0) THEN
							           FLOW_PARAM=FLOWA
						           ELSE
							           FLOW_PARAM=FLOWI
						           ENDIF
							   ENDIF
						   ENDIF
						ENDIF						
					ENDIF								

					IF (FLOW_PARAM>0.0) THEN
                         FRAC_BINGHAM_ACTIVE=FRAC_BINGHAM_ACTIVE+0.5
 					XCHANGE=DIAGONAL_WEIGHTING*FLOW_DIFFUSIVITY*(FLOW_PARAM**3/3.0 + &
					   yeild_parameter*FLOW_PARAM**2/2.0)/SHEARFAC**2                        
                    SMULT=SIGN(1.0,CFNE(I,J))
					FLOW_VOLUME_CHANGE(I,J)=FLOW_VOLUME_CHANGE(I,J)-SMULT*XCHANGE 
					FLOW_VOLUME_CHANGE(IE,JN)=FLOW_VOLUME_CHANGE(IE,JN)+SMULT*XCHANGE
					FLOW_VOLUME_FLUX(I,J)=FLOW_VOLUME_FLUX(I,J)+XCHANGE
					FLOW_VOLUME_FLUX(IE,JN)=FLOW_VOLUME_FLUX(IE,JN)+XCHANGE
                    FLGRAVTOADD=FLGRAVTOADD+ABS(XCHANGE*(ELEVATION(I,J)-ELEVATION(IE,JN)))
                    FLOWTOADD=FLOWTOADD+ABS(XCHANGE)
					ENDIF
					TI=MAX1(0.0,REGOLITH(IW,JN))
					TA=(TO+TI)/2.0
					IF (BINGHAM_MAXIMUM_THICKNESS>0.0) THEN
						TI=MIN1(TI,BINGHAM_MAXIMUM_THICKNESS)
						TA=MIN1(TA,BINGHAM_MAXIMUM_THICKNESS)
                    ENDIF					
					SHEARFAC=ABS(CFNW(I,J))*B_ALPHA
					SHEARI=SHEARFAC*TI
					SHEARO=SHEARFAC*TO
					SHEARA=SHEARFAC*TA
					FLOWI=SHEARI-yeild_parameter
					FLOWO=SHEARO-yeild_parameter
					FLOWA=SHEARA-yeild_parameter
					IF (CFNW(I,J)>=0.0) THEN
					    IF (SCHEME1) THEN
					       FLOW_PARAM=FLOWO
					    ELSE
					       IF (SCHEME2) THEN
						       FLOW_PARAM=FLOWO
						   ELSE
							   IF (SCHEME3) THEN
							        FLOW_PARAM=FLOWA
							   ELSE
						           IF (FLOWO>0.0) THEN
							           FLOW_PARAM=FLOWA
						           ELSE
							           FLOW_PARAM=FLOWO
						           ENDIF
							   ENDIF
						   ENDIF
						ENDIF
					ELSE
                        IF (SCHEME1) THEN
					       FLOW_PARAM=FLOWO
					    ELSE
					       IF (SCHEME2) THEN
						       FLOW_PARAM=FLOWI
						   ELSE
							   IF (SCHEME3) THEN
							        FLOW_PARAM=FLOWA
							   ELSE
						           IF (FLOWO>0.0) THEN
							           FLOW_PARAM=FLOWA
						           ELSE
							           FLOW_PARAM=FLOWI
						           ENDIF
							   ENDIF
						   ENDIF
						ENDIF						
					ENDIF					
					IF (FLOW_PARAM>0.0) THEN
                         FRAC_BINGHAM_ACTIVE=FRAC_BINGHAM_ACTIVE+0.5
  					XCHANGE=DIAGONAL_WEIGHTING*FLOW_DIFFUSIVITY*(FLOW_PARAM**3/3.0 + &
					   yeild_parameter*FLOW_PARAM**2/2.0)/SHEARFAC**2                       
                    SMULT=SIGN(1.0,CFNW(I,J))
					FLOW_VOLUME_CHANGE(I,J)=FLOW_VOLUME_CHANGE(I,J)-SMULT*XCHANGE 
					FLOW_VOLUME_CHANGE(IW,JN)=FLOW_VOLUME_CHANGE(IW,JN)+SMULT*XCHANGE
					FLOW_VOLUME_FLUX(I,J)=FLOW_VOLUME_FLUX(I,J)+XCHANGE
					FLOW_VOLUME_FLUX(IW,JN)=FLOW_VOLUME_FLUX(IW,JN)+XCHANGE
                    FLGRAVTOADD=FLGRAVTOADD+ABS(XCHANGE*(ELEVATION(I,J)-ELEVATION(IW,JN)))
                    FLOWTOADD=FLOWTOADD+ABS(XCHANGE)
					ENDIF				
				ENDDO
            ENDDO
            ! This calculates statistics about the fluxes and elevation changes
            DO J=1,MY
				DO I=1,MX
                    XREGAVERAGE=XREGAVERAGE+REGOLITH(I,J)
                    IF(REGOLITH(I,J)>XREGMAX) THEN
                        XREGMAX=REGOLITH(I,J)
                        I_REGMAX=I
                        J_REGMAX=J
                    ENDIF
                    IF (REGOLITH(I,J)<XREGMIN) THEN
                        XREGMIN=REGOLITH(I,J)
                        I_REGMIN=I
                        J_REGMIN=J
                    ENDIF
                AVG_FLOW_VOLUME_FLUX = AVG_FLOW_VOLUME_FLUX + ABS(FLOW_VOLUME_FLUX(I,J))
                AVG_FLOW_VOLUME_CHANGE = AVG_FLOW_VOLUME_FLUX + ABS(FLOW_VOLUME_CHANGE(I,J))          
                CUMULATIVE_FLOW_VOLUME_FLUX(I,J)=CUMULATIVE_FLOW_VOLUME_FLUX(I,J)+FLOW_VOLUME_FLUX(I,J)
                NUM_FLUX=NUM_FLUX+1.0           
                IF(ABS(FLOW_VOLUME_CHANGE(I,J))>MAX_FLOW_VOLUME_CHANGE)  THEN
                    MAX_FLOW_VOLUME_CHANGE = ABS(FLOW_VOLUME_CHANGE(I,J))
                    II_MAX = I
                    JJ_MAX = J
                    ENDIF
                IF(ABS(FLOW_VOLUME_FLUX(I,J))>MAX_FLOW_VOLUME_FLUX)  THEN
                    MAX_FLOW_VOLUME_FLUX = ABS(FLOW_VOLUME_FLUX(I,J))
                    II_FLUX_MAX = I
                    JJ_FLUX_MAX = J
                    ENDIF    
				ENDDO
            ENDDO
            xregaverage=xregaverage/(mx*my)
            write(*,100) II_MAX,JJ_MAX,MAX_FLOW_VOLUME_CHANGE,II_FLUX_MAX,JJ_FLUX_MAX,MAX_FLOW_VOLUME_FLUX
100         FORMAT('I=',I5,' J=',I5,'MAX_CHANGE=',G13.6,/,'I=',I5,' J=',I5,' MAX_FLUX=',G13.6)
            write(*,200) i_regmin,j_regmin,xregmin,i_regmax,j_regmax,xregmax,xregaverage
200 format('i=',i5,' j=',i5,' regmin=',g13.6,/,'i=',i5,' j=',i5,' regmax=',g13.6,' regaverage=',g13.6)
			return
	end ! SUBROUTINE DO_BINGHAM_MASS_FLOW
	!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DO_GLEN_LAW_FLOW()
        USE ERODE_GLOBALS
		USE MASS_FLOW_GLOBALS
        IMPLICIT NONE
        !      ********************************************************************
        !        This subroutine models slope mass wasting and returns the net elevation
        !         change in ERODE_SLOPE(i,j)
        !   MODIFIES:  CFW, CFN, CFNE, CFNW, ERODE_SLOPE, IS_ROCK_SURFACE, NOMINAL_ERODE_SLOPE
        !              SLOPETOADD, SLGRAVTOADD
        !      ********************************************************************
        REAL (4) :: RAPID_CREEP
        EXTERNAL RAPID_CREEP
        INTEGER :: I,J,I_NEW,J_NEW,LOCAL_DIRECTION,I_OLD,J_OLD,IE,IW,IPATH,JS,JN
        integer :: ii_max, jj_max, ii_flux_max, jj_flux_max,i_regmax,j_regmax,i_regmin,j_regmin
        REAL (4) :: GRADIENT,SUMSLOPE
        REAL (4) :: REGEMIN,REGEMAX,REGEAVG,REGENUM
        REAL (4) :: ROCEMIN,ROCEMAX,ROCEAVG,ROCENUM
        real (4) :: xregmin,xregmax,term1,term2,term3,term4,term5,alpha1,alpha2,xthick,alphaterm
        REAL (4) :: FLOW_THICKNESS,b_alpha,bsine,flow_param,xchange,shearnot,yeild_parameter
        real (4) :: smult,shearfac
          !!!!!!!!write(outhist,4833)
4833    format('glen called')
        !      ********************************************************************
        !       This is the main loop determining rates of Glaen law mass flow
        !      ********************************************************************
	    MAX_FLOW_VOLUME_FLUX = -1.0
        MAX_FLOW_VOLUME_CHANGE = -1.0    !  added by (orkan) March 30
		FLOW_VOLUME_CHANGE=0.0
        FLOW_VOLUME_FLUX = 0.0
        AVG_FLOW_VOLUME_CHANGE = 0.0
        AVG_FLOW_VOLUME_FLUX = 0.0
        FRAC_BINGHAM_ACTIVE = 0.0
        II_MAX=0
        JJ_MAX=0
        II_FLUX_MAX=0
        JJ_FLUX_MAX=0
        I_REGMIN=0
        I_REGMAX=0
        J_REGMIN=0
        J_REGMAX=0
        XREGMIN=1.0E25
        XREGMAX=-1.0E25
        ! This code is directly comparable to the Bingham mass flow, only it uses the Glen's Law
        !   rheology.  Because all cells with regolith (ice) have flow, the flow depth is
        !   calculated using the average of the regolith thickness of the local and adjacent cells.
		B_ALPHA=FLOW_DIFFUSIVITY*(FLOW_VOLUME_DENSITY*GRAVITY)**GLEN_LAW_ALPHA/((GLEN_LAW_ALPHA+1.0)) !*(GLEN_LAW_ALPHA+2.0))
		ALPHA1=GLEN_LAW_ALPHA+1.0
		ALPHA2=GLEN_LAW_ALPHA+2.0
        ALPHATERM=(1.0-1.0/ALPHA2)
            DO  J = 1, MY
					JN=J-1
					IF (JN<1) THEN
						IF(IS_Y_PERIODIC) THEN
							JN=MY
						ELSE
							JN=2
						ENDIF
					ENDIF
					JS=J+1
					IF (JS>MY) THEN
						IF (IS_Y_PERIODIC) THEN
							JS=1
						ELSE
							JS=MY-1
						ENDIF
					ENDIF
                DO  I = 1, MX
					IW=I-1
					IF (IW<1) THEN
						IF(IS_X_PERIODIC) THEN
							IW=MX
						ELSE
							IW=2
						ENDIF
					ENDIF
					IE=I+1
					IF (IE>MX) THEN
						IF (IS_X_PERIODIC) THEN
							IE=1
						ELSE
							IE=MX-1
						ENDIF
					ENDIF
                    !      ********************************************************************
                    !       Determine driving topographic gradients.  Using sin(theta) rather than tan(theta)
                    !      ********************************************************************
                    GRADIENT=(ELEVATION(I,J)-ELEVATION(IW,J))/CELL_SIZE
					BSINE=SQRT(GRADIENT**2/(GRADIENT**2+1.0))
					CFW(I,J)=BSINE
					IF (GRADIENT<0.0) CFW(I,J)=-CFW(I,J)
	                GRADIENT=(ELEVATION(I,J)-ELEVATION(IE,JN))/(1.414*CELL_SIZE)
					BSINE=SQRT(GRADIENT**2/(GRADIENT**2+1.0))
					CFNE(I,J)=BSINE
					IF (GRADIENT<0.0) CFNE(I,J)=-CFNE(I,J)
                    GRADIENT=(ELEVATION(I,J)-ELEVATION(I,JN))/CELL_SIZE
					BSINE=SQRT(GRADIENT**2/(GRADIENT**2+1.0))
                    CFN(I,J)=BSINE
					IF (GRADIENT<0.0) CFN(I,J)=-CFN(I,J)
                    GRADIENT=(ELEVATION(I,J)-ELEVATION(IW,JN))/(1.414*CELL_SIZE)
					BSINE=SQRT(GRADIENT**2/(GRADIENT**2+1.0))
					CFNW(I,J)=BSINE
					IF (GRADIENT<0.0) CFNW(I,J)=-CFNW(I,J)
				ENDDO
			ENDDO
			DO  J = 1, MY
					JN=J-1
					IF (JN<1) THEN
						IF(IS_Y_PERIODIC) THEN
							JN=MY
						ELSE
							JN=2
						ENDIF
					ENDIF
					JS=J+1
					IF (JS>MY) THEN
						IF (IS_Y_PERIODIC) THEN
							JS=1
						ELSE
							JS=MY-1
						ENDIF
					ENDIF
                DO  I = 1, MX
					IW=I-1
					IF (IW<1) THEN
						IF(IS_X_PERIODIC) THEN
							IW=MX
						ELSE
							IW=2
						ENDIF
					ENDIF
					IE=I+1
					IF (IE>MX) THEN
						IF (IS_X_PERIODIC) THEN
							IE=1
						ELSE
							IE=MX-1
						ENDIF
					ENDIF
					FLOW_THICKNESS=(MAX1(0.0,REGOLITH(I,J))+MAX1(0.0,REGOLITH(I,JN)))/2.0
	                IF (MAXIMUM_FLOW_DEPTH_EROSION>0.0) THEN
                        FLOW_THICKNESS=MIN1(FLOW_THICKNESS,MAXIMUM_FLOW_DEPTH_EROSION)
                    ENDIF
					XTHICK=FLOW_THICKNESS-CRITICAL_FLOW_THICKNESSS
					IF (XTHICK>0.0) THEN
						SHEARFAC=CROSS_WEIGHTING*B_ALPHA*ABS(CFN(I,J))**GLEN_LAW_ALPHA
						IF (CRITICAL_FLOW_THICKNESSS>0.0) THEN
							XTHICK=FLOW_THICKNESS-CRITICAL_FLOW_THICKNESSS
                            XCHANGE=SHEARFAC*ALPHATERM*(FLOW_THICKNESS**ALPHA2-XTHICK**ALPHA2)
						ELSE
						    XCHANGE=SHEARFAC*(FLOW_THICKNESS**ALPHA2)*ALPHATERM
						ENDIF
						SMULT=SIGN(1.0,CFN(I,J))
						FLOW_VOLUME_CHANGE(I,J)=FLOW_VOLUME_CHANGE(I,J)-SMULT*XCHANGE 
						FLOW_VOLUME_CHANGE(I,JN)=FLOW_VOLUME_CHANGE(I,JN)+SMULT*XCHANGE
						FLOW_VOLUME_FLUX(I,J)=FLOW_VOLUME_FLUX(I,J)+XCHANGE
						FLOW_VOLUME_FLUX(I,JN)=FLOW_VOLUME_FLUX(I,JN)+XCHANGE
					ENDIF
					   
					FLOW_THICKNESS=(MAX1(0.0,REGOLITH(I,J))+MAX1(0.0,REGOLITH(IW,J)))/2.0
	                IF (MAXIMUM_FLOW_DEPTH_EROSION>0.0) THEN
                        FLOW_THICKNESS=MIN1(FLOW_THICKNESS,MAXIMUM_FLOW_DEPTH_EROSION)
                    ENDIF
					XTHICK=FLOW_THICKNESS-CRITICAL_FLOW_THICKNESSS
					IF (XTHICK>0.0) THEN
						SHEARFAC=CROSS_WEIGHTING*B_ALPHA*ABS(CFW(I,J))**GLEN_LAW_ALPHA
						IF (CRITICAL_FLOW_THICKNESSS>0.0) THEN
							XTHICK=FLOW_THICKNESS-CRITICAL_FLOW_THICKNESSS
                            XCHANGE=SHEARFAC*ALPHATERM*(FLOW_THICKNESS**ALPHA2-XTHICK**ALPHA2)
						ELSE
						    XCHANGE=SHEARFAC*(FLOW_THICKNESS**ALPHA2)*ALPHATERM
						ENDIF
						SMULT=SIGN(1.0,CFW(I,J))
						FLOW_VOLUME_CHANGE(I,J)=FLOW_VOLUME_CHANGE(I,J)-SMULT*XCHANGE 
						FLOW_VOLUME_CHANGE(IW,J)=FLOW_VOLUME_CHANGE(IW,J)+SMULT*XCHANGE
						FLOW_VOLUME_FLUX(I,J)=FLOW_VOLUME_FLUX(I,J)+XCHANGE
						FLOW_VOLUME_FLUX(IW,J)=FLOW_VOLUME_FLUX(IW,J)+XCHANGE
					ENDIF					
                    FLOW_THICKNESS=(MAX1(0.0,REGOLITH(I,J))+MAX1(0.0,REGOLITH(IE,JN)))/2.0               
	                IF (MAXIMUM_FLOW_DEPTH_EROSION>0.0) THEN
                        FLOW_THICKNESS=MIN1(FLOW_THICKNESS,MAXIMUM_FLOW_DEPTH_EROSION)
                    ENDIF
					XTHICK=FLOW_THICKNESS-CRITICAL_FLOW_THICKNESSS
					IF (XTHICK>0.0) THEN
						SHEARFAC=DIAGONAL_WEIGHTING*B_ALPHA*ABS(CFNE(I,J))**GLEN_LAW_ALPHA
                       IF (CRITICAL_FLOW_THICKNESSS>0.0) THEN
							XTHICK=FLOW_THICKNESS-CRITICAL_FLOW_THICKNESSS
                            XCHANGE=SHEARFAC*ALPHATERM*(FLOW_THICKNESS**ALPHA2-XTHICK**ALPHA2)
						ELSE
						    XCHANGE=SHEARFAC*(FLOW_THICKNESS**ALPHA2)*ALPHATERM
						ENDIF 
						SMULT=SIGN(1.0,CFNE(I,J))
						FLOW_VOLUME_CHANGE(I,J)=FLOW_VOLUME_CHANGE(I,J)-SMULT*XCHANGE 
						FLOW_VOLUME_CHANGE(IE,JN)=FLOW_VOLUME_CHANGE(IE,JN)+SMULT*XCHANGE
						FLOW_VOLUME_FLUX(I,J)=FLOW_VOLUME_FLUX(I,J)+XCHANGE
						FLOW_VOLUME_FLUX(IE,JN)=FLOW_VOLUME_FLUX(IE,JN)+XCHANGE
					ENDIF					
  					FLOW_THICKNESS=(MAX1(0.0,REGOLITH(I,J))+MAX1(0.0,REGOLITH(IW,JN)))/2.0                  
	                IF (MAXIMUM_FLOW_DEPTH_EROSION>0.0) THEN
                        FLOW_THICKNESS=MIN1(FLOW_THICKNESS,MAXIMUM_FLOW_DEPTH_EROSION)
                    ENDIF
					XTHICK=FLOW_THICKNESS-CRITICAL_FLOW_THICKNESSS
					IF (XTHICK>0.0) THEN
						SHEARFAC=DIAGONAL_WEIGHTING*B_ALPHA*ABS(CFNW(I,J))**GLEN_LAW_ALPHA
                        IF (CRITICAL_FLOW_THICKNESSS>0.0) THEN
							XTHICK=FLOW_THICKNESS-CRITICAL_FLOW_THICKNESSS
                            XCHANGE=SHEARFAC*ALPHATERM*(FLOW_THICKNESS**ALPHA2-XTHICK**ALPHA2)
						ELSE
						    XCHANGE=SHEARFAC*(FLOW_THICKNESS**ALPHA2)*ALPHATERM
						ENDIF
						SMULT=SIGN(1.0,CFNW(I,J))
						FLOW_VOLUME_CHANGE(I,J)=FLOW_VOLUME_CHANGE(I,J)-SMULT*XCHANGE 
						FLOW_VOLUME_CHANGE(IW,JN)=FLOW_VOLUME_CHANGE(IW,JN)+SMULT*XCHANGE
						FLOW_VOLUME_FLUX(I,J)=FLOW_VOLUME_FLUX(I,J)+XCHANGE
						FLOW_VOLUME_FLUX(IW,JN)=FLOW_VOLUME_FLUX(IW,JN)+XCHANGE
					ENDIF						
				ENDDO
            ENDDO
            DO J=1,MY
				DO I=1,MX
                    IF(REGOLITH(I,J)>XREGMAX) THEN
                        XREGMAX=REGOLITH(I,J)
                        I_REGMAX=I
                        J_REGMAX=J
                    ENDIF
                    IF (REGOLITH(I,J)<XREGMIN) THEN
                        XREGMIN=REGOLITH(I,J)
                        I_REGMIN=I
                        J_REGMIN=J
                    ENDIF
                AVG_FLOW_VOLUME_FLUX = AVG_FLOW_VOLUME_FLUX + ABS(FLOW_VOLUME_FLUX(I,J))
                AVG_FLOW_VOLUME_CHANGE = AVG_FLOW_VOLUME_FLUX + ABS(FLOW_VOLUME_CHANGE(I,J))
                CUMULATIVE_FLOW_VOLUME_FLUX(I,J)=CUMULATIVE_FLOW_VOLUME_FLUX(I,J)+FLOW_VOLUME_FLUX(I,J)
                NUM_FLUX=NUM_FLUX+1.0           
                IF(ABS(FLOW_VOLUME_CHANGE(I,J))>MAX_FLOW_VOLUME_CHANGE)  THEN
                    MAX_FLOW_VOLUME_CHANGE = ABS(FLOW_VOLUME_CHANGE(I,J))
                    II_MAX = I
                    JJ_MAX = J
                    ENDIF
                IF(ABS(FLOW_VOLUME_FLUX(I,J))>MAX_FLOW_VOLUME_FLUX)  THEN
                    MAX_FLOW_VOLUME_FLUX = ABS(FLOW_VOLUME_FLUX(I,J))
                    II_FLUX_MAX = I
                    JJ_FLUX_MAX = J
                    ENDIF    
				ENDDO
            ENDDO
            !write(*,100) II_MAX,JJ_MAX,MAX_FLOW_VOLUME_CHANGE,II_FLUX_MAX,JJ_FLUX_MAX,MAX_FLOW_VOLUME_FLUX
100         FORMAT('I=',I5,' J=',I5,'MAX_CHANGE=',G13.6,/,'I=',I5,' J=',I5,' MAX_FLUX=',G13.6)
            !write(*,200) i_regmin,j_regmin,xregmin,i_regmax,j_regmax,xregmax
200 format('i=',i5,' j=',i5,' regmin=',g13.6,/,'i=',i5,' j=',i5,' regmax=',g13.6)
			RETURN
	END ! SUBROUTINE DO_GLEN_LAW_FLOW
					
					