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
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	SUBROUTINE READ_MASS_WASTING_PARAMETERS
	USE ERODE_GLOBALS
	IMPLICIT NONE
	CHARACTER (80) :: TEXTLINE
	INTEGER :: MASS_WASTING_PARAMS, ISLOPEUSE, GRITICAL_GRADIENT_USE, ROERINGUSE 
	INTEGER :: IDO_DEPTHCREEP, DO_ROUTE 
	MASS_WASTING_PARAMS=93
        !        *******************************************************************
        !        *****************  mass wasting parameters  ************************
        !        *******************************************************************
		OPEN(MASS_WASTING_PARAMS,FILE='MASS_WASTING_PARAMETERS.PRM',ACTION='READ')
        READ(MASS_WASTING_PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
22      FORMAT(A)        
        READ(MASS_WASTING_PARAMS,*) ISLOPEUSE
        IF (ISLOPEUSE > 0) THEN
            DO_MODEL_SLOPES=.TRUE.
        ELSE
            DO_MODEL_SLOPES=.FALSE.
        ENDIF
        READ(MASS_WASTING_PARAMS,*) SLOPE_DIFFUSIVITY
        READ(MASS_WASTING_PARAMS,*) ALVCREEPFAC
        READ(MASS_WASTING_PARAMS,*) CRITICAL_GRADIENT_USE
        READ(MASS_WASTING_PARAMS,*) ROERINGUSE
        READ(MASS_WASTING_PARAMS,*) IDO_DEPTHCREEP
        READ(MASS_WASTING_PARAMS,*) CREEP_RATE_HALF_DEPTH
        IF (CRITICAL_GRADIENT_USE > 0) THEN
            USE_CRITICAL_SLOPE_CRADIENT =.TRUE.
        ELSE
            USE_CRITICAL_SLOPE_CRADIENT = .FALSE.
        ENDIF
        IF(IDO_DEPTHCREEP>0) THEN
            DEPTH_DEPENDENT_CREEP = .TRUE.
        ELSE
            DEPTH_DEPENDENT_CREEP = .FALSE.
        ENDIF
        IF (CREEP_RATE_HALF_DEPTH>0.0) THEN
            CREEP_DEPTH_CONSTANT = -LOG(0.5)/CREEP_RATE_HALF_DEPTH
        ELSE
            CREEP_DEPTH_CONSTANT = 0.0
        ENDIF
        READ(MASS_WASTING_PARAMS,*) CRITICAL_SLOPE_GRADIENT
        READ(MASS_WASTING_PARAMS,* )SLOPE_GRADIENT_EXPONENT
        READ(MASS_WASTING_PARAMS,*) MAXIMUM_DIFFUSIVITY_INCREASE
        READ(MASS_WASTING_PARAMS,*) SLOPE_FAILURE_DIFFUSIVITY
        READ(MASS_WASTING_PARAMS,*) DO_ROUTE
        READ(MASS_WASTING_PARAMS,*) AVALANCHE_RATE_CONSTANT
        READ(MASS_WASTING_PARAMS,*) AVALANCHE_SLOPE_EXPONENT
        READ(MASS_WASTING_PARAMS,*) AVALANCHE_FLUX_EXPONENT
        READ(MASS_WASTING_PARAMS,*) AVALANCHE_CRITICAL_VALUE
        READ(MASS_WASTING_PARAMS,*) MINIMUM_BEDROCK_GRADIENT
        READ(MASS_WASTING_PARAMS,*) MINIMUM_ROUTING_GRADIENT
	    CLOSE (MASS_WASTING_PARAMS)
        IF (DO_ROUTE>0) THEN
            ROUTE_REGOLITH_OVER_ROCK=.TRUE.
        ELSE
            ROUTE_REGOLITH_OVER_ROCK=.FALSE.
        ENDIF
        IF (ROERINGUSE > 0) THEN
            SLOPE_GRADIENT_EXPONENT=2.0
            SLOPE_FAILURE_DIFFUSIVITY=1.0
            USE_ROERING_MASS_WASTING=.TRUE.
        ELSE
            USE_ROERING_MASS_WASTING=.FALSE.
        ENDIF
        CRITICAL_GRADIENT_TERM = 1.0/CRITICAL_SLOPE_GRADIENT**SLOPE_GRADIENT_EXPONENT
	
        !        *******************************************************************
        !        *****************  mass wasting parameters  ************************
        !        *******************************************************************
        IF (DO_MODEL_SLOPES) THEN
            WRITE(OUTHIST,508) SLOPE_DIFFUSIVITY
            508 FORMAT(' *************** MASS WASTING PARAMETERS ***************'&
            ,/,' SLOPE EROSION IS DIRECTLY MODELED',/, &
            ' MASS MOVEMENT SLOPE SCALING=',G12.5)
            WRITE(OUTHIST,569) ALVCREEPFAC
            569   FORMAT(' RELATIVE RATE OF CREEP ON ALLUVIAL COVER =',G12.5)
            IF (USE_ROERING_MASS_WASTING) THEN
                WRITE(OUTHIST,11520)
                11520 FORMAT(' MASS WASTING USES ROERING PARAMETERIZATION')
            ELSE
                WRITE(OUTHIST,11521)
                11521 FORMAT(' HOWARD MASS WASTING PARAMETERIZATION USED')
            ENDIF
            IF (USE_CRITICAL_SLOPE_CRADIENT) WRITE(OUTHIST,513)CRITICAL_SLOPE_GRADIENT, &
            SLOPE_GRADIENT_EXPONENT,SLOPE_FAILURE_DIFFUSIVITY,MAXIMUM_DIFFUSIVITY_INCREASE
            513         FORMAT(' CRITICAL SLOPE GRADIENT=',G12.5,&
            ' CRITICAL GRADIENT POWER=',G12.5,/, &
            ' CRITICAL GRADIENT MULTIPLIER - SLOPE_FAILURE_DIFFUSIVITY=',G12.5,/, &
            'FACTOR GOVERNING MAXIMUM ENHANCEMENT OF CREEP RATE NEAR FAILURE=',G12.5)
            IF (DEPTH_DEPENDENT_CREEP) THEN
                WRITE(OUTHIST,5496) CREEP_RATE_HALF_DEPTH, CREEP_DEPTH_CONSTANT
5496            FORMAT('CREEP RATE IS DEPTH-DEPENDENT, DEPTH TO HALF-RATE=',G12.5,/, &
                     'CREEP_DEPTH_CONSTANT=',G12.5)
            ENDIF
            IF (DO_AVALANCHE) THEN
                WRITE(OUTHIST, 2495) AVALANCHE_RATE_CONSTANT, AVALANCHE_SLOPE_EXPONENT,AVALANCHE_FLUX_EXPONENT &
                    ,AVALANCHE_CRITICAL_VALUE
2495            FORMAT('AVALANCHE RATE CONSTANT=',G12.5,' AVALANCHE SLOPE EXPONENT=',G12.5 &
                            ,/,' AVALANCHE FLUX EXPONENT =',G12.5' AVALANCHE_CRITICAL_VALUE',G12.5)
            ENDIF
            WRITE(OUTHIST, 3495) MINIMUM_ROUTING_GRADIENT
3495 FORMAT('MINIMUM GRADIENT FOR WHICH WEATHERED BEDROCK IS ROUTED DOWNSLOPE=',G12.5)            
        ELSE
            WRITE(OUTHIST,438)
            438   FORMAT(' *******  SLOPE EROSION IS NOT DIRECTLY MODELED *******')
        ENDIF	
	RETURN
    END !SUBROUTINE READ_MASS_WASTING_PARAMETERS
	!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DO_MASS_WASTING()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        !      ********************************************************************
        !        This subroutine models slope mass wasting and returns the net elevation
        !         change in ERODE_SLOPE(i,j)
        !   MODIFIES:  CFW, CFN, CFNE, CFNW, ERODE_SLOPE, IS_ROCK_SURFACE, NOMINAL_ERODE_SLOPE
        !              SLOPETOADD, SLGRAVTOADD
        !      ********************************************************************
        REAL(4) :: RAPID_CREEP
        EXTERNAL RAPID_CREEP
        INTEGER :: I,J,I_NEW,J_NEW,LOCAL_DIRECTION,I_OLD,J_OLD,IE,IW,IPATH,JS,JN
        REAL(4) :: GRADIENT,SUMSLOPE
        LOGICAL (1) :: REGROUTED(IMMX,JMMX)
        INTEGER :: CHANGE_STATE(IMMX,JMMX) !TESTING STATE CHANGE 1/27/19 (0=no change, 1=to bedrock,2=to regolith)
        REAL(4) :: MAX_REGOLITH_SLOPE(IMMX,JMMX)    !! Maximum downslope gradient on regolith cells, not including adjacent bedrock cells
        REAL(4) :: REGEMIN,REGEMAX,REGEAVG,REGENUM
        REAL(4) :: ROCEMIN,ROCEMAX,ROCEAVG,ROCENUM
        REAL(4) :: DMULT,LOCAL_AVALANCHE_VALUE,AVALANCHE_FRACTION
        REAL(4) :: maxflux, maxthresh,maxgrd !test
        REAL(4) :: upslope_contribute,CREEP_DEPTH_FACTOR,EXP_ARG,MINARG,MAXARG !DIAGNOSTIC
        LOGICAL :: IREG, WRITEDETAIL,REVERSE,TESTSTOP
        !!!!!!!!write(outhist,4833)
4833    format('mass called')
        TESTSTOP=.FALSE.
        MINARG=1.0e25
        MAXARG=-MINARG
        maxflux = 0.0
        maxthresh = 0.0
        maxgrd = 0.0
        WRITEDETAIL=.FALSE.
        ROCEMIN=1.0E+25
        REGEMIN=ROCEMIN
        ROCEMAX=-ROCEMIN
        REGEMAX=-ROCEMIN
        REGEAVG=0.0
        REGENUM=0.0
        ROCEAVG=0.0
        ROCENUM=0.0
        SLOPETOADD=0.0
        SLGRAVTOADD=0.0
        CHANGE_REGOLITH_STATE=0 !TEST 1/27/19
        MAX_REGOLITH_SLOPE=0.0
        IF (DO_AVALANCHE) THEN
            AVALANCHE_FRACTION=0.0
            AVALANCHE_FLUX=0.0
        ENDIF
        !      ********************************************************************
        !       This is the main loop determining rates of mass wasting
        !      ********************************************************************
        REGROUTED=.FALSE.
            DO  J = 1, MY
			 	JS=J+1
				JN=J-1
				IF (J == 1) THEN
					IF (IS_Y_PERIODIC) THEN
                       JN=MY
                    ELSE
                       JN=2
                    ENDIF                      
				ENDIF
				IF (J== MY) THEN
				    IF (IS_Y_PERIODIC) THEN
                       JS= 1
                    ELSE
                       JS = MY-1
                    ENDIF
			    ENDIF
                DO  I = 1, MX
				    IE=I+1
					IW=I-1
				    IF (I == 1) THEN
					    IF (IS_X_PERIODIC) THEN
                           IW=MX
                        ELSE
                           IW=2
                        ENDIF
					ENDIF
					IF (I== MX) THEN
						IF (IS_X_PERIODIC) THEN
                           IE= 1
                        ELSE
                           IE = MX-1
                        ENDIF
					ENDIF

					   
                    !      ********************************************************************
                    !       first determine mass wasting rate to the west, cfw(i,j)
                    !       initially use the absolute gradient
                    !      ********************************************************************
                    GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IW,J))/CELL_SIZE
                    !      ********************************************************************
                    !       the mass wasting rate dependency on gradient is a function of the
                    !         input parameter powertouse
                    !      ********************************************************************
                    CFW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                    !      ********************************************************************
                    !       if there is provision for a critical slope gradient then adjust
                    !         the mass wasting rate in accordance with critical_gradient_term, gradient,
                    !         slope_failure_diffusivity, and slope_gradient_exponent.  
                    !       Also make sure that the rate enhancement
                    !         does not exceed maximum_diffusivity_increase - this is done in subroutine failsafe
                    !      ********************************************************************
                    IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                        IF (USE_ROERING_MASS_WASTING) THEN
                            CFW(I,J)=CFW(I,J)*RAPID_CREEP(GRADIENT)
                        ELSE
                            CFW(I,J)=CFW(I,J)+RAPID_CREEP(GRADIENT)
                        ENDIF
                    ENDIF
                    !      ********************************************************************
                    !       Change the sign of cfw if the gradient is negative
                    !      ********************************************************************
                    IF (ELEVATION(I,J) < ELEVATION(IW,J)) THEN
					     CFW(I,J)=-CFW(I,J)
						 REVERSE=.TRUE.
                         IF (IS_SEDIMENT_COVERED(IW,J)) CFW(I,J)=CFW(I,J)*ALVCREEPFAC
					ELSE
					     REVERSE=.FALSE.
                         IF (IS_SEDIMENT_COVERED(I,J)) CFW(I,J)=CFW(I,J)*ALVCREEPFAC
					ENDIF
					IF (DEPTH_DEPENDENT_CREEP) THEN
					    IF (REVERSE) THEN
						    IF (REGOLITH(IW,J)/=0.0) THEN
                                EXP_ARG=CREEP_DEPTH_CONSTANT*ABS(REGOLITH(IW,J))
                                IF (EXP_ARG>10.0) THEN
                                    CREEP_DEPTH_FACTOR=1.0
                                ELSE
                                    CREEP_DEPTH_FACTOR=1.0-EXP(-EXP_ARG)
                                ENDIF
                            ELSE
                                CREEP_DEPTH_FACTOR=0.0
                            ENDIF
						ELSE
                            IF (REGOLITH(I,J)/=0.0) THEN
                                EXP_ARG=CREEP_DEPTH_CONSTANT*ABS(REGOLITH(I,J))
                                IF (EXP_ARG>10.0) THEN
                                    CREEP_DEPTH_FACTOR=1.0
                                ELSE
                                    CREEP_DEPTH_FACTOR=1.0-EXP(-EXP_ARG)
                                ENDIF
                            ELSE
                                CREEP_DEPTH_FACTOR=0.0
                            ENDIF
						ENDIF
                        CFW(I,J)=CFW(I,J)*CREEP_DEPTH_FACTOR
					ENDIF
                    !      ********************************************************************
                    !       Do the same for the other directions, including cfn, and if nine-
                    !          point divergence is used, cfne and cfnw -- also do it slightly
                    !         differently for matrix edges - taking into account that we may have periodic
                    !         boundaries  or  no-mass-flow boundaries 
                    !      ********************************************************************
                    GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IE,JN))/CELL_SIZE
                    CFNE(I,J)=SLOPE_DIFFUSIVITY*GRADIENT             
                    IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                        GRADIENT=GRADIENT*0.7071
                        IF (USE_ROERING_MASS_WASTING) THEN
                            CFNE(I,J)=CFNE(I,J)*RAPID_CREEP(GRADIENT)
                        ELSE
                            CFNE(I,J)=CFNE(I,J)+RAPID_CREEP(GRADIENT)
                        ENDIF
                    ENDIF
					IF (ELEVATION(I,J) < ELEVATION(IE,JN)) THEN
					     CFNE(I,J)=-CFNE(I,J)
						 REVERSE=.TRUE.
                         IF (IS_SEDIMENT_COVERED(IE,JN)) CFNE(I,J)=CFNE(I,J)*ALVCREEPFAC
					ELSE
					     REVERSE=.FALSE.
                         IF (IS_SEDIMENT_COVERED(I,J)) CFNE(I,J)=CFNE(I,J)*ALVCREEPFAC
					ENDIF
					IF (DEPTH_DEPENDENT_CREEP) THEN
					    IF (REVERSE) THEN
						    IF (REGOLITH(IE,JN)/=0.0) THEN
                                EXP_ARG=CREEP_DEPTH_CONSTANT*ABS(REGOLITH(IE,JN))
                                IF (EXP_ARG>10.0) THEN
                                    CREEP_DEPTH_FACTOR=1.0
                                ELSE
                                    CREEP_DEPTH_FACTOR=1.0-EXP(-EXP_ARG)
                                ENDIF
                            ELSE
                                CREEP_DEPTH_FACTOR=0.0
                            ENDIF
						ELSE
                            IF (REGOLITH(I,J)/=0.0) THEN
                                EXP_ARG=CREEP_DEPTH_CONSTANT*ABS(REGOLITH(I,J))
                                IF (EXP_ARG>10.0) THEN
                                    CREEP_DEPTH_FACTOR=1.0
                                ELSE
                                    CREEP_DEPTH_FACTOR=1.0-EXP(-EXP_ARG)
                                ENDIF
                            ELSE
                                CREEP_DEPTH_FACTOR=0.0
                            ENDIF
						ENDIF
                        CFNE(I,J)=CFNE(I,J)*CREEP_DEPTH_FACTOR
					ENDIF
                    GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I,JN))/CELL_SIZE
                    CFN(I,J)=SLOPE_DIFFUSIVITY*GRADIENT       
                    IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                        IF (USE_ROERING_MASS_WASTING) THEN
                            CFN(I,J)=CFN(I,J)*RAPID_CREEP(GRADIENT)
                        ELSE
                            CFN(I,J)=CFN(I,J)+RAPID_CREEP(GRADIENT)
                        ENDIF
                    ENDIF
					IF (ELEVATION(I,J) < ELEVATION(I,JN)) THEN
					     CFN(I,J)=-CFN(I,J)
						 REVERSE=.TRUE.
                         IF (IS_SEDIMENT_COVERED(I,JN)) CFN(I,J)=CFN(I,J)*ALVCREEPFAC
					ELSE
					     REVERSE=.FALSE.
                         IF (IS_SEDIMENT_COVERED(I,J)) CFN(I,J)=CFN(I,J)*ALVCREEPFAC
					ENDIF
					IF (DEPTH_DEPENDENT_CREEP) THEN
					    IF (REVERSE) THEN
						    IF (REGOLITH(I,JN)/=0.0) THEN
                                EXP_ARG=CREEP_DEPTH_CONSTANT*ABS(REGOLITH(I,JN))
                                IF (EXP_ARG>10.0) THEN
                                    CREEP_DEPTH_FACTOR=1.0
                                ELSE
                                    CREEP_DEPTH_FACTOR=1.0-EXP(-EXP_ARG)
                                ENDIF
                            ELSE
                                CREEP_DEPTH_FACTOR=0.0
                            ENDIF
						ELSE
                            IF (REGOLITH(I,J)/=0.0) THEN
                                EXP_ARG=CREEP_DEPTH_CONSTANT*ABS(REGOLITH(I,J))
                                IF (EXP_ARG>10.0) THEN
                                    CREEP_DEPTH_FACTOR=1.0
                                ELSE
                                    CREEP_DEPTH_FACTOR=1.0-EXP(-EXP_ARG)
                                ENDIF
                            ELSE
                                CREEP_DEPTH_FACTOR=0.0
                            ENDIF
						ENDIF
                        CFN(I,J)=CFN(I,J)*CREEP_DEPTH_FACTOR
					ENDIF
                    GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IW,JN))/CELL_SIZE
                    CFNW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                    IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                        GRADIENT=GRADIENT*0.7071
                        IF (USE_ROERING_MASS_WASTING) THEN
                            CFNW(I,J)=CFNW(I,J)*RAPID_CREEP(GRADIENT)
                        ELSE
                            CFNW(I,J)=CFNW(I,J)+RAPID_CREEP(GRADIENT)
                        ENDIF
                    ENDIF
					IF (ELEVATION(I,J) < ELEVATION(IW,JN)) THEN
					     CFNW(I,J)=-CFNW(I,J)
						 REVERSE=.TRUE.
                         IF (IS_SEDIMENT_COVERED(IW,JN)) CFNW(I,J)=CFNW(I,J)*ALVCREEPFAC
					ELSE
					     REVERSE=.FALSE.
                         IF (IS_SEDIMENT_COVERED(I,J)) CFNW(I,J)=CFNW(I,J)*ALVCREEPFAC
					ENDIF
					IF (DEPTH_DEPENDENT_CREEP) THEN
					    IF (REVERSE) THEN
						    IF (REGOLITH(IW,JN)/=0.0) THEN
                                EXP_ARG=CREEP_DEPTH_CONSTANT*ABS(REGOLITH(IW,JN))
                                IF (EXP_ARG>10.0) THEN
                                    CREEP_DEPTH_FACTOR=1.0
                                ELSE
                                    CREEP_DEPTH_FACTOR=1.0-EXP(-EXP_ARG)
                                ENDIF
                            ELSE
                                CREEP_DEPTH_FACTOR=0.0
                            ENDIF
						ELSE
                            IF (REGOLITH(I,J)/=0.0) THEN
                                EXP_ARG=CREEP_DEPTH_CONSTANT*ABS(REGOLITH(I,J))
                                IF (EXP_ARG>10.0) THEN
                                    CREEP_DEPTH_FACTOR=1.0
                                ELSE
                                    CREEP_DEPTH_FACTOR=1.0-EXP(-EXP_ARG)
                                ENDIF
                            ELSE
                                CREEP_DEPTH_FACTOR=0.0
                            ENDIF
						ENDIF
                        CFNW(I,J)=CFNW(I,J)*CREEP_DEPTH_FACTOR
					ENDIF
                ENDDO
            ENDDO                  
            IF (ROUTE_REGOLITH_OVER_ROCK)  THEN 
                !      ********************************************************************
                !       This loop
                !       routes mass wasting material downslope on steep bedrock slopes
                !       We first determine the maximum downslope gradient on regolith-covered
                !         cells, not including any adjacent bedrock cells.
                !      ********************************************************************
                DO J=MYY,1,-1
                    DO I=1,MX
                        IW=I-1
                        IE=I+1
                        IF (IS_X_PERIODIC) THEN
                            IF (IW < 1) IW=MX
                            IF (IE > MX) IE=1
                        ENDIF
						if (iw.lt.1) iw=1
						if (ie.gt.mx) ie=mx
                        JN=J-1
                        JS=J+1
                        IF (IS_Y_PERIODIC) THEN
                            IF (JN < 1) JN=MY
                            IF (JS > MY) JS=1
                        ENDIF
						if (jn.lt.1) jn=1
						if (js.gt.my) js=my
                        IF (.NOT.IS_ROCK_SURFACE(I,J)) THEN
                             IF (.NOT.IS_ROCK_SURFACE(IW,J)) THEN
                                 GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IW,J))/CELL_SIZE
                                 IF (GRADIENT.GT.MAX_REGOLITH_SLOPE(I,J)) MAX_REGOLITH_SLOPE(I,J)=GRADIENT
                             ENDIF
                             IF (.NOT.IS_ROCK_SURFACE(IE,J)) THEN
                                 GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IE,J))/CELL_SIZE
                                 IF (GRADIENT.GT.MAX_REGOLITH_SLOPE(I,J)) MAX_REGOLITH_SLOPE(I,J)=GRADIENT
                             ENDIF
                             IF (.NOT.IS_ROCK_SURFACE(I,JS)) THEN
                                 GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I,JS))/CELL_SIZE
                                 IF (GRADIENT.GT.MAX_REGOLITH_SLOPE(I,J)) MAX_REGOLITH_SLOPE(I,J)=GRADIENT
                             ENDIF
                             IF (.NOT.IS_ROCK_SURFACE(I,JN)) THEN
                                 GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I,JN))/CELL_SIZE
                                 IF (GRADIENT.GT.MAX_REGOLITH_SLOPE(I,J)) MAX_REGOLITH_SLOPE(I,J)=GRADIENT
                             ENDIF
                             IF (.NOT.IS_ROCK_SURFACE(IW,JN)) THEN
                                 GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IW,JN))/CELL_SIZE
                                 GRADIENT=GRADIENT*0.7071
                                 IF (GRADIENT.GT.MAX_REGOLITH_SLOPE(I,J)) MAX_REGOLITH_SLOPE(I,J)=GRADIENT
                             ENDIF
                             IF (.NOT.IS_ROCK_SURFACE(IW,JS)) THEN
                                 GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IW,JS))/CELL_SIZE
                                 GRADIENT=GRADIENT*0.7071
                                 IF (GRADIENT.GT.MAX_REGOLITH_SLOPE(I,J)) MAX_REGOLITH_SLOPE(I,J)=GRADIENT
                             ENDIF
                             IF (.NOT.IS_ROCK_SURFACE(IE,JN)) THEN
                                 GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IE,JN))/CELL_SIZE
                                 GRADIENT=GRADIENT*0.7071
                                 IF (GRADIENT.GT.MAX_REGOLITH_SLOPE(I,J)) MAX_REGOLITH_SLOPE(I,J)=GRADIENT
                             ENDIF
                             IF (.NOT.IS_ROCK_SURFACE(IE,JS)) THEN
                                 GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IE,JS))/CELL_SIZE
                                 GRADIENT=GRADIENT*0.7071
                                 IF (GRADIENT.GT.MAX_REGOLITH_SLOPE(I,J)) MAX_REGOLITH_SLOPE(I,J)=GRADIENT
                             ENDIF
                        ENDIF
					enddo
				enddo
                        
                L2086: DO  J=MYY,1,-1
                    M2086: DO  I=1,MX
                        !      ********************************************************************
                        !       We're only looking at steep bedrock slopes
                        !      ********************************************************************
                        IF (.NOT.IS_ROCK_SURFACE(I,J)) CYCLE M2086 !    goto 2086
                        IF (D8_GRADIENT(I,J)<MINIMUM_ROUTING_GRADIENT) CYCLE M2086 !SKIP ROUTING ON LOW-GRADIENT BEDROCK SLOPES !test
                        IF (USE_EROSION_MASK) THEN
                            IF(EROSION_MASK(I,J)) CYCLE M2086
                        ENDIF                        
                        !      ********************************************************************
                        !       If the direction of mass transport is towards the point i,j then
                        !       add in the contributions from the neighboring points, but not if the
                        !       surrounding point is bedrock
                        !      ********************************************************************
                        upslope_contribute=0.0
                        REGROUTED(I,J)=.TRUE.
                        IPATH=0
                        IW=I-1
                        IE=I+1
                        IF (IS_X_PERIODIC) THEN
                            IF (IW < 1) IW=MX
                            IF (IE > MX) IE=1
                        ENDIF
                        JN=J-1
                        JS=J+1
                        IF (IS_Y_PERIODIC) THEN
                            IF (JN < 1) JN=MY
                            IF (JS > MY) JS=1
                        ENDIF
                        IF ((J > 1).OR.IS_Y_PERIODIC) THEN
                            IF ((CFN(I,J) < 0.0).AND.(.NOT.IS_ROCK_SURFACE(I,JN))) &
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFN(I,J)
                        ENDIF
                        IF ((I > 1).OR.IS_X_PERIODIC) THEN
                            IF ((CFW(I,J) < 0.0).AND.(.NOT.IS_ROCK_SURFACE(IW,J))) &
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFW(I,J)
                        ENDIF
                        IF ((IS_Y_PERIODIC).OR.(J < MYY)) THEN
                            IF ((CFN(I,JS) > 0.0).AND.(.NOT.IS_ROCK_SURFACE(I,JS))) &
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFN(I,JS)
                        ENDIF
                        IF ((I < MX).OR.IS_X_PERIODIC) THEN
                            IF ((CFW(IE,J) > 0.0).AND.(.NOT.IS_ROCK_SURFACE(IE,J)))  &
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFW(IE,J)
                        ENDIF
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)*CROSS_WEIGHTING
                        IF ((J > 1).OR.IS_Y_PERIODIC) THEN
                            IF ((I < MX).OR.IS_X_PERIODIC) THEN
                                IF ((CFNE(I,J) < 0.0).AND.(.NOT.IS_ROCK_SURFACE(IE,JN))) &
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNE(I,J)*DIAGONAL_WEIGHTING
                            ENDIF
                            IF ((I > 1).OR.IS_X_PERIODIC) THEN
                                IF ((CFNW(I,J) < 0.0).AND.(.NOT.IS_ROCK_SURFACE(IW,JN))) &
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNW(I,J)*DIAGONAL_WEIGHTING
                            ENDIF
                        ENDIF
                        IF ((IS_Y_PERIODIC).OR.(J < MYY)) THEN
                            IF ((I > 1).OR.IS_X_PERIODIC) THEN
                                IF ((CFNE(IW,JS) > 0.0).AND.(.NOT.IS_ROCK_SURFACE(IW,JS))) &
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNE(IW,JS)*DIAGONAL_WEIGHTING
                            ENDIF
                            IF ((I < MX).OR.IS_X_PERIODIC) THEN
                                IF ((CFNW(IE,JS) > 0.0).AND.(.NOT.IS_ROCK_SURFACE(IE,JS))) &
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNW(IE,JS)*DIAGONAL_WEIGHTING
                            ENDIF
                        ENDIF
                        upslope_contribute=erode_slope(i,j)
                        !      ********************************************************************
                        !       Here we add in the material weathered locally as stored as a negative
                        !           number in REGOLITH(I,J) 
                        !      *******************************************************************
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-REGOLITH(I,J)*VOLUME_CHANGE_COEFFICIENT
                        IF (DO_AVALANCHE) THEN
                           IF(OLD_AVALANCHE_FLUX(I,J)>0.0) THEN
                            LOCAL_AVALANCHE_VALUE = OLD_AVALANCHE_FLUX(I,J)**AVALANCHE_FLUX_EXPONENT &
                                 * D8_GRADIENT(I,J)**AVALANCHE_SLOPE_EXPONENT-AVALANCHE_CRITICAL_VALUE 
                            IF (LOCAL_AVALANCHE_VALUE>0.0) THEN
                                 AVALANCHE_FRACTION=AVALANCHE_FRACTION+1.0
                                 ERODE_SLOPE(I,J) = ERODE_SLOPE(I,J) + LOCAL_AVALANCHE_VALUE*AVALANCHE_RATE_CONSTANT
                                 AVALANCHE_FLUX(I,J)=AVALANCHE_FLUX(I,J)+ERODE_SLOPE(I,J)
                            ENDIF
                          ENDIF
                        ENDIF
                        !      ********************************************************************
                        !       We now need to route all mass wasting contributions from surrounding
                        !        points downslope to the nearest non-bedrock location.  We assume that
                        !        debris transport on bedrock slopes is "instantaneous"
                        !      ********************************************************************
                        I_OLD=I
                        J_OLD=J
                        L2087: DO
                            LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                            !      ********************************************************************
                            !       If we are in a depression we need to deposit, so branch
                            !      ********************************************************************
                            IF (LOCAL_DIRECTION < 2) THEN
                                IREG=.FALSE.
							    I_NEW=I_OLD
								J_NEW=J_OLD
								IS_ROCK_SURFACE(I_NEW,J_NEW)=.FALSE.
                                EXIT L2087
                            ENDIF
                            J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                            I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
				
                            IF (IS_X_PERIODIC) THEN
                                IF (I_NEW < 1) I_NEW=MX
                                IF (I_NEW > MX) I_NEW=1
                            ENDIF
                            IPATH=IPATH+1
                            IF (IPATH > RMMX) THEN
                                ABORT_SIMULATION=.TRUE.
                            ENDIF
                            !      ********************************************************************
                            !       Forget it if we reach a non-eroding or border location
                            !      ********************************************************************
                            IF (IS_Y_PERIODIC) THEN
                                IF (J_NEW < 1) J_NEW=MY
                                IF (J_NEW > MY) J_NEW=1
                            ELSE
                                IF (J_NEW == MY) THEN
                                    ERODE_SLOPE(I,J)=REGOLITH(I,J)*VOLUME_CHANGE_COEFFICIENT
                                    CYCLE M2086
                                ENDIF
                            ENDIF
                            IF (DO_AVALANCHE.and.is_rock_surface(i_new,j_new)) THEN
                                AVALANCHE_FLUX(I_NEW,J_NEW)=AVALANCHE_FLUX(I_NEW,J_NEW)+ERODE_SLOPE(I,J)
                            ENDIF
                            IF ((LOCAL_DIRECTION < 2).or.(d8_gradient(I_old,j_old).lt. &
                                (MINIMUM_ROUTING_GRADIENT*max_regolith_slope(i_new,j_new)))) THEN
                                IREG=.FALSE.
                                EXIT L2087
                            ENDIF			              
                            !      ********************************************************************
                            !       Branch if we have reached a regolith-mantled slope
                            !       IREG false means the local bedrock cell is to be changed to
                            !          a regolith covered cell
                            !      ********************************************************************
                            IF (DO_AVALANCHE.AND.TESTSTOP) THEN
                                IF (D8_GRADIENT(I_NEW,J_NEW)<MINIMUM_ROUTING_GRADIENT) THEN
                                    IREG=.FALSE.
                                    EXIT L2087
                                ENDIF
                            ENDIF
                            IREG=.FALSE.
                            IF ((.NOT.IS_ROCK_SURFACE(I_NEW,J_NEW)).OR.IS_SEDIMENT_COVERED(I_NEW,J_NEW) &
                                     .OR.(D8_GRADIENT(I_NEW,J_NEW)<MINIMUM_ROUTING_GRADIENT)) THEN		 																								   
       							    if(d8_gradient(I_old,j_old).lt. &
                                    (MINIMUM_ROUTING_GRADIENT*max_regolith_slope(i_new,j_new))) then
								ireg=.false.
								else
                                IREG=.TRUE.
								endif
                                EXIT L2087
                            ENDIF
                            !      ********************************************************************
                            !       Otherwise continue marching downslope
                            !      ********************************************************************
                            I_OLD=I_NEW
                            J_OLD=J_NEW
                            IF (.NOT.ABORT_SIMULATION) CYCLE L2087
                            IF (ABORT_SIMULATION) THEN
                                IABORTMIN=I_OLD-3
                                IABORTMAX=I_OLD+3
                                JABORTMIN=J_OLD-3
                                JABORTMAX=J_OLD+3
                                WRITE(*,814) I_OLD,J_OLD,I,J
                                WRITE(OUTHIST,814) I_OLD,J_OLD,I,J
                                814         FORMAT(' *** INFINITE CYCLE AT ',2I6,' ***' &
                                ,/,' IN REGOLITH ROUTING STARTING FROM: ',2I6)
                                IF ((.NOT.IS_Y_PERIODIC).AND.(J_NEW == MY)) THEN
                                    CYCLE M2086
                                ELSE
                                    EXIT L2086
                                ENDIF
                            ENDIF
                        ENDDO L2087
                        !        ********************************************************************
                        !         We have routed eroded regolith downslope until a regolith covered
                        !              slope has been reached and if the local bedrock cell
                        !              slope gradient is less than that of the cell to which the
                        !              eroded regolith is delivered, the bedrock cell is converted
                        !              to regolith
                        !        ********************************************************************
                        IF (.NOT.IREG) THEN
                            I_NEW=I_OLD
                            J_NEW=J_OLD
                            CHANGE_REGOLITH_STATE(I_NEW,J_NEW)= 2 !TEST 1/27/19
                        ELSE
                            CHANGE_REGOLITH_STATE(I_NEW,J_NEW)= 2 !TEST 1/27/19
                        ENDIF
                        !      ********************************************************************
                        !       We're at a regolith mantled slope, so we need to deposit our sediment,
                        !       which includes material mass wasted from sourrounding points and
                        !       the local contribution of weathered material
                        !      ********************************************************************
                        ERODE_SLOPE(I_NEW,J_NEW)=ERODE_SLOPE(I_NEW,J_NEW)+ERODE_SLOPE(I,J)
                        IF (DO_AVALANCHE) THEN
                                AVALANCHE_FLUX(I_NEW,J_NEW)=0.0
                        ENDIF
                        SLOPETOADD=SLOPETOADD+ERODE_SLOPE(I,J)*ABS(DFLOAT((I-I_NEW)**2) &
                        + DFLOAT((J-J_NEW)**2))
                        SLGRAVTOADD=SLGRAVTOADD+ERODE_SLOPE(I,J)*(ELEVATION(I,J)-ELEVATION(I_NEW,J_NEW))
                        !      ********************************************************************
                        !       Make sure that our local rate of erosion is just 
                        !       the local production function for now
                        !      ********************************************************************
                        ERODE_SLOPE(I,J)=REGOLITH(I,J)*VOLUME_CHANGE_COEFFICIENT
                        IF (DO_AVALANCHE) THEN
                          IF (old_avalanche_flux(i,j)>0.0) THEN
                            LOCAL_AVALANCHE_VALUE = OLD_AVALANCHE_FLUX(I,J)**AVALANCHE_FLUX_EXPONENT &
                                 * D8_GRADIENT(I,J)**AVALANCHE_SLOPE_EXPONENT-AVALANCHE_CRITICAL_VALUE 
                            IF (LOCAL_AVALANCHE_VALUE>0.0) THEN
                                 ERODE_SLOPE(I,J) = REGOLITH(I,J) - LOCAL_AVALANCHE_VALUE*AVALANCHE_RATE_CONSTANT
                                 AVALANCHE_FLUX(I,J)=AVALANCHE_FLUX(I,J)+ERODE_SLOPE(I,J)
                                 CUMULATIVE_AVALANCHE_EROSION(I,J)=CUMULATIVE_AVALANCHE_EROSION(I,J)+ &
                                     LOCAL_AVALANCHE_VALUE*AVALANCHE_RATE_CONSTANT
                            ENDIF
                          ENDIF
                        ENDIF                 
                    ENDDO M2086
                ENDDO L2086
                IF (DO_AVALANCHE) OLD_AVALANCHE_FLUX=AVALANCHE_FLUX
            ENDIF
            !      ********************************************************************
            !       This is the main loop for regolith mantled slopes and low gradient bedrock slopes
            !      ********************************************************************
            L90: DO  J=MYY,1,-1
                M90: DO  I=1,MX
                    IF (USE_EROSION_MASK) THEN
                        IF(EROSION_MASK(I,J)) CYCLE M90
                    ENDIF                    
                    !      ********************************************************************
                    !       Include contributions from surrounding points if they are regolith
                    !         mantled and losses to surrounding points if they are either
                    !         regolith or bedrock
                    !      *******************************************************************
                    IE=I+1
                    IW=I-1
                    IF (IS_X_PERIODIC) THEN
                        IF (IE > MX) IE=1
                        IF (IW < 1) IW=MX
                    ENDIF
                    JN=J-1
                    JS=J+1
                    IF (IS_Y_PERIODIC) THEN
                        IF (JN < 1) JN=MY
                        IF (JS > MY) JS=1
                    ENDIF
                    !      ********************************************************************
                    !       Skip this if it is non-eroding border point or if it is bedrock
                    !      ********************************************************************
                    IF (.NOT.REGROUTED(I,J)) THEN
                        IF ((IE <= MX).OR.IS_X_PERIODIC) THEN
                            IF (CFW(IE,J) <= 0.0) THEN
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFW(IE,J)*CROSS_WEIGHTING
                            ELSE
                                IF (.NOT.REGROUTED(IE,J)) THEN
                                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFW(IE,J)*CROSS_WEIGHTING
                                ENDIF
                            ENDIF
                        ENDIF
                        IF (CFN(I,JS) <= 0.0) THEN
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFN(I,JS)*CROSS_WEIGHTING
                        ELSE
                            IF (.NOT.REGROUTED(I,JS))  &
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFN(I,JS)*CROSS_WEIGHTING
                        ENDIF
                        IF (CFW(I,J) >= 0.0) THEN
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFW(I,J)*CROSS_WEIGHTING
                        ELSE
                            IF ((I > 1).OR.IS_X_PERIODIC) THEN
                                IF (.NOT.REGROUTED(IW,J)) THEN
                                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFW(I,J)*CROSS_WEIGHTING
                                ENDIF
                            ENDIF
                        ENDIF
                        IF (CFN(I,J) >= 0.0) THEN
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFN(I,J)*CROSS_WEIGHTING
                        ELSE
                            IF ((J > 1).OR.(IS_Y_PERIODIC)) THEN
                                IF (.NOT.REGROUTED(I,JN)) THEN
                                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFN(I,J)*CROSS_WEIGHTING
                                ENDIF
                            ENDIF
                        ENDIF
                        IF (CFNE(I,J) > 0.0) THEN
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNE(I,J)*DIAGONAL_WEIGHTING
                        ELSE
                            IF ((J > 1).OR.(IS_Y_PERIODIC)) THEN
                                IF ((I < MX).OR.IS_X_PERIODIC) THEN
                                    IF (.NOT.REGROUTED(IE,JN)) THEN
                                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNE(I,J)*DIAGONAL_WEIGHTING
                                    ENDIF
                                ENDIF
                            ENDIF
                        ENDIF
                        IF (CFNW(I,J) > 0.0) THEN
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNW(I,J)*DIAGONAL_WEIGHTING
                        ELSE
                            IF ((J > 1).OR.(IS_Y_PERIODIC)) THEN
                                IF ((I > 1).OR.IS_X_PERIODIC) THEN
                                    IF (.NOT.REGROUTED(IW,JN)) THEN
                                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNW(I,J)*DIAGONAL_WEIGHTING
                                    ENDIF
                                ENDIF
                            ENDIF
                        ENDIF
                        IF ((I > 1).OR.IS_X_PERIODIC) THEN
                            IF (CFNE(IW,JS) <= 0.0) THEN
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNE(IW,JS)*DIAGONAL_WEIGHTING
                            ELSE
                                IF (.NOT.REGROUTED(IW,JS)) THEN
                                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNE(IW,JS)*DIAGONAL_WEIGHTING
                                ENDIF
                            ENDIF
                        ENDIF
                        IF ((I < MX).OR.IS_X_PERIODIC) THEN
                            IF (CFNW(IE,JS) <= 0.0) THEN
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNW(IE,JS)*DIAGONAL_WEIGHTING
                            ELSE
                                IF (.NOT.REGROUTED(IE,JS)) THEN
                                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNW(IE,JS)*DIAGONAL_WEIGHTING
                                ENDIF
                            ENDIF
                        ENDIF
                        DMULT=CELL_SIZE
                        IF (USE_ROERING_MASS_WASTING) THEN
                            SLOPETOADD=SLOPETOADD+  &
                            (SLOPE_DIFFUSIVITY*D8_GRADIENT(I,J)*RAPID_CREEP(D8_GRADIENT(I,J)))*DMULT
                            SLGRAVTOADD=SLGRAVTOADD+  &
                            (SLOPE_DIFFUSIVITY*D8_GRADIENT(I,J)*RAPID_CREEP(D8_GRADIENT(I,J)))*DMULT*D8_GRADIENT(I,J)
                        ELSE
                            SLOPETOADD=SLOPETOADD+  &
                            (SLOPE_DIFFUSIVITY*D8_GRADIENT(I,J)+RAPID_CREEP(D8_GRADIENT(I,J)))*DMULT
                            SLGRAVTOADD=SLGRAVTOADD+  &
                            (SLOPE_DIFFUSIVITY*D8_GRADIENT(I,J)+RAPID_CREEP(D8_GRADIENT(I,J)))*DMULT*D8_GRADIENT(I,J)
                        ENDIF
                    ENDIF
                    NOMINAL_ERODE_SLOPE(I,J)=(CFN(I,JS)-CFN(I,J)-CFW(I,J))*CROSS_WEIGHTING  &
                    +(-CFNE(I,J)-CFNW(I,J))*DIAGONAL_WEIGHTING
                    IF ((I < MX).OR.IS_X_PERIODIC) THEN
                        NOMINAL_ERODE_SLOPE(I,J)=NOMINAL_ERODE_SLOPE(I,J)+CFW(IE,J)*CROSS_WEIGHTING &
                        + CFNW(IE,JS)*DIAGONAL_WEIGHTING
                    ENDIF
                    IF ((I > 1).OR.IS_X_PERIODIC) THEN
                        NOMINAL_ERODE_SLOPE(I,J)=NOMINAL_ERODE_SLOPE(I,J)+CFNE(IW,JS)*DIAGONAL_WEIGHTING
                    ENDIF
                ENDDO M90
            ENDDO L90
        SUMSLOPE=0.0
        DO  J=1,MY
            DO  I=1,MX
                SUMSLOPE=SUMSLOPE+ERODE_SLOPE(I,J)           
            ENDDO
        ENDDO
        IF (WRITEDETAIL) WRITE(OUTHIST,4234) SUMSLOPE
        4234  FORMAT ('NET MASS WASTING ELEV CHANGE=',G12.5)
        DO J=1,MY
            DO I=1,MX
                IF (IS_ROCK_SURFACE(I,J))THEN
                    IF (ERODE_SLOPE(I,J) > ROCEMAX) ROCEMAX=ERODE_SLOPE(I,J)
                    IF (ERODE_SLOPE(I,J) < ROCEMIN) ROCEMIN=ERODE_SLOPE(I,J)
                    ROCEAVG=ROCEAVG+ERODE_SLOPE(I,J)
                    ROCENUM=ROCENUM+1.0
                ELSE
                    IF (ERODE_SLOPE(I,J) > REGEMAX) REGEMAX=ERODE_SLOPE(I,J)
                    IF (ERODE_SLOPE(I,J) < REGEMIN) REGEMIN=ERODE_SLOPE(I,J)
                    REGEAVG=REGEAVG+ERODE_SLOPE(I,J)
                    REGENUM=REGENUM+1.0
                ENDIF
            ENDDO
        ENDDO
        IF (ROCENUM > 0.0) ROCEAVG=ROCEAVG/ROCENUM
        IF (REGENUM > 0.0) REGEAVG=REGEAVG/REGENUM
        IF (WRITEDETAIL) WRITE(OUTHIST,14234) ROCEMIN,ROCEAVG,ROCEMAX,ROCENUM,REGEMIN, &
        REGEAVG,REGEMAX,REGENUM
14234   FORMAT(' ROCK SLOPE EROSION: MIN=',G12.5,' AVG=',G12.5,' MAX=', &
        G12.5,' N=',G12.5,/, &
        ' REGOLITH SLOPE EROSION: MIN=',G12.5,' AVG=',G12.5,' MAX=', &
        G12.5,' N=',G12.5)
        !write(*,2045) maxflux,maxthresh,maxgrd
2045    format('maxflux=',g12.5,' maxthresh=',g12.5,' maxgrd=',g12.5)
        IF (DO_AVALANCHE) THEN
             AVALANCHE_FRACTION=AVALANCHE_FRACTION/(MX*MY)
             WRITE(*,2047) AVALANCHE_FRACTION
2047         FORMAT('AF=',G12.5)
        ENDIF
        RETURN
    END !  SUBROUTINE DO_MASS_WASTING()
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL(4) FUNCTION RAPID_CREEP(GRADIENT)
        USE ERODE_GLOBALS
        !      ********************************************************************
        !        This is the look-up function for slope mass-wasting rate if there is a
        !         critical gradient for slope failure
        !      ********************************************************************
        IMPLICIT NONE
        LOGICAL USE_LOOKUP
        REAL(4), INTENT(IN)  :: GRADIENT
        INTEGER ICAT
       !USE_LOOKUP=.FALSE.  ! COMMENT OUT ONE OF THESE STATEMENTS
        USE_LOOKUP=.TRUE.						 
        IF (GRADIENT <= 0.0) THEN
            RAPID_CREEP=0.0
            RETURN
        ENDIF
        IF (GRADIENT >= GRADMAX) THEN
            RAPID_CREEP=FAILMAX
        ELSE
            IF (USE_LOOKUP) THEN
                ICAT=INT((100.0*GRADIENT)/GRADMAX)+1
                RAPID_CREEP= &
                RATECAT(ICAT)+(GRADIENT-GRADCAT(ICAT))*DELRATE(ICAT)
            ELSE
                RAPID_CREEP= &
                1.0/(1.0-(GRADIENT/GRADMAX)**SLOPE_GRADIENT_EXPONENT)
            ENDIF
        ENDIF
        RETURN
    END !REAL(4) FUNCTION RAPID_CREEP(GRADIENT) 
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL(4) FUNCTION ROCK_MASS_WASTING(GRADIENT)
        USE ERODE_GLOBALS
        !      ********************************************************************
        !        This is the look-up function for bedrock mass-wasting rate if there is a
        !         critical gradient for slope failure
        !      ********************************************************************
        IMPLICIT NONE
        REAL(4), INTENT(IN) :: GRADIENT
        INTEGER ICAT
        LOGICAL USE_LOOKUP
        USE_LOOKUP=.FALSE.
        IF (GRADIENT <= 0.0) THEN
           ROCK_MASS_WASTING=0.0
           RETURN
        ENDIF
        IF (GRADIENT >= RGRADMAX) THEN
            ROCK_MASS_WASTING=FAILMAX
        ELSE
            IF (USE_LOOKUP) THEN
                ICAT=INT((100.0*GRADIENT)/RGRADMAX)+1
                ROCK_MASS_WASTING= &
                RRATECAT(ICAT)+(GRADIENT-RGRADCAT(ICAT))*RDELRATE(ICAT)
            ELSE
                ROCK_MASS_WASTING= &
                    1.0/(1.0-(GRADIENT/RGRADMAX)**SLOPE_GRADIENT_EXPONENT)
            ENDIF
        ENDIF
        RETURN
    END ! REAL(4) FUNCTION ROCK_MASS_WASTING(GRADIENT)
    
