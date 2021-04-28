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
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	SUBROUTINE READ_WEATHERING_PARAMETERS
	USE ERODE_GLOBALS
    USE AREA_GLOBALS
	IMPLICIT NONE
	    CHARACTER (80) :: TEXTLINE
	    INTEGER WEATHERING_PARAMS, ISEEPAGEWEATHER,WEATHER2USE,SPATIAL_WEATHERING_USE
        REAL(4):: LOW_WEATHERING_SCALE,HIGH_WEATHERING_SCALE
		WEATHERING_PARAMS=93
        !        *******************************************************************
        !        ****************  regolith weathering parameters  ******************
        !        *******************************************************************
		OPEN(WEATHERING_PARAMS,FILE='WEATHERING_PARAMETERS.PRM',ACTION='READ')
        READ(WEATHERING_PARAMS,22) TEXTLINE
22      FORMAT(A)        
        WRITE(*,22) TEXTLINE
        READ(WEATHERING_PARAMS,*) ROCK_WEATHERING_RATE
        READ(WEATHERING_PARAMS,*) WEATHER_DECAY_RATE
        READ(WEATHERING_PARAMS,*) INITIAL_REGOLITH_THICKNESS
        READ(WEATHERING_PARAMS,*) VOLUME_CHANGE_COEFFICIENT
        READ(WEATHERING_PARAMS,*) WEATHER2USE
        READ(WEATHERING_PARAMS,*) WEATHERING_TERM_2
        READ(WEATHERING_PARAMS,*) WEATHERING_DECAY_2
        READ(WEATHERING_PARAMS,*) ISEEPAGEWEATHER
        READ(WEATHERING_PARAMS,*) SEEPAGE_WEATHERING_SCALING
        READ(WEATHERING_PARAMS,*) SEEPAGE_WEATHERING_EXPONENT
        IF (ISEEPAGEWEATHER > 0) THEN
            SEEPAGE_WEATHERING=.TRUE.
        ELSE
            SEEPAGE_WEATHERING=.FALSE.
        ENDIF
        IF (WEATHER2USE > 0) THEN
            TWO_TERM_WEATHERING = .TRUE.
        ELSE
            TWO_TERM_WEATHERING = .FALSE.
        ENDIF
        READ(WEATHERING_PARAMS,*) CRITICAL_BEDROCK_GRADIENT
        READ(WEATHERING_PARAMS,*) WEATHER_MULT
        READ(WEATHERING_PARAMS,*) WEATHER_DIVERGENCE
        READ(WEATHERING_PARAMS,*) READREGOLITH
        !IF (.NOT.NEW_SIMULATION) READREGOLITH=1
        READ(WEATHERING_PARAMS,*) SPATIAL_WEATHERING_USE
        READ(WEATHERING_PARAMS,*) LOW_WEATHERING_SCALE
        READ(WEATHERING_PARAMS,*) HIGH_WEATHERING_SCALE
        IF (DO_SPATIAL_VARIATION.AND.(SPATIAL_WEATHERING_USE>0)) THEN
            USE_SPATIAL_WEATHERING=.TRUE.
            WEATHERING_AVAL=(HIGH_WEATHERING_SCALE-LOW_WEATHERING_SCALE)/(SPATIAL_MAX-SPATIAL_MIN)
            WEATHERING_BVAL=LOW_WEATHERING_SCALE-WEATHERING_AVAL*SPATIAL_MIN
        ELSE
            USE_SPATIAL_WEATHERING=.FALSE.
        ENDIF
        CLOSE(WEATHERING_PARAMS)		
       !        *******************************************************************
        !        ****************  regolith weathering parameters  ******************
        !        *******************************************************************
        WRITE(OUTHIST,788) ROCK_WEATHERING_RATE,WEATHER_DECAY_RATE,INITIAL_REGOLITH_THICKNESS
        788   FORMAT(' ******************* REGOLITH WEATHERING PARAMETERS ****' &
        ,/,' MAXIMUM WEATHERING RATE=',G12.5,' DEPTH TO HALF RATE=', &
        G12.5,' INITIAL REGOLITH THICKNESS=',G12.5)
        IF (TWO_TERM_WEATHERING) THEN
            WRITE(OUTHIST,779) WEATHERING_TERM_2,WEATHERING_DECAY_2
            779    FORMAT(' TWO EXPONENTIAL WEATHERING RATE USED',/,&
            ' RATE CONSTANT 2=',G12.5,' DECAY CONSTANT 2=',G12.5)
        ENDIF
        WRITE(OUTHIST,1788) VOLUME_CHANGE_COEFFICIENT
1788    FORMAT('VOLUME OF WEATHERED REGOLITH RELATIVE TO UNWEATHERED ROCK=',G12.5)
        WRITE(OUTHIST,787) CRITICAL_BEDROCK_GRADIENT,WEATHER_DIVERGENCE
        787   FORMAT(' GRADIENT VALUE GIVING DOUBLE WEATERING RATE=',G12.5,  &
        /,' DIVERGENCE VALUE GIVING DOUBLE WEATHERING RATE=', &
        G12.5)
        IF (WEATHER_DIVERGENCE > 0.0) THEN
            WEATHER_DIVERGENCE  = LOG(2.0E0) / WEATHER_DIVERGENCE
        ENDIF
        WRITE(OUTHIST,786) CRITICAL_BEDROCK_GRADIENT,WEATHER_DIVERGENCE,READREGOLITH
        786   FORMAT(' SCALED WEATHER GRADIENT=',G12.5,' SCALED DIVERGENCE', &
        ' VALUE=',G12.5,/,' READ REGOLITH THICKNESS FROM FILE', &
        ' (0=NO, 1=YES): ',I5)
        IF (SEEPAGE_WEATHERING) THEN
            WRITE(OUTHIST,1516) SEEPAGE_WEATHERING_SCALING,SEEPAGE_WEATHERING_EXPONENT
            1516 FORMAT (' ROCK WEATHERING BY EMERGENT SEEPAGE IS MODELED' &
            ,/,' SEEPAGE WEATHERING RATE FACTOR=',G12.5,' SEEPAGE WEATHERING EXPONENT=',G12.5)
        ENDIF
        WRITE(OUTHIST,2516) WEATHER_MULT
2516    FORMAT('SCALING FACTOR FOR GRADIENT-DEPENDENT MASS WASTING OF STEEP BEDROCK=',G12.5)
        IF (READREGOLITH>0) THEN
            WRITE(OUTHIST,2576)
2576        FORMAT('REGOLITH IS READ FROM "INREG.DAT"')
        ENDIF
        IF (DO_SPATIAL_VARIATION.AND.(SPATIAL_WEATHERING_USE>0)) THEN
            WRITE(OUTHIST,2586) LOW_WEATHERING_SCALE,HIGH_WEATHERING_SCALE
2586        FORMAT('SPATIALY-VARYING WEATHERING BEING USED',/, &
               'WITH MINIMUM RELATIVE RATE=',G12.5,' AND MAXIMUM RELATIVE RATE=',G12.5)
        ENDIF
	RETURN
	END !SUBROUTINE READ_WEATHERING_PARAMETERS
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
    !       ********************************************************************
    !        calculate topgraphic divergence
    !    ACCESSES: MX,MY,IS_X_PERIODIC,IS_Y_PERIODIC,ELEVATION
    !    MODIFIES: DIVERGENCE
    !       ********************************************************************
    SUBROUTINE CALCULATE_DIVERGENCE()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,IE,IW,JN,JS
        DO  J=1,MY
            DO  I=1,MX
                IW=I-1
                IE=I+1
                IF (IS_X_PERIODIC) THEN
                    IF (IW < 1) IW=MX
                    IF (IE > MX) IE=1
                ELSE
                    IF (IW < 1) IW=2
                    IF (IE > MX) IE=MX-1
                ENDIF
                JN=J-1
                JS=J+1
                IF (IS_Y_PERIODIC) THEN
                    IF (JN < 1) JN=MY
                    IF (JS > MY) JS=1
                ELSE
                    IF (JN < 1) JN=2
                    IF (JS > MY) JS=MY-1
                ENDIF
                !      ********************************************************************
                !      Divergence based on 9 point formula with reflection at top
                !         and bottom boundaries and periodic left&right boundaries
                !      ********************************************************************
                DIVERGENCE(I,J)=(ELEVATION(IW,JN)+ELEVATION(IE,JN)+ELEVATION(IW,JS)+ELEVATION(IE,JS)+ &
                4.0*(ELEVATION(IE,J)+ELEVATION(IW,J)+ELEVATION(I,JN)+ELEVATION(I,JS))-20.0*ELEVATION(I,J))/ &
                (6.0*CELL_AREA)
            ENDDO
        ENDDO
        RETURN
    END ! SUBROUTINE CALCULATE_DIVERGENCE
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
     !       ********************************************************************
    !        calculate topgraphic divergence
    !    This routine can utilize weighting filters of size 3x3,5x5,7x7,or 9x9
    !      depending upon the size of surrounding terrain used in the calculation.
    !    ACCESSES: MX,MY,IS_X_PERIODIC,IS_Y_PERIODIC,ELEVATION,WEIGHTS_3,WEIGHTS_5
    !              WEIGHTS_7,WEIGHTS_9
    !    MODIFIES: DIVERGENCE
    !       ********************************************************************   
    SUBROUTINE NEW_CALCULATE_DIVERGENCE(NN)
    USE ERODE_GLOBALS
    IMPLICIT NONE
    INTEGER I,J,M,N,NN,NNLOW,NNHIGH,II,JJ,NNN,MMM
    REAL*8 DIVSUM,DIVSQ,DIVAVG,DIVSD
    NNLOW=-NN / 2
    NNHIGH= NN / 2 
   ! WRITE(10,455) NN, NNLOW,NNHIGH
455 FORMAT('NN=',I5,' NNLOW=',I5,' NNHIGH=',I5)
    DIVERGENCE=0.0
    DIVSUM=0.0
    DIVSQ=0.0
    DO I=1,MX
        DO J=I,MY
            DO N=NNLOW,NNHIGH
                NNN=N+NN-NN/2
                DO M=NNLOW,NNHIGH
                    MMM=M+NN-NN/2
                    II=I+N
                    JJ=J+M
                    IF (II<1) THEN
                        IF (IS_X_PERIODIC) THEN
                            II=MX+N+I
                        ELSE
                            II=2-N+I
                        ENDIF
                    ENDIF
                    IF (II>MX) THEN
                        IF (IS_X_PERIODIC) THEN
                            II=N+I-MX
                        ELSE
                            II=2*MX-I-N
                        ENDIF
                    ENDIF
                    IF (JJ<1) THEN
                        IF (IS_Y_PERIODIC) THEN
                            JJ=MY+M+J
                        ELSE
                            JJ=2-M+J
                        ENDIF
                    ENDIF
                    IF (JJ>MY) THEN
                        IF (IS_Y_PERIODIC) THEN
                            JJ=M+J-MY
                        ELSE
                            JJ=2*MY-J-M
                        ENDIF
                    ENDIF
                    IF (NN==3) THEN
                        DIVERGENCE(I,J)=DIVERGENCE(I,J)+WEIGHTS_3(NNN,MMM)*ELEVATION(II,JJ)
                    ELSE
                        IF (NN==5) THEN
                            DIVERGENCE(I,J)=DIVERGENCE(I,J)+WEIGHTS_5(NNN,MMM)*ELEVATION(II,JJ)
                        ELSE
                            IF (NN==7) THEN
                                DIVERGENCE(I,J)=DIVERGENCE(I,J)+WEIGHTS_7(NNN,MMM)*ELEVATION(II,JJ)
                            ELSE
                                IF (NN==9) THEN
                                    DIVERGENCE(I,J)=DIVERGENCE(I,J)+WEIGHTS_9(NNN,MMM)*ELEVATION(II,JJ)
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
            DIVERGENCE(I,J)=DIVERGENCE(I,J)/CELL_AREA
        ENDDO
    ENDDO
  !  WRITE(10,200) NN,DIVAVG,DIVSD
200 FORMAT('NN=',I3,' DIVAVG=',G12.5,' DIVSD=',G12.5)
    RETURN
    END  !       SUBROUTINE NEW_CALCULATE_DIVERGENCE
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DO_WEATHERING()
        USE ERODE_GLOBALS
        USE ACCRETION_GLOBALS
        USE GROUNDWATER_VARIABLES
		USE MASS_FLOW_GLOBALS
        USE AREA_GLOBALS
        IMPLICIT NONE
        !      ********************************************************************
        !       This subroutine determines the weathering rate on bare bedrock
        !         and weathering rate and regolith thickness on regolith covered
        !         slopes
        !       For regolith-mantled slopes the matrix regolith(i,j) is the regolith thickness
        !       For bedrock slopes regolith(i,j) is less than zero and is the rate of weathering of rock into regolith
        !   ACCESSES: MX, MY, MYY, SUBMERGED, ROCK_WEATHERING_RATE, WEATHERING_TERM_2
        !             SEEPAGE_WEATHERING_SCALING,GROUNDWATER_FLUX,SEEPAGE_WEATHERING_EXPONENT
        !             CRITICAL_BEDROCK_GRADIENT,D8_GRADIENT, WEATHER_MULT,TWO_TERM_WEATHERING
        !             RESISTANT_SURFACE_LAYER,RAD_ERODE,RAD_CONST,RELATIVE_RESISTANCE 
        !             TIME_INCREMENT,IS_DEPRESSION,WEATHER_DECAY_RATE,SEEPAGE_WEATHERING
        !             D8_GRADIENT,USE_SOLAR_EROSION,DIVERGENCE,
        !   MODIFIES: REGOLITH,IS_ROCK_SURFACE 
        !   CALLS: ROCK_MASS_WASTING, FIND_DEPRESSION
        !      ********************************************************************
        REAL(4) :: DIVERGE,MAXRATE,WTERM,FAILTERM,ROCK_MASS_WASTING,MAXCOSFACTOR,RAD_ERODE
        REAL(4) :: COSFACTOR,AVGWTERM,MAXWTERM,XAVGSEEP,MAXSEEP,XSEEP,XTOTAL
        REAL(4) :: FRACTBEDROCK,MINREGOLITH,MAXREGOLITH,AVGREGOLITH, REGCHANGE
		real (4) :: FLOW_BASE_ELEVATION,MAX_ACTIVE_DEPTH
		REAL (4) :: SCOURED_CHANGE, STERM, ACTIVE_DEPTH, SINE_TERM, FLOW_THICKNESS 
        EXTERNAL ROCK_MASS_WASTING, RAD_ERODE
        LOGICAL :: IS_DEPRESSION,WRITEDETAIL,SKIPBEDROCK
        INTEGER :: I,J
         !!!!!!!!write(outhist,4833)
4833    format('weathering called')       
        WRITEDETAIL=.FALSE.
        FRACTBEDROCK=0.0
        MINREGOLITH=1.0E+25
        MAXREGOLITH=-1.0E+25
        AVGREGOLITH=0.0
        XTOTAL=0.0
        MAXSEEP=-1.0E+25
        MAXWTERM=-1.0E+25
        XAVGSEEP=0.0
        MAXSEEP=0.0
        AVGWTERM=0.0
        MAXWTERM=0.0
        MAXCOSFACTOR=10.0
        TOTALRADERODE=0.0
        L110: DO  J=1,MYY
            M110: DO  I=1,MX
                SKIPBEDROCK=.FALSE.
                COSFACTOR=AMIN1(MAXCOSFACTOR,SQRT(1.0+D8_GRADIENT(I,J)**2))
                IF (SUBMERGED(I,J)) CYCLE M110
                IF (USE_EROSION_MASK) THEN
                  IF(EROSION_MASK(I,J)) CYCLE M110
                ENDIF                
                IF (IS_ROCK_SURFACE(I,J)) THEN
                    IF (TWO_TERM_WEATHERING) THEN
                        MAXRATE = COSFACTOR*(ROCK_WEATHERING_RATE - WEATHERING_TERM_2)
                    ELSE
                        MAXRATE = COSFACTOR*ROCK_WEATHERING_RATE
                    ENDIF
                    IF (DO_SPATIAL_VARIATION.AND.USE_SPATIAL_WEATHERING) THEN
                        MAXRATE=MAXRATE*(WEATHERING_AVAL*SPATIAL_VALUE(I,J)+WEATHERING_BVAL)
                    ENDIF
                    IF (SEEPAGE_WEATHERING) THEN
                        IF (SEEPAGE_AVERAGING) THEN
                            XSEEP=SEEPAGE_WEATHERING_SCALING*(FILTERED_GROUNDWATER_FLUX(I,J))**SEEPAGE_WEATHERING_EXPONENT
                        ELSE
                            XSEEP=SEEPAGE_WEATHERING_SCALING*(GROUNDWATER_FLUX(I,J))**SEEPAGE_WEATHERING_EXPONENT
                        ENDIF
                        IF (XSEEP > 0.0) THEN
                             REGOLITH(I,J)=-(XSEEP/RELATIVE_RESISTANCE(I,J))
                        ENDIF
                        XAVGSEEP=XAVGSEEP+XSEEP
                        XTOTAL=XTOTAL+1.0
                        IF (XSEEP > MAXSEEP) MAXSEEP=XSEEP
                        MAXRATE=MAXRATE+XSEEP
                    ENDIF
                    !     *********************************************************************
                    !      Following is for bare bedrock slopes
                    !      The first part calculates a divergence-dependent weathering term
                    !      (if weather_divergence is greater than zero
                    !     *********************************************************************
                    IF (.NOT.SKIPBEDROCK) THEN
                        IF (WEATHER_DIVERGENCE > 0.0) THEN
                            DIVERGE=EXP(-WEATHER_DIVERGENCE*DIVERGENCE(I,J))
                        ELSE
                        !     *********************************************************************
                        !      Just set divergence to unity if there is no divergence dependency
                        !     *********************************************************************
                            DIVERGE=1.0
                        ENDIF
                        IF (CRITICAL_BEDROCK_GRADIENT > 0.0) THEN
                            FAILTERM=WEATHER_MULT*ROCK_MASS_WASTING(D8_GRADIENT(I,J))
                        ELSE
                            FAILTERM=0.0
                        ENDIF
                         !     *********************************************************************
                         !      Weathering rate is equal to a flat surface value possibly enhanced
                         !         by slope gradient (if critical_bedrock_gradient>0) and by surface convexity
                         !         (if weather_divergence>0)
                         !     *********************************************************************
                        IF (USE_SOLAR_EROSION) THEN
                            REGOLITH(I,J)= (-RAD_CONST*(RAD_ERODE(I,J)+FAILTERM))
                            TOTALRADERODE=TOTALRADERODE+RAD_CONST*RAD_ERODE(I,J)
                        ELSE
                            REGOLITH(I,J)=REGOLITH(I,J) -MAXRATE*(DIVERGE+FAILTERM)     &
                            /RELATIVE_RESISTANCE(I,J)
                            IF (RESISTANT_SURFACE_LAYER) REGOLITH(I,J)=REGOLITH(I,J)/RELATIVE_RESISTANCE(I,J)
                        ENDIF
                    ENDIF
                    CUMULATIVE_WEATHERING(I,J)=CUMULATIVE_WEATHERING(I,J)-REGOLITH(I,J)*TIME_INCREMENT    
                ELSE
                    !     *********************************************************************
                    !      For regolith-mantled slopes calculate increase in regolith thickness
                    !         Rate law assumes negative exponential dependency on regolith
                    !         thickness with rock_weathering_rate being the maximum rate at zero
                    !         thickness.  Other functional dependencies could be used, including
                    !         a maximum rate at a finite thickness.
                    !     *********************************************************************
                    WTERM=COSFACTOR*ROCK_WEATHERING_RATE
                    IF (DO_SPATIAL_VARIATION.AND.USE_SPATIAL_WEATHERING) THEN
                        WTERM=WTERM*(WEATHERING_AVAL*SPATIAL_VALUE(I,J)+WEATHERING_BVAL)
                    ENDIF
                   
                    STERM=0.0
                    !    **********************************************************************
                    !     This section calculates any bedrock erosion by mass flow and adds
                    !        The eroded material to the mobile regolith
                    !    **********************************************************************
                    IF (USE_FLOW_VOLUME) THEN
                        IF (MAXIMUM_FLOW_DEPTH_EROSION>0.0) THEN
                            FLOW_THICKNESS=MIN1(REGOLITH(I,J),MAXIMUM_FLOW_DEPTH_EROSION)
                        ENDIF
                                         !  added by ORKAN 12/19/2014
                        IF ((FLOW_VOLUME_FLUX(I,J)>0.0).and.(MASS_FLOW_EROSION_RATE>0.0)) THEN
                           SINE_TERM=sqrt(D8_GRADIENT(I,J)**2/(D8_GRADIENT(I,J)**2+1.0))
                           if (this_is_bingham_flow) then
                                ACTIVE_DEPTH=FLOW_THICKNESS-DELTA_HB/SINE_TERM
                           else
                               active_depth=flow_thickness
                           endif
                           IF ((CRITICAL_SCOUR_THICKNESS>0.0).AND.(ACTIVE_DEPTH>CRITICAL_SCOUR_THICKNESS)) THEN
                          IF (DEPTH_EXPONENT>0.0) THEN                             
                            STERM=MASS_FLOW_EROSION_RATE*(ACTIVE_DEPTH-CRITICAL_SCOUR_THICKNESS)
                          ELSE
                              IF (FLUX_EXPONENT>0.0) THEN
                                  STERM=MASS_FLOW_EROSION_RATE*FLOW_VOLUME_FLUX(I,J)
                              ELSE
                                  STERM=MASS_FLOW_EROSION_RATE
                              ENDIF
                          ENDIF
                          ENDIF
                          STERM_NUMBER=STERM_NUMBER+1.0
                          IF (STERM>STERM_MAX) STERM_MAX=STERM
                          IF (STERM<STERM_MIN) STERM_MIN=STERM
                          STERM_AVERAGE=STERM_AVERAGE+STERM
                          STERM=MASS_FLOW_EROSION_RATE*(STERM - MASS_FLOW_CRITICAL_VALUE)
                          IF (COSFACTOR>MAXCOSFACTOR) MAXCOSFACTOR=COSFACTOR
                          IF (ACTIVE_DEPTH>MAX_ACTIVE_DEPTH) MAX_ACTIVE_DEPTH=ACTIVE_DEPTH
                          STERM=MAX1(0.0,STERM)
                        ENDIF
                        ! End of mass flow erosion section
                    IF (REGOLITH(I,J) <= 0.0) THEN					
                        REGOLITH(I,J)=0.0
                        IF (DO_FLUVIAL_DETACHMENT.OR.DO_SEDIMENT_TRANSPORT) THEN
                            IF ((D8_GRADIENT(I,J) >= 0.0).AND.((ELEVATION(I,J)-SEDIMENT_BASE(I,J)) <= 0.0 )) THEN
                                IS_ROCK_SURFACE(I,J)=.TRUE.
                                IS_SEDIMENT_COVERED(I,J)=.FALSE.
                            ELSE
                                IS_ROCK_SURFACE(I,J)=.FALSE.
                            ENDIF
                        ELSE
                            IF (D8_GRADIENT(I,J) >= 0.0) THEN
                                IS_ROCK_SURFACE(I,J)=.TRUE.
                                IS_SEDIMENT_COVERED(I,J)=.FALSE.
                            ELSE
                                IS_ROCK_SURFACE(I,J)=.FALSE.
                            ENDIF
                        ENDIF
                    ELSE
                        WTERM = WTERM*EXP(-WEATHER_DECAY_RATE*REGOLITH(I,J))
                    ENDIF
                    IF (USE_FLOW_VOLUME) THEN                                                  !  added by ORKAN 12/19/2014
                        STERM = STERM*ABS(FLOW_VOLUME_FLUX(I,J))                                  !  added by ORKAN 12/19/2014  - this model might change
                        ENDIF                                                                   !  added by ORKAN 12/19/2014
                    ENDIF 					
                    IF (SEEPAGE_WEATHERING) THEN
                        IF (SEEPAGE_AVERAGING) THEN
                            XSEEP=SEEPAGE_WEATHERING_SCALING*(FILTERED_GROUNDWATER_FLUX(I,J))**SEEPAGE_WEATHERING_EXPONENT
                        ELSE
                            XSEEP=SEEPAGE_WEATHERING_SCALING*(GROUNDWATER_FLUX(I,J))**SEEPAGE_WEATHERING_EXPONENT
                        ENDIF
                        XSEEP=XSEEP*EXP(-WEATHER_DECAY_RATE*REGOLITH(I,J))
                        IF (XSEEP < 0.0) XSEEP=0.0
                        XAVGSEEP=XAVGSEEP+XSEEP
                        AVGWTERM=AVGWTERM+WTERM
                        XTOTAL=XTOTAL+1.0
                        IF (WTERM > MAXWTERM) MAXWTERM=WTERM
                        IF (XSEEP > MAXSEEP) MAXSEEP=XSEEP
                        WTERM=WTERM+XSEEP
                    ENDIF
                    IF (TWO_TERM_WEATHERING) THEN
                        WTERM = WTERM - WEATHERING_TERM_2*EXP(-WEATHERING_DECAY_2*REGOLITH(I,J))
                    ENDIF
                    IF (RESISTANT_SURFACE_LAYER) WTERM=WTERM/RELATIVE_RESISTANCE(I,J)
                    REGCHANGE=TIME_INCREMENT * WTERM*VOLUME_CHANGE_COEFFICIENT/RELATIVE_RESISTANCE(I,J)
                    REGOLITH(I,J)=REGOLITH(I,J) + REGCHANGE
                    CUMULATIVE_WEATHERING(I,J)=CUMULATIVE_WEATHERING(I,J)+REGCHANGE
                ENDIF
                IF (IS_ROCK_SURFACE(I,J)) FRACTBEDROCK=FRACTBEDROCK+1.0
                IF (REGOLITH(I,J) > MAXREGOLITH) MAXREGOLITH=REGOLITH(I,J)
                IF (REGOLITH(I,J) < MINREGOLITH) MINREGOLITH=REGOLITH(I,J)
                AVGREGOLITH=AVGREGOLITH+REGOLITH(I,J)
            ENDDO M110
        ENDDO L110
        IF (USE_SOLAR_EROSION) THEN
            DO J=1,MY
                DO I=1,MX
                    CALL FIND_DEPRESSION(I,J,IS_DEPRESSION)
                    IF (IS_DEPRESSION) THEN
                        IS_ROCK_SURFACE(I,J)=.FALSE.
                    ENDIF
                ENDDO
            ENDDO
        ENDIF
        FRACTBEDROCK=FRACTBEDROCK/(MX*MY)
        AVGREGOLITH=AVGREGOLITH/(MX*MY)
        IF (WRITEDETAIL) WRITE(*,782) FRACTBEDROCK,MINREGOLITH,AVGREGOLITH,MAXREGOLITH,COSFACTOR
        782   FORMAT(' FRACTROCK=',G12.5,' MINREG=',G12.5,' AVGREG=',G12.5,/, &
        ' MAXREG=',G12.5,' COSFACTOR=',G12.5)
        IF (XTOTAL > 0.0) THEN
            XAVGSEEP=XAVGSEEP/XTOTAL
            AVGWTERM=AVGWTERM/XTOTAL
        ENDIF
        IF (WRITEDETAIL) THEN
            IF (SEEPAGE_WEATHERING) THEN
                WRITE(*,298) AVGWTERM,MAXWTERM,XAVGSEEP,MAXSEEP
                298 FORMAT ('AVG & MAX WTERM=',2G12.5,'AVG & MAX XSEEP=',2G12.5)
            ENDIF
            IF (USE_SOLAR_EROSION) THEN
                WRITE(*,844) TOTALRADERODE
                844 FORMAT(' TOTALRADERODE=',G12.5)
            ENDIF
        ENDIF
        IF (USE_FLOW_VOLUME.AND.(STERM_NUMBER>0.0)) THEN
            STERM_AVERAGE=STERM_AVERAGE/STERM_NUMBER
            WRITE(OUTHIST,1186) STERM_MIN,STERM_AVERAGE,STERM_MAX,STERM_NUMBER, &
               MAXCOSFACTOR,MAX_ACTIVE_DEPTH,MINREGOLITH,AVGREGOLITH,MAXREGOLITH 
               1186           FORMAT('STERM_MIN=',G12.5,' STERM_AVERAGE=',G12.5,' STERM_MAX=',G12.5,' STERM_NUMBER=',G12.5 &
                           ,/,' MAX_COSFACTOR=',G12.5,' MAX_ACTIVE_DEPTH=',G12.5 &
                           ,/,' MIN REGOLITH=',G12.5,' AVG REGOLITH=',G12.5,' MAX REGOLITH=',G12.5)
        ENDIF		
        RETURN
    END ! SUBROUTINE DO_WEATHERING
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       ********************************************************************
    !        Find depressions
    !    ACCESSES:  IS_X_PERIODIC, IS_Y_PERIODIC, ELEVATION, MX, MY
    !       ********************************************************************
    SUBROUTINE FIND_DEPRESSION(I,J,IS_DEPRESSION)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,J
        INTEGER :: IL,IH,JL,JH,II,JJ,IC,JC
        REAL(4) :: ECOMP
        LOGICAL, INTENT(OUT) :: IS_DEPRESSION
        ECOMP=ELEVATION(I,J)
        IL=I-1
        IH=I+1
        JL=J-1
        JH=J+1
        IS_DEPRESSION=.TRUE.
        DO JJ=JL,JH
                JC=JJ
                IF (IS_Y_PERIODIC) THEN
                    IF (JC < 1) JC=MY
                    IF (JC > MY) JC=1
                ELSE
                    IF (JC < 1) CYCLE
                    IF (JC > MY) CYCLE
                ENDIF
            DO II=IL,IH
                IC=II
                IF (IS_X_PERIODIC) THEN
                    IF (IC < 1) IC=MX
                    IF (IC > MX) IC=1
                ELSE
                    IF (IC < 1) CYCLE
                    IF (IC > MY) CYCLE
                ENDIF
                IF ((IC == I).AND.(JC == J)) CYCLE
                IF (ELEVATION(IC,JC) <= ECOMP) IS_DEPRESSION=.FALSE.
            ENDDO
        ENDDO
        RETURN
    END ! SUBROUTINE FIND_DEPRESSION
