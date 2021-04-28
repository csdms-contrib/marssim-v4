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
	SUBROUTINE READ_EOLIAN_PARAMETERS
	USE ERODE_GLOBALS
	USE EOLIAN_GLOBALS
    USE ACCRETION_GLOBALS
	IMPLICIT NONE
	    INTEGER :: EOLIAN_PARAMS, IEXPOSUREMETHOD, IUSENORMAL, ICHOSE, DO_ONLY_EOLIAN_DEPOSITION, I
        REAL(4) :: TEMP2,TEMP3
		CHARACTER (80) :: TEXTLINE
		EOLIAN_PARAMS=93
		OPEN(EOLIAN_PARAMS,FILE='EOLIAN_PARAMETERS.PRM',ACTION='READ')
        !       *********************************************************************
        !
        !       Read in parameters for eolian erosion-deposition
        !
        !       *********************************************************************
        READ(EOLIAN_PARAMS,22) TEXTLINE
  22    FORMAT(A)
        WRITE(*,22) TEXTLINE
        READ(EOLIAN_PARAMS,*) EOLIAN_EVENT_PROBABILITY
        IF (FLUVIAL_AND_SLOPE_MODELING) EOLIAN_EVENT_PROBABILITY=EOLIAN_EVENT_PROBABILITY/DEFAULT_TIME_INCREMENT
        READ(EOLIAN_PARAMS,*) EOLIAN_TIME_INCREMENT
        READ(EOLIAN_PARAMS,*) IEXPOSUREMETHOD
        READ(EOLIAN_PARAMS,*) IUSENORMAL
        USE_TOTAL_EXPOSURE=.FALSE.
        USE_DIVERGENCE_EXPOSURE=.FALSE.
        IF (IEXPOSUREMETHOD==1) THEN
            USE_TOTAL_EXPOSURE=.TRUE.
            WRITE(OUTHIST,7201)
        ELSE
            IF (IEXPOSUREMETHOD == 2) THEN
               USE_DIVERGENCE_EXPOSURE=.TRUE.
               WRITE(OUTHIST,7202)
            ELSE
                WRITE(OUTHIST,7203)
            ENDIF
        ENDIF
7201    FORMAT('USING TOTAL EXPOSURE METHOD')
7202    FORMAT('USING DIVERGENCE EXPOSURE METHOD')
7203    FORMAT('USING DEFAULT EXPOSURE METHOD')
        IF (IUSENORMAL > 0) THEN
            DEFAULT_EOLIAN_PROCESS=.TRUE.
        ELSE
            DEFAULT_EOLIAN_PROCESS=.FALSE.
        ENDIF
        READ(EOLIAN_PARAMS,*) MINIMUM_EOLIAN_DEPOSIT_RATE
        READ(EOLIAN_PARAMS,*) MAXIMUM_EOLIAN_DEPOSIT_RATE
        READ(EOLIAN_PARAMS,*) ICHOSE
        READ(EOLIAN_PARAMS,*) DO_ONLY_EOLIAN_DEPOSITION 
        EOLIAN_CONSTANT_1=(MAXIMUM_EOLIAN_DEPOSIT_RATE-MINIMUM_EOLIAN_DEPOSIT_RATE)/2.0
        EOLIAN_CONSTANT_2=(MAXIMUM_EOLIAN_DEPOSIT_RATE+MINIMUM_EOLIAN_DEPOSIT_RATE)/2.0
        IF (ICHOSE == 1) THEN
            READ(EOLIAN_PARAMS,*) EXPOSURE_10_PERCENT
            READ(EOLIAN_PARAMS,*) EXPOSURE_90_PERCENT
            EXPOSURE_50_PERCENT=(EXPOSURE_10_PERCENT+EXPOSURE_90_PERCENT)/2.0
            EOLIAN_CONSTANT_3=(MAXIMUM_EOLIAN_DEPOSIT_RATE-EOLIAN_CONSTANT_2)/((EXPOSURE_90_PERCENT-EXPOSURE_50_PERCENT)* &
            SQRT(EOLIAN_CONSTANT_1**2+(MAXIMUM_EOLIAN_DEPOSIT_RATE-EOLIAN_CONSTANT_2)**2))
        ENDIF
        IF (ICHOSE == 2) THEN
            READ(EOLIAN_PARAMS,*) EXPOSURE_50_PERCENT
            READ(EOLIAN_PARAMS,*) EXPOSURE_90_PERCENT
            EOLIAN_CONSTANT_3=(MAXIMUM_EOLIAN_DEPOSIT_RATE-EOLIAN_CONSTANT_2)/((EXPOSURE_90_PERCENT-EXPOSURE_50_PERCENT)*  &
            SQRT(EOLIAN_CONSTANT_1**2+(MAXIMUM_EOLIAN_DEPOSIT_RATE-EOLIAN_CONSTANT_2)**2))
        ENDIF
        IF (ICHOSE == 3) THEN
            IF (MINIMUM_EOLIAN_DEPOSIT_RATE < 0.0) THEN
                READ(EOLIAN_PARAMS,*) ZERO_PERCENT_EXPOSURE
                READ(EOLIAN_PARAMS,*) EXPOSURE_90_PERCENT
                TEMP2=EOLIAN_CONSTANT_2/SQRT(EOLIAN_CONSTANT_1**2-EOLIAN_CONSTANT_2**2)
                TEMP3=(MAXIMUM_EOLIAN_DEPOSIT_RATE-EOLIAN_CONSTANT_2)/SQRT(EOLIAN_CONSTANT_1**2+ &
                (MAXIMUM_EOLIAN_DEPOSIT_RATE-EOLIAN_CONSTANT_2)**2)
                EXPOSURE_50_PERCENT=(ZERO_PERCENT_EXPOSURE+EXPOSURE_90_PERCENT*TEMP2/TEMP3)/(1.0+TEMP2/TEMP3)
                EOLIAN_CONSTANT_3=TEMP3/(EXPOSURE_90_PERCENT-EXPOSURE_50_PERCENT)
            ELSE
                WRITE(*,350)
                350           FORMAT(' RL MUST BE LESS THAN ZERO')
                STOP
            ENDIF
        ENDIF
        IF (ICHOSE == 4) THEN
            IF (MINIMUM_EOLIAN_DEPOSIT_RATE < 0.0) THEN
                READ(EOLIAN_PARAMS,*)RATE0
                READ(EOLIAN_PARAMS,*) EXPOSURE_50_PERCENT
                IF (RATE0 >= ((MINIMUM_EOLIAN_DEPOSIT_RATE+MAXIMUM_EOLIAN_DEPOSIT_RATE)/2.0)) THEN
                    WRITE(*,361)
                    361               FORMAT('R0 MUST BE LESS THAT MEAN RATE')
                    STOP
                ENDIF
                IF (EXPOSURE_50_PERCENT <= 0) THEN
                    WRITE(*,362)
                    362               FORMAT(' E50 MUST BE GREATER THAN ZERO')
                    STOP
                ENDIF
                EOLIAN_CONSTANT_3=-(RATE0-EOLIAN_CONSTANT_2)/(EXPOSURE_50_PERCENT*  &
                SQRT(EOLIAN_CONSTANT_1**2-(RATE0-EOLIAN_CONSTANT_2)**2))
            ELSE
                 READ(EOLIAN_PARAMS,*) RATE0
                 READ(EOLIAN_PARAMS,*) EXPOSURE_50_PERCENT
                EOLIAN_CONSTANT_1 = 0.0
                EOLIAN_CONSTANT_3 = 0.0
            ENDIF
        ENDIF
        IF (ICHOSE == 0) THEN
             READ(EOLIAN_PARAMS,*) RATE0
             READ(EOLIAN_PARAMS,*) EXPOSURE_50_PERCENT
                EOLIAN_CONSTANT_1 = 0.0
                EOLIAN_CONSTANT_3 = 0.0
        ENDIF
        IF (DO_ONLY_EOLIAN_DEPOSITION>0) THEN
            ONLY_EOLIAN_DEPOSITION=.TRUE.
        ELSE
            ONLY_EOLIAN_DEPOSITION=.FALSE.
        ENDIF
        !       ********************************************************************
        !       Distance weighting factors used for either simulated eolian erosion/deposition or accretion/ablation
        !       ********************************************************************
        READ(EOLIAN_PARAMS,*) DISTANCE_DECAY_FACTOR
        READ(EOLIAN_PARAMS,*) WEIGHTING_CALCULATION_DISTANCE
        READ(EOLIAN_PARAMS,*) MAX_OUT_DISTANCE
        READ(EOLIAN_PARAMS,*) WRITE_INITIAL_EXPOSURE
		CLOSE(EOLIAN_PARAMS)      
        WEIGHTING_DECAY_FACTOR=LOG(0.5)/DISTANCE_DECAY_FACTOR
        IF (DISTANCE_DECAY_FACTOR /= 0.0) THEN
            MAXIMUM_WEIGHT_DISTANCE=MIN(WEIGHTING_CALCULATION_DISTANCE,INT(LOG(0.01)/(WEIGHTING_DECAY_FACTOR)))
        ELSE
            MAXIMUM_WEIGHT_DISTANCE=WEIGHTING_CALCULATION_DISTANCE
        ENDIF
        WRITE(OUTHIST,3563) MAXIMUM_WEIGHT_DISTANCE,WEIGHTING_DECAY_FACTOR
        3563  FORMAT(' WEIGHTING RANGE=',I6,' WEIGHTING DECAY=',G12.5)
        CALL SETUP_DISTANCE_WEIGHTING() 
        !        *********************************************************************
        !         Write out eolian parameters
        !        *********************************************************************
        IF (MODEL_EOLIAN_CHANGES) THEN
            WRITE(OUTHIST,2094)
            2094 FORMAT('*********************EOLIAN PARAMETERS***************************')
            WRITE(OUTHIST,1539) EOLIAN_EVENT_PROBABILITY,EOLIAN_TIME_INCREMENT
            1539 FORMAT(' EOLIAN EVENT PROBABILITY=',G12.5,' EOLIAN TIME INCREMENT=',G12.5)
            IF (USE_TOTAL_EXPOSURE) THEN
                WRITE(OUTHIST,1540)
                1540 FORMAT(' ALL CELLS WITHIN WINDOW USED FOR WEIGHTING')
            ELSE
                WRITE(OUTHIST,1541)
                1541 FORMAT(' ONLY VISIBLE CELLS WITHIN WINDOW USED FOR WEIGHTING')
            ENDIF
            IF (DEFAULT_EOLIAN_PROCESS) THEN
                WRITE(OUTHIST,1542)
                1542 FORMAT(' DEPOSITION AND EROSION ARE NORMAL TO THE LAND SURFACE')
            ELSE
                WRITE(OUTHIST,1543)
                1543 FORMAT(' DEPOSITION AND EROSION OCCUR IN VERTICAL DIRECTION')
            ENDIF
            WRITE(OUTHIST,2100) MINIMUM_EOLIAN_DEPOSIT_RATE,MAXIMUM_EOLIAN_DEPOSIT_RATE
            2100  FORMAT(' INPUT PARAMETERS: MINIMUM_EOLIAN_DEPOSIT_RATE=',G12.5, &
            ' MAXIMUM_EOLIAN_DEPOSIT_RATE=',G12.5)
            SELECT CASE (ICHOSE)
                CASE (1)
                    WRITE(OUTHIST,1544) EXPOSURE_10_PERCENT, EXPOSURE_90_PERCENT
                    1544 FORMAT(' 10TH AND 90TH PERCENT EXPOSURE INDEX VALUES: ', G12.5,' ',G12.5)
                CASE (2)
                    WRITE(OUTHIST,1545) EXPOSURE_50_PERCENT, EXPOSURE_90_PERCENT
                    1545 FORMAT(' 50TH AND 90TH PERCENT EXPOSURE INDEX VALUES: ', G12.5,' ',G12.5)
                CASE (3)
                    WRITE(OUTHIST,1546) ZERO_PERCENT_EXPOSURE,EXPOSURE_90_PERCENT
                    1546 FORMAT(' 0TH AND 90TH PERCENT EXPOSURE INDEX VALUES: ',G12.5,' ',G12.5)
                CASE(4)
                    WRITE(OUTHIST,1547) RATE0, EXPOSURE_50_PERCENT
                    1547 FORMAT(' RATE AT 0.0 EXPOSURE AND 90TH PERCENT EXPOSURE INDEX VALUES: ',G12.5,' ',G12.5)
                CASE DEFAULT
                    WRITE(OUTHIST,1548)
                    1548 FORMAT(' INVALID CHOICE OF INPUT PARAMETERS')
                END SELECT
            IF (ONLY_EOLIAN_DEPOSITION) THEN
               WRITE(OUTHIST,2548)
2548           FORMAT('NO EOLIAN EROSION IS MODELED, ONLY DEPOSITION')
            ENDIF
            IF (WRITE_INITIAL_EXPOSURE>0) THEN
                WRITE(OUTHIST,2549)
2549            FORMAT('A FILE OF INITIAL EXPOSURE VALUES IS WRITTEN')
            ENDIF             
        ENDIF
        IF (MODEL_EOLIAN_CHANGES.OR.MODEL_ACCRETION_AND_ABLATION.OR.USE_SOLAR_EROSION) THEN
            WRITE(OUTHIST,1549) DISTANCE_DECAY_FACTOR,WEIGHTING_CALCULATION_DISTANCE
            1549 FORMAT(' DECAY FACTOR FOR DISTANCE WEIGHTING=',G12.5,' MAXIMUM CALCULATION DISTANCE=',I5)        
            WRITE(OUTHIST,2095) MAXIMUM_WEIGHT_DISTANCE,WEIGHTING_DECAY_FACTOR
            2095    FORMAT(' MAXIMUM RANGE=',I5, &
            ' WEIGHTING_DECAY_FACTOR=',G12.5)
            DO  I=1,350
                IF (WEIGHTING_DECAY_FACTOR /= 0.0) THEN
                    WEIGHTX(I)=EXP(WEIGHTING_DECAY_FACTOR*I)
                    WEIGHTD(I)=EXP(WEIGHTING_DECAY_FACTOR*I*1.414)/(1.414)
                ELSE
                    WEIGHTX(I)=1.0
                    WEIGHTD(I)=1.0/(1.414)
                ENDIF
            ENDDO
        ENDIF	
	RETURN
	END !SUBROUTINE READ_EOLIAN_PARAMETERS
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      *********************************************************************
    !       Do eolian deposion and erosion based upon an "exposure" index which is
    !       positive if the location is on a relative high, zero on a level plain,
    !       and negative in depressions or valleys
    !       For examples of use and model presentation see:
    !          forsberg-taylor, n. k., howard, a. d., and craddock, r. a., 2004,
    !          j. geophys. res., 109, e05002, doi:10.1029/2004je002242
    !       and
    !          luo, w., and howard, a. d., 2008, j. geophysical research, 113,
    !          e05002,doi:10.1029/2007je002981
    !       Two variant subroutines to define relative exposure are exposure and totalexposure.
    !         In the exposure subroutine only cells visible to the local cell are used in computing
    !           the weighted exposure
    !         In the totalexposure subroutine all cells within the calculation window are used to
    !           compute exposure
    !       this routine is not mass-conservative.  It assumes that deposition is primarily by
    !       suspension from imported sediment and eroded sediment is exported as suspended
    !       material.  Therefore it is not suitable for modeling eolian deposition and
    !       erosion by saltation.
    !       If default_eolian_process is true, deposition occurs normal to the surface otherwise
    !         it is vertically-directed (without the inverse cosine correction)
    !       If eolian erosion occurs and fluvial and slope modeling is also used, then
    !        the rate of erosion is dependent upon whether bedrock or regolith occurs at the surface
    !   MODIFIES: ERODE_SLOPE, ELEVATION, CUMULATIVE_EOLIAN_CHANGE, IS_ROCK_SURFACE
    !             SEDIMENT_BASE
    !   CALLS: TOTAL_EXPOSURE, EXPOSURE
    !      *********************************************************************
    SUBROUTINE DO_EOLIAN_CHANGE()
        USE ERODE_GLOBALS
        USE EOLIAN_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,IW,IE,JN,JS
        REAL(4) :: EXPOSE,ABSCHANGE,RATEMULT
        REAL(4) :: TEMP1,EOLRATE,MAXVAL,MINVAL
        REAL(4) :: MINEOLIAN,MAXEOLIAN,ETEMP,MINEXPOSERATE,MAXEXPOSERATE
        REAL(4) :: MINEXPOSE,MAXEXPOSE,AVGEXPOSE,EXPOSESSQ
        REAL(4), ALLOCATABLE, DIMENSION(:,:) :: OLDVAL
        MAXVAL=-1.0E+25
        MINVAL=-MAXVAL
        MINEXPOSE=1.0E+25
        MAXEXPOSE=-MINEXPOSE
        AVGEXPOSE=0.0
        EXPOSESSQ=0.0
        MINEXPOSERATE=1.0E+25
        MAXEXPOSERATE=-MINEXPOSERATE
        MINEOLIAN=MINVAL
        MAXEOLIAN=MAXVAL
        L310: DO  J=1,MY
            M310: DO  I=1,MX
                IF (USE_TOTAL_EXPOSURE) THEN
                    CALL TOTAL_EXPOSURE(I,J,EXPOSE)
                ELSE
                    IF (USE_DIVERGENCE_EXPOSURE) THEN
                        CALL DIVERGENCE5X5(I,J,EXPOSE)
                    ELSE
                        CALL EXPOSURE(I,J,EXPOSE)
                    ENDIF
                ENDIF
                AVGEXPOSE=AVGEXPOSE+EXPOSE
                EXPOSESSQ=EXPOSESSQ+EXPOSE**2
                TEMP1=EOLIAN_CONSTANT_3*(EXPOSE-EXPOSURE_50_PERCENT)
                EOLRATE=EOLIAN_CONSTANT_2+EOLIAN_CONSTANT_1*TEMP1/SQRT(1.0+TEMP1**2)
                IF (ONLY_EOLIAN_DEPOSITION) THEN
                    IF (EOLRATE<0.0) EOLRATE=0.0
                ENDIF
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                    IF (IS_ROCK_SURFACE(I,J)) THEN
                        IF (EOLRATE < 0.0) THEN
                            ERODE_SLOPE(I,J)=EOLRATE/REGOLITH_CRITICAL_SHEAR_FACTOR
                        ELSE
                            ERODE_SLOPE(I,J)=EOLRATE
                        ENDIF
                    ELSE
                        IF (DEFAULT_EOLIAN_PROCESS) THEN
                            ERODE_SLOPE(I,J)=SQRT(D8_GRADIENT(I,J)**2+1.0)*EOLRATE
                        ELSE
                            ERODE_SLOPE(I,J)=EOLRATE
                        ENDIF
                    ENDIF
                ELSE
                    IF (DEFAULT_EOLIAN_PROCESS) THEN
                        ERODE_SLOPE(I,J)=SQRT(D8_GRADIENT(I,J)**2+1.0)*EOLRATE
                    ELSE
                        ERODE_SLOPE(I,J)=EOLRATE
                    ENDIF
                ENDIF
                IF (.NOT.FLUVIAL_AND_SLOPE_MODELING) THEN
                    IF ((ERODE_SLOPE(I,J) < 0.0).AND.(ELEVATION(I,J) < INITIAL_ELEVATION(I,J)))&
                    THEN
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)/REGOLITH_CRITICAL_SHEAR_FACTOR
                    ENDIF
            ENDIF

                IF (EXPOSE > MAXEXPOSE) THEN
                    MAXEXPOSE=EXPOSE
                    MAXEXPOSERATE=ERODE_SLOPE(I,J)
                ENDIF
                IF (EXPOSE < MINEXPOSE) THEN
                    MINEXPOSE=EXPOSE
                    MINEXPOSERATE=ERODE_SLOPE(I,J)
                ENDIF
                IF (ERODE_SLOPE(I,J) < MINEOLIAN) THEN
                    MINEOLIAN=ERODE_SLOPE(I,J)
                ENDIF
                IF (ERODE_SLOPE(I,J) > MAXEOLIAN) THEN
                    MAXEOLIAN=ERODE_SLOPE(I,J)
                ENDIF
            ENDDO M310
        ENDDO L310
        MAXIMUM_ELEVATION_CHANGE=0.0
        DO  J=1,MY
            DO  I=1,MX
                ABSCHANGE=ABS(ERODE_SLOPE(I,J))
                IF (ABSCHANGE > MAXIMUM_ELEVATION_CHANGE) MAXIMUM_ELEVATION_CHANGE=ABSCHANGE
            ENDDO
        ENDDO
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            RATEMULT=TIME_INCREMENT
        ELSE
            RATEMULT=EOLIAN_TIME_INCREMENT
        ENDIF
        EXPOSESSQ=SQRT((EXPOSESSQ-AVGEXPOSE**2/(MX*MY))/&
        (MX*MY-1))
        AVGEXPOSE=AVGEXPOSE/(MX*MY)
        WRITE(OUTHIST,11310) MINEOLIAN,MAXEOLIAN,EOLIAN_TIME_INCREMENT
        !WRITE(*,11310) MINEOLIAN,MAXEOLIAN,EOLIAN_TIME_INCREMENT
        11310         FORMAT(' MINEOLIAN=',G12.5,' MAXEOLIAN=',G12.5, &
        ' EOLIAN_TIME_INCREMENT=',G12.5)
        WRITE(OUTHIST,21310) MINEXPOSE,MINEXPOSERATE, &
        MAXEXPOSE,MAXEXPOSERATE,  &
        AVGEXPOSE,EXPOSESSQ
        !WRITE(*,21310) MINEXPOSE,MINEXPOSERATE, &
!        MAXEXPOSE,MAXEXPOSERATE, &
!        AVGEXPOSE,EXPOSESSQ
        21310         FORMAT('MINEXPOSE=',G12.5,' MINEXPOSERATE=',G12.5,/,&
        ' MAXEXPOSE=',G12.5,' MAXEXPOSERATE=',G12.5,        &
        ' AVG=',G12.5,' SD=',G12.5)
        DO  J=1,MY
            DO  I=1,MX 
                ETEMP=ERODE_SLOPE(I,J)*RATEMULT
                IF ((.NOT.IS_ROCK_SURFACE(I,J)).AND.(FLUVIAL_AND_SLOPE_MODELING)) THEN
                    IF (ETEMP.LT.0.0) THEN
                        ETEMP=MAX(-REGOLITH(I,J),ETEMP)
                    ENDIF
                    REGOLITH(I,J)=REGOLITH(I,J)+ETEMP
                    IF (REGOLITH(I,J).LE.0.0) IS_ROCK_SURFACE(I,J)=.TRUE.
                ENDIF
                IF ((IS_ROCK_SURFACE(I,J)).AND.(FLUVIAL_AND_SLOPE_MODELING)) THEN
                    IF (ETEMP.GT.0.0) THEN
                        REGOLITH(I,J)=ETEMP
                        IS_ROCK_SURFACE(I,J)=.FALSE.
                    ENDIF
                ENDIF
                ELEVATION(I,J)=ELEVATION(I,J)+ETEMP
                CUMULATIVE_EOLIAN_CHANGE(I,J)=CUMULATIVE_EOLIAN_CHANGE(I,J)+ETEMP
                ACCRETE_WORK=ACCRETE_WORK+ETEMP
                ERODE_SLOPE(I,J)=0.0
            ENDDO
        ENDDO
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            DO  J=1,MY
                DO  I=1,MX
                    IF (ELEVATION(I,J) > SEDIMENT_BASE(I,J)) THEN
                        IS_ROCK_SURFACE(I,J)=.FALSE.
                    ELSE
                        IS_ROCK_SURFACE(I,J)=.TRUE.
                        SEDIMENT_BASE(I,J)=ELEVATION(I,J)
                    ENDIF
                ENDDO
            ENDDO
        ENDIF
        RETURN
    END ! SUBROUTINE DO_EOLIAN_CHANGE()
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE EXPOSURE(I,J,EXPOSE)
        USE ERODE_GLOBALS
        USE EOLIAN_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,J
        INTEGER :: K,ILOC,JLOC
        REAL(4) :: ECOMP,ELOC,WLOC,LSLOPE,WEIGHTXUM,WEIGHTDUM
        REAL(4) :: SLOPE1SUM,SLOPE2SUM,SLOPE3SUM,SLOPE4SUM
        REAL(4) :: SLOPE5SUM,SLOPE6SUM,SLOPE7SUM,SLOPE8SUM, MAXSLOPE
        REAL(4), INTENT(OUT) :: EXPOSE
        LOGICAL DOCALC
        DOCALC=.TRUE.
        ECOMP=ELEVATION(I,J)
        ILOC=I+1
        IF (ILOC > MX) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=1
            ELSE
                SLOPE1SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,J)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE1SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L100: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I+K
                IF (ILOC > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC-MX
                    ELSE
                        EXIT L100
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,J)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTX(K)
                    SLOPE1SUM=SLOPE1SUM+LSLOPE*WLOC
                    WEIGHTXUM=WEIGHTXUM+WLOC
                ENDIF
            ENDDO L100
            SLOPE1SUM=SLOPE1SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        ILOC=I-1
        IF (ILOC < 1) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=MX
            ELSE
                SLOPE2SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,J)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE2SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L110: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I-K
                IF (ILOC < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC+MX
                    ELSE
                        EXIT L110
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,J)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTX(K)
                    SLOPE2SUM=SLOPE2SUM+LSLOPE*WLOC
                    WEIGHTXUM=WEIGHTXUM+WLOC
                ENDIF
            ENDDO L110
            SLOPE2SUM=SLOPE2SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J+1
        IF (JLOC > MY) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=1
            ELSE
                SLOPE3SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(I,JLOC)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE3SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L120: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J+K
                IF (JLOC > MY) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC-MY
                    ELSE
                        EXIT L120
                    ENDIF
                ENDIF
                ELOC=ELEVATION(I,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTX(K)
                    SLOPE3SUM=SLOPE3SUM+LSLOPE*WLOC
                    WEIGHTXUM=WEIGHTXUM+WLOC
                ENDIF
            ENDDO L120
            SLOPE3SUM=SLOPE3SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J-1
        IF (JLOC < 1) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=MY
            ELSE
                SLOPE4SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(I,JLOC)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE4SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L130: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J-K
                IF (JLOC < 1) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC+MY
                    ELSE
                        EXIT L130
                    ENDIF
                ENDIF
                ELOC=ELEVATION(I,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTX(K)
                    SLOPE4SUM=SLOPE4SUM+LSLOPE*WLOC
                    WEIGHTXUM=WEIGHTXUM+WLOC
                ENDIF
            ENDDO L130
            SLOPE4SUM=SLOPE4SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        ILOC=I+1
        JLOC=J+1
        IF (ILOC > MX) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=1
            ELSE
                SLOPE5SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (JLOC > MY) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=1
            ELSE
                SLOPE5SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=(ELOC-ECOMP)
            MAXSLOPE=LSLOPE
            SLOPE5SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L140: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I+K
                JLOC=J+K
                IF (ILOC > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC-MX
                    ELSE
                        EXIT L140
                    ENDIF
                ENDIF
                IF (JLOC > MY) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC-MY
                    ELSE
                        EXIT L140
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTD(K)
                    SLOPE5SUM=SLOPE5SUM+LSLOPE*WLOC
                    WEIGHTDUM=WEIGHTDUM+WLOC
                ENDIF
            ENDDO L140
            SLOPE5SUM=SLOPE5SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J-1
        ILOC=I-1
        IF (ILOC < 1) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=MX
            ELSE
                SLOPE6SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (JLOC < 1) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=MY
            ELSE
                SLOPE6SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE6SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L150: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I-K
                JLOC=J-K
                IF (ILOC < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC+MX
                    ELSE
                        EXIT L150
                    ENDIF
                ENDIF
                IF (JLOC < 1) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC+MY
                    ELSE
                        EXIT L150
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTD(K)
                    SLOPE6SUM=SLOPE6SUM+LSLOPE*WLOC
                    WEIGHTDUM=WEIGHTDUM+WLOC
                ENDIF
            ENDDO L150
            SLOPE6SUM=SLOPE6SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J+1
        ILOC=I-1
        IF (JLOC > MY) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=1
            ELSE
                SLOPE7SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (ILOC < 1) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=MX
            ELSE
                SLOPE7SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE7SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L160: DO  K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J+K
                ILOC=I-K
                IF (JLOC > MY) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC-MY
                    ELSE
                        EXIT L160
                    ENDIF
                ENDIF
                IF (ILOC < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC+MX
                    ELSE
                        EXIT L160
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTD(K)
                    SLOPE7SUM=SLOPE7SUM+LSLOPE*WLOC
                    WEIGHTDUM=WEIGHTDUM+WLOC
                ENDIF
            ENDDO L160
            SLOPE7SUM=SLOPE7SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J-1
        ILOC=I+1
        IF (JLOC < 1) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=MY
            ELSE
                SLOPE8SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (ILOC > MX) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=1
            ELSE
                SLOPE8SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE8SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L170: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J-K
                ILOC=I+K
                IF (JLOC < 1) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC+MY
                    ELSE
                        EXIT L170
                    ENDIF
                ENDIF
                IF (ILOC > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC-MX
                    ELSE
                        EXIT L170
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTD(K)
                    SLOPE8SUM=SLOPE8SUM+LSLOPE*WLOC
                    WEIGHTDUM=WEIGHTDUM+WLOC
                ENDIF
            ENDDO L170
            SLOPE8SUM=SLOPE8SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        EXPOSE=(SLOPE1SUM+SLOPE2SUM+SLOPE3SUM+SLOPE4SUM+SLOPE5SUM+ &
        SLOPE6SUM+SLOPE7SUM+SLOPE8SUM)/(8.0*CELL_SIZE)
        RETURN
    END ! SUBROUTINE EXPOSURE(I,J,EXPOSE)

    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE TOTAL_EXPOSURE(I,J,EXPOSE)
        USE ERODE_GLOBALS
        USE EOLIAN_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,J
        INTEGER :: K,ILOC,JLOC
        REAL(4) :: ECOMP,ELOC,WLOC,LSLOPE,WEIGHTXUM,WEIGHTDUM
        REAL(4) :: SLOPE1SUM,SLOPE2SUM,SLOPE3SUM,SLOPE4SUM
        REAL(4) :: SLOPE5SUM,SLOPE6SUM,SLOPE7SUM,SLOPE8SUM
        REAL(4), INTENT(OUT) :: EXPOSE
        LOGICAL DOCALC
        DOCALC=.TRUE.
        ECOMP=ELEVATION(I,J)
        ILOC=I+1
        IF (ILOC > MX) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=1
            ELSE
                SLOPE1SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,J)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            SLOPE1SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L100: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I+K
                IF (ILOC > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC-MX
                    ELSE
                        EXIT L100
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,J)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTX(K)
                SLOPE1SUM=SLOPE1SUM+LSLOPE*WLOC
                WEIGHTXUM=WEIGHTXUM+WLOC
            ENDDO L100
            SLOPE1SUM=SLOPE1SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        ILOC=I-1
        IF (ILOC < 1) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=MX
            ELSE
                SLOPE2SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,J)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            SLOPE2SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L110: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I-K
                IF (ILOC < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC+MX
                    ELSE
                        EXIT L110
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,J)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTX(K)
                SLOPE2SUM=SLOPE2SUM+LSLOPE*WLOC
                WEIGHTXUM=WEIGHTXUM+WLOC
            ENDDO L110
            SLOPE2SUM=SLOPE2SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J+1
        IF (JLOC > MY) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=1
            ELSE
                SLOPE3SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(I,JLOC)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            SLOPE3SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L120: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J+K
                IF (JLOC > MY) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC-MY
                    ELSE
                        EXIT L120
                    ENDIF
                ENDIF
                ELOC=ELEVATION(I,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTX(K)
                SLOPE3SUM=SLOPE3SUM+LSLOPE*WLOC
                WEIGHTXUM=WEIGHTXUM+WLOC
            ENDDO L120
            SLOPE3SUM=SLOPE3SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J-1
        IF (JLOC < 1) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=MY
            ELSE
                SLOPE4SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(I,JLOC)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            SLOPE4SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L130: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J-K
                IF (JLOC < 1) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC+MY
                    ELSE
                        EXIT L130
                    ENDIF
                ENDIF
                ELOC=ELEVATION(I,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTX(K)
                SLOPE4SUM=SLOPE4SUM+LSLOPE*WLOC
                WEIGHTXUM=WEIGHTXUM+WLOC
            ENDDO L130
            SLOPE4SUM=SLOPE4SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        ILOC=I+1
        JLOC=J+1
        IF (ILOC > MX) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=1
            ELSE
                SLOPE5SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (JLOC > MY) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=1
            ELSE
                SLOPE5SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=(ELOC-ECOMP)
            SLOPE5SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L140: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I+K
                JLOC=J+K
                IF (ILOC > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC-MX
                    ELSE
                        EXIT L140
                    ENDIF
                ENDIF
                IF (JLOC > MY) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC-MY
                    ELSE
                        EXIT L140
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTD(K)
                SLOPE5SUM=SLOPE5SUM+LSLOPE*WLOC
                WEIGHTDUM=WEIGHTDUM+WLOC
            ENDDO L140
            SLOPE5SUM=SLOPE5SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J-1
        ILOC=I-1
        IF (ILOC < 1) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=MX
            ELSE
                SLOPE6SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (JLOC < 1) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=MY
            ELSE
                SLOPE6SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=ELOC-ECOMP
            SLOPE6SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L150: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I-K
                JLOC=J-K
                IF (ILOC < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC+MX
                    ELSE
                        EXIT L150
                    ENDIF
                ENDIF
                IF (JLOC < 1) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC+MY
                    ELSE
                        EXIT L150
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTD(K)
                SLOPE6SUM=SLOPE6SUM+LSLOPE*WLOC
                WEIGHTDUM=WEIGHTDUM+WLOC
            ENDDO L150
            SLOPE6SUM=SLOPE6SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J+1
        ILOC=I-1
        IF (JLOC > MY) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=1
            ELSE
                SLOPE7SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (ILOC < 1) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=MX
            ELSE
                SLOPE7SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=ELOC-ECOMP
            SLOPE7SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L160: DO  K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J+K
                ILOC=I-K
                IF (JLOC > MY) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC-MY
                    ELSE
                        EXIT L160
                    ENDIF
                ENDIF
                IF (ILOC < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC+MX
                    ELSE
                        EXIT L160
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTD(K)
                SLOPE7SUM=SLOPE7SUM+LSLOPE*WLOC
                WEIGHTDUM=WEIGHTDUM+WLOC
            ENDDO L160
            SLOPE7SUM=SLOPE7SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J-1
        ILOC=I+1
        IF (JLOC < 1) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=MY
            ELSE
                SLOPE8SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (ILOC > MX) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=1
            ELSE
                SLOPE8SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=ELOC-ECOMP
            SLOPE8SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L170: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J-K
                ILOC=I+K
                IF (JLOC < 1) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC+MY
                    ELSE
                        EXIT L170
                    ENDIF
                ENDIF
                IF (ILOC > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC-MX
                    ELSE
                        EXIT L170
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTD(K)
                SLOPE8SUM=SLOPE8SUM+LSLOPE*WLOC
                WEIGHTDUM=WEIGHTDUM+WLOC
            ENDDO L170
            SLOPE8SUM=SLOPE8SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        EXPOSE=(SLOPE1SUM+SLOPE2SUM+SLOPE3SUM+SLOPE4SUM+SLOPE5SUM+ &
        SLOPE6SUM+SLOPE7SUM+SLOPE8SUM)/(8.0*CELL_SIZE)
        RETURN
    END ! SUBROUTINE TOTAL_EXPOSURE(I,J,EXPOSE)
