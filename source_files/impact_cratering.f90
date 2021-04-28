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
	SUBROUTINE READ_CRATER_PARAMETERS
	USE ERODE_GLOBALS
	USE CRATER_GLOBALS
	IMPLICIT NONE
	CHARACTER (80) :: TEXTLINE
	INTEGER :: CRATER_PARAMS, ISOFTCRATER, USE_CRATER_EDGE_ABORT, DO_CORRECT_BIAS, DO_CENTRAL_PEAK, IFOLD
    INTEGER :: NRANDOM,IERROR, REAL_CRATERS_USE
    REAL(4) :: DIAM1,DIAM2
    INTEGER, DIMENSION(:), ALLOCATABLE :: RANDSEED
	    CRATER_PARAMS=93
		OPEN(CRATER_PARAMS,FILE='CRATERING_PARAMETERS.PRM',ACTION='READ')
	        !       **********************************************************************
        !
        !       Part 4:  Read in cratering parameters
        !
        !       **********************************************************************
        READ(CRATER_PARAMS,22) TEXTLINE
22      FORMAT(A)
        WRITE(*,22) TEXTLINE
        READ(CRATER_PARAMS,*) IMPACT_PROBABILITY
        READ(CRATER_PARAMS,*) IFOLD
        READ(CRATER_PARAMS,*) ISOFTCRATER
        READ(CRATER_PARAMS,*) USE_CRATER_EDGE_ABORT
        READ(CRATER_PARAMS,*) MINIMUM_HARD_DIAMETER
		READ(CRATER_PARAMS,*) DO_CORRECT_BIAS
        READ(CRATER_PARAMS,*) DO_CENTRAL_PEAK
        IF (ISOFTCRATER > 0) THEN
            IS_REGOLITH_CRATER=.TRUE.
        ELSE
            IS_REGOLITH_CRATER=.FALSE.
        ENDIF
        IF (USE_CRATER_EDGE_ABORT > 0) THEN
            CRATER_EDGE_ABORT = .TRUE.
        ELSE
            CRATER_EDGE_ABORT = .FALSE.
        ENDIF
	    IF (DO_CORRECT_BIAS>0) THEN
		    CORRECT_BIAS=.TRUE.
	    ELSE
		    CORRECT_BIAS=.FALSE.
        ENDIF
    	IF (DO_CENTRAL_PEAK>0) THEN
		    MAKE_CENTRAL_PEAK=.TRUE.
	    ELSE
		    MAKE_CENTRAL_PEAK=.FALSE.
	    ENDIF
		IF (CORRECT_BIAS) ALLOCATE (BWEIGHT(MX,MY),STAT=IERROR)
		IF (IERROR > 0) CALL ALLOCERROR()
        READ(CRATER_PARAMS,*) LARGE_CRATER_DEPTH_SCALE
        READ(CRATER_PARAMS,*) LARGE_CRATER_DEPTH_EXPONENT
        READ(CRATER_PARAMS,*) LARGE_CRATER_RIM_SCALE
        READ(CRATER_PARAMS,*) LARGE_CRATER_RIM_EXPONENT
        READ(CRATER_PARAMS,*) TRANSITION_DIAMETER
        READ(CRATER_PARAMS,*) SMALL_CRATER_DEPTH_SCALE
        READ(CRATER_PARAMS,*) SMALL_CRATER_DEPTH_EXPONENT
        READ(CRATER_PARAMS,*) SMALL_CRATER_RIM_SCALE
        READ(CRATER_PARAMS,*) SMALL_CRATER_RIM_EXPONENT
        READ(CRATER_PARAMS,*) LARGE_CRATER_SHAPE_SCALE
        READ(CRATER_PARAMS,*) LARGE_CRATER_SHAPE_EXPONENT
        READ(CRATER_PARAMS,*) SMALL_CRATER_SHAPE_SCALE
        READ(CRATER_PARAMS,*) SMALL_CRATER_SHAPE_EXPONENT
        READ(CRATER_PARAMS,*) CRATER_FREQUENCY_EXPONENT
        READ(CRATER_PARAMS,*) FREQUENCY_CUTOFF_SCALING
        READ(CRATER_PARAMS,*) SMALLEST_POSSIBLE_CRATER
        READ(CRATER_PARAMS,*) SMALLEST_MODELED_CRATER
        READ(CRATER_PARAMS,*) LARGEST_MODELED_CRATER
        READ(CRATER_PARAMS,*) EJECTA_THICKNESS_VARIABILITY
        READ(CRATER_PARAMS,*) NOISESD
        READ(CRATER_PARAMS,*) INHERITANCE_PARAMETER
        READ(CRATER_PARAMS,*) MAXIMUM_RIM_GRADIENT
        READ(CRATER_PARAMS,*) EJECTA_FRACTION_RETAINED
        READ(CRATER_PARAMS,*) INHERIT_EXPONENT
		READ(CRATER_PARAMS,*) PEAK_HEIGHT_MULT
        READ(CRATER_PARAMS,*) PEAK_HEIGHT_EXPONENT
        READ(CRATER_PARAMS,*) PEAK_DIAMETER_MULT
        READ(CRATER_PARAMS,*) PEAK_DIAMETER_EXPONENT
        READ(CRATER_PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE		
        READ(CRATER_PARAMS,*) REAL_CRATERS_USE
        READ(CRATER_PARAMS,*) RADIUS_MAX_INHERIT
        READ(CRATER_PARAMS,*) RADIUS_MAX_USE
        READ(CRATER_PARAMS,*) MINIMUM_REAL_CRATER_DIAMETER
        READ(CRATER_PARAMS,*) MAXIMUM_REAL_CRATER_DIAMETER
        READ(CRATER_PARAMS,*) COSINE_POWER		
        USE_REAL_CRATERS = .FALSE.
        IF ((REAL_CRATERS_USE > 0).AND.MODEL_IMPACT_CRATERING) THEN
			USE_REAL_CRATERS = .TRUE.
            OPEN(105,FILE='REAL_CRATERS.TXT',ACTION='READ')
            READ(105,207) CRATER_DATABASE_LOCATION
207 FORMAT(A160)
            CRATER_DATABASE_LOCATION = TRIM(CRATER_DATABASE_LOCATION)
            NUMBER_REAL_CRATERS = 1
            DO
                READ(105,207,END=205) CRATERFILENAMES(NUMBER_REAL_CRATERS)
                READ(105,*,END=205) DIAM1,DIAM2
                CRATERDIAM(NUMBER_REAL_CRATERS)=1000.0*(DIAM1+DIAM2)/2.0
                CRATERFILENAMES(NUMBER_REAL_CRATERS)=TRIM(CRATERFILENAMES(NUMBER_REAL_CRATERS))
                NUMBER_REAL_CRATERS=NUMBER_REAL_CRATERS+1
            ENDDO
205         CONTINUE
            CLOSE(105)
            IF (NUMBER_REAL_CRATERS <= 1) USE_REAL_CRATERS=.FALSE.

        ENDIF
        IF (SMALLEST_MODELED_CRATER < SMALLEST_POSSIBLE_CRATER) SMALLEST_MODELED_CRATER=SMALLEST_POSSIBLE_CRATER
		CLOSE(CRATER_PARAMS)
        ALLOCATE(CRATER_RESET(MX,MY))
        CRATER_RESET=.FALSE.
		        !       **********************************************************************
        !
        !        Write cratering parameters
        !
        !       **********************************************************************
        IF (MODEL_IMPACT_CRATERING.OR.DO_EVENTS) THEN
            WRITE(OUTHIST,1550) IMPACT_PROBABILITY
            1550 FORMAT('***********DOING IMPACT CRATERING*************',/, &
            ' PROBABILITY OF IMPACT EVENT=',G12.5)
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                IF(IS_REGOLITH_CRATER) THEN
                    WRITE(OUTHIST,1551)
                    1551 FORMAT(' CRATERS AND EJECTA ARE MODELED AS REGOLITH')
                ELSE
                    WRITE(OUTHIST,1552)
1552                FORMAT(' CRATERS AND EJECTA ARE MODELED AS BEDROCK')
                    WRITE(OUTHIST,2552) MINIMUM_HARD_DIAMETER
2552                FORMAT('CRATERS LESS THAN DIAMETER=',G12.5,' ARE MODELED AS REGOLITH')
                ENDIF
            ENDIF
            IF ((IS_X_PERIODIC).AND.(IS_Y_PERIODIC).AND.(IFOLD > 0)) &
            THEN
                DO_EJECTA_WRAPAROUND=.TRUE.
                WRITE(OUTHIST,4958)
                4958       FORMAT(' EJECTA MODELING FOLDS OVER PERIODIC LATERAL BOUNDS ')
            ELSE
                DO_EJECTA_WRAPAROUND=.FALSE.
            ENDIF
            WRITE(OUTHIST,4841) LARGE_CRATER_DEPTH_SCALE,LARGE_CRATER_DEPTH_EXPONENT, &
            LARGE_CRATER_RIM_SCALE,LARGE_CRATER_RIM_EXPONENT,TRANSITION_DIAMETER
            4841   FORMAT(' LARGE_CRATER_DEPTH_SCALE=',G12.5,' LARGE_CRATER_DEPTH_EXPONENT=', &
            G12.5,/,' LARGE_CRATER_RIM_SCALE=',G12.5, &
            ' LARGE_CRATER_RIM_EXPONENT=',G12.5,' SIMPLE-COMPLEX TRANSITION DIAMETER=',&
            G12.5)
           ! WRITE(OUTHIST,4821)NIMPACTS,DIFFINTERVAL,NOISEINTERVAL
            !4821   FORMAT(' NIMPACTS=',I8,/,&
            !' DIFFINTERVAL=',I8,' NOISEINTERVAL=',I8)
            WRITE(OUTHIST,4822)CRATER_FREQUENCY_EXPONENT,FREQUENCY_CUTOFF_SCALING, &
            SMALLEST_MODELED_CRATER,LARGEST_MODELED_CRATER
            4822   FORMAT(' CRATER_FREQUENCY_EXPONENT=',G12.5,' FREQUENCY_CUTOFF_SCALING=',G12.5,/, &
            ' SMALLEST_MODELED_CRATER=',G12.5, &
            ' LARGEST_MODELED_CRATER=',G12.5)
            WRITE(OUTHIST,4823)EJECTA_THICKNESS_VARIABILITY,NOISESD
            4823   FORMAT(' EJECTA_THICKNESS_VARIABILITY=',G12.5,' NOISESD=',G12.5)
            WRITE(OUTHIST,4824)INHERITANCE_PARAMETER
            4824   FORMAT(' INHERITANCE_PARAMETER=',G12.5)
            WRITE(OUTHIST,4825) EJECTA_FRACTION_RETAINED
4825        FORMAT(' FRACTION OF EXCAVATED BOWL RETAINED AS EJECTA (NOMINALLY 1.0)',G12.5)
            IF (CRATER_EDGE_ABORT) THEN
                WRITE(OUTHIST,1234)
1234            FORMAT('IF SIMULATION USES FIXED BOUNDARIES, CRATERS BEYOND EDGES ARE NOT MODELED')
            ENDIF
            IF (CORRECT_BIAS) THEN
                WRITE(OUTHIST,2234) 
2234            FORMAT('CRATER MORPHOLOGY IS ADJUSTED FOR ZERO NET ELEVATION CHANGE')
            ELSE
                WRITE(OUTHIST,3234)
3234            FORMAT('CRATER MORPHOLOGY IS NOT ADJUSTED FOR ZERO NET ELEVATION CHANGE')
            ENDIF
            IF (MAKE_CENTRAL_PEAK) THEN
                WRITE(OUTHIST, 1235) PEAK_HEIGHT_MULT,PEAK_HEIGHT_EXPONENT,PEAK_DIAMETER_MULT,PEAK_DIAMETER_EXPONENT                
1235            FORMAT('CRATER CENTRAL PEAKS ARE MODELED, PARAMETERS GOVERNING PEAK HEIGHT:',/, &
                '   PEAK_HEIGHT_MULT=',G12.5,' PEAK_HEIGHT_EXPONENT=',G12.5,/, &
                'PARAMETERS GOVERNING PEAK DIAMETER:',/, &
                '   PEAK_DIAMETER_MULT=',G12.5,' PEAK_DIMATER_EXPONENT=',G12.5)
            ELSE
                WRITE(OUTHIST,1236)
1236            FORMAT('CENTRAL PEAKS ARE NOT MODELED')
            ENDIF
            IF (USE_REAL_CRATERS) THEN
                WRITE(OUTHIST,1237) MINIMUM_REAL_CRATER_DIAMETER,MAXIMUM_REAL_CRATER_DIAMETER
1237            FORMAT('CRATERS BETWEEN DIAMETER ',G12.5,' AND ',G12.5, &
                     ' ARE MODELED USING ACTUAL FRESH MARTIAN CRATER DTMS')
                WRITE(OUTHIST,1238) RADIUS_MAX_INHERIT,RADIUS_MAX_USE,COSINE_POWER
1238            FORMAT('PARAMETERS GOVERNING SUPERPOSITION OF REAL CRATERS ON SIMULATION DOMAIN:',/, &
                     'RADIUS_MAX_INHERIT=',G12.5,' RADIUS_MAX_USE=',G12.5,' COSINE_POWER=',G12.5)
            ENDIF
            IF (USE_REAL_CRATERS) THEN
                WRITE(OUTHIST,1206) RADIUS_MAX_INHERIT,RADIUS_MAX_USE,MINIMUM_REAL_CRATER_DIAMETER, &
                    MAXIMUM_REAL_CRATER_DIAMETER,COSINE_POWER
                WRITE(OUTHIST,1205) NUMBER_REAL_CRATERS
1206        FORMAT(' MAXIMUM INHERIT RADIUS=',G12.5,' MAXIMUM RADIUS USED=',G12.5,/, &
            ' MINIMUM DIAMETER OF REAL CRATERS=',G12.5,' MAXIMUM DIAMETER OF REAL CRATERS=',G12.5,' COSINE POWER=',G12.5)            
1205        FORMAT('NUMBER OF REAL CRATERS=',I6)
            ENDIF
            IF (EJECTA_THICKNESS_VARIABILITY > 0.0) THEN
                RANDOM_EJECTA_THICKNESS=.TRUE.
                !       **********************************************************************
                !
                !        ******* Renormalize ejecta_thickness_variability for a lognormal distribution of
                !                deposit depths with mean unity
                !
               !       **********************************************************************
                   EJECTA_THICKNESS_VARIABILITY=EJECTA_THICKNESS_VARIABILITY/SQRT(EXP(1.0E0)*(EXP(1.0E0)-1.0E0))
            ELSE
                RANDOM_EJECTA_THICKNESS=.FALSE.
            ENDIF
            IF (NOISESD > 0.0) THEN
                MICRONOISE=.TRUE.
            ELSE
                MICRONOISE=.FALSE.
            ENDIF
            ITOTHITS=0
            ALINV=1.0/CRATER_FREQUENCY_EXPONENT
            !       **********************************************************************
            !
            !       ******* Determine the maximum range and values of xmin,xmax,ymin,ymax
            !               using the largest crater size
            !
            !       **********************************************************************
            CCON=1.0/(1.0/LARGEST_MODELED_CRATER**CRATER_FREQUENCY_EXPONENT-1.0/SMALLEST_MODELED_CRATER**CRATER_FREQUENCY_EXPONENT)
            BCON=-CCON/SMALLEST_MODELED_CRATER**CRATER_FREQUENCY_EXPONENT
            DIAMETER=LARGEST_MODELED_CRATER
            RADIUS=DIAMETER/2.0
            CRATER_DEPTH=LARGE_CRATER_DEPTH_SCALE*DIAMETER**LARGE_CRATER_DEPTH_EXPONENT
            RIM_HEIGHT=LARGE_CRATER_RIM_SCALE*DIAMETER**LARGE_CRATER_RIM_EXPONENT
            INTERIOR_SHAPE_EXPONENT=LARGE_CRATER_SHAPE_SCALE*DIAMETER**LARGE_CRATER_SHAPE_EXPONENT
            EXTERIOR_SHAPE_EXPONENT=2.0-(RIM_HEIGHT/(EJECTA_FRACTION_RETAINED*((RIM_HEIGHT-CRATER_DEPTH)/2.0+ &
                    CRATER_DEPTH/(INTERIOR_SHAPE_EXPONENT+2))))
            MAXRANGE=AMIN1(MX*CELL_SIZE,RADIUS*(1.0/0.01)**(1.0/(EXTERIOR_SHAPE_EXPONENT-1.0)))
            !       **********************************************************************
            !
            !        ******* The total playing field size is scaled so that all significant
            !                depostion from a 50 km crater off the target field will be
            !                accounted for
            !
            !       **********************************************************************
            CALL FIND_MODIFICATION_RANGE()
        ENDIF
        CALL RANDOM_SEED(SIZE=NRANDOM)
        ALLOCATE(RANDSEED(NRANDOM),STAT=IERROR)
		IF (IERROR > 0) CALL ALLOCERROR()
        RANDSEED=ISEED
        CALL RANDOM_SEED(PUT=RANDSEED)
        DEALLOCATE (RANDSEED)
        !       **********************************************************************
	RETURN
	END !SUBROUTINE READ_CRATER_PARAMETERS

    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DO_IMPACT_CRATERING(IS_RANDOM,CRATER_DIAMETER,X_LOCATION,Y_LOCATION)
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        IMPLICIT NONE
        LOGICAL :: IS_RANDOM, DO_ABORT
        LOGICAL :: CENTRAL_SEDIMENT_COVER
        REAL*4 CRATER_DIAMETER,X_LOCATION,Y_LOCATION
        REAL*4 :: SEDIMENT_COVER_HEIGHt
        INTEGER ICNTR,JCNTR
        !!!!!!!!write(outhist,4833)
4833    format('crater called')        
        CENTRAL_SEDIMENT_COVER = .TRUE.
        SEDIMENT_COVER_HEIGHT=0.05
        !     **********************************************************************
        !       This routine creates an impact crater in a random location with a
        !         size determined by a production function and geometry governed by
        !         a number of input parameters.  the model is documeted in:
        !          forsberg-taylor, n. k., howard, a. d., and craddock, r. a., 2004,
        !            j. geophys. res., 109, e05002, doi:10.1029/2004je002242
        !         and
        !          howard, a. d, 2007, geomorphology, 91, 332-363
        !  CALLS: GET_CRATER_SIZE, FIND_IMPACT_SITE, FIND_REFERENCE_ELEVATION, CREATE_CRATER
        !     **********************************************************************
            COUNT2=0
            COUNT10=0
        !     **********************************************************************
        !
        !     ******** Generate a random crater size from the assumed population
        !              distribution
        !
        !     **********************************************************************
            CALL GET_CRATER_SIZE(IS_RANDOM,CRATER_DIAMETER)
        !     **********************************************************************
        !
        !     ******** Where in the total playing field did it hit?
        !
        !     **********************************************************************
            DO_ABORT = .FALSE.
            CALL FIND_IMPACT_SITE(IS_RANDOM,X_LOCATION,Y_LOCATION,DO_ABORT,CRATER_DIAMETER,ICNTR,JCNTR)
            IF (DO_ABORT) RETURN
        !     **********************************************************************
        !
        !     ******** Determine the average elevation of the impactor footprint
        !
        !     **********************************************************************
        CALL FIND_REFERENCE_ELEVATION()
        !     **********************************************************************
        !
        !     ******** Do the excavation and deposition
        !
        !     **********************************************************************
        CALL CREATE_CRATER(ICNTR,JCNTR)
        IF (FLUVIAL_AND_SLOPE_MODELING.AND.DO_SEDIMENT_TRANSPORT.AND.DO_SEDIMENT_ROUTING &
            .AND.CENTRAL_SEDIMENT_COVER.AND.(.NOT.COMPLETE_RUNOFF)) THEN
        ENDIF
        IF (IS_RANDOM) TOTAL_RANDOM_CRATERS = TOTAL_RANDOM_CRATERS+1
        RETURN
    END ! SUBROUTINE DO_IMPACT_CRATERING
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE GET_CRATER_SIZE(IS_RANDOM,CRATER_DIAMETER)
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        IMPLICIT NONE
        LOGICAL :: IS_RANDOM
        REAL (4) :: TEMP
        REAL(4) :: TEMP1,DIAMRATIO
        REAL(4) :: RRAND, CRATER_DIAMETER
		REAL(4) :: RX,T1,T2,VOLI,VOLE1
        EXTERNAL RRAND
        !     **********************************************************************
        !
        !     ******* Determine a random impactor diameter from a negative
        !             power function cumulative diameter size.  Also determine
        !             parameters governing excavation and deposition shape
        !             functions and the maximum range of the impactor footprint
        !  MODIFIES: DIAMETER, RADIUS, CRATER_DEPTH, RIM_HEIGHT,INTERIOR_SHAPE_EXPONENT
        !            EXTERIOR_SHAPE_EXPONENT, DIAMRATIO, MAXRANGE
        !
        !     **********************************************************************
        IF (IS_RANDOM) THEN
            L110: DO
                TEMP=RRAND()
                IF (TEMP == 0.0) CYCLE L110
                DIAMETER=(CCON/(TEMP-BCON))**ALINV
                TEMP=RRAND()
                TEMP1=EXP(FREQUENCY_CUTOFF_SCALING*(SMALLEST_POSSIBLE_CRATER-DIAMETER))
                IF (TEMP < TEMP1) CYCLE L110
                EXIT L110
            ENDDO L110
        ELSE
            DIAMETER=CRATER_DIAMETER
        ENDIF
        IF (IS_RANDOM) THEN
            WRITE(*,7000) DIAMETER
            WRITE(OUTHIST,7000) DIAMETER
            WRITE(OUTCRATER,7000) DIAMETER
7000        FORMAT('RANDOM CRATERING, DIAM=',G12.5)
        ELSE
            WRITE(*,7001) DIAMETER
            WRITE(OUTHIST,7001) DIAMETER
            WRITE(OUTCRATER,7001) DIAMETER
7001        FORMAT('EVENT-BASED CRATERING, DIAM=',G12.5)
        ENDIF
        RADIUS=DIAMETER/2.0
        DIAMRATIO=DIAMETER/TRANSITION_DIAMETER
        IF (DIAMRATIO >= 1.0) THEN
            CRATER_DEPTH=LARGE_CRATER_DEPTH_SCALE*DIAMETER**LARGE_CRATER_DEPTH_EXPONENT
            RIM_HEIGHT=LARGE_CRATER_RIM_SCALE*DIAMETER**LARGE_CRATER_RIM_EXPONENT
            INTERIOR_SHAPE_EXPONENT=LARGE_CRATER_SHAPE_SCALE*DIAMETER**LARGE_CRATER_SHAPE_EXPONENT
        ELSE
            CRATER_DEPTH=SMALL_CRATER_DEPTH_SCALE*DIAMETER**SMALL_CRATER_DEPTH_EXPONENT
            RIM_HEIGHT=SMALL_CRATER_RIM_SCALE*DIAMETER**SMALL_CRATER_RIM_EXPONENT
            INTERIOR_SHAPE_EXPONENT=SMALL_CRATER_SHAPE_SCALE*DIAMETER**SMALL_CRATER_SHAPE_EXPONENT
        ENDIF
        RX=RADIUS*((CRATER_DEPTH-RIM_HEIGHT)/CRATER_DEPTH)**(1.0/INTERIOR_SHAPE_EXPONENT)
      IF ((MAKE_CENTRAL_PEAK).AND.(DIAMETER.GT.7000.0).AND.(DIAMETER.LE.100000.0)) THEN
          CENTRAL_PEAK_HEIGHT=PEAK_HEIGHT_MULT*DIAMETER**PEAK_HEIGHT_EXPONENT
          CENTRAL_PEAK_DIAMETER=PEAK_DIAMETER_MULT*DIAMETER**PEAK_DIAMETER_EXPONENT
          CENTRAL_PEAK_VOLUME=3.14159*(CENTRAL_PEAK_DIAMETER/2.0)**2*CENTRAL_PEAK_HEIGHT/3.0
!          WRITE(*,2984) CENTRAL_PEAK_HEIGHT,CENTRAL_PEAK_DIAMETER,CENTRAL_PEAK_VOLUME
2984 FORMAT('CPH=',G12.5,' CPD=',G12.5,' CPV=',G12.5)
      ELSE
          CENTRAL_PEAK_VOLUME=0.0
          CENTRAL_PEAK_HEIGHT=0.0
          CENTRAL_PEAK_DIAMETER=0.0
      ENDIF
      T1=(CRATER_DEPTH-RIM_HEIGHT)*RX**2/2.0
      T2=(CRATER_DEPTH*RX**(INTERIOR_SHAPE_EXPONENT+2.0))/((INTERIOR_SHAPE_EXPONENT+2.0) &
	           *RADIUS**INTERIOR_SHAPE_EXPONENT)
      VOLI=3.14159*2.0*(T1-T2)
      VOLE1=CRATER_DEPTH*(RADIUS**(INTERIOR_SHAPE_EXPONENT+2.0) &
	     -RX**(INTERIOR_SHAPE_EXPONENT+2.0))/((INTERIOR_SHAPE_EXPONENT+2.0)*RADIUS**INTERIOR_SHAPE_EXPONENT)
      EXTERIOR_SHAPE_EXPONENT=2.0+2.0*3.14159*RIM_HEIGHT*RADIUS**2/(VOLI-VOLE1-CENTRAL_PEAK_VOLUME)
!      WRITE(*,8514) RX,RADIUS,T1,T2,VOLI,VOLE1,M2
8514 FORMAT('RX=',G12.5,' R=',G12.5,' T1=',G12.5,' T2=',G12.5,'VOLI=',G12.5, &
          ' VOLE1=',G12.5,' M2=',G12.5)
      MAXRANGE=AMIN1(mx*CELL_SIZE,RADIUS*(0.001)**(1.0/(-EXTERIOR_SHAPE_EXPONENT)))
        RETURN
    END !SUBROUTINE GET_CRATER_SIZE
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE FIND_MODIFICATION_RANGE()
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
    !  MODIFIES: XMIN, XMAX, YMIN, YMAX
        IMPLICIT NONE
        IF (IS_X_PERIODIC) THEN
            XMIN=0.0
            XMAX=MX*CELL_SIZE
        ELSE
            XMIN=-MAXRANGE
            XMAX=MX*CELL_SIZE+MAXRANGE
        ENDIF
        IF (IS_Y_PERIODIC) THEN
            YMIN=0.0
            YMAX=MY*CELL_SIZE
        ELSE
            YMIN=-MAXRANGE
            YMAX=MY*CELL_SIZE+MAXRANGE
        ENDIF
        RETURN
    END ! SUBROUTINE FIND_MODIFICATION_RANGE
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE FIND_IMPACT_SITE(IS_RANDOM,X_LOCATION,Y_LOCATION, DO_ABORT,CRATER_DIAMETER,ICNTR,JCNTR)
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        IMPLICIT NONE
        LOGICAL :: IS_RANDOM, DO_ABORT
        REAL(4) :: X_LOCATION,Y_LOCATION, CRATER_DIAMETER
        REAL (4) :: XTEMP,YTEMP
        INTEGER :: ICNTR,JCNTR,IRANGE,I1,I2,I3
        REAL(4) :: RRAND
        EXTERNAL RRAND
        !     **********************************************************************
        !
        !     ****** Find a random location within the virtual target domain defined
        !            by xmin,xmax,ymin,ymax and determine the matrix indices of the
        !            footprint
        !  MODIFIES: XCENTER, YCENTER, ICNTR, JCNTR, IMIN, IMAX, JMIN, JMAX
        !            IMINR, IMAXR, JMINR, JMAXR
        !     **********************************************************************
        IF (IS_RANDOM) THEN
            XTEMP=RRAND()
            YTEMP=RRAND()
            XCENTER=XMIN+(XMAX-XMIN)*XTEMP
            YCENTER=YMIN+(YMAX-YMIN)*YTEMP
        ELSE
            XCENTER=X_LOCATION
            YCENTER=Y_LOCATION
        ENDIF
        ICNTR=XCENTER/CELL_SIZE+1
        JCNTR=YCENTER/CELL_SIZE+1
        IMIN=MAX0(1,INT(1+(XCENTER-MAXRANGE)/CELL_SIZE))
        IMAX=MIN0(MX,INT(1+(XCENTER+MAXRANGE)/CELL_SIZE))
        JMIN=MAX0(1,INT(1+(YCENTER-MAXRANGE)/CELL_SIZE))
        JMAX=MIN0(MY,INT(1+(YCENTER+MAXRANGE)/CELL_SIZE))
        !     **********************************************************************
        !
        !     ***** If periodic left and right boundaries are used, iminr and imaxr
        !           are the larger domain that allow boundary overlap
        !
        !     **********************************************************************
        IF (IS_X_PERIODIC) THEN
            IF (DO_EJECTA_WRAPAROUND) THEN
                IMINR=MAX0(-MX+1,INT(1+(XCENTER-MAXRANGE)/CELL_SIZE))
                IMAXR=MIN0(2*MX,INT(1+(XCENTER+MAXRANGE)/CELL_SIZE))
            ELSE
                IMINR=MAX0(ICNTR-MX/2+1,INT(1+(XCENTER-MAXRANGE)/CELL_SIZE))
                IMAXR=MIN0(ICNTR+MX/2,INT(1+(XCENTER+MAXRANGE)/CELL_SIZE))
            ENDIF
        ELSE
            IMINR=IMIN
            IMAXR=IMAX
        ENDIF
        IF (IS_Y_PERIODIC) THEN
            IF (DO_EJECTA_WRAPAROUND) THEN
                JMINR=MAX0(-MY+1,INT(1+(YCENTER-MAXRANGE)/CELL_SIZE))
                JMAXR=MIN0(2*MY,INT(1+(YCENTER+MAXRANGE)/CELL_SIZE))
            ELSE
                JMINR=MAX0(JCNTR-MY/2+1,INT(1+(YCENTER-MAXRANGE)/CELL_SIZE))
                JMAXR=MIN0(JCNTR+MY/2,INT(1+(YCENTER+MAXRANGE)/CELL_SIZE))
            ENDIF
        ELSE
            JMINR=JMIN
            JMAXR=JMAX
        ENDIF
        IF     &
        ((IMINR == -128).OR.(IMAXR == 256).OR.(JMINR == -128).OR. &
        (JMAXR == 256)) THEN
        ENDIF
        IF (CRATER_EDGE_ABORT.AND.DO_FLOW_BOUNDARIES) THEN
            IRANGE = INT(RADIUS/(CELL_SIZE)) + 4
            ! WRITE(OUTHIST,9984) RADIUS, CELL_SIZE, IRANGE, ICNTR,JCNTR
9984         FORMAT('R=',G13.6,' DX=',G13.6,' IRANGE=',I7,' ICNTR=',I7,' JCNTR=',I7)
            IF ((ICNTR<1).OR.(ICNTR>MX).OR.(JCNTR.LT.1).OR.(JCNTR.GT.MY)) DO_ABORT =.TRUE.
            IF ((ICNTR+IRANGE).GT.MX) DO_ABORT = .TRUE.
            IF ((ICNTR-IRANGE).LT.1) DO_ABORT = .TRUE.
            IF ((JCNTR+IRANGE).GT.MY) DO_ABORT = .TRUE.
            IF ((JCNTR-IRANGE).LT.1) DO_ABORT = .TRUE.
            IF (DO_ABORT) THEN
                WRITE(*,9986)
9986            FORMAT('***Out of Range***')
            ENDIF
        ENDIF
        RETURN
    END ! SUBROUTINE FIND_IMPACT_SITE
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE CREATE_CRATER(ICNTR,JCNTR)
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,ITEMP,JTEMP,ICNTR,JCNTR
        REAL(4) :: XCOMP,YCOMP,DCOMP,TEMP1,TEMP2,DEPOSITTHICK
        REAL(4) :: EUSE,LOGNORMAL_RANDOM_DEVIATE,TOTALVOL,BIASADD,TOTPOINTS
        REAL(4) :: LOCALGRAV,LOCALWORK,BWEIGHTSUM,BWEIGHTNUMB, MAXWEIGHTDIST
        REAL(4) :: NETCHANGE,NETPOINTS,DELELEV,INTERIOR_CHANGE,EXTERIOR_CHANGE
        LOGICAL (1) :: NOTDONE(IMMX,JMMX),NEWSEDIMENT_BASE(IMMX,JMMX),IS_SOFT
        !     **********************************************************************
        !
        !      ******  This subroutine does the actual work of excavating and
        !              depositing.  Small craters off of the target area are
        !              bypassed
        !  MODIFIES: ERODE_SLOPE, CFNE, REGOLITH, IS_ROCK_SURFACE, IS_SEDIMENT_COVERED
        !            ELEVATION, CUMULATIVE_CRATERING_CHANGE, SEDIMENT_BASE
        !            CRATERWORK, LOCALWORK, CRATERGRAV, LOCALGRAV
        !  CALLS: LOGNORMAL_RANDOM_DEVIATE
        !
        !     **********************************************************************
        IF ((IMAX >= IMIN).AND.(JMAX >= JMIN)) THEN
            TOTALVOL=0.0
            TOTPOINTS=0.0
            NETCHANGE=0.0
            NETPOINTS=0.0
            LOCALGRAV=0.0
            LOCALWORK=0.0
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                WRITE(OUTCRATER,120) DIAMETER &
                ,XCENTER,YCENTER,ITERATION,PRESENT_TIME
                120   FORMAT(3(G12.5,' '),I9,' ',G12.5)
            ELSE
                WRITE(OUTCRATER,1120) DIAMETER &
                ,XCENTER,YCENTER,ITERATION
                1120   FORMAT(3(G12.5,' '),I9)
            ENDIF
            WRITE (OUTCRATER,100) MAXRANGE, CRATER_DEPTH, RIM_HEIGHT,INTERIOR_SHAPE_EXPONENT, EXTERIOR_SHAPE_EXPONENT
100     FORMAT(' MAXRANGE=',G12.5,' CRATER_DEPTH=',G12.5,' RIM_HEIGHT=',G12.5, &
        ' INTERIOR_SHAPE_EXPONENT=',G12.5,' EXTERIOR_SHAPE_EXPONENT=',G12.5)
            !     **********************************************************************
            !
            !      ****** If we're here we model the impact.  Update counters.
            !
            !     **********************************************************************
            ITOTHITS=ITOTHITS+1
            IF (DIAMETER > 2.0) COUNT2=COUNT2+1
            IF (DIAMETER > 10.0) COUNT10=COUNT10+1
            NOISECALL=.TRUE.
		    BWEIGHT=0.0
            BWEIGHTSUM=0.0
            BWEIGHTNUMB=0.0
            MAXWEIGHTDIST=ABS(IMAXR-IMINR)*CELL_SIZE/2.0
            L2100: DO  JTEMP=JMINR,JMAXR
                M2100: DO ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
	                XCOMP=(ITEMP-1)*CELL_SIZE
                    YCOMP=(JTEMP-1)*CELL_SIZE
                    DCOMP=SQRT((XCOMP-XCENTER)**2+(YCOMP-YCENTER)**2)
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    IF (DCOMP.LE.MAXWEIGHTDIST) THEN
                       BWEIGHT(I,J)=1.0-DCOMP/MAXWEIGHTDIST
                       BWEIGHTSUM=BWEIGHTSUM+BWEIGHT(I,J)
                       BWEIGHTNUMB=BWEIGHTNUMB+1.0
                    ENDIF
                    ERODE_SLOPE(I,J)=0.0
                    NOTDONE(I,J)=.TRUE.
                    NEWSEDIMENT_BASE(I,J)=.FALSE.
                ENDDO M2100
            ENDDO L2100
            BWEIGHTSUM=BWEIGHTSUM/BWEIGHTNUMB
            !     **********************************************************************
            !
            !     ****** Cycle across the footprint area on the target matrix
            !
            !     **********************************************************************
            L100: DO JTEMP=JMINR,JMAXR
                M100: DO ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    XCOMP=(ITEMP-1)*CELL_SIZE
                    YCOMP=(JTEMP-1)*CELL_SIZE
                    DCOMP=SQRT((XCOMP-XCENTER)**2+(YCOMP-YCENTER)**2)
                    IF (DCOMP<(1.5*RADIUS)) CRATER_RESET(I,J)=.TRUE.
                    !     **********************************************************************
                    !
                    !     *****   dcomp is the radial  distance of the present location on
                    !             the footprint from the crater center.  We model the interior
                    !             and exterior of the crater separately
                    !
                    !     **********************************************************************
                    IF (DCOMP <= RADIUS) THEN
                        !     **********************************************************************
                        !
                        !     *****   Model the crater interior
                        !
                        !     **********************************************************************
                        CRATER_RESET(I,J)=.TRUE.
                        IF (USE_CRATER_EVENT) CRATER_EVENT(I,J)=RADIUS
                        DEPOSITTHICK=CRATER_DEPTH*(DCOMP/RADIUS)**INTERIOR_SHAPE_EXPONENT
                        !     **********************************************************************
                        !
                        !     *****   Random variability of crater elevation is presently allowed
                        !             only in areas of net deposition.  Perhaps it also ought to
                        !             be allowed on the lower crater rim
                        !
                        !     **********************************************************************
                        IF (RANDOM_EJECTA_THICKNESS) THEN
                            DEPOSITTHICK=DEPOSITTHICK &
                            *LOGNORMAL_RANDOM_DEVIATE(EJECTA_THICKNESS_VARIABILITY)
                        ENDIF
                        DEPOSITTHICK=DEPOSITTHICK+(RIM_HEIGHT-CRATER_DEPTH)
                        !     **********************************************************************
                        !
                        !     *****   The amount of inheritance of the original topography increases
                        !             from zero at the crater center (where the average elevation of
                        !             the footprint is used) to inheritance_parameter at the rim
                        !
                        !     **********************************************************************
                        EUSE=INHERITANCE_PARAMETER*(DCOMP/RADIUS)**INHERIT_EXPONENT
                        INTERIOR_CHANGE= (1.0-EUSE)*(ELEVBASE-ELEVATION(I,J))+ &
                        DEPOSITTHICK
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+  INTERIOR_CHANGE
                        CUMULATIVE_CRATER_EXCAVATION(I,J)=CUMULATIVE_CRATER_EXCAVATION(I,J)+INTERIOR_CHANGE
                        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                            IS_SEDIMENT_COVERED(I,J)=.FALSE.
                            NEWSEDIMENT_BASE(I,J)=.TRUE.
                            CFNE(I,J)=0.0
                            IF (IS_REGOLITH_CRATER.OR.(DIAMETER<MINIMUM_HARD_DIAMETER) )  THEN
                                IS_ROCK_SURFACE(I,J)=.FALSE.
                                REGOLITH(I,J)=INITIAL_REGOLITH_THICKNESS
                                IS_SOFT = .TRUE.

                            ELSE
                                IS_SOFT = .FALSE.
                                IS_ROCK_SURFACE(I,J)=.TRUE.
                                REGOLITH(I,J)=-ROCK_WEATHERING_RATE
                            ENDIF
                        ENDIF
                    ELSE
                        !     **********************************************************************
                        !
                        !     ******  Model the crater exterior.  The fractional inheritance of the
                        !             original topography goes from inheritance_parameter at the rim to 100%
                        !             at the far margins of deposition
                        !
                        !     **********************************************************************
                        IF (USE_CRATER_EVENT) THEN
                            IF (DCOMP<=(1.3*RADIUS)) CRATER_EVENT(I,J)=RADIUS
                        ENDIF
                        TEMP1=(RADIUS/DCOMP)**EXTERIOR_SHAPE_EXPONENT
                        TEMP2=AMIN1((1.0-INHERITANCE_PARAMETER),TEMP1)
                        DEPOSITTHICK=RIM_HEIGHT*TEMP1
                        !     **********************************************************************
                        !
                        !     ******  Determine the random variability of depositional amount.
                        !
                        !     **********************************************************************
                        IF (RANDOM_EJECTA_THICKNESS) THEN
                            DEPOSITTHICK=DEPOSITTHICK*LOGNORMAL_RANDOM_DEVIATE(EJECTA_THICKNESS_VARIABILITY)
                        ENDIF
                        EXTERIOR_CHANGE=(ELEVBASE-ELEVATION(I,J))*TEMP2  &
                        +DEPOSITTHICK
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+EXTERIOR_CHANGE
                        CUMULATIVE_EJECTA_DEPOSITION(I,J)=CUMULATIVE_EJECTA_DEPOSITION(I,J)+EXTERIOR_CHANGE
                        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                            IF (ABS(ERODE_SLOPE(I,J)) > 100.0*ROCK_WEATHERING_RATE*TIME_INCREMENT) THEN
                                IS_SEDIMENT_COVERED(I,J)=.FALSE.
                                NEWSEDIMENT_BASE(I,J)=.TRUE.
                                IF (IS_REGOLITH_CRATER.OR.(DIAMETER<MINIMUM_HARD_DIAMETER) )  THEN
                                    IS_ROCK_SURFACE(I,J)=.FALSE.
                                    REGOLITH(I,J)=INITIAL_REGOLITH_THICKNESS
                                    IS_SOFT = .TRUE.
                                ELSE
                                    IS_ROCK_SURFACE(I,J)=.TRUE.
                                    REGOLITH(I,J)=-ROCK_WEATHERING_RATE
                                    IS_SOFT = .FALSE.
                                ENDIF
                                CFNE(I,J)=0.0
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO M100
            ENDDO L100
            IF (USE_CRATER_EVENT) CRATER_EVENT(ICNTR,JCNTR)=-RADIUS
            IF (IS_SOFT) THEN
               WRITE(*,779)
779                             FORMAT ('***SOFT CRATER***')
            ELSE
               WRITE(*,778)
778                             FORMAT ('***HARD CRATER***')
            ENDIF
            L3100: DO JTEMP=JMINR,JMAXR
                M3100: DO ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    IF (NOTDONE(I,J)) THEN
                        TOTALVOL=TOTALVOL+ERODE_SLOPE(I,J)
                        TOTPOINTS=TOTPOINTS+1.0
                        NOTDONE(I,J)=.FALSE.
                    ENDIF
                ENDDO M3100
            ENDDO L3100
            IF (IS_X_PERIODIC.AND.IS_X_PERIODIC) THEN
                BIASADD=-TOTALVOL/TOTPOINTS
            !            write(*,3110) totalvol,totpoints,biasadd
3110        FORMAT(' OV=',G12.5,' T=',G12.5,' BA=',G12.5)
            ELSE
                BIASADD=0.0
            ENDIF
            L4100: DO  JTEMP=JMINR,JMAXR
                M4100: DO ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    NOTDONE(I,J)=.TRUE.
                ENDDO M4100
            ENDDO L4100
            L1100: DO  JTEMP=JMINR,JMAXR
                M1100: DO  ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    XCOMP=(ITEMP-1)*CELL_SIZE
                    YCOMP=(JTEMP-1)*CELL_SIZE
                    DCOMP=SQRT((XCOMP-XCENTER)**2+(YCOMP-YCENTER)**2)
                    IF (NOTDONE(I,J)) THEN
                        DELELEV=+ERODE_SLOPE(I,J)+BIASADD*BWEIGHT(I,J)/BWEIGHTSUM
                        ELEVATION(I,J)=ELEVATION(I,J)+DELELEV
                        CRATERWORK=CRATERWORK+DCOMP*DELELEV
                        LOCALWORK=LOCALWORK+DCOMP*DELELEV
                        CRATERGRAV=CRATERGRAV-0.5*ABS(DELELEV)
                        LOCALGRAV=LOCALGRAV-0.5*ABS(DELELEV)
                        CUMULATIVE_CRATERING_CHANGE(I,J)=CUMULATIVE_CRATERING_CHANGE(I,J) &
                        + DELELEV
                        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                            IF (NEWSEDIMENT_BASE(I,J)) SEDIMENT_BASE(I,J)=ELEVATION(I,J)
                        ENDIF
                        NETCHANGE=NETCHANGE+ERODE_SLOPE(I,J)+BIASADD
                        NETPOINTS=NETPOINTS+1.0
                        NOTDONE(I,J)=.FALSE.
                        IF (MODEL_GROUNDWATER) REFERENCE_ELEVATION(I,J)=ELEVATION(I,J)
                    ENDIF
                ENDDO M1100
            ENDDO L1100
            NETCHANGE=NETCHANGE/NETPOINTS
            IF ((MAKE_CENTRAL_PEAK).AND.(DIAMETER.GT.7000.0).AND.(DIAMETER.LE.100000.0)) &
			       CALL ADD_CENTRAL_PEAK()
        ELSE
            !WRITE(*,130)
            130   FORMAT(' MISS')
        ENDIF
        !WRITE(*,7552) LOCALWORK,LOCALGRAV
        7552  FORMAT(' LOCALWORK=',G12.5,' LOCALGRAV=',G12.5)
        RETURN
    END ! SUBROUTINE CREATE_CRATER
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE ADD_CENTRAL_PEAK()
   USE ERODE_GLOBALS
   USE CRATER_GLOBALS
   IMPLICIT NONE
   INTEGER I,J,ITEMP,JTEMP
   REAL*4 DCOMP,XCOMP,YCOMP, LOGNORMAL_RANDOM_DEVIATE
         DO ITEMP=IMINR,IMAXR
          DO JTEMP=JMINR,JMAXR
              I=ITEMP
              J=JTEMP
              XCOMP=(ITEMP-1)*CELL_SIZE
              YCOMP=(JTEMP-1)*CELL_SIZE
              DCOMP=SQRT((XCOMP-XCENTER)**2+(YCOMP-YCENTER)**2)

! **********************************************************************
!
! ******   IF THE LATERAL BOUNDARIES ARE PERIODIC, THEN ALLOW THE
!          BOUNDARY OVERLAP
!
! **********************************************************************
              IF (IS_X_PERIODIC) THEN
                  IF (ITEMP.LT.1) I=ITEMP+MX
                  IF (ITEMP.GT.MX) I=ITEMP-MX
              ENDIF
              IF (IS_Y_PERIODIC) THEN
                  IF (JTEMP.LT.1) J=JTEMP+MY
                  IF (JTEMP.GT.MY) J=JTEMP-MY
              ENDIF
              IF (DCOMP.LT.(CENTRAL_PEAK_DIAMETER/2.0)) THEN
                  ELEVATION(I,J)=ELEVATION(I,J)+(CENTRAL_PEAK_HEIGHT*(1.0-2.0*DCOMP/CENTRAL_PEAK_DIAMETER)) &
                      * LOGNORMAL_RANDOM_DEVIATE(EJECTA_THICKNESS_VARIABILITY)
              ENDIF
		  ENDDO
		ENDDO
        RETURN
        END ! SUBROUTINE ADD_CENTRAL_PEAK

    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE FIND_REFERENCE_ELEVATION()
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        IMPLICIT NONE
        REAL(4) :: TEMP,WEIGHT,XCOMP,YCOMP,DCOMP
        INTEGER :: I,J,ITEMP,JTEMP
        !     **********************************************************************
        !
        !     ***** Finds the average elevation of the footprint. The only trick
        !           is that the locations lying within the new crater interior are
        !           given greatest weighting and that points outside the crater rim
        !           are given a weighting in determining the average elevation in
        !           proportion to the amount of deposition that will occur
        !  MODIFIES: ELEVBASE
        !
        !     **********************************************************************
        ELEVBASE=0.0
        WEIGHT=0.0
        IF ((IMAX >= IMIN).AND.(JMAX >= JMIN)) THEN
            L100: DO  JTEMP=JMINR,JMAXR
                M100: DO ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    XCOMP=(I-1)*CELL_SIZE
                    YCOMP=(J-1)*CELL_SIZE
                    DCOMP=SQRT((XCOMP-XCENTER)**2+(YCOMP-YCENTER)**2)
                    IF (DCOMP <= RADIUS) THEN
                        ELEVBASE=ELEVBASE+ELEVATION(I,J)
                        WEIGHT=WEIGHT+1.0
                    ELSE
                        TEMP=(RADIUS/DCOMP)**EXTERIOR_SHAPE_EXPONENT
                        ELEVBASE=ELEVBASE+ELEVATION(I,J)*TEMP
                        WEIGHT=WEIGHT+TEMP
                    ENDIF
                ENDDO M100
            ENDDO L100
            ELEVBASE=ELEVBASE/WEIGHT
        ENDIF
        RETURN
    END !  SUBROUTINE FIND_REFERENCE_ELEVATION
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE CREATE_SEDIMENT_COVER(SEDIMENT_COVER_HEIGHT)
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        IMPLICIT NONE
        REAL* 8 :: SEDIMENT_COVER_HEIGHT, ELEVATION_REFERENCE,DCOMP,XCOMP,YCOMP
        INTEGER :: I,J,ITEMP,JTEMP,ICNTR,JCNTR
        IF ((IMAX >= IMIN).AND.(JMAX >= JMIN)) THEN
        ICNTR=XCENTER/CELL_SIZE+1
        JCNTR=YCENTER/CELL_SIZE+1
        ELEVATION_REFERENCE=ELEVATION(ICNTR,JCNTR)
             L2100: DO  JTEMP=JMINR,JMAXR
                M2100: DO ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                ENDDO M2100
             ENDDO L2100
            !     **********************************************************************
            !
            !     ****** Cycle across the footprint area on the target matrix
            !
            !     **********************************************************************
            L100: DO JTEMP=JMINR,JMAXR
                M100: DO ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    XCOMP=(ITEMP-1)*CELL_SIZE
                    YCOMP=(JTEMP-1)*CELL_SIZE
                    DCOMP=SQRT((XCOMP-XCENTER)**2+(YCOMP-YCENTER)**2)
                    !     **********************************************************************
                    !
                    !     *****   dcomp is the radial  distance of the present location on
                    !             the footprint from the crater center.  We model the interior
                    !             and exterior of the crater separately
                    !
                    !     **********************************************************************
                    IF (DCOMP <= RADIUS) THEN
                        !     **********************************************************************
                        !
                        !     *****   Only look at the crater interior
                        !
                        !     **********************************************************************
                       IF ((ELEVATION(I,J)-ELEVATION_REFERENCE)<=(SEDIMENT_COVER_HEIGHT*CRATER_DEPTH)) THEN
                           IS_SEDIMENT_COVERED(I,J) = .TRUE.
                           IS_ROCK_SURFACE(i,J) = .FALSE.
                           SEDIMENT_BASE(I,J)=ELEVATION(I,J)-1.0
                       ENDIF
                    ENDIF
                ENDDO M100
            ENDDO L100
            ENDIF
            RETURN
    END !  SUBROUTINE CREATE_SEDIMENT_COVER
!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       SUBROUTINE CREATE_REAL_CRATER(IS_RANDOM,CRATER_DIAMETER,X_LOCATION,Y_LOCATION)
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        IMPLICIT NONE
        LOGICAL :: IS_RANDOM, DO_ABORT,INORMAL,JNORMAL,done_once,do_this
        INTEGER :: I,ISELECT,MXC,MYC,J, STATUS,ICNTR,JCNTR,mpoints,npoints,newnumb
        INTEGER :: III,JJJ,ICRATCENT,JCRATCENT,IERROR,ITEMP,JTEMP,RNUMBER
        INTEGER :: II,JJ,IW,IE,JS,JN,IIRANGE
        REAL*4 :: ESUM,ENUM,LXMIN,LXMAX,LYMIN,LYMAX
        REAL*4 DDIF,AVAL,BVAL,DCOMP,ELEVCOMP,MEANELEV,T1,T2,CLOSEST_DIAM
        REAL*4 CRATER_DIAMETER,X_LOCATION,Y_LOCATION,RRAND,euse,depositthick,interior_change
        REAL*4 :: WEIGHT,XCOMP,YCOMP,TEMP,ZMIN,ZMAX,RADMAX,LOCAL_MAX_USE
        real*4 :: sumx,sumx2,sumy,sumy2,sumxy,sumzx,sumzy,sumz,aterm,bterm,cterm
        real*4 ::  term1,rint,rslope,rloc,rrim,eref,dweight,costerm,rinherit,term2
        real*4 :: router,rinner,interiorelev,newmean,ix,jy,scale_distance,xicent,yjcent
        REAL*4, ALLOCATABLE, DIMENSION(:,:) :: CRATERELEVS,newel
        EXTERNAL RRAND
        CHARACTER*160 :: FILETYPE,CRATERFILE
        IF (IS_RANDOM) THEN
            WRITE(*,401)
401 FORMAT('RETURNING - IS A RANDOM CRATER')
            RETURN
        ENDIF
        INORMAL=.TRUE.
        JNORMAL=.TRUE.
        done_once=.false.
        TEMP=RRAND()
        IF (TEMP.GE.0.25) THEN
            IF (TEMP.LT.0.5) THEN
                INORMAL=.FALSE.
            ELSE
                IF (TEMP.LT.0.75) THEN
                    JNORMAL=.FALSE.
                ELSE
                    INORMAL=.FALSE.
                    JNORMAL=.FALSE.
                ENDIF
            ENDIF
        ENDIF
        CLOSEST_DIAM=1.0E+25
        DO I=1, NUMBER_REAL_CRATERS
            DDIF=ABS(CRATER_DIAMETER-CRATERDIAM(I))
            IF (DDIF<CLOSEST_DIAM) THEN
                CLOSEST_DIAM=DDIF
                ISELECT=I
            ENDIF
        ENDDO
        WRITE(OUTCRATER,410) ISELECT,CRATER_DIAMETER,CRATERDIAM(ISELECT)
        WRITE(*,410) ISELECT,CRATER_DIAMETER,CRATERDIAM(ISELECT)
        WRITE(OUTHIST,410) ISELECT,CRATER_DIAMETER,CRATERDIAM(ISELECT)
410     FORMAT('REAL CRATER, ISELECT=',I6,' Crater diameter=',g12.5,' Modeled diameter=',g12.5)
        IF (CLOSEST_DIAM>15000.0) THEN
            WRITE(*,321) CRATER_DIAMETER,CRATERDIAM(ISELECT),CLOSEST_DIAM
            321         FORMAT('CLOSEST CRATER DIFFERENCE TOO LARGE',&
                        'CRATER SIZE=',G12.5,' SELECTED DIAMETER=',G12.5,'DDIF=',G13.6)
            RETURN
        ENDIF
        write(outcrater,739) iteration, present_time
        write(*,739) iteration, present_time
        write(outhist,739) iteration, present_time
739  format('Real crater at iteration=',i9,' time=',g13.6)
        RADIUS=CRATERDIAM(ISELECT)/2.0
        DIAMETER=CRATERDIAM(ISELECT)
        CRATERFILE=TRIM(CRATER_DATABASE_LOCATION)//TRIM(CRATERFILENAMES(ISELECT))
        WRITE(OUTCRATER,415) CRATER_DATABASE_LOCATION
        WRITE(OUTCRATER,415)TRIM(CRATERFILENAMES(ISELECT))
        WRITE(*,415) CRATERFILE
        WRITE(OUTCRATER,415) CRATERFILE
415     FORMAT(A160)
        OPEN(103,FILE=CRATERFILE,ERR=399,ACTION='read')
        READ(103,200) FILETYPE
200     FORMAT(A)
        READ(103,*) MXC,MYC
        ALLOCATE(CRATERELEVS(MXC,MYC),newel(mxc,myc),STAT=IERROR)
        READ(103,*) LXMIN,LXMAX
        READ(103,*) LYMIN,LYMAX
        READ(103,*) ZMIN,ZMAX
        WRITE(outcrater,411) MXC,MYC,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
411     FORMAT('MXC=',I6,' MYC=',I6,' XMIN=',G12.5,' XMAX=',G12.5 &
        ,' YMIN=',G12.5,' YMAX=',G12.5,' ZMIN=',G12.5,' ZMAX=',G12.5)
        IF (JNORMAL) THEN
            DO J=MYC,1,-1
                IF (INORMAL) THEN
                    READ(103,*) (CRATERELEVS(I,J),I=1,MXC)
                ELSE
                    READ(103,*) (CRATERELEVS(I,J),I=MXC,1,-1)
                ENDIF
            ENDDO
        ELSE
            DO J=1,MYC
                IF (INORMAL) THEN
                    READ(103,*) (CRATERELEVS(I,J),I=1,MXC)
                ELSE
                    READ(103,*) (CRATERELEVS(I,J),I=MXC,1,-1)
                ENDIF
            ENDDO
        ENDIF
        CLOSE(103)
        LOCAL_MAX_USE=MIN(RADIUS_MAX_USE,(1000.0*(LXMAX-LXMIN)/(2.0*RADIUS)))
        WRITE(outcrater,405) RADIUS,LOCAL_MAX_USE
405     FORMAT('CRATER RADIUS=',G13.6,' LOCAL MAX=',G13.6)
        newel=0.0
        router=local_max_use*radius/cell_size
        rinner=router-0.1*radius/cell_size
        rinherit=radius_max_inherit*radius/cell_size
        scale_distance=router
        rrim=radius/cell_size
        write(outcrater,767) rinner,router,rinherit,scale_distance
767     format(' rinner=',g12.5,' router=',g12.5,' rinherit=',g12.5,' scale_distance=',g12.5)
        ICRATCENT=MXC/2+1
        JCRATCENT=MYC/2+1
        ICNTR=X_LOCATION/CELL_SIZE+1
        JCNTR=Y_LOCATION/CELL_SIZE+1
        IF ((ICNTR<1).OR.(ICNTR>MX).OR.(JCNTR.LT.1).OR.(JCNTR.GT.MY)) THEN
                WRITE(*,9986)
9986            FORMAT('***Crater Out of Range***')
                RETURN
        ENDIF
        WRITE(outcrater,406) ICRATCENT,JCRATCENT,ICNTR,JCNTR
406     FORMAT('ICRATCENT=',I6,' JCRATCENT=',I6,' ICNTR=',I6,' JCNTR=',I6)
        IF ((ICNTR>MX).OR.(ICNTR<1).OR.(JCNTR>MY).OR.(JCNTR<1)) RETURN
        ELEVCOMP=0.0
        WEIGHT=0.0
        IMIN=AMAX0(ICNTR-MXC/2,1)
        IMAX=AMIN0(ICNTR+MXC/2,MX)
        JMIN=AMAX0(JCNTR-MYC/2,1)
        JMAX=AMIN0(JCNTR+MYC/2,MY)
        WRITE(outcrater,407) IMIN,IMAX,JMIN,JMAX
        AVAL=1.0/(RADIUS*(RADIUS_MAX_INHERIT-LOCAL_MAX_USE))
        BVAL=-AVAL*LOCAL_MAX_USE*RADIUS
407     FORMAT('IMIN=',I6,' IMAX=',I6,' JMIN=',I6,' JMAX=',I6)
        IF ((.NOT.DONE_ONCE).AND.(DIAMETER>50000.0)) THEN  !DIAGNOSTIC
!            OPEN (44,FILE='COS_DIAGNOSTICS.TXT',ACTION='WRITE')
!            WRITE(44,834) RADIUS, RADIUS_MAX_INHERIT,LOCAL_MAX_USE,RINHERIT,ROUTER
834         FORMAT('RADIUS=',G13.6,' RADIUS_MAX_INHERIT=',G12.5,' LOCAL_MAX_USE=',G12.5,' RINHERIT=',G13.6,' ROUTER=',G12.5)
        ENDIF
            L100: DO  JTEMP=JMIN,JMAX
                M100: DO ITEMP=IMIN,IMAX
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    XCOMP=(I-1)*CELL_SIZE
                    YCOMP=(J-1)*CELL_SIZE
                    DCOMP=SQRT((XCOMP-X_LOCATION)**2+(YCOMP-Y_LOCATION)**2)
                    IF (DCOMP <= RADIUS_MAX_INHERIT*RADIUS) THEN
                        IF (USE_CRATER_EVENT) CRATER_EVENT(I,J)=RADIUS
                        ELEVCOMP=ELEVCOMP+ELEVATION(I,J)
                         CRATER_RESET(I,J)=.TRUE.
                        WEIGHT=WEIGHT+1.0
                    ELSE
                        IF (DCOMP<=LOCAL_MAX_USE*RADIUS) THEN
                            IF (USE_CRATER_EVENT) THEN
                                IF (DCOMP<=(1.3*RADIUS)) CRATER_EVENT(I,J)=RADIUS
                            ENDIF
                            TEMP=COS(3.14159*((DCOMP-RADIUS_MAX_INHERIT*RADIUS)/ &
                            ((LOCAL_MAX_USE-RADIUS_MAX_INHERIT)*RADIUS))/2.0)**COSINE_POWER
                            IF ((.NOT.DONE_ONCE).AND.(DIAMETER>50000.0).AND.(J==JCNTR)) THEN
                            !WRITE(44,835) I,DCOMP,TEMP
835                         FORMAT(I7,' ',G12.5,' ',G12.5)
                            ENDIF
                            ELEVCOMP=ELEVCOMP+ELEVATION(I,J)*TEMP
                            WEIGHT=WEIGHT+TEMP
                        ENDIF
                    ENDIF
                ENDDO M100
            ENDDO L100
        ELEVCOMP=ELEVCOMP/WEIGHT
        RADMAX=(RADIUS*LOCAL_MAX_USE)
        MEANELEV=0.0
        RNUMBER=0
        L200: DO J=1,MYC
            M200: DO I=1,MXC
                DCOMP=SQRT(((I-ICRATCENT)*CELL_SIZE)**2+((J-JCRATCENT)*CELL_SIZE)**2)
                IF (DCOMP<=RADMAX) THEN
                    RNUMBER=RNUMBER+1
                    MEANELEV=MEANELEV+CRATERELEVS(I,J)
                ENDIF
            ENDDO M200
        ENDDO L200
        MEANELEV=MEANELEV/RNUMBER
        WRITE(outcrater,408) ELEVCOMP,MEANELEV
408     FORMAT('ELEVCOMP=',G13.6,' MEANELEV=',G13.6)
    npoints=0
    mpoints=0
    meanelev=0.0
    interiorelev=0.0
    sumx=0.0
    sumx2=0.0
    sumy=0.0
    sumy2=0.0
    sumxy=0.0
    sumz=0.0
    sumzx=0.0
    sumzy=0.0
    xicent=float(icratcent)
    yjcent=float(jcratcent)
    do j=1,myc
        jy=float(j)
        do i=1,mxc
            ix=float(i)
            rloc=sqrt((ix-xicent)**2+(jy-yjcent)**2)
            if((rloc>=rinner).and.(rloc<=router)) then
                npoints=npoints+1
                meanelev=meanelev+craterelevs(i,j)
                sumx=sumx+ix
                sumx2=sumx2+ix*ix
                sumy=sumy+jy
                sumy2=sumy2+jy*jy
                sumxy=sumxy+ix*jy
                sumz=sumz+craterelevs(i,j)
                sumzx=sumzx+craterelevs(i,j)*ix
                sumzy=sumzy+craterelevs(i,j)*jy
            endif
            if (rloc<scale_distance) then
                interiorelev=interiorelev+craterelevs(i,j)
                mpoints=mpoints+1
            endif
        enddo
    enddo
    if (npoints>0) then
        meanelev=meanelev/npoints
        aterm=(sumzx*sumy2-sumzy*sumxy)/(sumx2*sumy2-sumxy**2)
        bterm=(sumzy*sumx2-sumzx*sumxy)/(sumx2*sumy2-sumxy**2)
        cterm=(sumz-aterm*sumx-bterm*sumy)/npoints
    endif
    if (mpoints>0) then
        interiorelev=interiorelev/mpoints
    endif
    write(outcrater,757) npoints,mpoints,meanelev,aterm,bterm,cterm
757 format('npoints-',i7,' mpoints=',i7,' meanelev=',g13.6,' aterm=',g13.6,' bterm=',g13.6,' cterm=',g13.6)
       do j=1,myc
        jy=float(j)
        do i=1,mxc
            ix=float(i)
                rloc=sqrt((ix-xicent)**2+(jy-yjcent)**2)
                if (rloc<=rinherit) then
                    newel(i,j)=craterelevs(i,j)-meanelev
                    newmean=newmean+newel(i,j)
                    newnumb=newnumb+1.0
                else
                    if (rloc>(router)) then
                        newel(i,j)=0.0
                    else
                        eref=craterelevs(i,j)-(aterm*ix+jy*bterm+cterm)
                        rslope=1.0/(rrim*(1.0-2.0))
                        rint=1.0-rslope*rrim
                        costerm=COS(3.14159*((rloc*cell_size/RADIUS-RADIUS_MAX_INHERIT) &
                        /((LOCAL_MAX_USE-RADIUS_MAX_INHERIT)))/2.0)**COSINE_POWER
                       IF ((.NOT.DONE_ONCE).AND.(DIAMETER>50000.0).AND.(J==(MYC/2+1))) THEN
                            !WRITE(44,836) I,RLOC,COSTERM
836                         FORMAT(I7,' ',G12.5,' ',G12.5)
                       ENDIF
                        t1=rloc*rslope+rint
                        t2=1.0-t1
                        term1=t1*(craterelevs(i,j)-meanelev)*costerm
                        term2=(1.0-costerm)
                        newel(i,j)=term1
                        newmean=newmean+newel(i,j)
                        newnumb=newnumb+1.0
                    endif
                endif
        enddo
    enddo
    IF ((.NOT.DONE_ONCE).AND.(DIAMETER>50000.0)) THEN
                DONE_ONCE=.TRUE.
                CLOSE(44)
    ENDIF
    write(outcrater,747) newmean/newnumb,interiorelev,rslope,rint
747 format('newmean=',g13.6,' interiorelev=',g13.6,' rslope=',g13.6,' rint=',g13.6,/)
    do_this=.false.
! is this necessary?
    if (do_this) then
    IIRANGE=2
       do j=1,myc
        jy=float(j)
        do i=1,mxc
            ix=float(i)
                rloc=sqrt((ix-xicent)**2+(jy-yjcent)**2)
                if ((rloc>=rinherit).and.(rloc<router)) then
                    iw=max(i-IIRANGE,1)
                    ie=min(i+IIRANGE,mxc)
                    jn=max(j-IIRANGE,1)
                    js=min(j+IIRANGE,myc)
                    esum=0.0
                    enum=0.0
                    do jj=jn,js
                        do ii=iw,ie
                            esum=esum+newel(i,j)
                            enum=enum+1.0
                        enddo
                    enddo
                    craterelevs(i,j)=esum/enum
                endif
        enddo
        enddo
      do j=1,myc
        jy=float(j)
        do i=1,mxc
            ix=float(i)
                rloc=sqrt((ix-xicent)**2+(jy-yjcent)**2)
                if ((rloc>=rinherit).and.(rloc<router)) then
                   newel(i,j)=craterelevs(i,j)
                endif
        enddo
      enddo
      endif
        L300: DO J=JMIN,JMAX
            JJJ=J-JCNTR+JCRATCENT
            M300: DO I=IMIN,IMAX
                    III=I-ICNTR+ICRATCENT
                    XCOMP=(I-1)*CELL_SIZE
                    YCOMP=(J-1)*CELL_SIZE
                    DCOMP=SQRT((XCOMP-X_LOCATION)**2+(YCOMP-Y_LOCATION)**2)
                    if (dcomp>(router*cell_size)) cycle M300
                    elevation(i,j)=dcomp*(elevation(i,j)-elevcomp-newel(iii,jjj))/(router*cell_size)+elevcomp+newel(iii,jjj)
            ENDDO M300
        ENDDO L300
        DEALLOCATE(CRATERELEVS,newel)
         L500: DO J=JMIN,JMAX
           M500:  DO I=IMIN,IMAX
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                       IS_SEDIMENT_COVERED(I,J)=.FALSE.
                      SEDIMENT_BASE(I,J)=ELEVATION(I,J)
                      CFNE(I,J)=0.0
                       IF (IS_REGOLITH_CRATER.OR.(DIAMETER<MINIMUM_HARD_DIAMETER) )  THEN
                         IS_ROCK_SURFACE(I,J)=.FALSE.
                         REGOLITH(I,J)=INITIAL_REGOLITH_THICKNESS

                       ELSE

                         IS_ROCK_SURFACE(I,J)=.TRUE.
                         REGOLITH(I,J)=-ROCK_WEATHERING_RATE
                      ENDIF
                ENDIF
                IF (MODEL_GROUNDWATER) REFERENCE_ELEVATION(I,J)=ELEVATION(I,J)
           ENDDO M500
         ENDDO L500
          IF (USE_CRATER_EVENT) CRATER_EVENT(ICNTR,JCNTR)=-RADIUS
        IF (IS_REGOLITH_CRATER.OR.(DIAMETER<MINIMUM_HARD_DIAMETER) )  THEN
                                     WRITE(*,403)
403                                  FORMAT('SOFT CRATER')
        ELSE
                                       WRITE(*,404)
404  FORMAT('HARD CRATER')
        ENDIF
        RETURN
399     CONTINUE
        WRITE(*,402)
402     FORMAT('ERROR IN READING FILE, RETURNING')
        RETURN
        END !SUBROUTINE CREATE_REAL_CRATER 
