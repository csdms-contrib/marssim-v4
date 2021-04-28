    !     ####################################################################################################
    !     this is one of several source files for the marssim landform evolution model
    !     copyright (c) 2009 alan d. howard
    !     developer can be contacted by ah6p`virginia.edu and department of environmental sciences, p.o. box 400123,
    !                university of virginia, charlottesville, va 22904-4123
    !     this program is free software; you can redistribute it and/or modify it under the terms of the gnu general public license
    !       as published by the free software foundation; either version 2 of the
    !       license, or (at your option) any later version.
    !     this program is distributed in the hope that it will be useful, but without any warranty;
    !        without even the implied warranty of merchantability or fitness for a particular purpose. see the gnu
    !        general public license for more details.
    !      you should have received a copy of the gnu general public license along with this program; if not, write to
    !        the free software foundation, inc., 51 franklin street, fifth floor, boston, ma 02110-1301 usa.
    !        a web link:   http://www.gnu.org/licenses/gpl-2.0.txt
    !     ####################################################################################################
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	SUBROUTINE READ_ACCRETION_ABLATION_PARAMETERS
	USE ERODE_GLOBALS
	USE ACCRETION_GLOBALS
	IMPLICIT NONE
	CHARACTER (80) :: TEXTLINE
	INTEGER :: ACCRETION_PARAMS, EXPOSECREEPUSE, RADIATION_USE, TOPEXPOSEUSE, INVERSEUSE
	INTEGER :: EXPOSESMOOTH,USE_REGOLITH_ABLATION, IERROR,REGOLITH_DEPOSIT,NORMAL_DEPOSIT
	    ACCRETION_PARAMS=93
		OPEN(ACCRETION_PARAMS,FILE='ACCRETION_ABLATION_PARAMETERS.PRM',ACTION='READ')
        !       **********************************************************************
        !
        !       Read in accretion parameters
        !
        !       **********************************************************************
        READ(ACCRETION_PARAMS,22) TEXTLINE
22      FORMAT(A)        
        WRITE(*,22) TEXTLINE
        READ(ACCRETION_PARAMS,*) ACCRETION_RATE
        READ(ACCRETION_PARAMS,*) EXPOSECREEPUSE
        READ(ACCRETION_PARAMS,*) RADIATION_USE
        READ(ACCRETION_PARAMS,*) TOPEXPOSEUSE
        READ(ACCRETION_PARAMS,*) INVERSEUSE
        READ(ACCRETION_PARAMS,*) EXPOSESMOOTH
        IF (EXPOSECREEPUSE > 0) THEN
            EXPOSURE_DEPENDENT_CREEP=.TRUE.
            ALLOCATE(EXPOSE_RATE(MX,MY),SMOOTH_EXPOSE(MX,MY),STAT=IERROR)
			IF (IERROR > 0) CALL ALLOCERROR()
        ELSE
            EXPOSURE_DEPENDENT_CREEP=.FALSE.
        ENDIF
        IF (RADIATION_USE > 0) THEN
            USE_SOLAR_EROSION=.TRUE.
        ELSE
            USE_SOLAR_EROSION=.FALSE.
        ENDIF
        IF (TOPEXPOSEUSE > 0) THEN
            USE_TOP_EXPOSURE=.TRUE.
        ELSE
            USE_TOP_EXPOSURE=.FALSE.
        ENDIF
        IF (INVERSEUSE > 0) THEN
            USE_INVERSE_EXPOSURE=.TRUE.
        ELSE
            USE_INVERSE_EXPOSURE=.FALSE.
        ENDIF
        IF (EXPOSESMOOTH > 0) THEN
            USE_EXPOSURE_SMOOTHING=.TRUE.
        ELSE
            USE_EXPOSURE_SMOOTHING=.FALSE.
        ENDIF
        READ(ACCRETION_PARAMS,*) RAD_CONST
        READ(ACCRETION_PARAMS,*) RAD_DUST_FACTOR
        READ(ACCRETION_PARAMS,*) RAD_THRESH_CONVEXITY
        READ(ACCRETION_PARAMS,*) RAD_DEPOSIT_RATE
		READ(ACCRETION_PARAMS,*) DEPOSITFRACTION 
        READ(ACCRETION_PARAMS,*) USE_REGOLITH_ABLATION
        DO_REGOLITH_ABLATION = .FALSE.
        IF (USE_REGOLITH_ABLATION>0) THEN
            DO_REGOLITH_ABLATION=.TRUE.
        ENDIF
        READ(ACCRETION_PARAMS,*) NORMAL_DEPOSIT
        IF (NORMAL_DEPOSIT>0) THEN
            DO_NORMAL_ACCRETION=.TRUE.
        ELSE
            DO_NORMAL_ACCRETION=.FALSE.
        ENDIF
        READ(ACCRETION_PARAMS,*) REGOLITH_DEPOSIT
        IF (REGOLITH_DEPOSIT>0) THEN
            DO_ACCRETE_REGOLITH = .TRUE.
        ELSE
            DO_ACCRETE_REGOLITH = .FALSE.
        ENDIF
		CLOSE(ACCRETION_PARAMS)
		       !    *************************************************************************
        !         Write accretion-ablation parameters
        ! ****************************************************************************
        IF (MODEL_ACCRETION_AND_ABLATION) THEN
            WRITE(OUTHIST,1553) ACCRETION_RATE
            1553 FORMAT('****************MODELING ACCRETION AND ABLATION***********',/, &
            ' ACCRETION RATE=',G12.5)
        ENDIF
        IF (EXPOSURE_DEPENDENT_CREEP) THEN
            WRITE(OUTHIST,4987)
            4987 FORMAT ('**********MODELING EXPOSURE-DEPENDENT-CREEP**************',/, &
            ' USES EXPOSURE-DEPENDENT CREEP')
        ENDIF
        IF (USE_SOLAR_EROSION) THEN
            WRITE(OUTHIST,4986) RAD_CONST, RAD_DUST_FACTOR,RAD_THRESH_CONVEXITY,RAD_DEPOSIT_RATE
            4986 FORMAT ('**********************USES SOLAR RADIATION EROSION************',/, &
            'RADIATION_EROSION_SCALING=',G12.5, &
            ' DUST ENHANCEMENT FACTIOR=',G12.5,/,' CRITICAL CONVEXITY FOR ICE DEPOSITION=',G12.5, &
            ' ICE_DEPOSITION_RATE=',G12.5)
        ENDIF
        IF (MODEL_ACCRETION_AND_ABLATION.OR.EXPOSURE_DEPENDENT_CREEP.OR.USE_SOLAR_EROSION) THEN
            IF (USE_TOP_EXPOSURE) THEN
                WRITE(OUTHIST,4984)
                4984 FORMAT ('USES TOP EXPOSURE CALCULATION')
            ENDIF
            IF (USE_INVERSE_EXPOSURE) THEN
                WRITE(OUTHIST,4983)
                4983 FORMAT(' USES INVERSE MEASURE OF EXPOSURE')
            ENDIF
            IF (USE_EXPOSURE_SMOOTHING) THEN
                WRITE(OUTHIST,4982)
                4982 FORMAT(' USES SMOOTHING OF EXPOSURE INDEX')
            ENDIF
        ENDIF
		IF (USE_SOLAR_EROSION) THEN
			WRITE(OUTHIST,4988) DEPOSITFRACTION
4988        FORMAT(' FRACTION OF SUBMIMATED VOLATILES REDEPOSITED ON SURFACE=',G12.5)
        ENDIF
        IF (DO_REGOLITH_ABLATION) THEN
            WRITE(OUTHIST,5988)
5988        FORMAT('REGOLITH IS SUBJECT TO ABLATION')
        ELSE
            WRITE(OUTHIST,5989)
5989        FORMAT('REGOLITH IS NOT SUBJECT TO ABLATION')
        ENDIF
        IF (NORMAL_DEPOSIT>0) THEN
            WRITE(OUTHIST,5990)
5990        FORMAT('SURFACE DEPOSITION/ABLATION IS NORNAL TO SLOPE')
        ELSE
            WRITE(OUTHIST,5991)
5991        FORMAT('SURFACE DEPOSITION/ABLATION IS PERPENDICULAR')
        ENDIF
        IF (REGOLITH_DEPOSIT>0) THEN
            WRITE(OUTHIST,5992)
5992        FORMAT('SURFACE ACCRETION OCCURS AS REGOLITH')
        ELSE
            WRITE(OUTHIST,5993)
5993        FORMAT('SURFACE ACCRETION DOES NOT CHANGE STATE OF SURFACE MATERIALS')
        ENDIF
	RETURN
	END !SUBROUTINE READ_ACCRETION_ABLATION_PARAMETERS
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
    SUBROUTINE DO_ACCRETION_ABLATION
        USE ERODE_GLOBALS
        USE ACCRETION_GLOBALS
        IMPLICIT NONE
        LOGICAL :: NORMAL_DEPOSITION
        INTEGER :: I,J,IW,IE,JN,JS
        REAL(4) :: MAXCORRECT,CORRECT
        REAL(4) :: Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9
        REAL(4) :: GG,HH,DNM,ECHANGE
        !     ***************************************************************************
        !      This routine models surface-normal erosion or deposition, such as solution, ablation, or
        !         mineral or volaltile crust deposition.  the local gradient of the surface
        !         is typically determined by fitting a plane to the 3x3 neighbor cells.  if usemaxgrad is
        !         set to true then the d8 gradient is used instead.  for surface-normal erosion or
        !         ablation the vertical change is proportional to the inverse cosine of gradient.
        !      The variable default_eolian_process is used to determine whether the surface normal is determined
        !        by the downslope (d8) gradient or the average gradient through the surrounding 9 points.
        !  MODIFIES: ELEVATION
        !     ***************************************************************************
        NORMAL_DEPOSITION=.FALSE.
        MAXCORRECT=0.0
        DO J=1,MY
            DO I=1,MX
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
                    IF (JS > MY) JS=1
                ELSE
                    IF (JS > MY) JS=MY-1
                ENDIF
                IF (IS_Y_PERIODIC) THEN
                    IF (JN < 1) JN=MY
                ELSE
                    IF (JN < 1) JN=2
                ENDIF
                Z1=ELEVATION(IW,JN)
                Z2=ELEVATION(I,JN)
                Z3=ELEVATION(IE,JN)
                Z4=ELEVATION(IW,J)
                Z5=ELEVATION(I,J)
                Z6=ELEVATION(IE,J)
                Z7=ELEVATION(IW,JS)
                Z8=ELEVATION(I,JS)
                Z9=ELEVATION(IE,JS)
                GG=(-Z1+Z3-Z4+Z6-Z7+Z9)/(6.0*CELL_SIZE)
                HH=(-Z1-Z2-Z3+Z7+Z8+Z9)/(6.0*CELL_SIZE)
                DNM=GG*GG+HH*HH
                IF (DO_NORMAL_ACCRETION) THEN
                    IF (DEFAULT_EOLIAN_PROCESS) THEN
                        CORRECT=SQRT(1.0+D8_GRADIENT(I,J)**2)
                    ELSE
                        CORRECT=SQRT(1.0+DNM)
                    ENDIF
                ELSE
                    CORRECT=1.0
                ENDIF
                IF (CORRECT > MAXCORRECT) MAXCORRECT=CORRECT
                ECHANGE= ACCRETION_RATE*TIME_INCREMENT*VOLUME_CHANGE_COEFFICIENT  &
                * CORRECT
                ELEVATION(I,J)=ELEVATION(I,J)+ ECHANGE
                ACCRETE_WORK=ACCRETE_WORK+ECHANGE
                IF (DO_ACCRETE_REGOLITH) THEN
                    IF (ECHANGE>0.0) THEN
                        REGOLITH(I,J)=REGOLITH(I,J)+ECHANGE
                        IS_ROCK_SURFACE(I,J)=.FALSE.
                    ELSE
                        REGOLITH(I,J)=MAX(0.0,(REGOLITH(I,J)+ECHANGE))
                        IF (REGOLITH(I,J)<=0.0) THEN
                            IS_ROCK_SURFACE(I,J)=.TRUE.
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        WRITE(*,100) MAXCORRECT, TIME_INCREMENT
        100   FORMAT(' MAXIMUM CORRECTION=',G12.5, ' TIME INCREMENT=',G12.5)
        RETURN
    END ! SUBROUTINE DO_ACCRETION_ABLATION
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DO_EXPOSURE_DEPENDENT_CREEP()
        USE ERODE_GLOBALS
        USE EOLIAN_GLOBALS
        USE ACCRETION_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,IW,IE,JN,JS,K,L,II,JJ
        REAL(4) :: TEMP1,GRADIENT,EEEE
        REAL(4) :: MAXEXPOSE,MINEXPOSE
        !     **********************************************************************
        !       This subroutine erodes or deposits on a presumably regolith-covered surface
        !       in proportion to the 'exposure' through gradient-dependent mass wasting.
        !       There is the option to use either of two different measures of 'exposure' (esposure and topexpose)
        !       and to have positive exposure enhance or diminish creep rate.
        !       The exposure index can also be smoothed prior to creep.
        !  MODIFIES: ERODE_SLOPE, CFN, CFW, CFNE, CFNW, ELEVATION , EXPOSE_RATE, SMOOTH_EXPOSE
        !  CALLS: FIND_TOP_EXPOSURE_INDEX, EXPOSURE
        !     **********************************************************************

        MINEXPOSE=1.0E+24
        MAXEXPOSE=-MINEXPOSE
        DO J=1,MY
            DO I=1,MX
                ERODE_SLOPE(I,J)=0.0
                IF (USE_TOP_EXPOSURE) THEN
                    CALL FIND_TOP_EXPOSURE_INDEX(I,J,EXPOSE_RATE(I,J))
                ELSE
                    CALL EXPOSURE(I,J,EXPOSE_RATE(I,J))
                ENDIF
                TEMP1=EOLIAN_CONSTANT_3*(EXPOSE_RATE(I,J)-EXPOSURE_50_PERCENT)
                EXPOSE_RATE(I,J)=EOLIAN_CONSTANT_2+EOLIAN_CONSTANT_1*TEMP1/SQRT(1.0+TEMP1**2)
                IF (USE_INVERSE_EXPOSURE) EXPOSE_RATE(I,J)=1.0-EXPOSE_RATE(I,J)
                IF (EXPOSE_RATE(I,J) > MAXEXPOSE) MAXEXPOSE=EXPOSE_RATE(I,J)
                IF (EXPOSE_RATE(I,J) < MINEXPOSE) MINEXPOSE=EXPOSE_RATE(I,J)
            ENDDO
        ENDDO
        IF (USE_EXPOSURE_SMOOTHING) THEN
            MINEXPOSE=1.0E+24
            MAXEXPOSE=-MINEXPOSE
            DO L=1,MY
                DO K=1,MX
                    SMOOTH_EXPOSE(K,L)=0.0
                    DO J=-1,1
                        DO I=-1,1
                            II=I
                            JJ=J
                            IF (IS_X_PERIODIC) THEN
                                IF (I < 1) II=II+MX
                                IF (I > MX) II=II-MX
                            ELSE
                                IF (I < 1) II=II-I+1
                                IF (I > MX) II=2*MX-I
                            ENDIF
                            IF (IS_Y_PERIODIC) THEN
                                IF (J > MY) JJ=JJ+MY
                            ELSE
                                IF (J > MY) JJ=JJ-MY
                            ENDIF
                            IF (IS_Y_PERIODIC) THEN
                                IF (J < 1) JJ=JJ-J+1
                            ELSE
                                IF (J < 1) JJ=2*MY-J
                            ENDIF
                            SMOOTH_EXPOSE(K,L)=SMOOTH_EXPOSE(K,L) &
                            +EXPOSE_RATE(K+II,L+JJ)
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
            DO J=1,MY
                DO I=1,MX
                    EXPOSE_RATE(I,J)=SMOOTH_EXPOSE(I,J)/9.0
                    IF (EXPOSE_RATE(I,J) > MAXEXPOSE) MAXEXPOSE=EXPOSE_RATE(I,J)
                    IF (EXPOSE_RATE(I,J) < MINEXPOSE) MINEXPOSE=EXPOSE_RATE(I,J)
                ENDDO
            ENDDO
        ENDIF
        WRITE (*,222) MINEXPOSE,MAXEXPOSE
        222   FORMAT('EXPOSECREEP, MN= ',G12.5,' MX= ',G12.5)
        DO J=1,MY
            DO I=1,MX
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
                    IF (JS > MY) JS=1
                ELSE
                    IF (JS > MY) JS=MY-1
                ENDIF
                IF (IS_Y_PERIODIC) THEN
                    IF (JN < 1) JN=MY
                ELSE
                    IF (JN < 1) JN=2
                ENDIF
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IW,J))/CELL_SIZE
                EEEE=0.5*(EXPOSE_RATE(I,J)+EXPOSE_RATE(IW,J))
                IF (ELEVATION(I,J) < ELEVATION(IW,J)) THEN
                    CFW(I,J)=-SLOPE_DIFFUSIVITY*GRADIENT*EEEE
                ELSE
                    CFW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT*EEEE
                ENDIF
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IE,JN))/CELL_SIZE
                EEEE=0.5*(EXPOSE_RATE(I,J)+EXPOSE_RATE(IE,JN))
                IF (ELEVATION(I,J) < ELEVATION(IE,JN)) THEN
                    CFNE(I,J)=-SLOPE_DIFFUSIVITY* &
                    GRADIENT*EEEE
                ELSE
                    CFNE(I,J)=SLOPE_DIFFUSIVITY*GRADIENT*EEEE
                ENDIF
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I,JN))/CELL_SIZE
                EEEE=0.5*(EXPOSE_RATE(I,J)+EXPOSE_RATE(I,JN))
                IF (ELEVATION(I,J) < ELEVATION(I,JN)) THEN
                    CFN(I,J)=-SLOPE_DIFFUSIVITY*      &
                    GRADIENT*EEEE
                ELSE
                    CFN(I,J)=SLOPE_DIFFUSIVITY*GRADIENT*EEEE
                ENDIF
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IW,JN))/CELL_SIZE
                EEEE=0.5*(EXPOSE_RATE(I,J)+EXPOSE_RATE(IW,JN))
                IF (ELEVATION(I,J) < ELEVATION(IW,JN)) THEN
                    CFNW(I,J)=-SLOPE_DIFFUSIVITY*GRADIENT*EEEE
                ELSE
                    CFNW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT*EEEE
                ENDIF
            ENDDO
        ENDDO
        DO J=1,MY
            DO I=1,MX
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
                    IF (JS > MY) JS=1
                ELSE
                    IF (JS > MY) JS=MY-1
                ENDIF
                IF (IS_Y_PERIODIC) THEN
                    IF (JN < 1) JN=MY
                ELSE
                    IF (JN < 1) JN=2
                ENDIF
                IF ((IE <= MX).OR.IS_X_PERIODIC) THEN
                    IF (CFW(IE,J) <= 0.0) THEN
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFW(IE,J)*CROSS_WEIGHTING
                    ELSE
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFW(IE,J)*CROSS_WEIGHTING
                    ENDIF
                ENDIF
                IF (CFN(I,JS) <= 0.0) THEN
                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFN(I,JS)*CROSS_WEIGHTING
                ELSE
                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFN(I,JS)*CROSS_WEIGHTING
                ENDIF
                IF (CFW(I,J) >= 0.0) THEN
                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFW(I,J)*CROSS_WEIGHTING
                ELSE
                    IF ((I > 1).OR.IS_X_PERIODIC) THEN
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFW(I,J)*CROSS_WEIGHTING
                    ENDIF
                ENDIF
                IF (CFN(I,J) >= 0.0) THEN
                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFN(I,J)*CROSS_WEIGHTING
                ELSE
                    IF ((J > 1).OR.(IS_Y_PERIODIC)) THEN
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFN(I,J)*CROSS_WEIGHTING
                    ENDIF
                ENDIF
                IF (CFNE(I,J) > 0.0) THEN
                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNE(I,J)*DIAGONAL_WEIGHTING
                ELSE
                    IF ((J > 1).OR.(IS_Y_PERIODIC)) THEN
                        IF ((I < MX).OR.IS_X_PERIODIC) THEN
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNE(I,J)*DIAGONAL_WEIGHTING
                        ENDIF
                    ENDIF
                ENDIF
                IF (CFNW(I,J) > 0.0) THEN
                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNW(I,J)*DIAGONAL_WEIGHTING
                ELSE
                    IF ((J > 1).OR.(IS_Y_PERIODIC)) THEN
                        IF ((I > 1).OR.IS_X_PERIODIC) THEN
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNW(I,J)*DIAGONAL_WEIGHTING
                        ENDIF
                    ENDIF
                ENDIF
                IF ((I > 1).OR.IS_X_PERIODIC) THEN
                    IF (CFNE(IW,JS) <= 0.0) THEN
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNE(IW,JS)*DIAGONAL_WEIGHTING
                    ELSE
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNE(IW,JS)*DIAGONAL_WEIGHTING
                    ENDIF
                ENDIF
                IF ((I < MX).OR.IS_X_PERIODIC) THEN
                    IF (CFNW(IE,JS) <= 0.0) THEN
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNW(IE,JS)*DIAGONAL_WEIGHTING
                    ELSE
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNW(IE,JS)*DIAGONAL_WEIGHTING
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        DO J=1,MY
            DO I=1,MX
                ELEVATION(I,J)=ELEVATION(I,J)+ERODE_SLOPE(I,J)
            ENDDO
        ENDDO
        !IF (MOD(ITERATION,IMAGE_OUTPUT_INTERVAL) == 0) THEN
            !OPEN(73,FILE='ACCRETION.PRN')
            !J=MY/2
            !DO I=1,MX
                !WRITE(73,7333) I,ELEVATION(I,J),EXPOSE_RATE(I,J),ERODE_SLOPE(I,J) &
                !,CFN(I,J),CFW(I,J),CFNE(I,J),CFNW(I,J)
            !ENDDO
            !CLOSE(73)
        !ENDIF
        7333  FORMAT(I5,7(' ',G12.5))
    END ! SUBROUTINE DO_EXPOSURE_DEPENDENT_CREEP
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE FIND_TOP_EXPOSURE_INDEX(I,J,EXPOSE)
        USE ERODE_GLOBALS
        USE EOLIAN_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,J
        INTEGER :: NN,II,JJ,INDX,JNDX
        REAL(4) :: ECOMP,ELOC,WLOC,LSLOPE,WEIGHTSUM
        REAL(4) :: SLOPESUM,STEP_DISTANCE
        REAL(4), INTENT(OUT) :: EXPOSE
        !     ***********************************************************************
        !      Tthis calculates an 'exposure' factor of a given location based upon the the relative distance
        !         and elevation of surrounding cells. the 'exposure' is more negative to the degree
        !         that surrounding points are lower than the given location. this definition of 'exposure'
        !         is an alternative to that calculated by the exposure subroutine in the eolian module.
        !     ***********************************************************************
        NN=MIN(MAXIMUM_WEIGHT_DISTANCE,25)
        ECOMP=ELEVATION(I,J)
        SLOPESUM=0.0
        WEIGHTSUM=0.0
        L100: DO  JNDX=-NN,NN
            JJ=J+JNDX
            IF (JJ < 1) THEN
                IF (IS_Y_PERIODIC) THEN
                    JJ=JJ+MY
                ELSE
                    CYCLE L100
                ENDIF
            ENDIF
            IF (JJ > MY) THEN
                IF (IS_Y_PERIODIC) THEN
                    JJ=JJ-MY
                ELSE
                    CYCLE L100
                ENDIF
            ENDIF
            L200: DO  INDX=-NN,NN
                II=I+INDX
                IF (II < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        II=II+MX
                    ELSE
                        CYCLE L200
                    ENDIF
                ENDIF
                IF (II > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        II=II-MX
                    ELSE
                        CYCLE L200
                    ENDIF
                ENDIF
                IF ((II == I).AND.(JJ == J)) THEN
                ELSE
                    STEP_DISTANCE=DISTANCE_VALUE(ICNNT-INDX,JCNNT-JNDX)
                    ELOC=ELEVATION(II,JJ)
                    LSLOPE=(ELOC-ECOMP)/STEP_DISTANCE
                    WLOC=DISTANCE_WEIGHTING(ICNNT-INDX,JCNNT-JNDX)
                    SLOPESUM=SLOPESUM+LSLOPE*WLOC
                    WEIGHTSUM=WEIGHTSUM+WLOC
                ENDIF
            ENDDO L200
        ENDDO L100
        EXPOSE=SLOPESUM/(WEIGHTSUM*CELL_SIZE)
        RETURN
    END ! SUBROUTINE FIND_TOP_EXPOSURE_INDEX
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE SETUP_DISTANCE_WEIGHTING()
        USE ERODE_GLOBALS
        USE EOLIAN_GLOBALS
        !     ***********************************************************************
        !      This sets up the weighting for determining the 'exposure' of a given location
        !        the weighting decreases as a negative exponential of distance
        !  MODIFIES: DISTANCE_VALUE, DISTANCE_WEIGHTING
        !     ***********************************************************************

        IMPLICIT NONE
        INTEGER I,J
        ICNNT=26
        JCNNT=26
        DO J=1,51
            DO I=1,51
                DISTANCE_VALUE(I,J)=SQRT((DFLOAT(I-ICNNT))**2+ &
                (DFLOAT(J-JCNNT))**2)
                DISTANCE_WEIGHTING(I,J)=EXP(WEIGHTING_DECAY_FACTOR*DISTANCE_VALUE(I,J))
            ENDDO
        ENDDO
        RETURN
    END !SUBROUTINE SETUP_DISTANCE_WEIGHTING 
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (4) FUNCTION RAD_ERODE(I,J)
        USE ERODE_GLOBALS
        USE EOLIAN_GLOBALS
        USE ACCRETION_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,J
        INTEGER :: IH,IL,K,ICOMP,JCOMP,JH,JL,ILL,JLL,IHH,JHH
        REAL(4) :: DDIV,COMPSLOPE,COMPANGLE, DIRECTIONBIAS !!! EXPERIMENTAL
        REAL(4) :: ANGLE_FACTOR,ANGLE_DIFF,SUMWEIGHT,SUMCHANGE
        REAL(4) :: DISTANCE,RAD_CHANGE,RAD_WEIGHT,LOCALSLOPE,LOCALANGLE,LOCALMAXTAN,LOCALTAN
        INTEGER :: LOCAL_WEIGHT_DISTANCE
        LOGICAL :: USE_COSINE_SQUARE
        LOCAL_WEIGHT_DISTANCE=MIN(MAX_OUT_DISTANCE,MAXIMUM_WEIGHT_DISTANCE)
        !     ***********************************************************************
        !      This calculates the ablation of a surface based upon reradiated infrared light
        !        from surrounding cells.  The amount of reradiated energy depends upon whether it is
        !        a dark dust-covered location or a lighter 'bedrock' location.  Radiation contribution from surrounding cells
        !        are calculated up to a maximum distance (maximum_weight_distance) with a weighting proportional to the
        !        inverse square of the distance.  See howard and moore, 2008, geophysical research letters, 35, l03203,
        !        doi:10.1029/2007gl032618; moore, j.m. et al., 2017, icarus, 287, 320-333
        !      An ice-covered surface is indicated by IS_SEDIMENT_COVERED being TRUE and IS_ROCK_SURFACE is FALSE
        !      A bedrock surface is indicate be IS_ROCK_SURFACE being TRUE and IS_SEDIMENT_COVERED is FALSE
        !      A dust (sediment-covered) surface is indicated by IS_ROCK_SURFACE being FALSE and IS_SEDIMENT_COVERED being FALSE
        !  CALLS: FIND_RAD_CHANGE
        !     ***********************************************************************
        USE_COSINE_SQUARE=.false.
        DIRECTIONBIAS=1.0
        IF (.NOT.IS_ROCK_SURFACE(I,J)) THEN
            IF (.NOT.DO_REGOLITH_ABLATION) THEN
                RAD_ERODE=0.0
                RETURN
            ENDIF
        ENDIF
        SUMCHANGE=0.0
        SUMWEIGHT=0.0
        IL=I-1
        IH=I+1
        IF (IS_X_PERIODIC) THEN
            IF (IL < 1) IL=IL+MX
            IF (IH > MX) IH=IH-MX
            DDIV=2.0*CELL_SIZE
        ELSE
            IF (IL < 1) THEN
                IL=1
                DDIV=1.0*CELL_SIZE
            ENDIF
            IF (IH > MX) THEN
                IH=MX
                DDIV=1.0*CELL_SIZE
            ENDIF
        ENDIF
        LOCALSLOPE=(ELEVATION(IH,J)-ELEVATION(IL,J))/DDIV
        LOCALANGLE=ATAN(LOCALSLOPE)
        LOCALMAXTAN=-1000.0
        L200: DO K=1,LOCAL_WEIGHT_DISTANCE
            IF (K == 0) CYCLE L200
            ICOMP=I+K
            JCOMP=J
            IF (IS_X_PERIODIC) THEN
                IF (ICOMP < 1) ICOMP=ICOMP+MX
                IF (ICOMP > MX) ICOMP=ICOMP-MX
            ELSE
                IF (ICOMP < 1) CYCLE L200
                IF (ICOMP > MX) CYCLE L200
            ENDIF
            ILL=ICOMP-1
            IHH=ICOMP+1
            IF (IS_X_PERIODIC) THEN
                IF (ILL < 1) ILL=ILL+MX
                IF (IHH > MX) IHH=IHH-MX
                DDIV=2.0*CELL_SIZE
            ELSE
                IF (ILL < 1) THEN
                    ILL=1
                    DDIV=1.0*CELL_SIZE
                ENDIF
                IF (IHH > MX) THEN
                    IHH=MX
                    DDIV=1.0*CELL_SIZE
                ENDIF
            ENDIF
            LOCALTAN=(ELEVATION(ICOMP,J)-ELEVATION(I,J))/(ABS(K)*CELL_SIZE)
            if (LOCALTAN.GE.LOCALMAXTAN) THEN
               LOCALMAXTAN=LOCALTAN
            ELSE
               CYCLE L200
            ENDIF
            COMPSLOPE=(ELEVATION(IHH,J)-ELEVATION(ILL,J))/DDIV
            COMPANGLE=ATAN(COMPSLOPE)
            ANGLE_DIFF=COMPANGLE-LOCALANGLE
            IF (ANGLE_DIFF > 0.0) THEN
              
                ANGLE_FACTOR=(SIN(0.5*ANGLE_DIFF))
                IF (USE_COSINE_SQUARE) ANGLE_FACTOR=ANGLE_FACTOR**2  !!NEW
                
            ELSE
                ANGLE_FACTOR=0.0
            ENDIF
            DISTANCE=SQRT((ABS(K)*CELL_SIZE)**2+(ELEVATION(ICOMP,JCOMP)-ELEVATION(I,J))**2)   !!NEW
            CALL FIND_RAD_CHANGE(ICOMP,JCOMP,ANGLE_FACTOR,DISTANCE,RAD_CHANGE,RAD_WEIGHT)
            SUMCHANGE=SUMCHANGE+RAD_CHANGE/DIRECTIONBIAS !!!
            SUMWEIGHT=SUMWEIGHT+RAD_WEIGHT
        ENDDO L200
        LOCALMAXTAN=-1000.0
        L201: DO K=-1,-LOCAL_WEIGHT_DISTANCE,-1
            IF (K == 0) CYCLE L201
            ICOMP=I+K
            JCOMP=J
            IF (IS_X_PERIODIC) THEN
                IF (ICOMP < 1) ICOMP=ICOMP+MX
                IF (ICOMP > MX) ICOMP=ICOMP-MX
            ELSE
                IF (ICOMP < 1) CYCLE L201
                IF (ICOMP > MX) CYCLE L201
            ENDIF
            IL=ICOMP-1
            IH=ICOMP+1
            IF (IS_X_PERIODIC) THEN
                IF (IL < 1) IL=IL+MX
                IF (IH > MX) IH=IH-MX
                DDIV=2.0*CELL_SIZE
            ELSE
                IF (IL < 1) THEN
                    IL=1
                    DDIV=1.0*CELL_SIZE
                ENDIF
                IF (IH > MX) THEN
                    IH=MX
                    DDIV=1.0*CELL_SIZE
                ENDIF
            ENDIF
            LOCALTAN=(ELEVATION(ICOMP,J)-ELEVATION(I,J))/(ABS(K)*CELL_SIZE)
            if (LOCALTAN.GE.LOCALMAXTAN) THEN
               LOCALMAXTAN=LOCALTAN
            ELSE
               CYCLE L201
            ENDIF
            COMPSLOPE=(ELEVATION(IH,J)-ELEVATION(IL,J))/DDIV
            COMPANGLE=ATAN(COMPSLOPE)
            ANGLE_DIFF=LOCALANGLE-COMPANGLE
            IF (ANGLE_DIFF > 0.0) THEN
                ANGLE_FACTOR=(SIN(0.5*ANGLE_DIFF))
                IF (USE_COSINE_SQUARE) ANGLE_FACTOR=ANGLE_FACTOR**2  !!NEW
            ELSE
                ANGLE_FACTOR=0.0
            ENDIF
            DISTANCE=SQRT((ABS(K)*CELL_SIZE)**2+(ELEVATION(ICOMP,JCOMP)-ELEVATION(I,J))**2)   !!NEW
            CALL FIND_RAD_CHANGE(ICOMP,JCOMP,ANGLE_FACTOR,DISTANCE,RAD_CHANGE,RAD_WEIGHT)
            SUMCHANGE=SUMCHANGE+RAD_CHANGE/DIRECTIONBIAS !!!
            SUMWEIGHT=SUMWEIGHT+RAD_WEIGHT
        ENDDO L201
        JL=J-1
        JH=J+1
        IF (IS_Y_PERIODIC) THEN
            IF (JL < 1) JL=JL+MY
            IF (JH > MY) JH=JH-MY
            DDIV=2.0*CELL_SIZE
        ELSE
            IF (JL < 1) THEN
                JL=1
                DDIV=1.0*CELL_SIZE
            ENDIF
            IF (JH > MY) THEN
                JH=MY
                DDIV=1.0*CELL_SIZE
            ENDIF
        ENDIF
        LOCALMAXTAN=-1000.0
        LOCALSLOPE=(ELEVATION(I,JH)-ELEVATION(I,JL))/DDIV
        LOCALANGLE=ATAN(LOCALSLOPE)
        L210: DO K=1,LOCAL_WEIGHT_DISTANCE
            IF (K == 0) CYCLE L210
            ICOMP=I
            JCOMP=J+K
            IF (IS_X_PERIODIC) THEN
                IF (JCOMP < 1) JCOMP=JCOMP+MY
                IF (JCOMP > MY) JCOMP=JCOMP-MY
            ELSE
                IF (JCOMP < 1) CYCLE L210
                IF (JCOMP > MY) CYCLE L210
            ENDIF
            JL=JCOMP-1
            JH=JCOMP+1
            IF (IS_X_PERIODIC) THEN
                IF (JL < 1) JL=JL+MY
                IF (JH > MY) JH=JH-MY
                DDIV=2.0*CELL_SIZE
            ELSE
                IF (JL < 1) THEN
                    JL=1
                    DDIV=1.0*CELL_SIZE
                ENDIF
                IF (JH > MY) THEN
                    JH=MY
                    DDIV=1.0*CELL_SIZE
                ENDIF
            ENDIF
            LOCALTAN=(ELEVATION(I,JCOMP)-ELEVATION(I,J))/(ABS(K)*CELL_SIZE)
            if (LOCALTAN.GE.LOCALMAXTAN) THEN
               LOCALMAXTAN=LOCALTAN
            ELSE
               CYCLE L210
            ENDIF         
            COMPSLOPE=(ELEVATION(I,JH)-ELEVATION(I,JL))/DDIV
            COMPANGLE=ATAN(COMPSLOPE)
            ANGLE_DIFF=COMPANGLE-LOCALANGLE
            IF (ANGLE_DIFF > 0.0) THEN
                ANGLE_FACTOR=(SIN(0.5*ANGLE_DIFF))
                IF (USE_COSINE_SQUARE) ANGLE_FACTOR=ANGLE_FACTOR**2  !!NEW
            ELSE
                ANGLE_FACTOR=0.0
            ENDIF
            DISTANCE=SQRT((ABS(K)*CELL_SIZE)**2+(ELEVATION(ICOMP,JCOMP)-ELEVATION(I,J))**2)   !!NEW
            CALL FIND_RAD_CHANGE(ICOMP,JCOMP,ANGLE_FACTOR,DISTANCE,RAD_CHANGE,RAD_WEIGHT)
            SUMCHANGE=SUMCHANGE+RAD_CHANGE/DIRECTIONBIAS !!!
            SUMWEIGHT=SUMWEIGHT+RAD_WEIGHT
        ENDDO L210
        LOCALMAXTAN=-1000.0
        L211: DO K=-1,-LOCAL_WEIGHT_DISTANCE,-1
            IF (K == 0) CYCLE L211
            ICOMP=I
            JCOMP=J+K
            IF (IS_X_PERIODIC) THEN
                IF (JCOMP < 1) JCOMP=JCOMP+MY
                IF (JCOMP > MY) JCOMP=JCOMP-MY
            ELSE
                IF (JCOMP < 1) CYCLE L211
                IF (JCOMP > MY) CYCLE L211
            ENDIF
            JL=JCOMP-1
            JH=JCOMP+1
            IF (IS_X_PERIODIC) THEN
                IF (JL < 1) JL=JL+MY
                IF (JH > MY) JH=JH-MY
                DDIV=2.0*CELL_SIZE
            ELSE
                IF (JL < 1) THEN
                    JL=1
                    DDIV=1.0*CELL_SIZE
                ENDIF
                IF (JH > MY) THEN
                    JH=MY
                    DDIV=1.0*CELL_SIZE
                ENDIF
            ENDIF
            LOCALTAN=(ELEVATION(I,JCOMP)-ELEVATION(I,J))/(ABS(K)*CELL_SIZE)
            if (LOCALTAN.GE.LOCALMAXTAN) THEN
               LOCALMAXTAN=LOCALTAN
            ELSE
               CYCLE L211
            ENDIF         
            COMPSLOPE=(ELEVATION(I,JH)-ELEVATION(I,JL))/DDIV
            COMPANGLE=ATAN(COMPSLOPE)
            ANGLE_DIFF=LOCALANGLE-COMPANGLE
            IF (ANGLE_DIFF > 0.0) THEN
                ANGLE_FACTOR=(SIN(0.5*ANGLE_DIFF))
                IF (USE_COSINE_SQUARE) ANGLE_FACTOR=ANGLE_FACTOR**2  !!NEW
            ELSE
                ANGLE_FACTOR=0.0
            ENDIF
            DISTANCE=SQRT((ABS(K)*CELL_SIZE)**2+(ELEVATION(ICOMP,JCOMP)-ELEVATION(I,J))**2)   !!NEW
            CALL FIND_RAD_CHANGE(ICOMP,JCOMP,ANGLE_FACTOR,DISTANCE,RAD_CHANGE,RAD_WEIGHT)
            SUMCHANGE=SUMCHANGE+RAD_CHANGE/DIRECTIONBIAS !!!
            SUMWEIGHT=SUMWEIGHT+RAD_WEIGHT
        ENDDO L211      
        IL=I-1
        IH=I+1
        DDIV=2.0*CELL_SIZE
        IF (IS_X_PERIODIC) THEN
            IF (IL < 1) IL=IL+MX
            IF (IH > MX) IH=IH-MX
        ELSE
            IF (IL < 1) THEN
                IL=1
                DDIV=1.0*CELL_SIZE
            ENDIF
            IF (IH > MX) THEN
                IH=MX
                DDIV=1.0*CELL_SIZE
            ENDIF
        ENDIF
        JL=J-1
        JH=J+1
        IF (IS_Y_PERIODIC) THEN
            IF (JL < 1) JL=JL+MY
            IF (JH > MY) JH=JH-MY
        ELSE
            IF (JL < 1) THEN
                JL=1
                DDIV=1.0*CELL_SIZE
            ENDIF
            IF (JH > MY) THEN
                JH=MY
                DDIV=1.0*CELL_SIZE
            ENDIF
        ENDIF
        LOCALMAXTAN=-1000.0
        LOCALSLOPE=(ELEVATION(IH,JL)-ELEVATION(IL,JH))/(1.414*DDIV)
        LOCALANGLE=ATAN(LOCALSLOPE)
        L230: DO K=1,LOCAL_WEIGHT_DISTANCE
            IF (K == 0) CYCLE L230
            DDIV=2.0*CELL_SIZE
            ICOMP=I+K
            JCOMP=J-K
            IF (IS_X_PERIODIC) THEN
                IF (JCOMP < 1) JCOMP=JCOMP+MY
                IF (JCOMP > MY) JCOMP=JCOMP-MY
            ELSE
                IF (JCOMP < 1) CYCLE L230
                IF (JCOMP > MY) CYCLE L230
            ENDIF
            IF (IS_X_PERIODIC) THEN
                IF (ICOMP < 1) ICOMP=ICOMP+MX
                IF (ICOMP > MX) ICOMP=ICOMP-MX
            ELSE
                IF (ICOMP < 1) CYCLE L230
                IF (ICOMP > MX) CYCLE L230
            ENDIF
            JL=JCOMP-1
            JH=JCOMP+1
            IF (IS_X_PERIODIC) THEN
                IF (JL < 1) JL=JL+MY
                IF (JH > MY) JH=JH-MY
                DDIV=2.0*CELL_SIZE
            ELSE
                IF (JL < 1) THEN
                    JL=1
                    DDIV=1.0*CELL_SIZE
                ENDIF
                IF (JH > MY) THEN
                    JH=MY
                    DDIV=1.0*CELL_SIZE
                ENDIF
            ENDIF
            IL=ICOMP-1
            IH=ICOMP+1
            IF (IS_X_PERIODIC) THEN
                IF (IL < 1) IL=IL+MX
                IF (IH > MX) IH=IH-MX
                DDIV=2.0*CELL_SIZE
            ELSE
                IF (IL < 1) THEN
                    IL=1
                    DDIV=1.0*CELL_SIZE
                ENDIF
                IF (IH > MX) THEN
                    IH=MX
                    DDIV=1.0*CELL_SIZE
                ENDIF
            ENDIF
            LOCALTAN=(ELEVATION(ICOMP,JCOMP)-ELEVATION(I,J))/(ABS(K)*CELL_SIZE*1.414)
            if (LOCALTAN.GE.LOCALMAXTAN) THEN
               LOCALMAXTAN=LOCALTAN
            ELSE
               CYCLE L230
            ENDIF
            
            COMPSLOPE=(ELEVATION(IH,JL)-ELEVATION(IL,JH))/(1.414*DDIV)
            COMPANGLE=ATAN(COMPSLOPE)
            ANGLE_DIFF=COMPANGLE-LOCALANGLE
            IF (ANGLE_DIFF > 0.0) THEN
                ANGLE_FACTOR=(SIN(0.5*ANGLE_DIFF))
                IF (USE_COSINE_SQUARE) ANGLE_FACTOR=ANGLE_FACTOR**2  !!NEW
            ELSE
                ANGLE_FACTOR=0.0
            ENDIF
            DISTANCE=SQRT((ABS(K)*CELL_SIZE*1.414)**2+(ELEVATION(ICOMP,JCOMP)-ELEVATION(I,J))**2)   !!NEW            
            CALL FIND_RAD_CHANGE(ICOMP,JCOMP,ANGLE_FACTOR,DISTANCE,RAD_CHANGE,RAD_WEIGHT)
            SUMCHANGE=SUMCHANGE+RAD_CHANGE*DIRECTIONBIAS !!!
            SUMWEIGHT=SUMWEIGHT+RAD_WEIGHT
        ENDDO L230
        LOCALMAXTAN=-1000.0
        L231: DO K=-1,-LOCAL_WEIGHT_DISTANCE,-1
            IF (K == 0) CYCLE L231
            DDIV=2.0*CELL_SIZE
            ICOMP=I+K
            JCOMP=J-K
            IF (IS_X_PERIODIC) THEN
                IF (JCOMP < 1) JCOMP=JCOMP+MY
                IF (JCOMP > MY) JCOMP=JCOMP-MY
            ELSE
                IF (JCOMP < 1) CYCLE L231
                IF (JCOMP > MY) CYCLE L231
            ENDIF
            IF (IS_X_PERIODIC) THEN
                IF (ICOMP < 1) ICOMP=ICOMP+MX
                IF (ICOMP > MX) ICOMP=ICOMP-MX
            ELSE
                IF (ICOMP < 1) CYCLE L231
                IF (ICOMP > MX) CYCLE L231
            ENDIF
            JL=JCOMP-1
            JH=JCOMP+1
            IF (IS_X_PERIODIC) THEN
                IF (JL < 1) JL=JL+MY
                IF (JH > MY) JH=JH-MY
                DDIV=2.0*CELL_SIZE
            ELSE
                IF (JL < 1) THEN
                    JL=1
                    DDIV=1.0*CELL_SIZE
                ENDIF
                IF (JH > MY) THEN
                    JH=MY
                    DDIV=1.0*CELL_SIZE
                ENDIF
            ENDIF
            IL=ICOMP-1
            IH=ICOMP+1
            IF (IS_X_PERIODIC) THEN
                IF (IL < 1) IL=IL+MX
                IF (IH > MX) IH=IH-MX
                DDIV=2.0*CELL_SIZE
            ELSE
                IF (IL < 1) THEN
                    IL=1
                    DDIV=1.0*CELL_SIZE
                ENDIF
                IF (IH > MX) THEN
                    IH=MX
                    DDIV=1.0*CELL_SIZE
                ENDIF
            ENDIF
            LOCALTAN=(ELEVATION(ICOMP,JCOMP)-ELEVATION(I,J))/(ABS(K)*CELL_SIZE*1.414)
            if (LOCALTAN.GE.LOCALMAXTAN) THEN
               LOCALMAXTAN=LOCALTAN
            ELSE
               CYCLE L231
            ENDIF     
            COMPSLOPE=(ELEVATION(IH,JL)-ELEVATION(IL,JH))/(1.414*DDIV)
            COMPANGLE=ATAN(COMPSLOPE)
            ANGLE_DIFF=LOCALANGLE-COMPANGLE
            IF (ANGLE_DIFF > 0.0) THEN
                ANGLE_FACTOR=(SIN(0.5*ANGLE_DIFF))
                IF (USE_COSINE_SQUARE) ANGLE_FACTOR=ANGLE_FACTOR**2  !!NEW
            ELSE
                ANGLE_FACTOR=0.0
            ENDIF
            DISTANCE=SQRT((ABS(K)*CELL_SIZE*1.414)**2+(ELEVATION(ICOMP,JCOMP)-ELEVATION(I,J))**2)   !!NEW            
            CALL FIND_RAD_CHANGE(ICOMP,JCOMP,ANGLE_FACTOR,DISTANCE,RAD_CHANGE,RAD_WEIGHT)
            SUMCHANGE=SUMCHANGE+RAD_CHANGE/DIRECTIONBIAS !!!
            SUMWEIGHT=SUMWEIGHT+RAD_WEIGHT
        ENDDO L231
        IL=I-1
        IH=I+1
        DDIV=2.0*CELL_SIZE
        IF (IS_X_PERIODIC) THEN
            IF (IL < 1) IL=IL+MX
            IF (IH > MX) IH=IH-MX
        ELSE
            IF (IL < 1) THEN
                IL=1
                DDIV=1.0*CELL_SIZE
            ENDIF
            IF (IH > MX) THEN
                IH=MX
                DDIV=1.0*CELL_SIZE
            ENDIF
        ENDIF
        JL=J-1
        JH=J+1
        IF (IS_Y_PERIODIC) THEN
            IF (JL < 1) JL=JL+MY
            IF (JH > MY) JH=JH-MY
        ELSE
            IF (JL < 1) THEN
                JL=1
                DDIV=1.0*CELL_SIZE
            ENDIF
            IF (JH > MY) THEN
                JH=MY
                DDIV=1.0*CELL_SIZE
            ENDIF
        ENDIF
        LOCALMAXTAN=-1000.0
        LOCALSLOPE=(ELEVATION(IL,JL)-ELEVATION(IH,JH))/(1.414*DDIV)
        LOCALANGLE=ATAN(LOCALSLOPE)
        L240: DO K=1,LOCAL_WEIGHT_DISTANCE
            IF (K == 0) CYCLE L240
            DDIV=2.0*CELL_SIZE
            ICOMP=I-K
            JCOMP=J-K
            IF (IS_X_PERIODIC) THEN
                IF (JCOMP < 1) JCOMP=JCOMP+MY
                IF (JCOMP > MY) JCOMP=JCOMP-MY
            ELSE
                IF (JCOMP < 1) CYCLE L240
                IF (JCOMP > MY) CYCLE L240
            ENDIF
            IF (IS_X_PERIODIC) THEN
                IF (ICOMP < 1) ICOMP=ICOMP+MX
                IF (ICOMP > MX) ICOMP=ICOMP-MX
            ELSE
                IF (ICOMP < 1) CYCLE L240
                IF (ICOMP > MX) CYCLE L240
            ENDIF
            JL=JCOMP-1
            JH=JCOMP+1
            IF (IS_X_PERIODIC) THEN
                IF (JL < 1) JL=JL+MY
                IF (JH > MY) JH=JH-MY
                DDIV=2.0*CELL_SIZE
            ELSE
                IF (JL < 1) THEN
                    JL=1
                    DDIV=1.0*CELL_SIZE
                ENDIF
                IF (JH > MY) THEN
                    JH=MY
                    DDIV=1.0*CELL_SIZE
                ENDIF
            ENDIF
            IL=ICOMP-1
            IH=ICOMP+1
            IF (IS_X_PERIODIC) THEN
                IF (IL < 1) IL=IL+MX
                IF (IH > MX) IH=IH-MX
                DDIV=2.0*CELL_SIZE
            ELSE
                IF (IL < 1) THEN
                    IL=1
                    DDIV=1.0*CELL_SIZE
                ENDIF
                IF (IH > MX) THEN
                    IH=MX
                    DDIV=1.0*CELL_SIZE
                ENDIF
            ENDIF
            LOCALTAN=(ELEVATION(ICOMP,JCOMP)-ELEVATION(I,J))/(ABS(K)*CELL_SIZE*1.414)
            if (LOCALTAN.GE.LOCALMAXTAN) THEN
               LOCALMAXTAN=LOCALTAN               
            ELSE
               CYCLE L240
            ENDIF  
            COMPSLOPE=(ELEVATION(IL,JL)-ELEVATION(IH,JH))/(1.414*DDIV)
            COMPANGLE=ATAN(COMPSLOPE)
            ANGLE_DIFF=COMPANGLE-LOCALANGLE
            IF (ANGLE_DIFF > 0.0) THEN
                ANGLE_FACTOR=(SIN(0.5*ANGLE_DIFF))
                IF (USE_COSINE_SQUARE) ANGLE_FACTOR=ANGLE_FACTOR**2  !!NEW
            ELSE
                ANGLE_FACTOR=0.0
            ENDIF
            DISTANCE=SQRT((ABS(K)*CELL_SIZE*1.414)**2+(ELEVATION(ICOMP,JCOMP)-ELEVATION(I,J))**2)   !!NEW           
            CALL FIND_RAD_CHANGE(ICOMP,JCOMP,ANGLE_FACTOR,DISTANCE,RAD_CHANGE,RAD_WEIGHT)
            SUMCHANGE=SUMCHANGE+RAD_CHANGE/DIRECTIONBIAS !!!
            SUMWEIGHT=SUMWEIGHT+RAD_WEIGHT
        ENDDO L240
        LOCALMAXTAN=-1000.0
        L241: DO K=-1,-LOCAL_WEIGHT_DISTANCE,-1
            IF (K == 0) CYCLE L241
            DDIV=2.0*CELL_SIZE
            ICOMP=I-K
            JCOMP=J-K
            IF (IS_X_PERIODIC) THEN
                IF (JCOMP < 1) JCOMP=JCOMP+MY
                IF (JCOMP > MY) JCOMP=JCOMP-MY
            ELSE
                IF (JCOMP < 1) CYCLE L241
                IF (JCOMP > MY) CYCLE L241
            ENDIF
            IF (IS_X_PERIODIC) THEN
                IF (ICOMP < 1) ICOMP=ICOMP+MX
                IF (ICOMP > MX) ICOMP=ICOMP-MX
            ELSE
                IF (ICOMP < 1) CYCLE L241
                IF (ICOMP > MX) CYCLE L241
            ENDIF
            JL=JCOMP-1
            JH=JCOMP+1
            IF (IS_X_PERIODIC) THEN
                IF (JL < 1) JL=JL+MY
                IF (JH > MY) JH=JH-MY
                DDIV=2.0*CELL_SIZE
            ELSE
                IF (JL < 1) THEN
                    JL=1
                    DDIV=1.0*CELL_SIZE
                ENDIF
                IF (JH > MY) THEN
                    JH=MY
                    DDIV=1.0*CELL_SIZE
                ENDIF
            ENDIF
            IL=ICOMP-1
            IH=ICOMP+1
            IF (IS_X_PERIODIC) THEN
                IF (IL < 1) IL=IL+MX
                IF (IH > MX) IH=IH-MX
                DDIV=2.0*CELL_SIZE
            ELSE
                IF (IL < 1) THEN
                    IL=1
                    DDIV=1.0*CELL_SIZE
                ENDIF
                IF (IH > MX) THEN
                    IH=MX
                    DDIV=1.0*CELL_SIZE
                ENDIF
            ENDIF
            LOCALTAN=(ELEVATION(ICOMP,JCOMP)-ELEVATION(I,J))/(ABS(K)*CELL_SIZE*1.414)
            if (LOCALTAN.GE.LOCALMAXTAN) THEN
               LOCALMAXTAN=LOCALTAN                
            ELSE
               CYCLE L241
            ENDIF
            COMPSLOPE=(ELEVATION(IL,JL)-ELEVATION(IH,JH))/(1.414*DDIV)
            COMPANGLE=ATAN(COMPSLOPE)
            ANGLE_DIFF=LOCALANGLE-COMPANGLE
            IF (ANGLE_DIFF > 0.0) THEN
                ANGLE_FACTOR=(SIN(0.5*ANGLE_DIFF))
                IF (USE_COSINE_SQUARE) ANGLE_FACTOR=ANGLE_FACTOR**2  !!NEW
            ELSE
                ANGLE_FACTOR=0.0
            ENDIF
            DISTANCE=SQRT((ABS(K)*CELL_SIZE*1.414)**2+(ELEVATION(ICOMP,JCOMP)-ELEVATION(I,J))**2)   !!NEW
            CALL FIND_RAD_CHANGE(ICOMP,JCOMP,ANGLE_FACTOR,DISTANCE,RAD_CHANGE,RAD_WEIGHT)
            SUMCHANGE=SUMCHANGE+RAD_CHANGE/DIRECTIONBIAS !!!
            SUMWEIGHT=SUMWEIGHT+RAD_WEIGHT
        ENDDO L241   
        RAD_ERODE=SUMCHANGE
        RETURN
    END !   REAL (4) FUNCTION RAD_ERODE(I,J)
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DEPOSIT_ICE()
        USE ERODE_GLOBALS
        USE ACCRETION_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,IMAX,JMAX,IMIN,JMIN
        REAL(4) :: EXPOSEAVG,EXPOSEMAX,EXPOSEMIN,EXPOSURE,EXPOSEFRACT
        REAL(4) :: COSFACTOR,MAXCOSFACT,TOTDEPOSITICE,SCALEFACTOR
        !     **********************************************************************
        !      This deposits ice on a planetary surface based upon the amount of reradiated infrared energy
        !        reaching the point from surrounding cells.  If that reradiated energy is greater than a
        !        critical value no deposition occurs.  The routine calls topexpose to determine reradiation
        !        intensity to the given point.
        !      An ice-covered surface is indicated by IS_SEDIMENT_COVERED being TRUE and IS_ROCK_SURFACE is FALSE
        !      A bedrock surface is indicate be IS_ROCK_SURFACE being TRUE
        !      A dust (sediment-covered) surface is indicated by IS_ROCK_SURFACE being FALSE and IS_SEDIMENT_COVERED being FALSE
        !  MODIFIES:  ELEVATION, IS_SEDIMENT_COVERED, SEDIMENT_BASE,IS_ROCK_SURFACE, CFNW 
        !  CALLS: CALCULATE_DIVERGENCE, FIND_TOP_EXPOSURE_INDEX
        !     **********************************************************************

        IF (RAD_DEPOSIT_RATE <= 0.0) RETURN
        EXPOSEAVG=0.0
        EXPOSEMAX=-1.0E25
        EXPOSEMIN=-EXPOSEMAX
        EXPOSEFRACT=0.0
        CFNW=0.0
        TOTDEPOSITICE=0.0
        IF (ITERATION == 1) THEN
            CALL CALCULATE_DIVERGENCE()
            DO J=1,MY
                DO I=1,MX
                    CALL FIND_TOP_EXPOSURE_INDEX(I,J,EXPOSURE)
                    IF (EXPOSURE > EXPOSEMAX) THEN
                        IMAX=I
                        JMAX=J
                        EXPOSEMAX=EXPOSURE
                    ENDIF
                    IF (EXPOSURE < EXPOSEMIN) THEN
                        IMIN=I
                        JMIN=J
                        EXPOSEMIN=EXPOSURE
                    ENDIF
                    EXPOSEAVG=EXPOSEAVG+EXPOSURE
                ENDDO
            ENDDO
            EXPOSEAVG=EXPOSEAVG/FLOAT(MX*MY)
            !WRITE(OUTHIST,110) IMAX,JMAX,EXPOSEMAX,DIVERGENCE(IMAX,JMAX)
            110   FORMAT(' IMAX=',I5,' JMAX=',I5,' MAX. EXPOSE=',G12.5,' DIVERGE=',G12.5)
            !WRITE(OUTHIST,120) IMIN,JMIN,EXPOSEMIN,DIVERGENCE(IMIN,JMIN)
            120   FORMAT(' IMIN=',I5,' JMIN=',I5,' MIN. EXPOSE=',G12.5,' DIVERGE=',G12.5)
            !WRITE(OUTHIST,140) EXPOSEAVG
            140   FORMAT('AVG. EXPOSURE=',G12.5)
        ENDIF
        MAXCOSFACT=10.0
        DO J=1,MY
            DO I=1,MX
                COSFACTOR=MIN(MAXCOSFACT,SQRT(1.0+D8_GRADIENT(I,J)**2))
                IF (ELEVATION(I,J) < SEDIMENT_BASE(I,J)) THEN
                    IS_SEDIMENT_COVERED(I,J)=.FALSE.
                    SEDIMENT_BASE(I,J)=ELEVATION(I,J)
                ENDIF
                CALL FIND_TOP_EXPOSURE_INDEX(I,J,EXPOSURE)
                IF (EXPOSURE < RAD_THRESH_CONVEXITY) THEN
                    EXPOSEFRACT=EXPOSEFRACT+1.0
                    IF (.NOT.IS_SEDIMENT_COVERED(I,J)) THEN
                        IS_SEDIMENT_COVERED(I,J)=.TRUE.
                        SEDIMENT_BASE(I,J)=ELEVATION(I,J)
                        IS_ROCK_SURFACE(I,J)=.FALSE.
                    ENDIF
                    CFNW(I,J)=COSFACTOR*RAD_DEPOSIT_RATE
                    TOTDEPOSITICE=TOTDEPOSITICE+CFNW(I,J)
                ENDIF
            ENDDO
        ENDDO
        SCALEFACTOR=TOTALRADERODE*DEPOSITFRACTION*TIME_INCREMENT/TOTDEPOSITICE
        WRITE(*,1099) SCALEFACTOR
        1099  FORMAT(' SCALEFACTOR=',G12.5)
        DO J=1,MY
            DO I=1,MX
                ELEVATION(I,J)=ELEVATION(I,J)+CFNW(I,J)*SCALEFACTOR
                ACCRETE_WORK=ACCRETE_WORK+CFNW(I,J)*SCALEFACTOR
            ENDDO
        ENDDO
        EXPOSEFRACT=EXPOSEFRACT/FLOAT(MX*MY)
        WRITE(*,150) EXPOSEFRACT
        150   FORMAT(' FRACT. ICE DEPOSIT=',G12.5)
        RETURN
    END ! SUBROUTINE DEPOSIT_ICE()
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE FIND_RAD_CHANGE(ICOMP,JCOMP,ANGLE_FACTOR,DISTANCE,RAD_CHANGE,RAD_WEIGHT)
        USE ERODE_GLOBALS
        USE ACCRETION_GLOBALS
        IMPLICIT NONE
        REAL (4), INTENT(IN) :: ANGLE_FACTOR,DISTANCE
        REAL (4),INTENT(INOUT) ::RAD_WEIGHT
        REAL (4), INTENT(OUT) :: RAD_CHANGE
        INTEGER, INTENT(IN) :: ICOMP,JCOMP
        !     ************************************************************************
        !      This determines the net radiation ablatin-deposition contribution from a given cell
        !        depending upon distance, relative orientation of the cell relative to the local cell,
        !        and whether the given cell is 'bedrock' or sediment-covered.
        !          see howard and moore, 2008, geophysical research letters, 35, l03203,
        !        doi:10.1029/2007gl032618
        !     ************************************************************************
        RAD_WEIGHT=1.0/DISTANCE  ! to compensate
        IF (ANGLE_FACTOR > 0.0) THEN
            IF (IS_ROCK_SURFACE(ICOMP,JCOMP).OR.DO_REGOLITH_ABLATION) THEN
                RAD_CHANGE=RAD_WEIGHT*ANGLE_FACTOR*VOLUME_CHANGE_COEFFICIENT &
                    * COS(ATAN(ABS(D8_GRADIENT(ICOMP,JCOMP)))) !EXPERIMENTAL
            ELSE
                RAD_CHANGE=RAD_DUST_FACTOR*RAD_WEIGHT*ANGLE_FACTOR
            ENDIF
        ELSE
            RAD_CHANGE=0.0
        ENDIF
        RETURN
    END ! SUBROUTINE DEPOSIT_ICE()
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
