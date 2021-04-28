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
	!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	SUBROUTINE READ_GROUNDWATER_PARAMETERS
	USE ERODE_GLOBALS
    USE GROUNDWATER_VARIABLES
	IMPLICIT NONE
	     CHARACTER (80) :: TEXTLINE
	    INTEGER :: GROUNDWATER_PARAMS, SEEPUSE, IEXPFLOW, IWATERLOWER, ISHOWCALC, QQUSE, IDIVAVG
        INTEGER :: IERROR
        REAL(4) :: TEMP1
		GROUNDWATER_PARAMS=93
		OPEN(GROUNDWATER_PARAMS,FILE='GROUNDWATER_PARAMETERS.PRM',ACTION='READ')
       !        *******************************************************************
        !        *********************** groundwater flow **************************
        !        *******************************************************************
        READ(GROUNDWATER_PARAMS,22) TEXTLINE
22      FORMAT(A)        
        WRITE(*,22) TEXTLINE
        READ (GROUNDWATER_PARAMS,*) SEEPUSE
        READ (GROUNDWATER_PARAMS,*) SEEPAGE_ITERATION_INTERVAL
        READ (GROUNDWATER_PARAMS,*) IWATERLOWER
        READ (GROUNDWATER_PARAMS,*) ISHOWCALC
        READ (GROUNDWATER_PARAMS,*) IEXPFLOW
        IF (SEEPUSE > 0) THEN
            MODEL_GROUNDWATER=.TRUE.
            ALLOCATE(REFERENCE_ELEVATION(MX,MY))
            REFERENCE_ELEVATION=INITIAL_ELEVATION
        ELSE
            MODEL_GROUNDWATER=.FALSE.
        ENDIF
        IF (MODEL_GROUNDWATER) THEN
            ALLOCATE(TRANSMISSIVITY_TERM(MX,MY),WATER_ELEVATION(MX,MY),STAT=IERROR)
			IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(FILTERED_GROUNDWATER_FLUX(MX,MY),STAT=IERROR)
			IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(GROUNDWATER_FLUX(MX,MY),STAT=IERROR)
			IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(FIXED_HEAD(MX,MY),STAT=IERROR)
			IF (IERROR > 0) CALL ALLOCERROR()
        ENDIF
        IF (IEXPFLOW > 0) THEN
            EXPONENTIAL_PERMEABILITY_DECAY=.TRUE.
        ELSE
            EXPONENTIAL_PERMEABILITY_DECAY=.FALSE.
        ENDIF
        IF (IWATERLOWER > 0) THEN
            PERMEABILITY_RESCALING=.TRUE.
        ELSE
            PERMEABILITY_RESCALING=.FALSE.
        ENDIF
        IF (ISHOWCALC > 0) THEN
            SHOW_GROUNDWATER_CALCULATIONS=.TRUE.
        ELSE
            SHOW_GROUNDWATER_CALCULATIONS=.FALSE.
        ENDIF
        READ(GROUNDWATER_PARAMS,*) YEARLY_RECHARGE
        READ (GROUNDWATER_PARAMS,*) VISCOSITY
        READ (GROUNDWATER_PARAMS,*) DARCIES
        READ (GROUNDWATER_PARAMS,*) GROUNDWATER_DEPTH_SCALE
        READ (GROUNDWATER_PARAMS,*) GROUNDWATER_FLOW_FRACTION
        READ (GROUNDWATER_PARAMS,*) INITIAL_GROUNDWATER_DEPTH
        READ (GROUNDWATER_PARAMS,*) EPOWER
        IF (EXPONENTIAL_PERMEABILITY_DECAY) THEN
            INITIAL_GROUNDWATER_DEPTH=GROUNDWATER_DEPTH_SCALE*INITIAL_GROUNDWATER_DEPTH
            GROUNDWATER_DEPTH_SCALE=-LOG(0.5E0)/(GROUNDWATER_DEPTH_SCALE**EPOWER)
        ENDIF
        GROUNDWATER_RECHARGE_RATE=YEARLY_RECHARGE/SECONDS_PER_YEAR
        METRIC_PERMEABILITY=DARCIES/1.0E+12
        HYDRAULIC_CONDUCTIVITY=METRIC_PERMEABILITY*WATER_DENSITY*GRAVITY/VISCOSITY
        GROUNDWATER_SCALE_FACTOR=HYDRAULIC_CONDUCTIVITY/(CELL_SIZE*CELL_SIZE*2.0)
        IF (EXPONENTIAL_PERMEABILITY_DECAY) GROUNDWATER_SCALE_FACTOR=GROUNDWATER_SCALE_FACTOR/GROUNDWATER_DEPTH_SCALE
        TEMP1=GROUNDWATER_RECHARGE_RATE/GROUNDWATER_SCALE_FACTOR
        IF (MODEL_GROUNDWATER) THEN
           WRITE(*,4215) HYDRAULIC_CONDUCTIVITY,GROUNDWATER_SCALE_FACTOR,TEMP1
           WRITE(OUTHIST,4215) HYDRAULIC_CONDUCTIVITY,GROUNDWATER_SCALE_FACTOR,TEMP1
4215       FORMAT(' HYDRAULIC CONDUCTIVITY=',G12.5,' GROUNDWATER_SCALE_FACTOR=',G12.5, &
           ' GROUNDWATER_RECHARGE_RATE/GROUNDWATER_SCALE_FACTOR=',G12.5)
        ENDIF
        READ(GROUNDWATER_PARAMS,*) MAXIMUM_GROUNDWATER_ITERATIONS
        READ (GROUNDWATER_PARAMS,*) MAXIMUM_GROUNDWATER_ERROR
        READ (GROUNDWATER_PARAMS,*) GROUNDWATER_RELAXATION_CONST
        IF (MODEL_GROUNDWATER.AND.(GROUNDWATER_FLOW_FRACTION > 0.5)) THEN
            SEDIMENT_DISCHARGE_FACTOR=SEDIMENT_DISCHARGE_FACTOR/EFFECTIVE_DISCHARGE_RATIO
            SEDIMENT_CONSTANT=SEDIMENT_CONSTANT/FLOW_FRACTION
            SKLAR_MULT=SKLAR_MULT/FLOW_FRACTION
        ENDIF
       ! WRITE(OUTHIST,32001) SEDIMENT_DISCHARGE_FACTOR,SEDIMENT_CONSTANT,SKLAR_MULT
       ! WRITE(*,32001) SEDIMENT_DISCHARGE_FACTOR,SEDIMENT_CONSTANT,SKLAR_MULT
        32001 FORMAT('CALCULATED SEDIMENT_DISCHARGE_FACTOR=',G12.5,' SEDIMENT_CONSTANT=',G12.5 &
        ,' SLARMULT=',G12.5)
        READ (GROUNDWATER_PARAMS,*) QQUSE
        READ (GROUNDWATER_PARAMS,*) IDIVAVG
		CLOSE (GROUNDWATER_PARAMS)
        IF (QQUSE > 0) THEN
            USE_GROUNDWATER_FLOW=.TRUE.
        ELSE
            USE_GROUNDWATER_FLOW=.FALSE.
        ENDIF
        IF (IDIVAVG > 0) THEN
            SEEPAGE_AVERAGING=.TRUE.
        ELSE
            SEEPAGE_AVERAGING=.FALSE.
        ENDIF
       !        *******************************************************************
        !        **************** groundwater flow**********************************
        !        *******************************************************************
        IF (MODEL_GROUNDWATER) THEN
             WRITE(OUTHIST,1525)
             1525 FORMAT(' **************** GROUNDWATER FLOW PARAMETERS ******************')
             WRITE(OUTHIST,1526) SEEPAGE_ITERATION_INTERVAL
             1526 FORMAT(' GROUNDWATER FLOW CALCULATED EVERY ',I6,' ITERATIONS')
             IF (EXPONENTIAL_PERMEABILITY_DECAY) THEN
                 WRITE(OUTHIST,1527)
                 1527 FORMAT(' PERMEABILITY DECREASES EXPONENTIALLY WITH DEPTH BELOW SURFACE')
             ELSE
                 WRITE(OUTHIST,1528)
                 1528 FORMAT(' PERMEABILITY IS CONSTANT WITH DEPTH')
             ENDIF
             IF (PERMEABILITY_RESCALING) THEN
                 WRITE(OUTHIST,1529)
                 1529 FORMAT(' PERMEABILITY IS REFERENCED TO CURRENT LAND SURFACE ELEVATION')
             ELSE
                 WRITE(OUTHIST,1530)
                 1530 FORMAT(' PERMEABILITY IS REFERENCED TO INITIAL LAND SURFACE ELEVATION')
             ENDIF
             IF (SHOW_GROUNDWATER_CALCULATIONS) THEN
                 WRITE(OUTHIST,1531)
                 1531 FORMAT(' DETAILS OF GROUNDWATER CALCULATIONS ARE PRINTED OUT')
             ENDIF
             WRITE(OUTHIST,1532) YEARLY_RECHARGE,VISCOSITY,DARCIES,GROUNDWATER_DEPTH_SCALE
             1532 FORMAT(' YEARLY RECHARGE=',G12.5,' GROUNDWATER VISCOSITY=',G12.5, &
                  ' DARCIES=',G12.5, &
             ' GROUNDWATER DEPTH SCALE=',G12.5)
             WRITE(OUTHIST,1533) GROUNDWATER_FLOW_FRACTION, INITIAL_GROUNDWATER_DEPTH, EPOWER
             1533 FORMAT(' FRACTION OF SURFACE WATER FLOW THAT IS GROUNDWATER=',G12.5,/, &
             ' INITIAL GROUNDWATER DEPTH BENEATH SURFACE=',G12.5,' RATE OF EXPOENTIAL DECAY=',G12.5)
             WRITE(*,4215) HYDRAULIC_CONDUCTIVITY,GROUNDWATER_SCALE_FACTOR,GROUNDWATER_RECHARGE_RATE
             WRITE(OUTHIST,14215) HYDRAULIC_CONDUCTIVITY,GROUNDWATER_SCALE_FACTOR,GROUNDWATER_RECHARGE_RATE
             14215  FORMAT(' HYDRAULIC CONDUCTIVITY=',G12.5,' GROUNDWATER_SCALE_FACTOR=',G12.5, &
             ' GROUNDWATER_RECHARGE_RATE=',G12.5)
             WRITE(OUTHIST,1534) MAXIMUM_GROUNDWATER_ITERATIONS, MAXIMUM_GROUNDWATER_ERROR,GROUNDWATER_RELAXATION_CONST
             1534 FORMAT(' MAXIMUM NUMBER OF ITERATIONS FOR FLOW CONVERGENCE=',I6,/, &
             ' MAXIMUM PERMITTED ERROR=',G12.5,' CALCULATION CONVERGENCE FACTOR=',G12.5)
             IF (GROUNDWATER_FLOW_FRACTION > 0.5) THEN
                 WRITE(OUTHIST,1535)
                 1535 FORMAT(' SEDIMENT TRANSPORT PARAMETERS RECALCULATED FOR GROUNDWATER FLOW')
                 WRITE(OUTHIST,32001) SEDIMENT_DISCHARGE_FACTOR,SEDIMENT_CONSTANT,SKLAR_MULT
                 WRITE(*,42001) SEDIMENT_DISCHARGE_FACTOR,SEDIMENT_CONSTANT,SKLAR_MULT
                 42001 FORMAT('CALCULATED SEDIMENT_DISCHARGE_FACTOR=',G12.5,' SEDIMENT_CONSTANT=',G12.5 &
                 ,' SLARMULT=',G12.5)
             ENDIF
             IF (USE_GROUNDWATER_FLOW) THEN
                 WRITE(OUTHIST,1536)
                 1536 FORMAT (' SEEPAGE MODELING USES GROUNDWATER FLUX (SUBSURFACE DISCHARGE)')
             ELSE
                 WRITE(OUTHIST,1537)
                 1537 FORMAT (' SEEPAGE MODELING USES GROUNDWATER FLUX DIVERGENCE (FLOW TO SURFACE)')
             ENDIF
             IF (SEEPAGE_AVERAGING) THEN
                 WRITE(OUTHIST,1538)
                 1538 FORMAT(' GROUNDWATER FLOW MEASUREMENTS ARE SPATIALLY AVERAGED')
             ENDIF
        ENDIF		
	RETURN
	END !SUBROUTINE READ_GROUNDWATER_PARAMETERS
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
    SUBROUTINE EXPONENTIAL_HYDR_COND_GRNDWTR()
        !      ********************************************************************
        !
        !      This solves for the water table and calculates flow discharge, etc.
        !
        !      This subroutine is used when the hydraulic conductivity declines exponentially with depth
        !
        !      The model details and the use of groundwater flow for groundwater seepage erosion and for groundwater
        !      weathering are presented in:
        !          luo, w., and howard, a. d., 2008, j. geophysical research, 113,
        !          e05002,doi:10.1029/2007je002981
        !  MODIFIES: FIXED_HEAD, WATER_ELEVATION, TRANSMISSIVITY_TERM, GROUNDWATER_FLUX
        !            FILTERED_GROUNDWATER_FLUX
        !
        !      ********************************************************************
        !     sweep interior nodes with finite difference form of poisson's eqn.
        USE ERODE_GLOBALS
        USE GROUNDWATER_VARIABLES
        IMPLICIT NONE
        INTEGER :: I,J,NUMIT
        REAL(4) :: NTOTAL
        REAL(4) :: AMAX,OLDVAL,ERR,QTOT
        REAL(4) :: EWMIN,EWMAX,VCMIN,VCMAX,QCMAX,QCAVG,NQQ,QCMIN
        REAL(4) :: QQ1,QQ2,QQ3,QQ4,QQ5,QQ6,QQ7,QQ8,QQ9,DX2
        REAL(4) :: CIP,CIM,CJP,CJM,ELOCAL
        REAL(4) :: QX,QY,MINEL,CELL_AREAX,FRACTFIXED,AVGHEFFECTIVE
        REAL(4) ::  WATERNEW(IMMX,JMMX)
        REAL(4) :: MAXWATERCHANGE,QQQ,CRITICAL_CHANGE
        INTEGER :: NFIXED,LASTNFIXED
        INTEGER :: IW,IE,JN,JS
        LOGICAL :: DOSTOP,CHECKRELAX
        !!!!!!!!write(outhist,4833)
4833    format('groundwater called')     
        MAXWATERCHANGE=2.0
        AVGHEFFECTIVE=0.0
        DX2=CELL_SIZE*2.0
        CRITICAL_CHANGE=0.02
        CELL_AREAX=2.0*CELL_SIZE*CELL_SIZE
        VCMIN=1.0E+25
        VCMAX=-VCMIN
        EWMIN=VCMIN
        EWMAX=VCMAX
        MINEL=1.0E+25
        NUMIT=0
        NFIXED=0
        LASTNFIXED=0
        NTOTAL=0.0
        DOSTOP=.FALSE.
        L1997:      DO  J=1,MY
            L1998:      DO  I=1,MX
                IF (ELEVATION(I,J) < MINEL) THEN
                    MINEL=ELEVATION(I,J)
                    ILOWEST=I
                    JLOWEST=J
                ENDIF
            ENDDO L1998
        ENDDO L1997
        !      ********************************************************************
        !
        !      Initially only the point of lowest elevation (or the lower boundary
        !         if it is fixed) is fixed_head elevation, and the water table elevation
        !         (water_elevation) is set to that value.  the initial transmissivity factor
        !         is calculated for this flat water table.  Alternatively the depth to the water table
        !         can be initialized to be a fixed depth beneath the surface, except at the lowest location
        !
        !      ********************************************************************
        if (first_groundwater_call) write(*,2299)
2299    format('first groundwater call')
        L700: DO  J=1,MY
            L701: DO  I=1,MX
                FIXED_HEAD(I,J)=.FALSE.
                IF (FIRST_GROUNDWATER_CALL) THEN
                       WATER_ELEVATION(I,J)=MAX(MINEL,(ELEVATION(I,J)-INITIAL_GROUNDWATER_DEPTH))
                       FIXED_HEAD(I,J)=.FALSE.
                ELSE
                    IF (WATER_ELEVATION(I,J)>=ELEVATION(I,J)) THEN
                        ! FIXED_HEAD(I,J)=.TRUE.
                        WATER_ELEVATION(I,J)=ELEVATION(I,J)
                    ELSE
                        !FIXED_HEAD(I,J)=.FALSE.
                    ENDIF
                ENDIF
                IF (PERMEABILITY_RESCALING.OR.(ELEVATION(I,J) > INITIAL_ELEVATION(I,J))) THEN
                    TRANSMISSIVITY_TERM(I,J)=EXP(-GROUNDWATER_DEPTH_SCALE*(ELEVATION(I,J)-WATER_ELEVATION(I,J))**EPOWER)
                ELSE
                    TRANSMISSIVITY_TERM(I,J)=EXP(-GROUNDWATER_DEPTH_SCALE*(INITIAL_ELEVATION(I,J)-WATER_ELEVATION(I,J))**EPOWER)
                ENDIF
                IF (ELEVATION(I,J) == MINEL) &
                FIXED_HEAD(I,J)=.TRUE.
                IF (WATER_ELEVATION(I,J) > EWMAX) EWMAX=WATER_ELEVATION(I,J)
                IF (WATER_ELEVATION(I,J) < EWMIN) EWMIN=WATER_ELEVATION(I,J)
                IF (FIXED_HEAD(I,J)) NFIXED=NFIXED+1
                NTOTAL=NTOTAL+1
            ENDDO L701
        ENDDO L700
        FIXED_HEAD(ILOWEST,JLOWEST)=.TRUE.
        FRACTFIXED=NFIXED/NTOTAL
        WRITE(OUTHIST,721) EWMIN,EWMAX,FRACTFIXED
        721   FORMAT(' EWMIN=',G12.5,' EWMAX=',G12.5,' % FIXED=',G12.5)
        NFIXED=0.0
        NTOTAL=0.0
        !      ********************************************************************
        !
        !      This is the main cycle for the iterative solution
        !
        !      ********************************************************************
        L045: DO NUMIT=1,MAXIMUM_GROUNDWATER_ITERATIONS
            AMAX=0.0
            !      ********************************************************************
            !
            !      Use standard finite difference sweep with overrelaxation
            !
            !      ********************************************************************
            L040:  DO  J=1,MYY
                L041: DO  I=1,MX
                    IF (FIXED_HEAD(I,J)) THEN
                        WATERNEW(I,J)=ELEVATION(I,J)
                        CYCLE L041
                    ENDIF
                    OLDVAL=WATER_ELEVATION(I,J)
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
                        If (JS > MY) JS=MY-1
                    ENDIF
                    ELOCAL=TRANSMISSIVITY_TERM(I,J)
                    CIP=(ELOCAL+TRANSMISSIVITY_TERM(IW,J))
                    CIM=(ELOCAL+TRANSMISSIVITY_TERM(IE,J))
                    CJP=(ELOCAL+TRANSMISSIVITY_TERM(I,JS))
                    CJM=(ELOCAL+TRANSMISSIVITY_TERM(I,JN))
                    WATERNEW(I,J)=(CIP*WATER_ELEVATION(IW,J)+CIM*WATER_ELEVATION(IE,J)+ &
                    CJP*WATER_ELEVATION(I,JS) &
                    +CJM*WATER_ELEVATION(I,JN)+GROUNDWATER_RECHARGE_RATE/GROUNDWATER_SCALE_FACTOR)/(CIP+CIM+CJP+CJM)
                    ERR=ABS(WATERNEW(I,J)-OLDVAL)
                    IF (ERR > AMAX) AMAX=ERR
                ENDDO L041
            ENDDO L040
            IF (AMAX > MAXWATERCHANGE) THEN
                L049: DO  J=1,MYY
                    L050:  DO  I=1,MX
                        IF (.NOT.FIXED_HEAD(I,J)) THEN
                            WATER_ELEVATION(I,J)=WATER_ELEVATION(I,J)  &
                            +MAXWATERCHANGE*(WATERNEW(I,J)-WATER_ELEVATION(I,J))/AMAX
                        ENDIF
                    ENDDO L050
                ENDDO L049
            ELSE
                WATER_ELEVATION=WATERNEW
            ENDIF
            !      ********************************************************************
            !
            !      If the predicted water table elevation is greater than the actual
            !         elevation, that point becomes a fixed_head head location
            !
            !      ********************************************************************
            NFIXED=0
            L043:  DO  J=1,MYY
                L044: DO  I=1,MX
                    IF (WATER_ELEVATION(I,J) >= ELEVATION(I,J)) THEN
                        NFIXED=NFIXED+1
                        FIXED_HEAD(I,J)=.TRUE.
                        WATER_ELEVATION(I,J)=ELEVATION(I,J)
                    ENDIF
                    !      ********************************************************************
                    !
                    !      Recalculate transmissivity factor
                    !
                    !      ********************************************************************
                ENDDO L044
            ENDDO L043
            IF (MOD(NUMIT,100) == 0) THEN
                IF (AMAX<=MAXIMUM_GROUNDWATER_ERROR) THEN
                    DOSTOP=.TRUE.
                    IF (SHOW_GROUNDWATER_CALCULATIONS) WRITE(*,2811) NFIXED,LASTNFIXED
                ENDIF
                IF (NFIXED == LASTNFIXED) THEN
                  !  DOSTOP=.TRUE.
                 !   IF (SHOW_GROUNDWATER_CALCULATIONS) WRITE(*,2811) NFIXED,LASTNFIXED
                    2811  FORMAT(' NFIXED=',I6,' LASTNFIXED=',I6)
                ENDIF
                LASTNFIXED=NFIXED
                CHECKRELAX=.TRUE.
            ELSE
                CHECKRELAX=.FALSE.
            ENDIF
            L731:          DO  J=1,MYY
                L732:          DO  I=1,MX
                    IF (PERMEABILITY_RESCALING.OR.(ELEVATION(I,J) > INITIAL_ELEVATION(I,J))) THEN
                        TRANSMISSIVITY_TERM(I,J)=  &
                        EXP(-GROUNDWATER_DEPTH_SCALE*(ELEVATION(I,J)-WATER_ELEVATION(I,J))**EPOWER)
                    ELSE
                        TRANSMISSIVITY_TERM(I,J)= &
                        EXP(-GROUNDWATER_DEPTH_SCALE*(INITIAL_ELEVATION(I,J)-WATER_ELEVATION(I,J))**EPOWER)
                    ENDIF
                ENDDO L732
            ENDDO L731
            IF (AMAX > MAXIMUM_GROUNDWATER_ERROR) THEN
                IF (SHOW_GROUNDWATER_CALCULATIONS) THEN
                    WRITE(*,1731) NUMIT,NFIXED, AMAX
                    1731     FORMAT('WITER=',I5,' FIXED=',I8,' AMAX=',G12.5)
                ENDIF
                NFIXED=0
                IF (CHECKRELAX) THEN
                    L832:        DO J=1,MYY
                        L833:        DO I=1,MX
                            IF (FIXED_HEAD(I,J)) THEN
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
                                    If (JS > MY) JS=MY-1
                                ENDIF
                                IF (WATER_ELEVATION(I,J) > MAX(WATER_ELEVATION(IW,J),WATER_ELEVATION(IE,J) &
                                ,WATER_ELEVATION(I,JN),WATER_ELEVATION(I,JS))) THEN
                                    QQQ=((WATER_ELEVATION(IW,J)-WATER_ELEVATION(I,J))   &
                                    *(TRANSMISSIVITY_TERM(I,J)+TRANSMISSIVITY_TERM(IW,J))  &
                                    -(WATER_ELEVATION(I,J)-WATER_ELEVATION(IE,J)) &
                                    *(TRANSMISSIVITY_TERM(I,J)+TRANSMISSIVITY_TERM(IE,J))) &
                                    +((WATER_ELEVATION(I,JS)-WATER_ELEVATION(I,J)) &
                                    *(TRANSMISSIVITY_TERM(I,J)+TRANSMISSIVITY_TERM(I,JS)) &
                                    -(WATER_ELEVATION(I,J)-WATER_ELEVATION(I,JN)) &
                                    *(TRANSMISSIVITY_TERM(I,J)+TRANSMISSIVITY_TERM(I,JN)))*GROUNDWATER_SCALE_FACTOR
                                    IF (QQQ < (-GROUNDWATER_RECHARGE_RATE))    THEN
                                        IF (SHOW_GROUNDWATER_CALCULATIONS) THEN 
                                            WRITE(*,2731) I,J,QQQ,ELEVATION(I,J),ELEVATION(IW,J),ELEVATION(IE,J) &
                                            ,ELEVATION(I,JN),ELEVATION(I,JS),WATER_ELEVATION(I,J),WATER_ELEVATION(IW,J) &
                                            ,WATER_ELEVATION(IE,J),WATER_ELEVATION(I,JN),WATER_ELEVATION(I,JS)
                                            2731  FORMAT('UNFIX,I=',I5,' J=',I5,' Q=',G12.5 &
                                            ,/,' E=',G12.5' EE=',G12.5,' EW=',G12.5,' EN=',G12.5,' ES=',G12.5 &
                                            ,/,' W=',G12.5' WE=',G12.5,' WW=',G12.5,' WN=',G12.5,' WS=',G12.5 &
                                            )
                                        ENDIF
                                        FIXED_HEAD(I,J)=.FALSE.
                                        WATER_ELEVATION(I,J)=MIN(ELEVATION(I,J),((WATER_ELEVATION(IW,J)+WATER_ELEVATION(IE,J) &
                                        +WATER_ELEVATION(I,JN)+WATER_ELEVATION(I,JS))/4.0))
                                    ENDIF
                                ENDIF
                            ENDIF
                        ENDDO L833
                    ENDDO L832
                ENDIF               
            ENDIF
            IF (DOSTOP) GOTO 45
        ENDDO L045
        45  CONTINUE
        !      ********************************************************************
        !
        !      Done.  Now calculate groundwater discharge and flow divergence
        !
        !      ********************************************************************
        EWMIN=1.0E+25
        EWMAX=-EWMIN
        QCMAX=-1.0E+25
        QCMIN=1.0E+25
        QCAVG=0.0
        NQQ=0.0
        L540:      DO  J=1,MYY
            L541:      DO  I=1,MX
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
                    If (JS > MY) JS=MY-1
                ENDIF
                IF (WATER_ELEVATION(I,J) > EWMAX) EWMAX=WATER_ELEVATION(I,J)
                IF (WATER_ELEVATION(I,J) < EWMIN) EWMIN=WATER_ELEVATION(I,J)
                QX=TRANSMISSIVITY_TERM(I,J)*(WATER_ELEVATION(IW,J)-WATER_ELEVATION(IE,J))/DX2
                QY=TRANSMISSIVITY_TERM(I,J)*(WATER_ELEVATION(I,JS)-WATER_ELEVATION(I,JN))/DX2
                !      ********************************************************************
                !
                !      If use_goundwater_flux is true the groundwater_flux is the groundwater flow rate, otherwise it
                !         is the flow divergence
                !
                !      ********************************************************************
                IF (USE_GROUNDWATER_FLOW) THEN
                    GROUNDWATER_FLUX(I,J)=HYDRAULIC_CONDUCTIVITY*SQRT(QX**2+QY**2)
                ELSE
                    GROUNDWATER_FLUX(I,J)=((WATER_ELEVATION(IW,J)-WATER_ELEVATION(I,J))  &
                    *(TRANSMISSIVITY_TERM(I,J)+TRANSMISSIVITY_TERM(IW,J))  &
                    -(WATER_ELEVATION(I,J)-WATER_ELEVATION(IE,J))   &
                    *(TRANSMISSIVITY_TERM(I,J)+TRANSMISSIVITY_TERM(IE,J))) &
                    +((WATER_ELEVATION(I,JS)-WATER_ELEVATION(I,J)) &
                    *(TRANSMISSIVITY_TERM(I,J)+TRANSMISSIVITY_TERM(I,JS))  &
                    -(WATER_ELEVATION(I,J)-WATER_ELEVATION(I,JN))  &
                    *(TRANSMISSIVITY_TERM(I,J)+TRANSMISSIVITY_TERM(I,JN)))
                    GROUNDWATER_FLUX(I,J)=HYDRAULIC_CONDUCTIVITY*GROUNDWATER_FLUX(I,J)/(CELL_AREAX*GROUNDWATER_DEPTH_SCALE)
                ENDIF
                IF (GROUNDWATER_FLUX(I,J) > QCMAX) QCMAX=GROUNDWATER_FLUX(I,J)
                IF (GROUNDWATER_FLUX(I,J) < QCMIN) QCMIN=GROUNDWATER_FLUX(I,J)
                FILTERED_GROUNDWATER_FLUX(I,J)=GROUNDWATER_FLUX(I,J)
            ENDDO L541
        ENDDO L540
        L531:      DO  J=1,MY
            L532:      DO  I=1,MX
                IF (FIXED_HEAD(I,J)) NFIXED=NFIXED+1
                NTOTAL=NTOTAL+1
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
                    If (JS > MY) JS=MY-1
                ENDIF
                QQ1=GROUNDWATER_FLUX(I,J)
                QQ2=GROUNDWATER_FLUX(IE,J)
                QQ3=GROUNDWATER_FLUX(IE,JN)
                QQ4=GROUNDWATER_FLUX(IE,JS)
                QQ5=GROUNDWATER_FLUX(I,JN)
                QQ6=GROUNDWATER_FLUX(I,JS)
                QQ7=GROUNDWATER_FLUX(IW,JN)
                QQ8=GROUNDWATER_FLUX(IW,J)
                QQ9=GROUNDWATER_FLUX(IW,JS)
                QTOT=0.0
                !      ********************************************************************
                !
                !      If divergence averaging is used, filtered_groundwater_flux is the average of the values
                !         of groundwater_flux for the local and surrounding 8 points.
                !         otherwise it is the maximum of the
                !         local and surrounding 8 points
                !
                !      ********************************************************************
                IF (SEEPAGE_AVERAGING) THEN
                    FILTERED_GROUNDWATER_FLUX(I,J)=0.0
                    IF (QQ1 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ1
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ2 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ2
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ3 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ3
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ4 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ4
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ5 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ5
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ6 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ6
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ7 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ7
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ8 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ8
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ9 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ9
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QTOT > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)/9.0
                    ENDIF
                ELSE
                    FILTERED_GROUNDWATER_FLUX(I,J)=AMAX1(QQ1,QQ2,QQ3,QQ4,QQ5,QQ6,QQ7,QQ8,QQ9)
                ENDIF
                IF (FILTERED_GROUNDWATER_FLUX(I,J) > 0.0) THEN
                    QCAVG=QCAVG+FILTERED_GROUNDWATER_FLUX(I,J)
                    NQQ=NQQ+1.0
                ENDIF
            ENDDO L532
        ENDDO L531
        DO J=1,MYY
            DO I=1,MX
                AVGHEFFECTIVE=AVGHEFFECTIVE+TRANSMISSIVITY_TERM(I,J)
            ENDDO
        ENDDO
        AVGHEFFECTIVE=AVGHEFFECTIVE/(MYY*MX*GROUNDWATER_DEPTH_SCALE)
        WRITE(*,922) AVGHEFFECTIVE
        WRITE(OUTHIST,922) AVGHEFFECTIVE
        922   FORMAT('AVERAGE EFFECTIVE AQUIFER DEPTH=',G12.5)
        FRACTFIXED=FLOAT(NFIXED)/NTOTAL
        IF (NQQ > 0.0) QCAVG=QCAVG/NQQ
        WRITE (OUTHIST,997) NUMIT,AMAX,EWMIN,EWMAX,QCMIN,QCAVG,QCMAX, &
        FRACTFIXED
        997 FORMAT('ITER=',I5,'AMAX=',G15.8,' HMIN=',G12.5, &
        ' HMAX=',G12.5,/,' CRITICAL_GROUNDWATER_FLUX=',G12.5,          &
        ' QAVG=',G12.5,' MAXIMUM_DISCHARGE=',G12.5,' % FIXED=',G12.5)
        RETURN
    END !SUBROUTINE EXPONENTIAL_HYDR_COND_GRNDWTR()
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE CONSTANT_HYDR_COND_GRNDWTR()
        !      ********************************************************************
        !
        !      This solves for the water table and calculates flow discharge, etc.
        !
        !      This is similar to subroutine solveexp except that the aquifer is assumed to have
        !       a finite depth but constant hydraulic conductivity
        !  MODIFIES: FIXED_HEAD, WATER_ELEVATION, TRANSMISSIVITY_TERM, GROUNDWATER_FLUX
        !            FILTERED_GROUNDWATER_FLUX
        !
        !      ********************************************************************
        !     sweep interior nodes with finite difference form of poisson's eqn.
        USE ERODE_GLOBALS
        USE GROUNDWATER_VARIABLES
        IMPLICIT NONE
        INTEGER :: I,J,NUMIT
        REAL(4) :: NFIXED,NTOTAL,LASTNFIXED
        REAL(4) :: AMAX,OLDVAL,ERR,QTOT
        REAL(4) :: EWMIN,EWMAX,VCMIN,VCMAX,QCMAX,QCAVG,NQQ,QCMIN
        REAL(4) :: QQ1,QQ2,QQ3,QQ4,QQ5,QQ6,QQ7,QQ8,QQ9,DX2
        REAL(4) :: CIP,CIM,CJP,CJM,ELOCAL,ATHICK
        REAL(4) :: QX,QY,MINEL,CELL_AREAX,FRACTFIXED
        REAL(4) ::  WATERNEW(IMMX,JMMX)
        REAL(4) :: MAXWATERCHANGE
        INTEGER :: IW,IE,JN,JS
         !!!!!!!!write(outhist,4833)
4833    format('groundwater called')       
        MAXWATERCHANGE=0.1*GROUNDWATER_DEPTH_SCALE
        DX2=CELL_SIZE*2.0
        CELL_AREAX=2.0*CELL_SIZE*CELL_SIZE
        VCMIN=1.0E+25
        VCMAX=-VCMIN
        EWMIN=VCMIN
        EWMAX=VCMAX
        MINEL=1.0E+25
        NUMIT=0
        LASTNFIXED=-100.0
        NFIXED=0.0
        NTOTAL=0.0
        L1997:      DO  J=1,MY
            L1998:      DO  I=1,MX
                IF (ELEVATION(I,J) < MINEL) THEN
                    MINEL=ELEVATION(I,J)
                    ILOWEST=I
                    JLOWEST=J
                ENDIF
            ENDDO L1998
        ENDDO L1997
        !      ********************************************************************
        !
        !      Initially only the point of lowest elevation (or the lower boundary
        !         if it is fixed) is fixed_head elevation, and the water table elevation
        !         (water_elevation) is set to that value.  The initial transmissivity factor
        !         is calculated for this flat water table
        !
        !      ********************************************************************
        L700:      DO  J=1,MY
            L701:      DO  I=1,MX
                ATHICK=ELEVATION(I,J)-INITIAL_ELEVATION(I,J)+GROUNDWATER_DEPTH_SCALE
                IF (ATHICK <= 0.0) THEN
                    FIXED_HEAD(I,J)=.TRUE.
                    WATER_ELEVATION(I,J)=ELEVATION(I,J)
                    TRANSMISSIVITY_TERM(I,J)=0.0
                ELSE
                    FIXED_HEAD(I,J)=.FALSE.
                    WATER_ELEVATION(I,J)=ELEVATION(I,J)-ATHICK*INITIAL_GROUNDWATER_DEPTH
                    TRANSMISSIVITY_TERM(I,J)=ATHICK*INITIAL_GROUNDWATER_DEPTH
                    IF (WATER_ELEVATION(I,J) < MINEL) THEN
                        WATER_ELEVATION(I,J)=MINEL
                        TRANSMISSIVITY_TERM(I,J)=MAX(0.00,MINEL-INITIAL_ELEVATION(I,J)+GROUNDWATER_DEPTH_SCALE)
                    ENDIF
                ENDIF
                IF ((I == ILOWEST).AND.(J == JLOWEST)) THEN
                    FIXED_HEAD(I,J)=.TRUE.
                ENDIF
                IF (WATER_ELEVATION(I,J) > EWMAX) EWMAX=WATER_ELEVATION(I,J)
                IF (WATER_ELEVATION(I,J) < EWMIN) EWMIN=WATER_ELEVATION(I,J)
                IF (FIXED_HEAD(I,J)) NFIXED=NFIXED+1
                NTOTAL=NTOTAL+1
            ENDDO L701
        ENDDO L700
        FIXED_HEAD(ILOWEST,JLOWEST)=.TRUE.
        FRACTFIXED=NFIXED/NTOTAL
        WRITE(OUTHIST,721) EWMIN,EWMAX,FRACTFIXED
        721   FORMAT(' EWMIN=',G12.5,' EWMAX=',G12.5,' % FIXED=',G12.5)
        NFIXED=0.0
        NTOTAL=0.0
        !      ********************************************************************
        !
        !      This is the main cycle for the iterative solution
        !
        !      ********************************************************************
        LNUMIT:      DO NUMIT=1,MAXIMUM_GROUNDWATER_ITERATIONS
            AMAX=0.0
            !      ********************************************************************
            !
            !      Use standard finite difference sweep with overrelaxation
            !
            !      ********************************************************************
            L040:      DO  J=1,MYY
                L041:      DO  I=1,MX
                    IF (FIXED_HEAD(I,J)) CYCLE
                    OLDVAL=WATER_ELEVATION(I,J)
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
                        If (JS > MY) JS=MY-1
                    ENDIF
                    ELOCAL=TRANSMISSIVITY_TERM(I,J)
                    CIP=(ELOCAL+TRANSMISSIVITY_TERM(IW,J))
                    CIM=(ELOCAL+TRANSMISSIVITY_TERM(IE,J))
                    CJP=(ELOCAL+TRANSMISSIVITY_TERM(I,JS))
                    CJM=(ELOCAL+TRANSMISSIVITY_TERM(I,JN))
                    WATERNEW(I,J)=(CIP*WATER_ELEVATION(IW,J)+CIM*WATER_ELEVATION(IE,J)+ &
                    CJP*WATER_ELEVATION(I,JS)    &
                    +CJM*WATER_ELEVATION(I,JN)+GROUNDWATER_RECHARGE_RATE/GROUNDWATER_SCALE_FACTOR)/(CIP+CIM+CJP+CJM)
                    ERR=ABS(WATERNEW(I,J)-OLDVAL)
                    IF (ERR > AMAX) AMAX=ERR
                ENDDO L041
            ENDDO L040
            IF (AMAX > MAXWATERCHANGE) THEN
                L049:      DO  J=1,MYY
                    L050:         DO  I=1,MX
                        WATER_ELEVATION(I,J)=WATER_ELEVATION(I,J) &
                        +MAXWATERCHANGE*(WATERNEW(I,J)-WATER_ELEVATION(I,J))/AMAX
                    ENDDO L050
                ENDDO L049
            ELSE
                WATER_ELEVATION=WATERNEW
            ENDIF
            !      ********************************************************************
            !
            !      If the predicted water table elevation is greater than the actual
            !         elevation, that point becomes a fixed_head head location
            !
            !      ********************************************************************
            NFIXED=0
            L043:      DO  J=1,MYY
                L044:      DO  I=1,MX
                    IF (WATER_ELEVATION(I,J) >= ELEVATION(I,J)) THEN
                        NFIXED=NFIXED+1
                        FIXED_HEAD(I,J)=.TRUE.
                        WATER_ELEVATION(I,J)=ELEVATION(I,J)
                    ENDIF
                    !      ********************************************************************
                    !
                    !      Recalculate transmissivity factor
                    !
                    !      ********************************************************************
                ENDDO L044
            ENDDO L043
            L731:          DO  J=1,MYY
                L732:          DO  I=1,MX
                    ATHICK=WATER_ELEVATION(I,J)-INITIAL_ELEVATION(I,J)+GROUNDWATER_DEPTH_SCALE
                    TRANSMISSIVITY_TERM(I,J)=MAX(0.00,ATHICK)
                    IF (ATHICK <= 0.00) FIXED_HEAD(I,J)=.TRUE.
                    IF (WATER_ELEVATION(I,J) < MINEL) THEN
                        WATER_ELEVATION(I,J)=MINEL
                        TRANSMISSIVITY_TERM(I,J)=MAX(0.00,MINEL-INITIAL_ELEVATION(I,J)+GROUNDWATER_DEPTH_SCALE)
                    ENDIF
                ENDDO L732
            ENDDO L731
            IF (AMAX > MAXIMUM_GROUNDWATER_ERROR) THEN
                IF (SHOW_GROUNDWATER_CALCULATIONS) THEN
                    WRITE(*,1731) NUMIT,NFIXED, AMAX
                    1731     FORMAT('WITER=',I5,' FIXED=',G12.5,' AMAX=',G12.5)
                ENDIF
                IF (MOD(NUMIT,100) == 0) THEN
                    WRITE(*,2345) NFIXED,LASTNFIXED
                    2345      FORMAT(' NF=',G12.5,' LNF=',G12.5)
                    IF (NFIXED == LASTNFIXED) THEN
                   !     GOTO 45
                    ENDIF
                    LASTNFIXED=NFIXED
                ENDIF
                NFIXED=0.0
                CYCLE
            ELSE
                GOTO 45
            ENDIF
        ENDDO LNUMIT
        45    CONTINUE
        !      ********************************************************************
        !
        !      Done.  Now calculate groundwater discharge and flow divergence
        !
        !      ********************************************************************
        EWMIN=1.0E+25
        EWMAX=-EWMIN
        QCMAX=-1.0E+25
        QCMIN=1.0E+25
        QCAVG=0.0
        NQQ=0.0
        L540:      DO  J=1,MYY
            L541:      DO  I=1,MX
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
                    If (JS > MY) JS=MY-1
                ENDIF
                IF (WATER_ELEVATION(I,J) > EWMAX) EWMAX=WATER_ELEVATION(I,J)
                IF (WATER_ELEVATION(I,J) < EWMIN) EWMIN=WATER_ELEVATION(I,J)
                QX=TRANSMISSIVITY_TERM(I,J)*(WATER_ELEVATION(IW,J)-WATER_ELEVATION(IE,J))/DX2
                QY=TRANSMISSIVITY_TERM(I,J)*(WATER_ELEVATION(I,JS)-WATER_ELEVATION(I,JN))/DX2
                !      ********************************************************************
                !
                !      If use_goundwater_flux is true the groundwater_flux is the groundwater flow rate, otherwise it
                !         is the flow divergence
                !
                !      ********************************************************************
                IF (USE_GROUNDWATER_FLOW) THEN
                    GROUNDWATER_FLUX(I,J)=HYDRAULIC_CONDUCTIVITY*SQRT(QX**2+QY**2)
                ELSE
                    GROUNDWATER_FLUX(I,J)=((WATER_ELEVATION(IW,J)-WATER_ELEVATION(I,J))  &
                    *(TRANSMISSIVITY_TERM(I,J)+TRANSMISSIVITY_TERM(IW,J))  &
                    -(WATER_ELEVATION(I,J)-WATER_ELEVATION(IE,J))  &
                    *(TRANSMISSIVITY_TERM(I,J)+TRANSMISSIVITY_TERM(IE,J))) &
                    +((WATER_ELEVATION(I,JS)-WATER_ELEVATION(I,J)) &
                    *(TRANSMISSIVITY_TERM(I,J)+TRANSMISSIVITY_TERM(I,JS))  &
                    -(WATER_ELEVATION(I,J)-WATER_ELEVATION(I,JN))  &
                    *(TRANSMISSIVITY_TERM(I,J)+TRANSMISSIVITY_TERM(I,JN)))
                    GROUNDWATER_FLUX(I,J)=HYDRAULIC_CONDUCTIVITY*GROUNDWATER_FLUX(I,J)/(CELL_AREAX)
                ENDIF
                IF (GROUNDWATER_FLUX(I,J) > QCMAX) QCMAX=GROUNDWATER_FLUX(I,J)
                IF (GROUNDWATER_FLUX(I,J) < QCMIN) QCMIN=GROUNDWATER_FLUX(I,J)
                FILTERED_GROUNDWATER_FLUX(I,J)=GROUNDWATER_FLUX(I,J)
            ENDDO L541
        ENDDO L540
        L531:      DO  J=1,MY
            L532:      DO  I=1,MX
                IF (FIXED_HEAD(I,J)) NFIXED=NFIXED+1
                NTOTAL=NTOTAL+1
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
                    If (JS > MY) JS=MY-1
                ENDIF
                QQ1=GROUNDWATER_FLUX(I,J)
                QQ2=GROUNDWATER_FLUX(IE,J)
                QQ3=GROUNDWATER_FLUX(IE,JN)
                QQ4=GROUNDWATER_FLUX(IE,JS)
                QQ5=GROUNDWATER_FLUX(I,JN)
                QQ6=GROUNDWATER_FLUX(I,JS)
                QQ7=GROUNDWATER_FLUX(IW,JN)
                QQ8=GROUNDWATER_FLUX(IW,J)
                QQ9=GROUNDWATER_FLUX(IW,JS)
                QTOT=0.0
                !      ********************************************************************
                !
                !      If divergence averaging is used, filtered_groundwater_flux is the average of the values
                !         of groundwater_flux for the local and surrounding 8 points.
                !         otherwise it is the maximum of the
                !         local and surrounding 8 points
                !
                !      ********************************************************************
                IF (SEEPAGE_AVERAGING) THEN
                    FILTERED_GROUNDWATER_FLUX(I,J)=0.0
                    IF (QQ1 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ1
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ2 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ2
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ3 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ3
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ4 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ4
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ5 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ5
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ6 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ6
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ7 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ7
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ8 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ8
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QQ9 > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)+QQ9
                        QTOT=QTOT+1.0
                    ENDIF
                    IF (QTOT > 0.0) THEN
                        FILTERED_GROUNDWATER_FLUX(I,J)=FILTERED_GROUNDWATER_FLUX(I,J)/9.0
                    ENDIF
                ELSE
                    FILTERED_GROUNDWATER_FLUX(I,J)=AMAX1(QQ1,QQ2,QQ3,QQ4,QQ5,QQ6,QQ7,QQ8,QQ9)
                ENDIF
                IF (FILTERED_GROUNDWATER_FLUX(I,J) > 0.0) THEN
                    QCAVG=QCAVG+FILTERED_GROUNDWATER_FLUX(I,J)
                    NQQ=NQQ+1.0
                ENDIF
            ENDDO L532
        ENDDO L531
        FRACTFIXED=NFIXED/NTOTAL
        IF (NQQ > 0.0) QCAVG=QCAVG/NQQ
        WRITE (*,997) NUMIT,AMAX,EWMIN,EWMAX,QCMIN,QCAVG,QCMAX,&
        FRACTFIXED,NFIXED,NTOTAL
        WRITE (OUTHIST,997) NUMIT,AMAX,EWMIN,EWMAX,QCMIN,QCAVG,QCMAX, &
        FRACTFIXED,NFIXED,NTOTAL
        997 FORMAT('ITER=',I5,'AMAX=',G15.8,' HMIN=',G12.5,  &
        ' HMAX=',G12.5,/,' CRITICAL_GROUNDWATER_FLUX=',G12.5, &
        ' QAVG=',G12.5,' MAXIMUM_DISCHARGE=',G12.5,' % FIXED=',G12.5 &
        ,/,' NFIXED=',G12.5,' NTOTAL=',G12.5)
        RETURN
    END ! SUBROUTINE CONSTANT_HYDR_COND_GRNDWTR()

