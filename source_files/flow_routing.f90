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
	SUBROUTINE READ_FLOW_PARAMETERS
	USE ERODE_GLOBALS
    USE AREA_GLOBALS
    USE LAKE_GLOBALS
    IMPLICIT NONE
    REAL(4) :: EEEMIN,TEMPX,LOW_RUNOFF_SCALE,HIGH_RUNOFF_SCALE  
	INTEGER :: FLOW_PARAMS, USEWET, IQCONSTANT, RANDDIRUSE, IMODEL_PELAGIC_DEPOSITION, DOINRIVER
    INTEGER :: IERROR, SEDIMENT_DEPENDENT_RUNOFF_USE, REGOLITH_DEPTH_RUNOFF_USE,SPATIAL_RUNOFF_USE
    INTEGER :: XMX,XMY,I,J,SLOPE_DEPENDENT_RUNOFF_USE
	CHARACTER (80) :: TEXTLINE
	    FLOW_PARAMS=93
		OPEN(FLOW_PARAMS,FILE='FLOW_PARAMETERS.PRM',ACTION='READ')
	    !        *******************************************************************
        !        ********************* flow routing and drainage area ***************
        !        *******************************************************************
        READ(FLOW_PARAMS,22) TEXTLINE
  22    FORMAT(A)
        WRITE(*,22) TEXTLINE
        READ(FLOW_PARAMS,*) USEWET
        READ(FLOW_PARAMS,*) IQCONSTANT
        IF (IQCONSTANT > 0) THEN
            RESCALE_DISCHARGES=.TRUE.
        ELSE
            RESCALE_DISCHARGES=.FALSE.
        ENDIF
        USE_DISCHARGE=.TRUE.
        IF (USEWET<0) USE_DISCHARGE=.FALSE.
        IF ((USEWET > 0).OR.MODEL_LAKE_EVAPORATION) THEN
            COMPLETE_RUNOFF=.TRUE.
        ELSE
            COMPLETE_RUNOFF=.FALSE.
        ENDIF
        READ(FLOW_PARAMS,*) RANDDIRUSE
        IF (RANDDIRUSE > 0) THEN
            USE_RANDOM_FLOWROUTING=.TRUE.
        ELSE
            USE_RANDOM_FLOWROUTING=.FALSE.
        ENDIF
        READ(FLOW_PARAMS,*) DISCHARGE_CONSTANT
        READ(FLOW_PARAMS,*) DISCHARGE_EXPONENT
        DISCHARGE_SCALE_FACTOR=DISCHARGE_CONSTANT**(1.0/DISCHARGE_EXPONENT)
        TEMPQCONSTANT=DISCHARGE_CONSTANT
        LASTQCONSTANT=DISCHARGE_CONSTANT
        READ(FLOW_PARAMS,*) AREAFACTOR
        READ(FLOW_PARAMS,*) RAINDEPTH
        READ(FLOW_PARAMS,*) RAINSTD
        READ(FLOW_PARAMS,*) USE_DIVERGENCE_DEPENDENT_RUNOFF
        READ(FLOW_PARAMS,*) DIVERGENCE_FOR_MEAN_RUNOFF
        READ(FLOW_PARAMS,*) HIGH_DIVERGENCE_RUNOFF
        READ(FLOW_PARAMS,*) LOW_DIVERGENCE_RUNOFF
        READ(FLOW_PARAMS,*) DIVERGENCE_SCALE_PARAMETER
        READ(FLOW_PARAMS,*) INTEGER_DIVERGENCE_RANGE ! Integer >=3
        IF (INTEGER_DIVERGENCE_RANGE<3) INTEGER_DIVERGENCE_RANGE=3
        if (integer_divergence_range>3) then
            write(*,11191) integer_divergence_range
11191       format('irange=',i5)
           if ((integer_divergence_range.eq.5).or.(integer_divergence_range.eq.7) &
               .or.(integer_divergence_range.eq.9)) then
           else
               write(outhist,6291)
               write(*,6291)
6291           format('integer divergence range must be 3, 5, 7, or 9, program stopping')
               stop
          endif
        endif
        READ(FLOW_PARAMS,*) SPATIAL_RUNOFF_USE
        READ(FLOW_PARAMS,*) LOW_RUNOFF_SCALE
        READ(FLOW_PARAMS,*) HIGH_RUNOFF_SCALE
        IF (SPATIAL_RUNOFF_USE>0) THEN
            USE_SPATIAL_RUNOFF=.TRUE.
            RUNOFF_AVAL=(HIGH_RUNOFF_SCALE-LOW_RUNOFF_SCALE)/(SPATIAL_MAX-SPATIAL_MIN)
            RUNOFF_BVAL=LOW_RUNOFF_SCALE-RUNOFF_AVAL*SPATIAL_MIN
        ELSE
            USE_SPATIAL_RUNOFF=.FALSE.
        ENDIF
        READ(FLOW_PARAMS,*) SEDIMENT_DEPENDENT_RUNOFF_USE
        READ(FLOW_PARAMS,*) SEDIMENT_RUNOFF_FACTOR 
        READ(FLOW_PARAMS,*) REGOLITH_DEPTH_RUNOFF_USE
        READ(FLOW_PARAMS,*) RUNOFF_SCALE_BEDROCK
        READ(FLOW_PARAMS,*) RUNOFF_SCALE_DEEP_REGOLITH
        READ(FLOW_PARAMS,*) RUNOFF_DEPTH_DECAY_RATE
        READ(FLOW_PARAMS,*) SLOPE_DEPENDENT_RUNOFF_USE
        READ(FLOW_PARAMS,*) MEDIAN_SLOPE
        MEDIAN_SLOPE=ALOG(MEDIAN_SLOPE)
        READ(FLOW_PARAMS,*) SLOPE_RUNOFF_SCALE_FACTOR
        READ(FLOW_PARAMS,*) NCALCEVAP
        READ(FLOW_PARAMS,*) NCHANGEEVAP
        READ(FLOW_PARAMS,*) X_RATIO
        READ(FLOW_PARAMS,*) EVAPORATION_STANDARD_DEVIATION
        READ(FLOW_PARAMS,*) IMODEL_PELAGIC_DEPOSITION
        READ(FLOW_PARAMS,*) WASHLOAD_FRACTION
        READ(FLOW_PARAMS,*) PELAGICCREEP
        READ(FLOW_PARAMS,*) MINIMUM_PELAGIC_LAKE_SIZE
        RUNOFF_DEPTH_DECAY_RATE=-LOG(0.5)/RUNOFF_DEPTH_DECAY_RATE
        EVAPORATION_MEAN=X_RATIO+1.0
        LOG_EVAP_MEAN=LOG(EVAPORATION_MEAN**2/SQRT(EVAPORATION_STANDARD_DEVIATION**2 &
            + EVAPORATION_MEAN**2))
        LOG_EVAP_SD=SQRT(LOG(EVAPORATION_STANDARD_DEVIATION**2/ &
            EVAPORATION_MEAN**2+1.0))
        IF (IMODEL_PELAGIC_DEPOSITION > 0) THEN
            MODEL_PELAGIC_DEPOSITION=.TRUE.
        ELSE
            MODEL_PELAGIC_DEPOSITION=.FALSE.
        ENDIF
        IF (SEDIMENT_DEPENDENT_RUNOFF_USE>0) THEN
            USE_SEDIMENT_DEPENDENT_RUNOFF=.TRUE.
        ELSE
            USE_SEDIMENT_DEPENDENT_RUNOFF=.FALSE.
        ENDIF
        IF (REGOLITH_DEPTH_RUNOFF_USE>0) THEN
            USE_REGOLITH_DEPTH_DEPENDENT_RUNOFF=.TRUE.
        ELSE
            USE_REGOLITH_DEPTH_DEPENDENT_RUNOFF=.FALSE.
        ENDIF
        IF (SLOPE_DEPENDENT_RUNOFF_USE>0) THEN
            USE_SLOPE_DEPENDENT_RUNOFF=.TRUE.
        ELSE
            USE_SLOPE_DEPENDENT_RUNOFF=.FALSE.
        ENDIF
        IF(FLUVIAL_AND_SLOPE_MODELING.AND.MODEL_PELAGIC_DEPOSITION) THEN
            ALLOCATE(PELAGICCHANGE(MX,MY),PELAGICEDIFF(MX,MY),STAT=IERROR)
			IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(PELAGIC_SEDIMENT_VOLUME(LMMX),STAT=IERROR)
			IF (IERROR > 0) CALL ALLOCERROR()
        ENDIF
        !       *****************************************************************
        !        Read the location of an incoming river, its discharge and
        !           sediment load
        !       *****************************************************************
        READ(FLOW_PARAMS,*) DOINRIVER
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            IS_INFLUENT_RIVER_LOCATION=.FALSE.
        ENDIF
        IF (DOINRIVER > 0) THEN
            HAVE_INFLUENT_RIVERS=.TRUE.
            OPEN(76,FILE='INRIVER.PRM')
            READ(76,*) NUMBER_OF_INFLUENT_RIVERS
            DO I=1,NUMBER_OF_INFLUENT_RIVERS
                READ(76,*) I_INFLUENT_RIVER(I),J_INFLUENT_RIVER(I),INFLUENT_RIVER_DISCHARGE(I),INFLUENT_RIVER_SEDIMENT_LOAD(I)
                IS_INFLUENT_RIVER_LOCATION(I_INFLUENT_RIVER(I),J_INFLUENT_RIVER(I))=.TRUE.
            ENDDO
            CLOSE(76)
        ELSE
            HAVE_INFLUENT_RIVERS=.FALSE.
        ENDIF
		CLOSE(FLOW_PARAMS)
        !       *****************************************************************
        !        Calculate parameters for determining specific runoff as a function
        !           of local divergence
        !       If USE_DIVERGENCE_DEPENDENT_RUNOFF equals zero specific runoff is not dependent on local divergence
        !
        !       HIGH_DIVERGENCE_RUNOFF and LOW_DIVERGENCE_RUNOFF should be scaled so that the mean runoff should be near unity.
        !       Note that these values are the highest and lowest absolute runoff values
        !       The value of topographic divergence with approximately the mean runoff should equal the desired value of the product
        !       of the actual topographic divergence at the mean value times the DIVERGECE_SCALE_PARAMETER, the latter which influences the range
        !       of divergence values over which the runoff goes from near the lowest value to the highest.
        !       If DIVERGENCE_SCALE_FACTOR <0.0 topographically concave regions (valleys and depressionsdivergence<0) 
        !          have highest runoff else
        !       if DIVERGENCE_SCALE_FACTOR >0.0 topographically convex regions (e.g. hills, divergence<0) have highest
        !          runoff values
        !       *****************************************************************
        IF (USE_DIVERGENCE_DEPENDENT_RUNOFF /= 0) THEN
            DIVERGENCE_DEPENDENT_RUNOFF=.TRUE.
        ELSE          
            DIVERGENCE_DEPENDENT_RUNOFF=.FALSE.
        ENDIF
       !        *******************************************************************
        !        ********************* flow routing and drainage area ***************
        !        *******************************************************************
        IF (.NOT.MODEL_LAKE_EVAPORATION) THEN
            IF (COMPLETE_RUNOFF) THEN
                WRITE(OUTHIST,547)
                547 FORMAT(' DEPRESSIONS ARE COMPLETELY FILLED AND SPILL' &
                    ' OR ARE MODELED USING LAKE EVAPORATION')
            ELSE
                WRITE(OUTHIST,548)
                548 FORMAT(' DEPRESSIONS ARE DRY')
            ENDIF
        ELSE
            WRITE(OUTHIST,2560)
2560         FORMAT ('MODELING EVAPORATION WITHIN LAKES')
            WRITE(OUTHIST,2549) X_RATIO,EVAPORATION_MEAN,NCHANGEEVAP, &
                EVAPORATION_STANDARD_DEVIATION
                2549 FORMAT('X_RATIO=',G12.5,' ,EQUIVALENT EVAPORATION MEAN=',G12.5,/, &
                  'EVAPORATION_RECALCULATED EVERY ',I5,' ITERATIONS',/, &
                  'STANDARD DEVIATION OF EVAPORATION=',G12.5)    
        ENDIF
        IF (USE_RANDOM_FLOWROUTING) THEN
            WRITE(OUTHIST,495)
            495   FORMAT(' RANDOM DOWNSTREAM FLOW ROUTING')
        ELSE
            WRITE(OUTHIST,496)
            496   FORMAT(' DETERMINISTIC D8 DOWNSTREAM FLOW ROUTING')
        ENDIF
        IF (RESCALE_DISCHARGES) THEN
            WRITE(OUTHIST,1496) DISCHARGE_CONSTANT,DISCHARGE_EXPONENT
1496        FORMAT(' DISCHARGES ARE RESCALED WITH RESPECT TO EFFECTIVE CONTRIBUTING AREA',/, &
                        'DISCHARGE CONSTANT=',G12.5,' DISCHARGE EXPONENT=',G12.5)
        ENDIF
        IF (AREAFACTOR.NE.0.0) THEN
            WRITE(OUTHIST,494) AREAFACTOR, RAINDEPTH,RAINSTD
494         FORMAT(' AREA CORRECTION FACTOR=',G12.5, 'RAINFALL MEAN DEPTH=',G12.5, &
              ' STANDARD DEVIATION OF RAINFALL=',G12.5)
        ENDIF
        IF (DIVERGENCE_DEPENDENT_RUNOFF) THEN
            WRITE(OUTHIST,783) DIVERGENCE_FOR_MEAN_RUNOFF,HIGH_DIVERGENCE_RUNOFF,LOW_DIVERGENCE_RUNOFF,  &
            RUNOFF_CONSTANT_1,RUNOFF_CONSTANT_2,RUNOFF_CONSTANT_3,DIVERGENCE_SCALE_PARAMETER
783         FORMAT(' PROGRAM USES VARIABLE WATER YIELD BASED ON SLOPE' &
            ,' DIVERGENCE',/,' DIVERGENCE FOR MEAN WATER YIELD=',G12.5, &
            ' DIVERGENCE FOR 90% YIELD=',G12.5,' LOW WATER YIELD='&
            ,G12.5,/, &
            ' DERIVATIVE WATER YIELD CONSTANTS 1-3=',3(G12.5,' '),/, &
               'DIVERGENCE SCALE PARAMETER=',G12.5)
        ENDIF
        IF (USE_SEDIMENT_DEPENDENT_RUNOFF) THEN
            WRITE(OUTHIST,1783) SEDIMENT_RUNOFF_FACTOR
1783        FORMAT(' RUNOFF ON SEDIMENT COVERED AREAS IS SCALED TO BE ',G12.5, &
                 ' TIMES THAT ON NON-SEDIMENT COVERED AREAS')
        ENDIF
        IF (USE_REGOLITH_DEPTH_DEPENDENT_RUNOFF) THEN
            WRITE(OUTHIST,2783) RUNOFF_SCALE_BEDROCK,RUNOFF_SCALE_DEEP_REGOLITH,RUNOFF_DEPTH_DECAY_RATE
2783        FORMAT('RUNOFF DEPENDS UPON REGOLITH DEPTH',/,'RUNOFF FROM BEDROCK SCALING=',G12.5,/, &
                 'RUNOFF FROM DEEP REGOLITH SCALING=',G12.5,' RUNOFF DEPTH SCALING FACTOR=',G12.5)
        ENDIF
        IF (USE_SLOPE_DEPENDENT_RUNOFF) THEN
            WRITE(OUTHIST,3783) EXP(MEDIAN_SLOPE),SLOPE_RUNOFF_SCALE_FACTOR
3783        FORMAT('RUNOFF DEPENDS ON LOCAL GRADIENT',/,'GRADIENT FOR MEDIAN RUNOFF VALUE=',G12.5, &
                 ' SCALING FACTOR FOR GRADIENT DEPENDENCY=',G12.5)
        ENDIF
        IF (DO_SPATIAL_VARIATION.AND.USE_SPATIAL_RUNOFF) THEN
            WRITE(OUTHIST,4783) LOW_RUNOFF_SCALE,HIGH_RUNOFF_SCALE
4783        FORMAT('RUNOFF VARIES SPATIALLY AS READ IN USING THE FILE "RELATIVE_DISCHARGES.DAT"',/, &
                 'SCALING FOR LOWEST DISCHARGE=',G12.5,' SCALING FOR HIGHEST DISCHARGES=',G12.5)
        ENDIF
        IF (MODEL_PELAGIC_DEPOSITION) THEN
           WRITE(OUTHIST,1504) WASHLOAD_FRACTION,PELAGICCREEP
           1504 FORMAT(' PELAGIC DEPOSITION IS MODELED',/, &
           '     FRACTION OF SEDIMENT THAT IS DELIVERED OFFSHORE VERSUS DELTAIC=',G12.5, &
           '     COEFFICIENT FOR PELAGIC SEDIMENT DIFFUSIVITY=',G12.5)
        ENDIF
        IF (HAVE_INFLUENT_RIVERS) THEN
           WRITE(OUTHIST,1505)
           1505 FORMAT(' MODELING INFLOW OF WATER AND SEDIMENT FROM ONE OR MORE RIVERS AT OUTER BOUNDARIES')
        ENDIF		
		RETURN
	END !SUBROUTINE READ_FLOW_PARAMETERS
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL(4) FUNCTION DISCHARGE_FROM_CELL(I,J)
        USE ERODE_GLOBALS
        USE GROUNDWATER_VARIABLES
        USE AREA_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,J
        INTEGER :: IW,IE,JS,JN
        REAL(4) :: TEMP1,TEMP2,TEMP3,TEMP4,RECHARGERATIO,LOCALDIVERGENCE
        LOGICAL SKIPDIVERGENCE

        IF ((MOD(ITERATION,DIVERGEINTERVAL) == 0).OR.(ITERATION<2)) THEN
        RECHARGERATIO=GROUNDWATER_RECHARGE_RATE*1.0
        !     **********************************************************************
        !      Determines the local specific runoff yield from unsubmerged land surfaces.  
        !         Can depend upon the local slope
        !         divergence, occurrence of groundwater efflux,
        !         and/or whether there is a bistable erosion rate.  This code
        !         can be modified to make runoff a function of a variety of local conditions, such as 
        !         elevation effects on precipitation amount, regolith thickness, exposure
        !         of bedrock, and slope steepness (perhaps greater on steep slopes and lesser
        !         from very flat surfaces), or sediment-versus-regolith cover.  Note: The 
        !         difference in runoff from dry land versus lake surfaces (SUBMERGED=.TRUE.) is 
        !         handled within the main code.
        !     **********************************************************************
        TEMP1=1.0

        IF (DIVERGENCE_DEPENDENT_RUNOFF) THEN
            if ((integer_divergence_range==5).or.(integer_divergence_range==7).or. &
                (integer_divergence_range==9)) then
                TEMP2 = DIVERGENCE_SCALE_PARAMETER * (divergence(i,j)-DIVERGENCE_FOR_MEAN_RUNOFF)/SQRT(2.0)
                TEMP1 = 0.5*(1.0+TEMP2/SQRT(1.0+TEMP2*2))
                TEMP1 = LOW_DIVERGENCE_RUNOFF+( HIGH_DIVERGENCE_RUNOFF-LOW_DIVERGENCE_RUNOFF)*TEMP1
            ELSE 
            SKIPDIVERGENCE=.FALSE.
                IW=I-INTEGER_DIVERGENCE_RANGE
                IE=I+INTEGER_DIVERGENCE_RANGE
                IF (IS_X_PERIODIC) THEN
                    IF (IW < 1) IW=MX+I-INTEGER_DIVERGENCE_RANGE
                    IF (IE > MX) IE=I-MX+INTEGER_DIVERGENCE_RANGE
                ELSE
                    IF (IW <1) SKIPDIVERGENCE=.TRUE.
                    IF (IE > MX) SKIPDIVERGENCE=.TRUE.
                ENDIF
                JN=J-1
                JS=J+1
                IF (IS_Y_PERIODIC) THEN
                    IF (JN < 1) JN=MY+I-INTEGER_DIVERGENCE_RANGE
                    IF (JS > MY) JS=I-MY+INTEGER_DIVERGENCE_RANGE
                ELSE 
                    IF (JN<1) SKIPDIVERGENCE=.TRUE.
                    IF (JS > my) SKIPDIVERGENCE=.TRUE.
                ENDIF
            IF (.NOT.SKIPDIVERGENCE) THEN
                LOCALDIVERGENCE=(ELEVATION(IW,JN)+ELEVATION(IE,JN)+ELEVATION(IW,JS)+ELEVATION(IE,JS)+ &
                    4.0*(ELEVATION(IE,J)+ELEVATION(IW,J)+ELEVATION(I,JN)+ELEVATION(I,JS))-20.0*ELEVATION(I,J))/ &
                    (6.0*CELL_AREA)
                TEMP2 = DIVERGENCE_SCALE_PARAMETER * (LOCALDIVERGENCE-DIVERGENCE_FOR_MEAN_RUNOFF)/SQRT(2.0)
                TEMP1 = 0.5*(1.0+TEMP2/SQRT(1.0+TEMP2*2))
                TEMP1 = LOW_DIVERGENCE_RUNOFF+( HIGH_DIVERGENCE_RUNOFF-LOW_DIVERGENCE_RUNOFF)*TEMP1
            ELSE 
                TEMP1= (HIGH_DIVERGENCE_RUNOFF+LOW_DIVERGENCE_RUNOFF)/2.0
            ENDIF
          endif
        ENDIF
        IF (USE_SEDIMENT_DEPENDENT_RUNOFF) THEN
            IF (IS_SEDIMENT_COVERED(I,J)) THEN
                TEMP1=TEMP1*SEDIMENT_RUNOFF_FACTOR
            ENDIF
        ENDIF
        IF (USE_SLOPE_DEPENDENT_RUNOFF) THEN
            IF (D8_GRADIENT(I,J)>0.0) THEN
               TEMP1=2.0/(1.0+EXP(-SLOPE_RUNOFF_SCALE_FACTOR*(ALOG(D8_GRADIENT(I,J))-MEDIAN_SLOPE)))
            ELSE
                TEMP1=0.0
            ENDIF
        ENDIF                
        IF(USE_REGOLITH_DEPTH_DEPENDENT_RUNOFF) THEN
            IF (REGOLITH(I,J)<=0.0) THEN
                TEMP1=TEMP1*RUNOFF_SCALE_BEDROCK
            ELSE
                TEMP1=TEMP1*RUNOFF_SCALE_DEEP_REGOLITH-(RUNOFF_SCALE_DEEP_REGOLITH-RUNOFF_SCALE_BEDROCK)  &
                    *EXP(-RUNOFF_DEPTH_DECAY_RATE*REGOLITH(I,J))
            ENDIF
        ENDIF
        IF (DO_SPATIAL_VARIATION.AND.USE_SPATIAL_RUNOFF) THEN
            TEMP1=TEMP1*(SPATIAL_VARIATION(I,J)*RUNOFF_AVAL+RUNOFF_BVAL)**(1.0/DISCHARGE_EXPONENT)
        ELSE
            TEMP1=DISCHARGE_SCALE_FACTOR*TEMP1
        ENDIF
        IF (ACCELERATED_EROSION(I,J)) THEN
            DISCHARGE_FROM_CELL=BISTABLE_RUNOFF_FACTOR*TEMP1
        ELSE
            DISCHARGE_FROM_CELL=TEMP1
        ENDIF
            
        100   FORMAT(' DISCHARGE_FROM_CELL=',G12.5)
        !     **********************************************************************
        !     *****   New groundwater code - add in groundwater flow component
        !             If groundwater_flow_fraction=1 all surface flow from groundwater plus
        !             recharge (rainfall) on saturated ground (fixed_head true)
        !     **********************************************************************
        IF (MODEL_GROUNDWATER) THEN
            TEMP2=AMAX1(0.0,GROUNDWATER_FLUX(I,J))
            IF (FIXED_HEAD(I,J)) THEN
                TEMP2=TEMP2+RECHARGERATIO
            ELSE
                TEMP2=0.0
            ENDIF
            TEMP3=DISCHARGE_FROM_CELL*(1.0-GROUNDWATER_FLOW_FRACTION)
            TEMP4=TEMP2*GROUNDWATER_FLOW_FRACTION
            DISCHARGE_FROM_CELL=TEMP3+TEMP4
            AVGRUNOFF=AVGRUNOFF+TEMP3
            AVGSEEP=AVGSEEP+TEMP4
            IF (TEMP2 > 0.0) THEN
                AVGSEEPFRACT=AVGSEEPFRACT+TEMP4/DISCHARGE_FROM_CELL
                NFRACTAVG=NFRACTAVG+1.0
            ENDIF
            NSEEPAVG=NSEEPAVG+1.0
        ENDIF
        DISCHARGE_FROM_CELL=DISCHARGE_FROM_CELL*CELL_AREA
        LOCAL_DISCHARGE(I,J)=DISCHARGE_FROM_CELL
        ELSE 
            DISCHARGE_FROM_CELL=LOCAL_DISCHARGE(I,J)
        ENDIF     
        RETURN
    END !  REAL(4) FUNCTION DISCHARGE_FROM_CELL
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE IS_IT_SUBMERGED()
        !      ************************************************************
        !        Determines if location i,j is submerged (submerged(i,j)=.true.)
        !           or not (submerged(i,j)=.false.)
        !   MODIFIES: SUBMERGED
        !      ***********************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,NNN
        SUBMERGED=.FALSE.
        DO  J=1,MYY
            DO  I=1,MX
                IF (USE_AN_OCEAN) THEN
                    IF (ELEVATION(I,J) <= OCEAN_ELEVATION) SUBMERGED(I,J)=.TRUE.
                ENDIF
                IF(.NOT.COMPLETE_RUNOFF) CYCLE
                NNN=BASIN_NUMBER(I,J)
                IF (MODEL_LAKE_EVAPORATION) THEN
                    IF ((COMPLETE_RUNOFF).AND.(ENCLOSED(NNN)).AND.  &
                    (ELEVATION(I,J) < LAKE_SURFACE_ELEVATION(NNN))) THEN
                        SUBMERGED(I,J)=.TRUE.
                    ENDIF
                ELSE
                    IF ((COMPLETE_RUNOFF).AND.(ENCLOSED(NNN)).AND.  &
                    (ELEVATION(I,J) < LAKE_OUTLET_ELEVATION(NNN))) THEN
                        SUBMERGED(I,J)=.TRUE.
                    ENDIF
                ENDIF
                IF (DO_FLOW_BOUNDARIES) THEN
                ELSE
                    IF (COMPLETE_RUNOFF.AND.(.NOT.IS_Y_PERIODIC)) THEN
                        IF (ELEVATION(I,J) <= ELEVATION(I,MY)) SUBMERGED(I,J)=.TRUE.
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        RETURN
    END ! SUBROUTINE IS_IT_SUBMERGED
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DRAINAGE_BASIN_AREA_FLOW()
        USE ERODE_GLOBALS, CHANNEL_EROSION_RATE => CFNE
        USE AREA_GLOBALS
        USE LAKE_GLOBALS
        USE GROUNDWATER_VARIABLES
        IMPLICIT NONE
        !      *******************************************************************
        !       Calculates drainage areas for each location via flow routing,
        !       including overflow from enclosed depressions
        !
        !       Also determines total sediment yield from channel erosion
        !
        !       Local vectors are defined in area_globals and lake_globals:
        !           (i=integer, l=logical, else real)
        !
        !           tobeadded(k)   - the cumulative drainage area of  enclosed drainage
        !                            basins upstream
        !       l   needtodo(k)    - working variable which indicates which enclosed
        !                            basins need to be considered for grouping with
        !                            neighboring basins
        !       l   outer(k) - true if the exit point from an enclosed basin is
        !                      lower than any exterior contiguous point outside the basin
        !           local_basin_discharge(k) - the total discharge from the enclosed basin
        !           qtoadd(k) - the cumulative contributed discharge from enclosed
        !                       drainage basins upstream
        !   MODIFIES: I_OUTFLOW, J_OUTFLOW, ENCLOSED, SEDIMENT_YIELD, FLOW_DIRECTION
        !             DISCHARGE, DRAINAGE_AREA, BASIN_NUMBER,LAKE_AREA,IS_EXIT_BASIN
        !             ROUTED_DRAINAGE_AREA, DOWNSTREAM_BASIN
        !             LAKE_OUTLET_ELEVATION, LAKE_SURFACE_ELEVATION,LAKE_AREA
        !             ERODETOADD,ERGRAVTOADD 
        !   CALLS:  DRAWLINE, DISCHARGE_FROM_CELL, BASIN_REPORT 
        !
        !      *******************************************************************
        INTEGER :: I,J,LOCAL_DIRECTION,I_OLD,J_OLD,I_NEW,J_NEW,K,L,KK,KKFROM
        INTEGER :: IFROM,JFROM,NUMBBASIN,NNN,KKK,III,JJJ
        INTEGER :: IOUTLET,JOUTLET,NCONTIG
        INTEGER :: II,JJ,ILL,JLL,IUU,JUU,ILOW,IHIGH,JLOW,JHIGH,KAROUND
        INTEGER :: KKKK,KEND,KKCOUNT,KL
        INTEGER :: KKSQUARE
        INTEGER :: NBOTTOMBASIN,NBOTTOMADD,CURRENT_BASIN 
        REAL(4) :: ADDAREA
        REAL(4) :: AAMAX,AAMIN,AAAVG,AASD
        LOGICAL :: CONTIGUOUS,ISFOUND,TOSEDCOVER,FLOW_UPPER_BOUNDARY
        LOGICAL :: TOLAKE,NEWBASIN,REDO,XCONTINUE,YCONTINUE
        REAL(4) :: XCOMP,SEDADD,RIVERADD
        REAL(4) :: WSUM,WNUM,DISCHARGE_FROM_CELL,QADD,ADDQ
        REAL(4) :: TBEDCONT,TSEDCONT,ELEVLOWEST
        REAL(4) :: THEAVG,THEMAX,THEMIN
        REAL(4) :: GLOBAL_MAXIMUM_DISCHARGE
        INTEGER :: KIO,KII,KIJ
        real(4) :: ttt1,ttt2,ttt3,ttt4
        integer :: kkk1,kkk2,kkk3,kkk4,kkk5,kkk6,kkk7,kkk8,kkk9
        logical :: testit
        ITERATION_MAXIMUM_SEDIMENT_YIELD=0.0
        AAMAX=-1.0E+25
        ELEVLOWEST=1.0E+25
        AAMIN=-AAMAX
        AAAVG=0.0
        AASD=0.0
        NBOTTOMBASIN=0
        NBOTTOMADD=0
        GLOBAL_MAXIMUM_DISCHARGE=-1.0E+25
        AVGRUNOFF=0.0
        AVGSEEP=0.0
        AVGSEEPFRACT=0.0
        NSEEPAVG=0.0
        
        !      *******************************************************************
        !       no_lakes_exist is index of presence of any enclosed basins
        !      *******************************************************************
        NO_LAKES_EXIST=.TRUE.
        ABORT_SIMULATION=.FALSE.
        !      *******************************************************************
        !       set initial values
        !            current_basin is the current basin number
        !      *******************************************************************
        CURRENT_BASIN=0
        TOBEADDED=0
        DOWNSTREAM_BASIN=0
        I_OUTFLOW=0
        J_OUTFLOW=0
        BASIN_DRAINAGE_AREA=0
        LOCAL_BASIN_DISCHARGE=0.0
        QTOADD=0.0
        ENCLOSED =.FALSE.
        IS_EXIT_BASIN=.FALSE.
        TBEDCONT=0.0
        TSEDCONT=0.0
        FLOW_UPPER_BOUNDARY=.FALSE.
        !      *********************************************************************
        !        Initialize matrices
        !         "sediment_yield" will hold the total sediment supplied by bed erosion of bedrock
        !         channels upstream from the given location, except if a depression
        !         intervenes in the flow path, which zeros the sediment yield.  Sediment_yield is
        !         defifined for both alluvial and bedrock areas
        !         "sediment_yield" is used for downstream sediment routing modeling of
        !         erosion/deposition
        !      *********************************************************************
        BASIN_NUMBER=0
        SEDIMENT_YIELD=0.0
        DISCHARGE=0.0
        DRAINAGE_AREA=0.0
        IF (IS_Y_PERIODIC.OR.REFLECTIVE_UPPER_BOUNDARY) THEN
            DO  J=1,MY
                DO  I=1,MX
                    IF (ELEVATION(I,J) < ELEVLOWEST) THEN
                        ELEVLOWEST=ELEVATION(I,J)
                        ILOWEST=I
                        JLOWEST=J
                    ENDIF
                ENDDO
            ENDDO
        ELSE
            ILOWEST=0
            JLOWEST=0
        ENDIF
        IF (IS_Y_PERIODIC) THEN
        ELSE
            DO  I=1,MX
                FLOW_DIRECTION(I,MY)=1
!                IF (REFLECTIVE_UPPER_BOUNDARY) FLOW_DIRECTION(I,1) = 1
            ENDDO
        ENDIF
        !      *******************************************************************
        !       Main loop for determination of drainage areas
        !      *******************************************************************
        L110: DO  J=1,MYY
            1105  FORMAT('J=',I5)
            M110: DO  I=1,MX
                1106  FORMAT('I=',I5)
                !      *******************************************************************
                !       tolake is index of whether basin is enclosed
                !       newbasin is index of whether a new basin has been detected
                !         (either a new exit or a new enclosed depression)
                !       riveradd is contribution (if any) from rivers entering from outside
                !         the matrix
                !      *******************************************************************
                TOLAKE=.FALSE.
                NEWBASIN=.FALSE.
                RIVERADD=0.0
                I_OLD=I
                J_OLD=J
                !      *******************************************************************
                !       Increment drainage area for locally contributed drainage and for
                !         incoming streams
                !      *******************************************************************
                IF (HAVE_INFLUENT_RIVERS.AND.IS_INFLUENT_RIVER_LOCATION(I,J)) THEN
                    DO KL=1,NUMBER_OF_INFLUENT_RIVERS
                        IF((I == I_INFLUENT_RIVER(KL)).AND.(J == J_INFLUENT_RIVER(KL)))THEN
                            RIVERADD=INFLUENT_RIVER_DISCHARGE(KL)
                        ENDIF
                    ENDDO
                ENDIF
                DRAINAGE_AREA(I,J)=DRAINAGE_AREA(I,J)+CELL_AREA+RIVERADD
                QADD=DISCHARGE_FROM_CELL(I,J)+RIVERADD
                DISCHARGE(I,J)=DISCHARGE(I,J)+QADD
                L120: DO
                    !      *******************************************************************
                    !       Determine local drainage direction and check for bad value
                    !       Implicit loop for traversing downstream
                    !      *******************************************************************
                    LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                    IF (DO_FLOW_BOUNDARIES) THEN
                        IF (LOCAL_DIRECTION == 1) THEN
                            IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) THEN
                                NEWBASIN=.TRUE.
                                TOLAKE=.FALSE.
                                NBOTTOMBASIN=NBOTTOMBASIN+1
                            ENDIF
                            EXIT L120
                        ENDIF
                    ENDIF
                    !      *******************************************************************
                    !       Check for enclosed basin [point has negative flow_direction(i,j)]
                    !      *******************************************************************
                    IF (LOCAL_DIRECTION < -1) THEN
                        IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) THEN
                        ELSE
                            TOLAKE=.TRUE.
                            NO_LAKES_EXIST=.FALSE.
                        ENDIF
                        !      *******************************************************************
                        !       If it is a depression then check if this is a new basin
                        !      *******************************************************************
                        IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) THEN
                            NEWBASIN=.TRUE.
                        ENDIF
                        !      *******************************************************************
                        !       Branch if this is an enclosed basin
                        !      *******************************************************************
                        EXIT L120
                    ENDIF
                    IF (REFLECTIVE_UPPER_BOUNDARY.AND.(J==1)) THEN
                        IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) THEN
                        ELSE
                            TOLAKE=.TRUE.
                            NO_LAKES_EXIST=.FALSE.
                        ENDIF
                        !      *******************************************************************
                        !       If it is a depression then check if this is a new basin
                        !      *******************************************************************
                        IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) THEN
                            NEWBASIN=.TRUE.
                        ENDIF
                        !      *******************************************************************
                        !       Branch if this is an enclosed basin
                        !      *******************************************************************
                        EXIT L120
                    ENDIF                    
                    !      *******************************************************************
                    !       If its not a depression then determine the next downstream location
                    !       Also change reference location to the next downstream point and
                    !          increment its drainage area and sediment yield
                    !       We increment in all areas except after depressions for sediment_yield but only in
                    !          bedrock channels for bederode
                    !      *******************************************************************
                    IF (LOCAL_DIRECTION > 0) THEN
                        J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                        I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                        IF (IS_X_PERIODIC) THEN
                            IF (I_NEW < 1) I_NEW = MX
                            IF (I_NEW > MX) I_NEW=1
                        ENDIF
                        IF (IS_Y_PERIODIC) THEN
                            IF (J_NEW < 1) J_NEW=MY
                            IF (J_NEW > MY) J_NEW=1
                        ENDIF
                        I_OLD=I_NEW
                        J_OLD=J_NEW
                        DRAINAGE_AREA(I_OLD,J_OLD)=DRAINAGE_AREA(I_OLD,J_OLD)+CELL_AREA+RIVERADD

                        DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)+QADD
                    ENDIF
                    !      *******************************************************************
                    !       Check for exit point if non-periodic lower boundary (j_old=my)
                    !            and also check for a new
                    !           basin, branching if at an exit
                    !      *******************************************************************
                    IF (DO_FLOW_BOUNDARIES) THEN
                        IF (FLOW_DIRECTION(I_OLD,J_OLD) == 1) THEN
                            IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) THEN
                                NEWBASIN=.TRUE.
                                TOLAKE=.FALSE.
                                NBOTTOMBASIN=NBOTTOMBASIN+1
                                IS_EXIT_BASIN(CURRENT_BASIN+1)=.TRUE.
                            ENDIF
                            EXIT L120
                        ENDIF
                    ELSE
                        IF (IS_Y_PERIODIC) THEN
                            IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) THEN
                                IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) NEWBASIN=.TRUE.
                                EXIT L120
                            ENDIF
                        ELSE
                            IF (J_OLD == MY) THEN
                                IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) THEN
                                    NEWBASIN=.TRUE.
                                    TOLAKE=.FALSE.
                                    NBOTTOMBASIN=NBOTTOMBASIN+1
                                    IS_EXIT_BASIN(CURRENT_BASIN+1)=.TRUE.
                                ENDIF
                                EXIT L120
                            ENDIF                    
                        ENDIF
                    ENDIF
                    !      *******************************************************************
                    !       If not aborted and not an exit point, continue downstream,
                    !       incrementing the drainage area and sediment yield
                    !      *******************************************************************
                    IF (ABORT_SIMULATION) EXIT L120
                ENDDO L120
                !      *******************************************************************
                !       If we are here then we have reached the end of a basin - either
                !       an outlet or depression
                !       First check for aborted run and exit
                !      *******************************************************************
                IF (ABORT_SIMULATION) THEN
                    IABORTMIN=I_OLD-3
                    IABORTMAX=I_OLD+3
                    JABORTMIN=J_OLD-3
                    JABORTMAX=J_OLD+3
                    WRITE(OUTHIST,814) I_OLD,J_OLD,BASIN_NUMBER(I_OLD,J_OLD),DRAINAGE_AREA(I_OLD,J_OLD)
                    814         FORMAT(' *** INFINITE CYCLE AT ',2I6,' ***'  &
                    ,/,' AT BASIN NUMBER: ',I5,' AR=',G12.5)
                    RETURN
                ENDIF
                !      *******************************************************************
                !       If this is a new basin then increment the basin counter and assign
                !       the current_basin to the basin endpoint
                !      *******************************************************************
                IF (NEWBASIN) THEN
                    CURRENT_BASIN=CURRENT_BASIN+1
                    BASIN_NUMBER(I_OLD,J_OLD)=CURRENT_BASIN
                    IF (TOLAKE) THEN
                        ENCLOSED(CURRENT_BASIN)=.TRUE.
                        !      *******************************************************************
                        !      Look for surrounding points with same elevation - if points have
                        !          the same elevation then we make sure they are included in the same
                        !          basin number
                        !      Initially look only at surrounding 8 points but widen the search
                        !         if a contiguous point at the same elevation is found
                        !      *****************************************************************
                        KAROUND=1
                        NCONTIG=0
                        IF ((FLOW_DIRECTION(I_OLD,J_OLD) < 1).OR.(D8_GRADIENT(I_OLD,J_OLD) <= 0.0)) THEN
                            L667: DO
                                ISFOUND=.FALSE.
                                IF (IS_X_PERIODIC) THEN
                                    ILOW=I_OLD-KAROUND
                                    IHIGH=I_OLD+KAROUND
                                    IF (ILOW < 1) ILOW=MX+ILOW
                                    IF (IHIGH > MX) IHIGH=IHIGH-MX
                                ELSE
                                    ILOW=MAX(1,I_OLD-KAROUND)
                                    IHIGH=MIN(MX,I_OLD+KAROUND)
                                ENDIF
                                IF (IS_Y_PERIODIC) THEN
                                    JLOW=J_OLD-KAROUND
                                    JHIGH=J_OLD+KAROUND
                                    IF (JLOW < 1) JLOW=MY+JLOW
                                    IF (JHIGH > MY) JHIGH=JHIGH-MY
                                ELSE
                                    JLOW=MAX(1,J_OLD-KAROUND)
                                    JHIGH=MIN(MYY,J_OLD+KAROUND)
                                ENDIF
                                II=ILOW
                                L665: DO
                                    JJ=JLOW
                                    L669: DO
                                        YCONTINUE=.TRUE.
                                        IF (((II == ILOW).OR.(II == IHIGH)).AND.  &
                                        ((JJ == JLOW).OR. &
                                        (JJ == JHIGH))) THEN
                                            !      *******************************************************************
                                            !       Ignore the present point and any exit points (jj=my)
                                            !      *******************************************************************
                                            IF ((II == I_OLD).AND.(JJ == J_OLD)) &
                                            YCONTINUE=.FALSE.
                                            IF (IS_Y_PERIODIC) THEN
                                                IF ((II == ILOWEST).AND.(JJ == JLOWEST)) &
                                                YCONTINUE=.FALSE. !    goto 666
                                            ELSE
                                                IF (JJ == MY) YCONTINUE=.FALSE. !    goto 666
                                            ENDIF
                                            IF (YCONTINUE) THEN
                                                IF(ELEVATION(II,JJ) == ELEVATION(I_OLD,J_OLD)) THEN
                                                    !      *******************************************************************
                                                    !       If we are looking at more than 8 points then things are a little
                                                    !       more subtle -- We need to make sure that the point at an equal
                                                    !       elevation is actually contiguous with an already identified basin
                                                    !       point - So we look around the identified point for one already
                                                    !       found
                                                    !      *******************************************************************
                                                    IF(KAROUND > 1) THEN
                                                        KKSQUARE=(KAROUND-1)**2
                                                        IF (IS_X_PERIODIC) THEN
                                                            ILL=II-1
                                                            IF (ILL < 1) ILL=ILL+MX
                                                            IUU=II+1
                                                            IF (IUU > MX) IUU=IUU-MX
                                                        ELSE
                                                            ILL=MAX(II-1,1)
                                                            IUU=MIN(II+1,MX)
                                                        ENDIF
                                                        IF (IS_Y_PERIODIC) THEN
                                                            JLL=JJ-1
                                                            IF (JLL < 1) JLL=JLL+MY
                                                            JUU=JJ+1
                                                            IF (JUU > MY) JUU=JUU-MY
                                                        ELSE
                                                            JLL=MAX(JJ-1,1)
                                                            JUU=MIN(JJ+1,MY-1)
                                                        ENDIF
                                                        CONTIGUOUS=.FALSE.
                                                        III=ILL
                                                        L664: DO
                                                            JJJ=JLL
                                                            L670: DO
                                                                XCONTINUE=.TRUE.
                                                                IF ((III == II).AND.(JJJ == JJ)) &
                                                                XCONTINUE=.FALSE.
                                                                IF ((((III-II)**2) <= KKSQUARE)  &
                                                                .AND.(((JJJ-JJ)**2) <= KKSQUARE)) &
                                                                XCONTINUE=.FALSE.
                                                                IF (IS_Y_PERIODIC) THEN
                                                                    IF ((II == ILOWEST).AND. &
                                                                    (JJ == JLOWEST)) &
                                                                    XCONTINUE=.FALSE.
                                                                ELSE
                                                                    IF (JJ == MY) &
                                                                    XCONTINUE=.FALSE.
                                                                ENDIF
                                                                IF (XCONTINUE) THEN
                                                                    IF (ELEVATION(III,JJJ) == ELEVATION(II,JJ)) THEN
                                                                        IF (  (((III-II)**2) <= 1).AND.(((JJJ-JJ)**2) <= 1)) THEN
                                                                            CONTIGUOUS=.TRUE.
                                                                        ELSE
                                                                            EXIT L667
                                                                        ENDIF
                                                                    ENDIF
                                                                ENDIF
                                                                IF (JJJ /= JUU) THEN
                                                                    JJJ=JJJ+1
                                                                    IF (IS_Y_PERIODIC) THEN
                                                                        IF (JJJ > MY) JJJ=1
                                                                    ENDIF
                                                                    CYCLE L670
                                                                ENDIF
                                                                EXIT L670
                                                            ENDDO L670
                                                            IF (III == IUU) EXIT L664
                                                            III=III+1
                                                            IF (IS_X_PERIODIC) THEN
                                                                IF (III > MX) III=1
                                                            ENDIF
                                                        ENDDO L664
                                                    ELSE
                                                        CONTIGUOUS=.TRUE.
                                                    ENDIF
                                                    !      *******************************************************************
                                                    !       If the equal-elevation point is contiguous it is given the same
                                                    !       basin number
                                                    !      *******************************************************************
                                                    IF (CONTIGUOUS) THEN
                                                        BASIN_NUMBER(II,JJ)=CURRENT_BASIN
                                                        NCONTIG=NCONTIG+1
                                                        IF ((.NOT.IS_Y_PERIODIC).AND. &
                                                        (JJ == MY-1)) THEN
                                                        ELSE
                                                            ISFOUND=.TRUE.
                                                        ENDIF
                                                    ENDIF
                                                ENDIF
                                            ENDIF
                                        ENDIF
                                        IF (JJ /= JHIGH) THEN
                                            JJ=JJ+1
                                            IF (IS_Y_PERIODIC) THEN
                                                IF (JJ > MY) JJ=1
                                            ENDIF
                                            CYCLE L669
                                        ENDIF
                                        EXIT L669
                                    ENDDO L669
                                    IF (II == IHIGH) EXIT L665
                                    II=II+1
                                    IF (IS_X_PERIODIC) THEN
                                        IF (II > MX) II=1
                                    ENDIF
                                ENDDO L665
                                KAROUND=KAROUND+1
                                IF (KAROUND > NMAX) THEN
                                    RETURN
                                ENDIF
                                !      *******************************************************************
                                !        We widen the search if a contiguous point at the same elevation
                                !        has been found to include points possibly contiguous to the new
                                !        point but not the original point
                                !      *******************************************************************
                                IF (.NOT.ISFOUND) EXIT L667
                            ENDDO L667
                        ENDIF
                        !      *******************************************************************
                        !        End of code relating to grouping surrounding points at the same
                        !            elevation
                        !      *******************************************************************
                    ENDIF
                ENDIF
                !      *******************************************************************
                !       We need to retrace our steps from the initial location, setting the
                !       basin number of all points along the path to that of the exit or sink point
                !      *******************************************************************
                NUMBBASIN=BASIN_NUMBER(I_OLD,J_OLD)
                I_OLD=I
                J_OLD=J
                L140: DO
                    !      *******************************************************************
                    !       End if we reach an already-identified basin location
                    !      *******************************************************************
                    IF (BASIN_NUMBER(I_OLD,J_OLD) /= 0) THEN
                        EXIT
                    ELSE
                        BASIN_NUMBER(I_OLD,J_OLD)=NUMBBASIN
                        LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                        J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                        I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                        IF (DO_FLOW_BOUNDARIES) THEN
                            IF (FLOW_DIRECTION(I_NEW,J_NEW) == 1) THEN
                                EXIT
                            ENDIF
                        ELSE
                            IF (IS_X_PERIODIC) THEN
                                IF (I_NEW < 1) I_NEW = MX
                                IF (I_NEW > MX) I_NEW=1
                            ENDIF
                            IF (IS_Y_PERIODIC) THEN
                                IF (J_NEW < 1) J_NEW=MY
                                IF (J_NEW > MY) J_NEW=1
                            ELSE
                                IF (J_NEW == MY) EXIT
                            ENDIF
                        ENDIF
                        I_OLD=I_NEW
                        J_OLD=J_NEW
                    ENDIF
                ENDDO L140
            ENDDO M110
        ENDDO L110
        !      *******************************************************************
        !       The initial drainage areas and basin numbers have now been determined
        !       Check to make sure that we have not exceeded the dimensioned
        !       number of basins
        !       Also we have determined the sediment yield from channel erosion
        !      *******************************************************************

        IF (CURRENT_BASIN > LMMX) THEN
            WRITE(OUTHIST,933) CURRENT_BASIN
            933     FORMAT(' TOO MANY BASINS - REDIMENSION.  CURRENT_BASIN=',I6)
            RETURN
        ENDIF
        !      *******************************************************************
        !       Look at boundary points and set values for those that have not
        !       received drainage from upstream to unity area and assign a basin number
        !      *******************************************************************
        IF (DO_FLOW_BOUNDARIES) THEN
            DO J=1,MY
                DO I=1,MX
                    IF (FLOW_DIRECTION(I,J) == 1) THEN
                        IF (DISCHARGE(I,J) == 0.0) THEN
                            DISCHARGE(I,J)=DISCHARGE_FROM_CELL(I,J)
                        ENDIF
                        IF (BASIN_NUMBER(I,J) == 0) THEN
                            CURRENT_BASIN=CURRENT_BASIN+1
                            BASIN_NUMBER(I,J)=CURRENT_BASIN
                            NBOTTOMADD=NBOTTOMADD+1
                            IS_EXIT_BASIN(CURRENT_BASIN)=.TRUE.
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
        ELSE
            IF (IS_Y_PERIODIC) THEN
                IF (DRAINAGE_AREA(ILOWEST,JLOWEST) == 0) THEN
                    DRAINAGE_AREA(ILOWEST,JLOWEST)=CELL_AREA
                    DISCHARGE(ILOWEST,JLOWEST)=DISCHARGE_FROM_CELL(ILOWEST,JLOWEST)
                ENDIF
                IF (BASIN_NUMBER(ILOWEST,JLOWEST) == 0) THEN
                    CURRENT_BASIN=CURRENT_BASIN+1
                    BASIN_NUMBER(ILOWEST,JLOWEST)=CURRENT_BASIN
                ENDIF
            ELSE
                DO  I=1,MX
                    IF (DRAINAGE_AREA(I,MY) == 0) THEN
                        DRAINAGE_AREA(I,MY)=CELL_AREA
                        DISCHARGE(I,MY)=DISCHARGE_FROM_CELL(I,MY)
                    ENDIF
                    IF (BASIN_NUMBER(I,MY) == 0) THEN
                        CURRENT_BASIN=CURRENT_BASIN+1
                        BASIN_NUMBER(I,MY)=CURRENT_BASIN
                        NBOTTOMADD=NBOTTOMADD+1
                        IS_EXIT_BASIN(CURRENT_BASIN)=.TRUE.
                    ENDIF
                    IF (.NOT.IS_EXIT_BASIN(BASIN_NUMBER(I,MY))) THEN
                        WRITE(*,8227) 
                        8227 FORMAT(' PROBLEM IN LOWER BOUNDARY BASIN IDENTIFICATION')
                        IS_EXIT_BASIN(BASIN_NUMBER(I,MY))=.TRUE.
                    ENDIF
                ENDDO
            ENDIF
        ENDIF
        TOTAL_BASINS=CURRENT_BASIN
        DO  I=2,MX-1
            DO  J=2,MYY

                IF (BASIN_NUMBER(I,J) == 0) THEN
                    DO  JJ=J-1,J+1
                        DO  II=I-1,I+1
                            WRITE (OUTHIST,4955) II,JJ,FLOW_DIRECTION(II,JJ),BASIN_NUMBER(II,JJ)
                            4955 FORMAT('NOBASIN, II=',I5,' JJ=',I5,' ID=',I4,' N=',I5)
                        ENDDO
                    ENDDO
                    WRITE (OUTHIST,4956)
                    4956  FORMAT(' ')
                ENDIF
            ENDDO
        ENDDO
        IF (.NOT.IS_Y_PERIODIC) THEN
            DO  K=1,TOTAL_BASINS
                KKK=0
                DO  J=1,MY
                    DO I=1,MX
                        IF (BASIN_NUMBER(I,J) == K) KKK=KKK+1
                    ENDDO
                ENDDO
                IF (KKK == 0) ENCLOSED(K)=.FALSE.
                IF ((.NOT.ENCLOSED(K)).AND.(.NOT.IS_EXIT_BASIN(K))) THEN
                    9221   FORMAT('BAD NOENCLOSE=',I7)
                ENDIF
            ENDDO
        ENDIF
        ROUTED_DISCHARGE=DISCHARGE
        !      *******************************************************************
        !       If there are enclosed basins (*and if runoff overflow is permitted*)
        !       then we identify those enclosed basins which need to be examined
        !       for grouping and flow routing to the outlets - skip all this
        !       otherwise
        !      *******************************************************************
        IF ((COMPLETE_RUNOFF).AND.(.NOT.NO_LAKES_EXIST)) THEN
            DO  K=1,TOTAL_BASINS
                IF (ENCLOSED(K)) THEN
                    NEEDTODO(K)=.TRUE.
                ELSE
                    NEEDTODO(K)=.FALSE.
                ENDIF
            ENDDO
            !      *******************************************************************
            !       We need to keep regrouping basins and identifying lake outlets and
            !       cumulative drainage areas in the case where there is enclosed basins
            !       and lake overflow - redo is an index of whether another level of
            !       basin consolidation is needed
            !      *******************************************************************
            L250: DO
                REDO=.FALSE.
                !      *******************************************************************
                !       If we still need to work in a basin (an outlet location has still
                !       not been determined) then we reset the outlet elevation and
                !       lake area
                !      *******************************************************************
                L150: DO  K=1,TOTAL_BASINS
                    IF (NEEDTODO(K)) THEN
                        LAKE_OUTLET_ELEVATION(K)=1.0E+25
                        !      *******************************************************************
                        !         Determine the present lake drainage area
                        !      *******************************************************************
                        BASIN_DRAINAGE_AREA(K)=0
                        LOCAL_BASIN_DISCHARGE(K)=0.0
                        L160: DO  J=1,MYY
                            M160: DO  I=1,MX
                                IF (BASIN_NUMBER(I,J) == K) THEN
                                    IF (FLOW_DIRECTION(I,J) < 0) THEN
                                        BASIN_DRAINAGE_AREA(K)=   &
                                        BASIN_DRAINAGE_AREA(K)+DRAINAGE_AREA(I,J)
                                        LOCAL_BASIN_DISCHARGE(K)=LOCAL_BASIN_DISCHARGE(K)+DISCHARGE(I,J)
                                    ENDIF
                                    !      *******************************************************************
                                    !       We need to look at all points on the periphery of the basin to
                                    !       determine if it is an outlet
                                    !      *******************************************************************
                                    L170: DO  L=2,9
                                        J_NEW=J+DOWNSTREAM(L,2)
                                        IF (IS_Y_PERIODIC) THEN
                                            IF (J_NEW < 1) J_NEW = MY
                                            IF (J_NEW > MY) J_NEW=1
                                        ELSE
                                            IF (J_NEW < 1) CYCLE L170
                                        ENDIF
                                        I_NEW=I+DOWNSTREAM(L,1)
                                        IF (IS_X_PERIODIC) THEN
                                            IF (I_NEW < 1) I_NEW = MX
                                            IF (I_NEW > MX) I_NEW=1
                                        ELSE
                                            IF ((I_NEW < 1).OR.(I_NEW > MX)) CYCLE L170
                                        ENDIF
                                        !      *******************************************************************
                                        !       Do the following if the present point is at a peripheral location
                                        !       Check to see if the periphery point or a contiguous point outside
                                        !         the basin is lower than any identified lake outlet
                                        !      *******************************************************************
                                        IF (BASIN_NUMBER(I_NEW,J_NEW) /= K) THEN
                                            XCOMP=MAX(ELEVATION(I_NEW,J_NEW),ELEVATION(I,J))
                                            !      *******************************************************************
                                            !       If this is a potential lake outlet then there are two cases:
                                            !         1.  The basin peripery is lower than any contiguous point outside
                                            !          the basin.  If this is the case (outer is true) the outlet is
                                            !             identified with the outside point and the outlet elevation is
                                            !             equal to that of the outside point.
                                            !         2.  The basin periphery is higher than the lowest contiguous point
                                            !          outside.  The outlet is the periphery point and the lake outlet
                                            !          elevation is its elevation.  In this case we need to record where
                                            !             the flow is going to and its direction (ifrom,jfrom,kfrom)
                                            !      *******************************************************************
                                            IF (XCOMP < LAKE_OUTLET_ELEVATION(K)) THEN
                                                IF (ELEVATION(I,J) <= ELEVATION(I_NEW,J_NEW)) THEN
                                                    OUTER(K)=.TRUE.
                                                    IOUTLET=I_NEW
                                                    JOUTLET=J_NEW
                                                    LAKE_OUTLET_ELEVATION(K)=ELEVATION(I_NEW,J_NEW)
                                                ELSE
                                                    OUTER(K)=.FALSE.
                                                    IOUTLET=I
                                                    JOUTLET=J
                                                    LAKE_OUTLET_ELEVATION(K)=ELEVATION(I,J)
                                                    IFROM=I_NEW
                                                    JFROM=J_NEW
                                                    KKFROM=L
                                                ENDIF
                                                !      *******************************************************************
                                                !       downstream_basin is the basin into which the present basin drains
                                                !      *******************************************************************
                                                DOWNSTREAM_BASIN(K)=BASIN_NUMBER(I_NEW,J_NEW)
                                            ENDIF
                                        ENDIF
                                    ENDDO L170
                                ENDIF
                            ENDDO M160
                        ENDDO L160
                        !      *******************************************************************
                        !       We have now identified the lowest outlet for the present basin
                        !       so we set its outlet location.  If the basin outlet is on the
                        !       periphery (not outer) then we need to redefine the drainage direction
                        !       and gradient of the peripery point for flow to the outlet rather
                        !       than into the present basin
                        !      *******************************************************************
                        I_OUTFLOW(K)=IOUTLET
                        J_OUTFLOW(K)=JOUTLET
                        IF (.NOT.OUTER(K)) THEN                      
                            if (FLOW_DIRECTION(IOUTLET,JOUTLET)/=1) FLOW_DIRECTION(IOUTLET,JOUTLET)=KKFROM
                            D8_GRADIENT(IOUTLET,JOUTLET)=(ELEVATION(IOUTLET,JOUTLET) &
                            -ELEVATION(IFROM,JFROM))/CELL_SIZE
                        ENDIF
                    ENDIF
                    !      *******************************************************************
                    !       We have now done this basin - we don't need to worry about it for
                    !       the time being.
                    !      *******************************************************************
                    NEEDTODO(K)=.FALSE.
                ENDDO L150
                !      *******************************************************************
                !       We have now completed one stage of basin agglomeration.  We need
                !       to check for the case where lake basins drain into each other
                !       if this is true then redo is set true and needtodo(k) is set true
                !       redo is false if each chain of enclosed basin outlets following downstream
                !       eventually leads to an external outlet
                !      *******************************************************************
                K=0
                KKK=0
                !      *******************************************************************
                !         Cycle through all basins
                !      *******************************************************************
                L260: DO
                    K=K+1
                    KKCOUNT=0
                    !      *******************************************************************
                    !         Exit from loop if we have looked at all basins
                    !      *******************************************************************
                    IF (K >= TOTAL_BASINS) EXIT L260
                    !      *******************************************************************
                    !         If we find an enclosed basin we need to check that tracing successive
                    !         outlets downstream eventually leads to an external outlet, otherwise
                    !         we need an additional stage of basin agglomeration
                    !      *******************************************************************
                    IF (ENCLOSED(K)) THEN
                        KK=K
                        IF ((K > LMMX).OR.(KK > LMMX).OR.(KKK > LMMX) &
                        .OR.(K == 0).OR.(KK == 0)) THEN
                            WRITE (*,4623) K,KK,KKK,DOWNSTREAM_BASIN(KK),  &
                            LMMX
                            4623  FORMAT(' ABORTING, K=',I7,' KK=',I7,' KKK=',I7, &
                            ' DOWNSTREAM_BASIN(KK)=',I7,' LMMX=',I7)
                            CALL BASIN_REPORT()
                            RETURN
                        ENDIF
                        L623: DO
                            KKK=DOWNSTREAM_BASIN(KK)
                            IF ((KKK > LMMX).OR.(KKK == 0)) THEN
                                WRITE (*,5623) KK,KKK
                                5623  FORMAT(' ABORTING, KK=',I7,' KKK=',I7)
                                CALL BASIN_REPORT()
                                RETURN
                            ENDIF

                            IF (DOWNSTREAM_BASIN(KKK) == K) EXIT L623
                            IF (.NOT.ENCLOSED(KKK)) CYCLE L260
                            IF (DO_FLOW_BOUNDARIES) THEN
                                IF ((I_OUTFLOW(KKK) == 1).OR.(I_OUTFLOW(KKK) == MX).OR. &
                                (J_OUTFLOW(KKK) == 1).OR.(J_OUTFLOW(KKK) == MY)) THEN
                                    CYCLE L260
                                ENDIF
                            ELSE
                                IF (IS_Y_PERIODIC) THEN
                                    IF ((I_OUTFLOW(KKK) == ILOWEST).AND.(J_OUTFLOW(KKK) == JLOWEST)) &
                                    CYCLE L260
                                ELSE
                                    IF (J_OUTFLOW(KKK) == MY) CYCLE L260 !    goto 260
                                ENDIF
                            ENDIF
                            KK=KKK
                            KKCOUNT=KKCOUNT+1
                            IF (KKCOUNT > 500) CYCLE L260
                        ENDDO L623
                        !      *******************************************************************
                        !       This is the case of mutually-draining basins
                        !      *******************************************************************
                        KK=K
                        KEND=KKK
                        REDO=.TRUE.
                        NEEDTODO(K)=.TRUE.
                        KKK=DOWNSTREAM_BASIN(KK)
                        !      *******************************************************************
                        !       If there are two mutually-draining basins, then they are arbitrarily
                        !       agglomerated by reassigning one basins points (kkk) to the other (k)
                        !      *******************************************************************
                        DO  J=1,MYY
                            DO  I=1,MX
                                IF (BASIN_NUMBER(I,J) == KKK) BASIN_NUMBER(I,J)=K
                            ENDDO
                        ENDDO
                        !      *******************************************************************
                        !       If other basins drain to the reassigned drainage basin, then their
                        !         exit basin have to be reassined to the agglomerated basin
                        !      *******************************************************************
                        L627: DO  KKKK=1,TOTAL_BASINS
                            IF (ENCLOSED(KKKK)) THEN
                                IF (KKKK == K) CYCLE L627
                                IF (KKKK == KKK) CYCLE L627
                                IF (DOWNSTREAM_BASIN(KKKK) == KKK) DOWNSTREAM_BASIN(KKKK)=K
                            ENDIF
                        ENDDO L627
                        !      *******************************************************************
                        !       We need to deactivate the various values of the deactivated
                        !          (agglomerated) basin (kkk)
                        !      *******************************************************************
                        ENCLOSED(KKK)=.FALSE.
                        BASIN_DRAINAGE_AREA(KKK)=0.0
                        LOCAL_BASIN_DISCHARGE(KKK)=0.0
                        I_OUTFLOW(KKK)=-1
                        J_OUTFLOW(KKK)=-1
                        KK=DOWNSTREAM_BASIN(KKK)
                        TOBEADDED(KKK)=0
                        DOWNSTREAM_BASIN(KKK)=0
                        !      *******************************************************************
                        !       Look at additional basins if we have not examined all of them
                        !      *******************************************************************
                        IF (KKK == KEND) CYCLE L260
                        KKK=KK
                    ENDIF
                    !      *******************************************************************
                    !       Go on to the next basin
                    !      *******************************************************************
                ENDDO L260
                !      *******************************************************************
                !       if we have not completely agglomerated drainage basins so that they
                !       flow to the exterior, then cycle back to do it again
                !      *******************************************************************
                IF (.NOT.REDO) THEN
                    EXIT L250
                ENDIF
            ENDDO L250
            !      *******************************************************************
            !       We now need to route flow downstream through all lakes, such that
            !       the drainage area downstream from the outlet includes the lake area
            !       and all drainage basins and upstream lakes draining into the given
            !       lake
            !      *******************************************************************
            L750: DO K=1,TOTAL_BASINS
                KK=K
                NNN=0
                !      *******************************************************************
                !         Check for enclosed drainage basins
                !      *******************************************************************
                L760: DO
                    IF (ENCLOSED(KK)) THEN
                        KKKK=1
                        IF (ENCLOSED(DOWNSTREAM_BASIN(KK))) THEN
                            KKK=1
                        ELSE
                            KKK=0
                        ENDIF
                    ELSE
                        KKKK=0
                    ENDIF
                    IF (ENCLOSED(KK)) THEN
                        !      *******************************************************************
                        !         If its outlet basin is also enclosed, then we need to continue
                        !         downstream.  tobeadded(k) is the cumulative added drainage from
                        !         enclosed basins upstream from the given basin
                        !      *******************************************************************
                        IF (ENCLOSED(DOWNSTREAM_BASIN(KK))) THEN
                            TOBEADDED(DOWNSTREAM_BASIN(KK))=  &
                            TOBEADDED(DOWNSTREAM_BASIN(KK))+BASIN_DRAINAGE_AREA(K)
                            QTOADD(DOWNSTREAM_BASIN(KK))=QTOADD(DOWNSTREAM_BASIN(KK)) &
                            +LOCAL_BASIN_DISCHARGE(K)
                            KK=DOWNSTREAM_BASIN(KK)
                            3444 FORMAT(' KK=',I7,' TO',I7,' ENC=',I5,' TO ENC=',I5,' N=',I7)
                            NNN=NNN+1
                            !      *******************************************************************
                            !       Check for infinite loops between basins (basins draining into one
                            !       another or through a loop of connecting basins).  Abort if this
                            !       occurs
                            !      *******************************************************************
                            IF (NNN > TOTAL_BASINS) THEN
                                WRITE(*,799) ITERATION,KK
                                WRITE(OUTHIST,799) ITERATION,KK
                                799                 FORMAT(' **** BASIN SCRAMBLE AT ITERATION',I9, &
                                ' AT BASIN',I5)
                                CALL BASIN_REPORT()
                                RETURN
                            ENDIF
                            CYCLE L760
                        ENDIF
                        GOTO 8443
                    ELSE
                        GOTO 8443
                    ENDIF
                ENDDO L760
                8443 CONTINUE
                IF (ENCLOSED(K)) THEN
                    IF(ENCLOSED(DOWNSTREAM_BASIN(K))) THEN
                        IF (ENCLOSED(DOWNSTREAM_BASIN(KK))) THEN
                            KKK=1
                        ELSE
                            KKK=0
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO L750
            !      *******************************************************************
            !       Now that we have calculated contributions from upstream enclosed basins
            !         we need to calculate the total lake drainage area of each basin
            !      *******************************************************************
            DO  K=1,TOTAL_BASINS
                IF (ENCLOSED(K)) THEN
                    BASIN_DRAINAGE_AREA(K)=BASIN_DRAINAGE_AREA(K)+TOBEADDED(K)
                    LOCAL_BASIN_DISCHARGE(K)=LOCAL_BASIN_DISCHARGE(K)+QTOADD(K)
                ELSE
                    IF(.NOT.IS_Y_PERIODIC) THEN
                        LAKE_OUTLET_ELEVATION(K)=ELEVATION(1,MY)
                    ENDIF
                ENDIF
            ENDDO
            !      *******************************************************************
            !       We now go downstream from each enclosed basin outlet and increment the
            !       drainage area by the amount of drainage flowing through the outlet
            !      *******************************************************************
            L180: DO K=1,TOTAL_BASINS
                IF (ENCLOSED(K)) THEN
                    I_OLD=I_OUTFLOW(K)
                    J_OLD=J_OUTFLOW(K)
                    IF ((I_OLD == I_OUTFLOW(DOWNSTREAM_BASIN(K))).AND. &
                    (J_OLD == J_OUTFLOW(DOWNSTREAM_BASIN(K)))) &
                    CYCLE L180 !    go to 180
                    IF ((I_OLD < 1).OR.(I_OLD > MX).OR. &
                    (J_OLD < 1).OR.(J_OLD > MY)) THEN
                        WRITE(OUTHIST,377) K,I_OLD,J_OLD
                        377 FORMAT(' ***BASIN INTEGRITY PROBLEM AT BASIN' &
                        ,I5,', I_OLD & J_OLD=', &
                        2I5)
                    ENDIF
                    ADDAREA=BASIN_DRAINAGE_AREA(K)
                    ADDQ=LOCAL_BASIN_DISCHARGE(K)
                    IF(OUTER(K)) THEN
                        DRAINAGE_AREA(I_OLD,J_OLD)=DRAINAGE_AREA(I_OLD,J_OLD)+ADDAREA
                        DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)+ADDQ
                        ROUTED_DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)
                    ELSE
                        DRAINAGE_AREA(I_OLD,J_OLD)=ADDAREA
                        DISCHARGE(I_OLD,J_OLD)=ADDQ
                        ROUTED_DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)
                    ENDIF
                    IF (DO_FLOW_BOUNDARIES) THEN
                        IF (FLOW_DIRECTION(I_OLD,J_OLD) == 1) THEN
                            CYCLE L180
                        ENDIF
                    ELSE
                        IF (IS_Y_PERIODIC) THEN
                            IF ((I_OLD == ILOWEST).AND. &
                            (J_OLD == JLOWEST)) CYCLE L180
                        ELSE
                            IF (J_OLD == MY) CYCLE L180
                        ENDIF
                    ENDIF
                    LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                    IF (LOCAL_DIRECTION < 0) CYCLE L180
                    L190: DO
                        J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                        I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                        IF (IS_X_PERIODIC) THEN
                            IF (I_NEW < 1) I_NEW = MX
                            IF (I_NEW > MX) I_NEW=1
                        ENDIF
                        IF (IS_Y_PERIODIC) THEN
                            IF (J_NEW < 1) J_NEW=MY
                            IF (J_NEW > MY) J_NEW=1
                        ENDIF
                        I_OLD=I_NEW
                        J_OLD=J_NEW
                        IF (BASIN_NUMBER(I_OLD,J_OLD) /= DOWNSTREAM_BASIN(K)) THEN
                            IF (OUTER(K)) THEN
                                KIJ=1
                            ELSE 
                                KIJ=0
                            ENDIF
                            IF (OUTER(BASIN_NUMBER(I_OLD,J_OLD))) THEN
                                KIO=1
                            ELSE
                                KIO=0
                            ENDIF
                            IF (OUTER(DOWNSTREAM_BASIN(K))) THEN
                                KII=1
                            ELSE
                                KII=0
                            ENDIF
                            WRITE(OUTHIST,823) K,DOWNSTREAM_BASIN(K),&
                            BASIN_NUMBER(I_OLD,J_OLD),I_OLD,J_OLD &
                            ,DRAINAGE_AREA(I_OLD,J_OLD),ADDAREA,KIJ,KIO,KII
                            823               FORMAT(' PROBLEM IN ADDAREA FROM BASIN',I5, &
                            ' TO BASIN' &
                            ,I5,/,' NOW IN BASIN',I5,' AT I=',I5,' J=',I5, &
                            ' AREA=',G12.5,' ADDAREA=',G12.5,' OUT_K=',I3,' OUT_OLD=',I3,' OUT_DOWN=',I3)
                            EXIT L190 !    goto 200
                            CALL BASIN_REPORT()
                            ABORT_SIMULATION=.TRUE.
                            RETURN
                        ENDIF
                        IF ((I_OLD == I_OUTFLOW(DOWNSTREAM_BASIN(K))).AND.(J_OLD == J_OUTFLOW(DOWNSTREAM_BASIN(K)))) &
                        EXIT L190
                        DRAINAGE_AREA(I_OLD,J_OLD)=DRAINAGE_AREA(I_OLD,J_OLD)+ADDAREA
                        DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)+ADDQ
                        ROUTED_DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)
                        LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                        IF (LOCAL_DIRECTION < 0) EXIT L190
                        IF (DO_FLOW_BOUNDARIES) THEN
                            IF (LOCAL_DIRECTION == 1) THEN
                                EXIT L190
                            ENDIF
                        ELSE
                            IF (IS_Y_PERIODIC) THEN
                                IF ((I_OLD == ILOWEST).AND. &
                                (J_OLD == JLOWEST)) EXIT L190
                            ELSE
                                IF (J_OLD == MY) EXIT L190
                            ENDIF
                        ENDIF
                    ENDDO L190
                ENDIF
            ENDDO L180
            !      *******************************************************************
            !       Finished with basin and drainage area determination
            !      *******************************************************************
        ENDIF
        ADDAREA=0
        ADDQ=0.0
        IF (.NOT.IS_Y_PERIODIC) THEN
            DO  I=1,MX
                ADDQ=ADDQ+DISCHARGE(I,MY)
                ADDAREA=ADDAREA+DRAINAGE_AREA(I,MY)
            ENDDO
            IF (ADDAREA < ((MX-1)*(MY-1)-50)) THEN
                WRITE (OUTHIST,889) ITERATION,ADDAREA
                889     FORMAT(' ITER=',I9,' TOTAL AREA = ',G12.5)
            ENDIF
        ENDIF
        DO J=1,MY
            DO I=1,MX
                IF (FLOW_DIRECTION(I,J) < 0) THEN
                    K=BASIN_NUMBER(I,J)
                    IF (ENCLOSED(K)) THEN
                        IF ((I_OUTFLOW(K) == I).AND.(J_OUTFLOW(K) == J)) THEN
                            KKK=0
                            DO JJ=1,MY
                                DO II=1,MX
                                    IF (BASIN_NUMBER(II,JJ) == K) KKK=KKK+1
                                ENDDO
                            ENDDO
                            WRITE(*,2890) I,J,K,KKK
                            2890  FORMAT('DEPRESSION IS OUTFLOW, I=',I5,' J=',I5, &
                            ' BN=',I6,' SIZE=',I6)
                        ENDIF
                    ELSE
                        !WRITE(*,2889) I,J,K
                        2889  FORMAT('NON-ENCLOSED DEPRESSION AT I=',I5,' J=',I5,' BN=',I6)
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        WSUM=0.0
        WNUM=0.0
        !      ********************************************************************
        !        If discharge is not directly proportional to drainage area, then
        !        correct for non-linearity
        !      ********************************************************************
        DO J=1,MY
            DO I=1,MX
                IF (DISCHARGE(I,J) > GLOBAL_MAXIMUM_DISCHARGE) GLOBAL_MAXIMUM_DISCHARGE=DISCHARGE(I,J)
            ENDDO
        ENDDO
        GLOBAL_MAXIMUM_DISCHARGE=GLOBAL_MAXIMUM_DISCHARGE/(CELL_AREA*DISCHARGE_SCALE_FACTOR)
        8549  FORMAT(' MAX. YIELD=',G12.5)
        IF (RESCALE_DISCHARGES)  THEN
            DISCHARGE=DISCHARGE**DISCHARGE_EXPONENT
        ENDIF
        DRAINAGE_AREA=DRAINAGE_AREA+AREAFACTOR
        IF (MODEL_GROUNDWATER) THEN
            IF (NSEEPAVG > 0.0) THEN
                AVGRUNOFF=AVGRUNOFF*CELL_AREA/NSEEPAVG
                AVGSEEP=AVGSEEP*CELL_AREA/NSEEPAVG
                IF (NFRACTAVG > 0.0) THEN
                    AVGSEEPFRACT=AVGSEEPFRACT/NFRACTAVG
                ENDIF
                WRITE(OUTHIST,2881) AVGRUNOFF,AVGSEEP,AVGSEEPFRACT,NSEEPAVG,NFRACTAVG
                2881  FORMAT(' AVG. RUNOFF=',G12.5,' AVG. SEEP FLOW=',G12.5,/, &
                ' AVG. FRACTION OF FLOW FROM SEEPAGE=',G12.5,' N=',G12.5,' N1=',G12.5)
            ENDIF
        ENDIF
        !     **********************************************************************
        !      Determine channel widths
        !     **********************************************************************
        CHANNEL_WIDTH=CHANNEL_WIDTH_CONSTANT*DISCHARGE**CHANNEL_WIDTH_EXPONENT
        !     **********************************************************************
        !      Determine sediment transport through each cell based upon elevation
        !         changes locally and upstream during last iteration.  Sediment is not
        !         routed through depressions
        !     **********************************************************************
        ERODETOADD=0.0
        ERGRAVTOADD=0.0
        !      *********************************************************************
        !       sedadd is the sediment load added to the stream flow by local bedrock
        !         erosion
        !       tsedcont sums up the sediment flux contributions to
        !         local sediment yield - used for debugging
        !         CHANNEL_EROSION_RATE holds the local bedrock channel erosion from the last iteration
        !      *********************************************************************
        !      *********************************************************************
        !        sediment_yield is accumulated downstream - it is assumed that all sediment in
        !         bedrock channels moves instantaneously downstream - i.e., moves much more
        !         rapidly than the timescale for vertical elevation changes
        !      *********************************************************************
        L710: DO  J=1,MYY
            M710: DO  I=1,MX
                I_OLD=I
                J_OLD=J
                IF (.NOT.IS_SEDIMENT_COVERED(I,J)) THEN
                    IF (CHANNEL_EROSION_RATE(I,J) > 0.0) THEN
                        SEDADD=CHANNEL_EROSION_RATE(I,J)*CELL_AREA
                    ELSE
                        SEDADD=0.0
                    ENDIF
                    !      *********************************************************************
                    !       If there is an incoming river from out of the matrix,
                    !         add in its sediment supply
                    !      *********************************************************************
                    IF (HAVE_INFLUENT_RIVERS.AND.IS_INFLUENT_RIVER_LOCATION(I,J)) THEN
                        DO KL=1,NUMBER_OF_INFLUENT_RIVERS
                            IF((I == I_INFLUENT_RIVER(KL)).AND.(J == J_INFLUENT_RIVER(KL)))THEN
                                SEDADD=SEDADD+INFLUENT_RIVER_SEDIMENT_LOAD(KL)
                            ENDIF
                        ENDDO
                    ENDIF
                    TSEDCONT=TSEDCONT+SEDADD
                    TOSEDCOVER=.FALSE.
                ELSE
                    CYCLE M710
                ENDIF
                SEDIMENT_YIELD(I,J)=SEDIMENT_YIELD(I,J)+SEDADD
                IF (SEDIMENT_YIELD(I,J)>ITERATION_MAXIMUM_SEDIMENT_YIELD) &
                    ITERATION_MAXIMUM_SEDIMENT_YIELD=SEDIMENT_YIELD(I,J)
                L720: DO
                    LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                    !     **********************************************************************
                    !      Branch if a drainage exit or depression
                    !     **********************************************************************
                    IF (LOCAL_DIRECTION < 2) THEN
                        !      *********************************************************************
                        !       Increments total work and transport variables
                        !      *********************************************************************
                        IF (DO_FLUVIAL_DETACHMENT) THEN
                            ERODETOADD=ERODETOADD+SEDADD*SQRT(REAL((I_OLD-I)**2+ &
                            (J_OLD-J)**2))
                            ERGRAVTOADD=ERGRAVTOADD+SEDADD*(ELEVATION(I,J)-ELEVATION(I_OLD,J_OLD))
                        ENDIF
                        EXIT L720
                    ENDIF
                    J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                    I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                    IF (IS_X_PERIODIC) THEN
                        IF (I_NEW < 1) I_NEW=MX
                        IF (I_NEW > MX) I_NEW=1
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (J_NEW > MY) J_NEW=1
                        IF (J_NEW < 1) J_NEW=MY
                    ENDIF
                    IF (IS_SEDIMENT_COVERED(I_NEW,J_NEW)) THEN
                        !      *********************************************************************
                        !       Increments total work and transport variables
                        !      *********************************************************************
                        IF (DO_FLUVIAL_DETACHMENT.AND.(.NOT.TOSEDCOVER)) THEN
                            ERODETOADD=ERODETOADD+SEDADD*SQRT(REAL((I_OLD-I)**2+ &
                            (J_OLD-J)**2))
                            ERGRAVTOADD=ERGRAVTOADD+SEDADD*(ELEVATION(I,J)-ELEVATION(I_OLD,J_OLD))
                        ENDIF
                        TOSEDCOVER=.TRUE.
                    ELSE
                        IF (TOSEDCOVER) EXIT L720 !    goto 710
                    ENDIF
                    I_OLD=I_NEW
                    J_OLD=J_NEW
                    !     **********************************************************************
                    !      Otherwise add in contribution from upstream
                    !     **********************************************************************
                    SEDIMENT_YIELD(I_OLD,J_OLD)=SEDIMENT_YIELD(I_OLD,J_OLD) + SEDADD
                   IF (SEDIMENT_YIELD(I_OLD,J_OLD)>ITERATION_MAXIMUM_SEDIMENT_YIELD)  &
                       ITERATION_MAXIMUM_SEDIMENT_YIELD=SEDIMENT_YIELD(I_OLD,J_OLD)
                    IF (DO_FLOW_BOUNDARIES) THEN
                        IF (FLOW_DIRECTION(I_OLD,J_OLD) == 1) THEN
                            EXIT L720
                        ENDIF
                    ELSE
                        IF (IS_Y_PERIODIC) THEN
                            IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) EXIT L720
                        ELSE
                            IF (J_OLD == MY) EXIT L720
                        ENDIF
                    ENDIF
                ENDDO L720
            ENDDO M710
        ENDDO L710
        IF (COMPLETE_RUNOFF) LAKE_SURFACE_ELEVATION=LAKE_OUTLET_ELEVATION
        CALL IS_IT_SUBMERGED()
        LAKE_AREA=0.0
        DO K=1,TOTAL_BASINS
            IF (ENCLOSED(K)) THEN
                DO J=1,MY
                    DO I=1,MX
                        IF (BASIN_NUMBER(I,J) == K) THEN
                            IF (ELEVATION(I,J) <= LAKE_SURFACE_ELEVATION(K)) THEN
                                LAKE_AREA(K)=LAKE_AREA(K)+CELL_AREA
                            ENDIF
                        ENDIF
                    ENDDO
                ENDDO
            ENDIF
        ENDDO

        !      *********************************************************************
        !       Print out some info for debugging purposes
        !      *********************************************************************
        TSEDCONT=TSEDCONT*TIME_INCREMENT
        IF (DO_SEDIMENT_TRANSPORT) THEN
            IF(MOD((ITERATION+1),OUTPUT_PRINT_INTERVAL) == 0) THEN
                WRITE(OUTHIST,764) TSEDCONT
                764           FORMAT(      &
                ' BED EROSION CONTRIBUTION TO SEDIMENT YIELD=',G12.5)
                WRITE(OUTHIST,436)
                436         FORMAT(/,' SEDIMENT YIELD')
                CALL SUMMARIZE_MATRIX_DATA(SEDIMENT_YIELD,THEAVG,THEMAX,THEMIN)
            ENDIF
        ENDIF
        IF (COMPLETE_RUNOFF) THEN
            L8554: DO KK=1,TOTAL_BASINS
                IF (ENCLOSED(KK)) THEN
                    K=DOWNSTREAM_BASIN(KK)
                    KKK=K
                    IF (ENCLOSED(K)) THEN
                        II=I_OUTFLOW(K)
                        JJ=J_OUTFLOW(K)
                        I_OLD=II
                        J_OLD=JJ
                        LOCAL_DIRECTION=FLOW_DIRECTION(I_OLD,J_OLD)
                        IF (SUBMERGED(I_OLD,J_OLD).OR.(LOCAL_DIRECTION < 2)) CYCLE L8554
                        L8556:    DO
                            J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                            I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                            IF (IS_X_PERIODIC) THEN
                                IF (I_NEW < 1) I_NEW=MX
                                IF (I_NEW > MX) I_NEW=1
                            ENDIF
                            IF (IS_Y_PERIODIC) THEN
                                IF (J_NEW > MY) J_NEW=1
                                IF (J_NEW < 1) J_NEW=MY
                            ENDIF
                            LOCAL_DIRECTION=FLOW_DIRECTION(I_NEW,J_NEW)
                            IF (SUBMERGED(I_NEW,J_NEW).OR.(LOCAL_DIRECTION < 2) &
                            .OR.(BASIN_NUMBER(I_NEW,J_NEW) /= BASIN_NUMBER(I_OLD,J_OLD) )) &
                            THEN
                                K=BASIN_NUMBER(I_OLD,J_OLD)
                                IF (ENCLOSED(K)) THEN
                                    IF  ((II >= 1).AND.(II <= MX) &
                                    .AND.(I_OUTFLOW(DOWNSTREAM_BASIN(K)) >= 1) &
                                    .AND.(J_OUTFLOW(K) >= 1).AND.(I_OUTFLOW(K) <= MX).AND. &
                                    (J_OUTFLOW(K) <= MY).AND.(JJ >= 1).AND.(JJ <= MY)) THEN
                                        CALL DRAWLINE(II,JJ,I_OUTFLOW(K), &
                                        J_OUTFLOW(K),ROUTED_DISCHARGE(I_NEW,J_NEW))
                                        !                write(outhist,9054) k,downstream_basin(k),i_old,j_old,i_outflow(downstream_basin(k))&
                                        !                ,j_outflow(downstream_basin(k)),routed_discharge(i_outflow(k),j_outflow(k))
                                        9054 FORMAT('DL, K=',I5,' TO=',I5,' I_OLD=',I5,' J_OLD=',I5,' ITO=',I5, &
                                        ' JTO=',I5,' QE=',G12.5)
                                    ENDIF
                                ENDIF
                                EXIT L8556
                            ENDIF
                            I_OLD=I_NEW
                            J_OLD=J_NEW
                            IF (DO_FLOW_BOUNDARIES) THEN
                                IF (FLOW_DIRECTION(I_OLD,J_OLD) == 1) THEN
                                    EXIT L8556
                                ENDIF
                            ELSE
                                IF (IS_Y_PERIODIC) THEN
                                    IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) EXIT L8556
                                    IF (J_OLD == MY) EXIT L8556
                                ENDIF
                            ENDIF
                        ENDDO L8556
                    ENDIF
                ENDIF
            ENDDO L8554
    ENDIF
            testit=.false.
        if (testit) then
            open(877,file="border.dat",action='write')
            j=1
            do i=1,mx
               kkk1=basin_number(i,j)
               if (submerged(i,j)) then
                   kkk2=1
               else
                   kkk2=0
               endif
               if (kkk1>0) then
                   if (is_exit_basin(kkk1)) then
                       kkk3=1
                   else
                       kkk3=0
                   endif
                   ttt1=lake_outlet_elevation(kkk1)
                   ttt2=basin_drainage_area(kkk1)
                   kkk4=downstream_basin(kkk1)
                   kkk5=i_outflow(kkk1)
                   kkk6=j_outflow(kkk1)
                   if (enclosed(kkk1)) then
                       kkk7 = 1
                   else
                       kkk7=0
                   endif
                   if (overflows(kkk1)) then
                       kkk8=1
                   else
                       kkk8=0
                   endif
               else
                   kkk3=0
                   kkk4=0
                   kkk5=0
                   kkk6=0
                   kkk7=0
                   kkk8=0
                   ttt1=0.0
                   ttt2=0.0
               endif
               write(877,977) i,kkk1,kkk2,kkk3,kkk4,kkk5,kkk6,kkk7,kkk8,ttt1,ttt2, &
                 flow_direction(i,j),elevation(i,j),d8_gradient(i,j),drainage_area(i,j),discharge(i,j)
977 format(9i5,2(g11.4,' '),/,'    ',5(g11.4,' '))               
            enddo
        endif
        999 RETURN
    END ! SUBROUTINE DRAINAGE_BASIN_AREA_FLOW
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE BASIN_REPORT()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        !      *******************************************************************
        !       Debugging report in case of improper agglomeration of drainage
        !       basins
        !      *******************************************************************
        INTEGER :: KKK,I,J,II,JJ
        REAL(4) :: TAREA,MAXAREA
        L934: DO KKK=1,TOTAL_BASINS
            TAREA=0
            MAXAREA=0
            DO  J=1,MY
                DO  I=1,MX
                    IF (BASIN_NUMBER(I,J) == KKK) THEN
                        TAREA=TAREA+CELL_AREA
                        IF (DRAINAGE_AREA(I,J) > MAXAREA) MAXAREA=DRAINAGE_AREA(I,J)
                    ENDIF
                ENDDO
            ENDDO
            IF(ENCLOSED(KKK)) THEN
                WRITE(OUTHIST,935)KKK,DOWNSTREAM_BASIN(KKK),I_OUTFLOW(KKK),J_OUTFLOW(KKK), &
                LAKE_OUTLET_ELEVATION(KKK),TAREA,MAXAREA,BASIN_DRAINAGE_AREA(KKK)
                935         FORMAT(I6,' ENCL., TO:',I5,' AT ',2I6,' ELEVATION=',E15.7,&
                ' AR=',G12.5, &
                ' MA=',G12.5,' LA=',G12.5)
                IF (DOWNSTREAM_BASIN(KKK) == 0) THEN
                    DO  II=1,MX
                        DO  JJ=1,MY
                            IF (BASIN_NUMBER(II,JJ) == KKK) THEN
                                WRITE(OUTHIST,4936) II,JJ
                                4936  FORMAT ('NULL AT I=',I5,' J=',I5)
                            ENDIF
                        ENDDO
                    ENDDO
                ENDIF
            ELSE
                IF(I_OUTFLOW(KKK) < 0) THEN
                    WRITE (OUTHIST,936) KKK
                    936             FORMAT(I6,' DECOMMISSIONED')
                ELSE
                    WRITE (OUTHIST,937) KKK,TAREA,MAXAREA
                    937             FORMAT(I6,' OUTLET, AR=',G12.5,' MA=',G12.5)
                ENDIF
            ENDIF
        ENDDO L934
        WRITE (OUTHIST, 942)
        942 FORMAT(' BASIN NUMBERS',/)
        DO  J=1,MY
            WRITE (OUTHIST,944)  J,(BASIN_NUMBER(I,J),I=1,MX)
            944     FORMAT(I5,10I5,12(/,'     ',10I5))
        ENDDO
        RETURN
    END ! SUBROUTINE BASIN_REPORT
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DRAINAGE_BASIN_LAKE_FLOW()
       USE ERODE_GLOBALS, CHANNEL_EROSION_RATE => CFNE
    !    USE ERODE_GLOBALS
        USE AREA_GLOBALS
        USE LAKE_GLOBALS
        IMPLICIT NONE
        !      *******************************************************************
        !       Calculates drainage areas for each location via flow routing,
        !       including overflow from enclosed depressions
        !
        !      *******************************************************************
        INTEGER :: I,J,LOCAL_DIRECTION,I_OLD,J_OLD,I_NEW,J_NEW,K,L,KK,KKFROM
        INTEGER :: IFROM,JFROM,NUMBBASIN,KKK,III,JJJ
        INTEGER :: IOUTLET,JOUTLET,NCONTIG,KKSQUARE
        INTEGER :: II,JJ,ILL,JLL,IUU,JUU,ILOW,IHIGH,JLOW,JHIGH,KAROUND
        INTEGER :: KKKK,KEND,KKCOUNT,IX,JY,NENC,NFLOWITERS
        INTEGER :: III1,III2,NDIVV,FLOWITERS
        INTEGER :: NENCLOSED,NOVERFLOWS,IOVERFLOWS,ICYCLES,I111,I222,I333
        INTEGER :: ILW,JLW,IHGH,JHGH
        INTEGER :: IWINDOW,JWINDOW
        INTEGER :: KL
        INTEGER :: NBOTTOMBASIN,NBOTTOMADD
        REAL(4) :: FRACTOVERFLOWS
        REAL(4) :: INSUPPLY,NEEDAREA,OLDAREA
        REAL(4) :: ADDAREA
        REAL(4) :: AAMAX,AAMIN,AAAVG,AASD,XXDIST
        REAL(4) :: ELEVLOWEST,MAXIMUMQ,TBEDCONT,TSEDCONT,RIVERADD
        REAL(4) :: THEAVG,THEMAX,THEMIN
        LOGICAL :: CONTIGUOUS,ISFOUND,TOSEDCOVER
        LOGICAL :: TOLAKE,NEWBASIN,REDO,XCONTINUE,YCONTINUE
        REAL(4) :: XCOMP,SEDADD,OLDBASIN_OUTFLUX
        REAL(4) :: DISCHARGE_FROM_CELL,QADD,ADDQ
        REAL(4) :: SUPPLYVOLUME,ROOTARG,RMAX
        REAL(4) :: ELOWER,EUPPER,AAAA,BBBB,LAREA,ELRANGE
        REAL(4) :: TTEMP
        INTEGER (4) :: CURRENT_BASIN
        LOGICAL :: RECYCLES,MAKEISBORDER,FIRSTAREA_DONE
        INTEGER :: KOVERFLOWS,KLINES
        !   MODIFIES: I_OUTFLOW, J_OUTFLOW, ENCLOSED, SEDIMENT_YIELD, FLOW_DIRECTION
        !             DISCHARGE, DRAINAGE_AREA, BASIN_NUMBER,LAKE_AREA,IS_EXIT_BASIN
        !             ROUTED_DRAINAGE_AREA, DOWNSTREAM_BASIN
        !             LAKE_OUTLET_ELEVATION, LAKE_SURFACE_ELEVATION,LAKE_AREA
        !             ERODETOADD,ERGRAVTOADD 
        !   CALLS:  DRAWLINE, DISCHARGE_FROM_CELL, BASIN_REPORT 
        FLOWITERS=20000
        ELEVLOWEST=1.0E+25
        NDIVV=20
        IWINDOW=50
        JWINDOW=50
        FIRSTAREA_DONE=.FALSE.
        JLW=1
        ILW=1
        IHGH=MX
        JHGH=MYY
        AAMAX=-1.0E+25
        ELOWEST=1.0E+25
        AAMIN=-AAMAX
        AAAVG=0.0
        AASD=0.0
        ICYCLES=0
        MAXIMUMQ=-1.0E+25
        RMAX=LMAX*CELL_SIZE*CELL_SIZE
        !      *******************************************************************
        !       no_lakes_exist is index of presence of any enclosed basins
        !      *******************************************************************
        NO_LAKES_EXIST=.TRUE.
        ABORT_SIMULATION=.FALSE.
        NFLOWITERS=0
        !      *******************************************************************
        !       Set initial values - most of these are vectors
        !            current_basin is the current basin number
        !      *******************************************************************
        ITERATION_MAXIMUM_SEDIMENT_YIELD=0.0
        CURRENT_BASIN=0
        DOWNSTREAM_BASIN=0
        I_OUTFLOW=0
        J_OUTFLOW=0
        LAKE_AREA=0.0
        LOCAL_BASIN_DISCHARGE=0.0
        BASIN_OUTFLUX=0.0
        OVERFLOWS=.FALSE.
        ENCLOSED =.FALSE.
        ISBORDER=.FALSE.
        TBEDCONT=0.0
        TSEDCONT=0.0
        SEDIMENT_YIELD=0.0
        ILOWEST=0
        JLOWEST=0
        NBOTTOMBASIN=0
        NCONTIG=0
        NBOTTOMADD=0
        IF (IS_Y_PERIODIC) THEN
		! testing 123
		  DO  J=1,MY
                DO  I=1,MX
                    IF (ELEVATION(I,J) < ELEVLOWEST) THEN
                        ELEVLOWEST=ELEVATION(I,J)
                        ILOWEST=I
                        JLOWEST=J
                    ENDIF
                ENDDO
            ENDDO
            FLOW_DIRECTION(ILOWEST,JLOWEST)=-1
		! end testing 123
            !               flow_direction(ilowest,jlowest)=-1
        ELSE
            DO  I=1,MX
                FLOW_DIRECTION(I,MY)=1
            ENDDO
        ENDIF
        BASIN_NUMBER=0
        DISCHARGE=0.0
        !      *******************************************************************
        !       Main loop for determination of drainage areas
        !      *******************************************************************
        L110: DO J=JLW,JHGH
            M110: DO I=ILW,IHGH
                !      *******************************************************************
                !       tolake is index of whether basin is enclosed
                !       newbasin is index of whether a new basin has been detected
                !         (either a new exit or a new enclosed depression)
                !      *******************************************************************
                TOLAKE=.FALSE.
                NEWBASIN=.FALSE.
                RIVERADD=0.0
                MAKEISBORDER=.FALSE.
                I_OLD=I
                J_OLD=J
                !      *******************************************************************
                !       Increment drainage area for locally contributed drainage
                !      *******************************************************************
                IF (HAVE_INFLUENT_RIVERS.AND.IS_INFLUENT_RIVER_LOCATION(I,J)) THEN
                    DO KL=1,NUMBER_OF_INFLUENT_RIVERS
                        IF((I == I_INFLUENT_RIVER(KL)).AND.(J == J_INFLUENT_RIVER(KL)))THEN
                            RIVERADD=INFLUENT_RIVER_DISCHARGE(KL)
                        ENDIF
                    ENDDO
                ENDIF
                QADD=DISCHARGE_FROM_CELL(I,J)+RIVERADD
                DISCHARGE(I,J)=DISCHARGE(I,J)+QADD
                IF (QADD <= 0.0) THEN
                  !  WRITE(*,7219) I,J
                    7219 FORMAT(' ZERO DISCHARGE, I=',I5,' J=',I5)
                ENDIF
                L120: DO
                    !      *******************************************************************
                    !       Determine local drainage direction and check for bad value
                    !       implicit loop for traversing downstream
                    !      *******************************************************************
                    LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                    !      *******************************************************************
                    !       Check for enclosed basin [point has negative flow_direction(i,j)]
                    !      *******************************************************************
                    IF (LOCAL_DIRECTION < -1) THEN
                        IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) THEN
                        ELSE
                            TOLAKE=.TRUE.
                            NO_LAKES_EXIST=.FALSE.
                        ENDIF
                        !      *******************************************************************
                        !       If it is a depression then check if this is a new basin
                        !      *******************************************************************
                        IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) THEN
                            NEWBASIN=.TRUE.
                        ENDIF
                        !      *******************************************************************
                        !       Branch if this is an enclosed basin or a point on the matrix edge
                        !         if this edge is not periodic
                        !      *******************************************************************
                        EXIT L120
                    ENDIF
                    !      *******************************************************************
                    !       If its not a depression then determine the next downstream location
                    !       also change reference location to the next downstream point and
                    !       increment its drainage area
                    !      *******************************************************************
                    J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                    I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                    IF (IS_X_PERIODIC) THEN
                        IF (I_NEW < 1) I_NEW = MX
                        IF (I_NEW > MX) I_NEW=1
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (J_NEW < 1) J_NEW=MY
                        IF (J_NEW > MY) J_NEW=1
                    ENDIF
                    I_OLD=I_NEW
                    J_OLD=J_NEW
                    DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)+QADD
                    !      *******************************************************************
                    !       Check for exit point if non-periodic lower boundary (j_old=my)
                    !            and also check for a new
                    !            basin, branching if at an exit
                    !      *******************************************************************
                    IF (DO_FLOW_BOUNDARIES) THEN
                        IF (FLOW_DIRECTION(I_OLD,J_OLD) == 1) THEN
                            IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) THEN
                                NEWBASIN=.TRUE.
                                TOLAKE=.FALSE.
                                NBOTTOMBASIN=NBOTTOMBASIN+1
                            ENDIF
                            EXIT L120
                        ENDIF
                    ELSE
                        IF (IS_Y_PERIODIC) THEN
                            IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) THEN
                                IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) NEWBASIN=.TRUE.
                                EXIT L120
                            ENDIF
                        ELSE
                            IF (J_OLD == MY) THEN
                                IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) THEN
                                    NEWBASIN=.TRUE.
                                    TOLAKE=.FALSE.
                                    NBOTTOMBASIN=NBOTTOMBASIN+1
                                    IS_EXIT_BASIN(CURRENT_BASIN+1)=.TRUE. !    temporary
                                ENDIF
                                EXIT L120
                            ENDIF
                        ENDIF
                    ENDIF
                    !      *******************************************************************
                    !       If not aborted and not an exit point, continue downstream,
                    !       incrementing the drainage area abd sediment yield
                    !      *******************************************************************
                    IF (ABORT_SIMULATION) EXIT L120
                ENDDO L120
                !      *******************************************************************
                !       If we are here then we have reached the end of a basin - either
                !       an outlet or depression
                !       First check for aborted run and exit
                !      *******************************************************************
                IF (ABORT_SIMULATION) THEN
                    IABORTMIN=MAX(1,I_OLD-3)
                    IABORTMAX=MIN(I_OLD+3,MX)
                    JABORTMIN=MAX(1,J_OLD-3)
                    JABORTMAX=J_OLD+3
                    WRITE(OUTHIST,814) I_OLD,J_OLD,BASIN_NUMBER(I_OLD,J_OLD)
                    WRITE(*,814) I_OLD,J_OLD,BASIN_NUMBER(I_OLD,J_OLD)
                    814         FORMAT(' *** INFINITE CYCLE AT ',2I6,' ***' &
                    ,/,' AT BASIN NUMBER: ',I5)
                    DO  IX=IABORTMIN,IABORTMAX
                        DO  JY=JABORTMIN,JABORTMAX
                            WRITE(OUTHIST,8222) IX,JY,FLOW_DIRECTION(IX,JY),BASIN_NUMBER(IX,JY) &
                            ,DISCHARGE(IX,JY),D8_GRADIENT(IX,JY),  &
                            ELEVATION(IX,JY)
                        ENDDO
                    ENDDO
                    8222   FORMAT(' IX=',I5,' JY=',I5,' ID=',I5,' BN=',I5,' Q=',G12.5, &
                    ' S=',G12.5,' E=',G12.5,' A=',G12.5)
                    RETURN !    goto 999
                ENDIF
                !      *******************************************************************
                !       If this is a new basin then increment the basin counter and assign
                !       the basinnumb to the basin endpoint
                !      *******************************************************************
                IF (NEWBASIN) THEN
                    CURRENT_BASIN=CURRENT_BASIN+1
                    BASIN_NUMBER(I_OLD,J_OLD)=CURRENT_BASIN
                    IF (TOLAKE) THEN
                        ENCLOSED(CURRENT_BASIN)=.TRUE.
                        !      *******************************************************************
                        !      Look for surrounding points with same elevation - if points have
                        !          the same elevation then we make sure they are included in the same
                        !          basin number
                        !      Initially look only at surrounding 8 points but widen the search
                        !         if a contiguous point at the same elevation is found
                        !      *******************************************************************
                        KAROUND=1
                        IF ((FLOW_DIRECTION(I_OLD,J_OLD) < 2).OR.(D8_GRADIENT(I_OLD,J_OLD) <= 0.0)) THEN
                            L667: DO
                                ISFOUND=.FALSE.
                                IF (IS_X_PERIODIC) THEN
                                    ILOW=I_OLD-KAROUND
                                    IHIGH=I_OLD+KAROUND
                                    IF (ILOW < 1) ILOW=MX+ILOW
                                    IF (IHIGH > MX) IHIGH=IHIGH-MX
                                ELSE
                                    ILOW=MAX(1,I_OLD-KAROUND)
                                    IHIGH=MIN(MX,I_OLD+KAROUND)
                                ENDIF
                                IF (IS_Y_PERIODIC) THEN
                                    JLOW=J_OLD-KAROUND
                                    JHIGH=J_OLD+KAROUND
                                    IF (JLOW < 1) JLOW=MY+JLOW
                                    IF (JHIGH > MY) JHIGH=JHIGH-MY
                                ELSE
                                    JLOW=MAX(1,J_OLD-KAROUND)
                                    JHIGH=MIN(MYY,J_OLD+KAROUND)
                                ENDIF
                                II=ILOW
                                L665: DO
                                    JJ=JLOW
                                    L669: DO
                                        YCONTINUE=.TRUE.
                                        IF ((II == ILOW).OR.(II == IHIGH).OR.    &
                                        (JJ == JLOW).OR.    &
                                        (JJ == JHIGH)) THEN
                                            !      *******************************************************************
                                            !       Ignore the present point and any exit points (jj=my)
                                            !      *******************************************************************
                                            IF ((II == I_OLD).AND.(JJ == J_OLD)) &
                                            YCONTINUE=.FALSE.
                                            IF (IS_Y_PERIODIC) THEN
                                                IF ((II == ILOWEST).AND.(JJ == JLOWEST)) &
                                                YCONTINUE=.FALSE.
                                            ELSE
                                                IF (JJ == MY) YCONTINUE=.FALSE.
                                            ENDIF
                                            IF (YCONTINUE) THEN
                                                IF(ELEVATION(II,JJ) == ELEVATION(I_OLD,J_OLD)) THEN
                                                    !      *******************************************************************
                                                    !       If we are looking at more than 8 points then things are a little
                                                    !       more subtle -- we need to make sure that the point at an equal
                                                    !       elevation is actually contiguous with an already identified basin
                                                    !       point - so we look around the identified point for one already
                                                    !       found
                                                    !      *******************************************************************
                                                    IF(KAROUND > 1) THEN
                                                        KKSQUARE=(KAROUND-1)**2
                                                        IF (IS_X_PERIODIC) THEN
                                                            ILL=II-1
                                                            IF (ILL < 1) ILL=ILL+MX
                                                            IUU=II+1
                                                            IF (IUU > MX) IUU=IUU-MX
                                                        ELSE
                                                            ILL=MAX(II-1,1)
                                                            IUU=MIN(II+1,MX)
                                                        ENDIF
                                                        IF (IS_Y_PERIODIC) THEN
                                                            JLL=JJ-1
                                                            IF (JLL < 1) JLL=JLL+MY
                                                            JUU=JJ+1
                                                            IF (JUU > MY) JUU=JUU-MY
                                                        ELSE
                                                            JLL=MAX(JJ-1,1)
                                                            JUU=MIN(JJ+1,MY)
                                                        ENDIF
                                                        CONTIGUOUS=.FALSE.
                                                        III=ILL
                                                        L664: DO
                                                            JJJ=JLL
                                                            L670: DO
                                                                XCONTINUE=.TRUE.
                                                                IF ((III == II).AND.(JJJ == JJ)) &
                                                                XCONTINUE=.FALSE.
                                                                IF ((((III-II)**2) <= KKSQUARE)  &
                                                                .AND.(((JJJ-JJ)**2) <= KKSQUARE)) &
                                                                XCONTINUE=.FALSE.
                                                                IF (IS_Y_PERIODIC) THEN
                                                                    IF ((II == ILOWEST).AND. &
                                                                    (JJ == JLOWEST)) &
                                                                    XCONTINUE=.FALSE.
                                                                ELSE
                                                                    IF (JJ == MY) &
                                                                    XCONTINUE=.FALSE.
                                                                ENDIF
                                                                IF (XCONTINUE) THEN
                                                                    IF (ELEVATION(III,JJJ) == ELEVATION(II,JJ)) THEN
                                                                        IF (  (((III-II)**2) <= 1).AND.(((JJJ-JJ)**2) <= 1)) THEN
                                                                            CONTIGUOUS=.TRUE.
                                                                        ELSE
                                                                            EXIT L667
                                                                        ENDIF
                                                                    ENDIF
                                                                ENDIF
                                                                IF (JJJ /= JUU) THEN
                                                                    JJJ=JJJ+1
                                                                    IF (IS_Y_PERIODIC) THEN
                                                                        IF (JJJ > MY) JJJ=1
                                                                    ENDIF
                                                                    CYCLE L670
                                                                ENDIF
                                                                EXIT L670
                                                            ENDDO L670
                                                            IF (III == IUU) EXIT L664
                                                            III=III+1
                                                            IF (IS_X_PERIODIC) THEN
                                                                IF (III > MX) III=1
                                                            ENDIF
                                                        ENDDO L664
                                                    ELSE
                                                        CONTIGUOUS=.TRUE.
                                                    ENDIF
                                                    !      *******************************************************************
                                                    !       If the equal-elevation point is contiguous it is given the same
                                                    !       basin number
                                                    !      *******************************************************************
                                                    IF (CONTIGUOUS) THEN
                                                        BASIN_NUMBER(II,JJ)=CURRENT_BASIN
                                                        NCONTIG=NCONTIG+1
                                                        IF ((.NOT.IS_Y_PERIODIC).AND. &
                                                        (JJ == MY-1)) THEN
                                                        ELSE
                                                            ISFOUND=.TRUE.
                                                        ENDIF
                                                    ENDIF
                                                ENDIF
                                            ENDIF
                                        ENDIF
                                        IF (JJ /= JHIGH) THEN
                                            JJ=JJ+1
                                            IF (IS_Y_PERIODIC) THEN
                                                IF (JJ > MY) JJ=1
                                            ENDIF
                                            CYCLE L669
                                        ENDIF
                                        EXIT L669
                                    ENDDO L669
                                    IF (II == IHIGH) EXIT L665
                                    II=II+1
                                    IF (IS_X_PERIODIC) THEN
                                        IF (II > MX) II=1
                                    ENDIF
                                ENDDO L665
                                KAROUND=KAROUND+1
                                IF (KAROUND > NMAX) THEN
                                    RETURN
                                ENDIF
                                !      *******************************************************************
                                !        We widen the search if a contiguous point at the same elevation
                                !        has been found to include points possibly contiguous to the new
                                !        point but not the original point
                                !      *******************************************************************
                                IF (.NOT.ISFOUND) EXIT
                            ENDDO L667
                        ENDIF
                        !      *******************************************************************
                        !        End of code relating to grouping surrounding points at the same
                        !            elevation
                        !      *******************************************************************
                    ENDIF
                ENDIF
                !      *******************************************************************
                !       We need to retrace our steps from the initial location, setting the
                !       basin number of all points along the path to that of the exit or sink point
                !      *******************************************************************
                NUMBBASIN=BASIN_NUMBER(I_OLD,J_OLD)
                I_OLD=I
                J_OLD=J
                L140: DO
                    !      *******************************************************************
                    !       End if we reach an already-identified basin location
                    !      *******************************************************************
                    IF (BASIN_NUMBER(I_OLD,J_OLD) /= 0) THEN
                        EXIT
                    ELSE
                        BASIN_NUMBER(I_OLD,J_OLD)=NUMBBASIN
                        LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                        IF (LOCAL_DIRECTION < 2) EXIT
                        J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                        I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                        IF (DO_FLOW_BOUNDARIES) THEN
                            IF (FLOW_DIRECTION(I_NEW,J_NEW) == 1) THEN
                                EXIT
                            ENDIF
                        ELSE
                            IF (IS_X_PERIODIC) THEN
                                IF (I_NEW < 1) I_NEW = MX
                                IF (I_NEW > MX) I_NEW=1
                            ENDIF
                            IF (IS_Y_PERIODIC) THEN
                                IF (J_NEW < 1) J_NEW=MY
                                IF (J_NEW > MY) J_NEW=1
                            ENDIF
                        ENDIF
                        I_OLD=I_NEW
                        J_OLD=J_NEW
                    ENDIF
                ENDDO L140
            ENDDO M110
        ENDDO L110
        !      *******************************************************************
        !       The initial drainage areas and basin numbers have now been determined
        !       check to make sure that we have not exceeded the dimensioned
        !       number of basins
        !       Also we have determined the sediment yield from channel erosion
        !      *******************************************************************
        IF (CURRENT_BASIN > LMMX) THEN
            WRITE(*,933) CURRENT_BASIN
            933     FORMAT(' TOO MANY BASINS - REDIMENSION.  CURRENT_BASIN=',I6)
            RETURN
        ENDIF
        !      *******************************************************************
        !       Look at boundary points and set values for those that have not
        !       received drainage from upstream to unity area and assign a basin number
        !       Do this first for any 'sinkholes' and then for all locations
        !      *******************************************************************
        IF (DO_FLOW_BOUNDARIES) THEN
            DO J=1,MY
                DO I=1,MX
                    IF (FLOW_DIRECTION(I,J) == 1) THEN
                        IF (DISCHARGE(I,J) == 0.0) THEN
                            DISCHARGE(I,J)=DISCHARGE_FROM_CELL(I,J)
                        ENDIF
                        IF (BASIN_NUMBER(I,J) == 0) THEN
                            CURRENT_BASIN=CURRENT_BASIN+1
                            BASIN_NUMBER(I,MY)=CURRENT_BASIN
                            NBOTTOMADD=NBOTTOMADD+1
                            IS_EXIT_BASIN(CURRENT_BASIN)=.TRUE.
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
        ELSE
            IF (IS_Y_PERIODIC) THEN
            ELSE
                DO  I=1,MX
                    IF (DISCHARGE(I,MY) == 0) THEN
                        DISCHARGE(I,MY)=DISCHARGE_FROM_CELL(I,MY)
                    ENDIF
                    IF (BASIN_NUMBER(I,MY) == 0) THEN
                        CURRENT_BASIN=CURRENT_BASIN+1
                        BASIN_NUMBER(I,MY)=CURRENT_BASIN
                        NBOTTOMADD=NBOTTOMADD+1
                        IS_EXIT_BASIN(CURRENT_BASIN)=.TRUE.
                    ENDIF
                ENDDO
            ENDIF
        ENDIF
        DO J=1,MY
            DO I=1,MX
                IF (BASIN_NUMBER(I,J) == 0) THEN
                    CURRENT_BASIN=CURRENT_BASIN+1
                    BASIN_NUMBER(I,J)=CURRENT_BASIN
                ENDIF
            ENDDO
        ENDDO
        TOTAL_BASINS=CURRENT_BASIN
        FIRSTAREA_DONE=.TRUE.
        LAKE_SURFACE_ELEVATION=-1.0E25
        ROUTED_DISCHARGE=DISCHARGE
        I333=0
        DO I=1,TOTAL_BASINS
            IF (ENCLOSED(I)) I333=I333+1
        ENDDO
        !    5994  continue
        I111=0
        IF (COMPLETE_RUNOFF) I111=1
        I222=0
        IF (NO_LAKES_EXIST) I222=1
        !      *******************************************************************
        !       If there are enclosed basins (*and if runoff overflow is permitted*)
        !       then we identify those enclosed basins which need to be examined
        !       for grouping and flow routing to the outlets - Skip all this
        !       otherwise .  We first set up initial vector values based upon whether
        !       the basin is enclosed or not
        !      *******************************************************************
        IF ((COMPLETE_RUNOFF).AND.(.NOT.NO_LAKES_EXIST)) THEN
            491   FORMAT (' AGGLOMERATING BASINS')
            NENC=0
            DO K=1,TOTAL_BASINS
                IF (ENCLOSED(K).AND.(.NOT.ISBORDER(K))) THEN
                    OVERFLOWS(K)=.FALSE.
                    NEXTCYCLE(K)=.FALSE.
                    NEEDTODO(K)=.TRUE.
                    NEWGEOMETRY(K)=.TRUE.
                    LOCAL_BASIN_DISCHARGE(K)=0.0
                    BASIN_OUTFLUX(K)=0.0
                    NENC=NENC+1
                    LAKE_SURFACE_ELEVATION(K)=-1.0E25
                ELSE
                    NEEDTODO(K)=.FALSE.
                    OVERFLOWS(K)=.FALSE.
                    NEWGEOMETRY(K)=.FALSE.
                    NEXTCYCLE(K)=.FALSE.
                    LOCAL_BASIN_DISCHARGE(K)=0.0
                    BASIN_OUTFLUX(K)=0.0
                ENDIF
            ENDDO
			IF (IS_Y_PERIODIC) THEN
				K=BASIN_NUMBER(ILOWEST,JLOWEST)
			ENDIF
            WRITE(*,937) TOTAL_BASINS,NENC
            937   FORMAT(' TOTAL_BASINS=',I8,' NO. ENCLOSED=',I8)
            !      *******************************************************************
            !       We need to keep regrouping basins and identifying lake outlets and
            !       cumulative drainage areas in the case where there is enclosed basins
            !       and lake overflow - redo is an index of whether another level of
            !       basin consolidation is needed after this one
            !      *******************************************************************
            L250: DO
                REDO=.FALSE.
                ICYCLES=ICYCLES+1
                !      *******************************************************************
                !       If we still need to work in a basin (an outlet location has still
                !       not been determined) then we reset the outlet elevation and
                !       lake discharge to just the discharge within the basin
                !      *******************************************************************
                DO K=1,TOTAL_BASINS
                    IF (NEEDTODO(K)) LOCAL_BASIN_DISCHARGE(K)=0.0
                    OLDOVERFLOWS(K)=BASIN_OUTFLUX(K)
                ENDDO
                DO J=JLW,JHGH
                    DO  I=ILW,IHGH
                        K=BASIN_NUMBER(I,J)
                        IF (NEEDTODO(K)) THEN
                            LOCAL_BASIN_DISCHARGE(K)=LOCAL_BASIN_DISCHARGE(K)+DISCHARGE_FROM_CELL(I,J)
                        ENDIF
                    ENDDO
                ENDDO
                DO K=1,TOTAL_BASINS
                    IF (NEWGEOMETRY(K)) THEN
                        LAKE_OUTLET_ELEVATION(K)=1.0E+25
                        !      *******************************************************************
                        !         Here we reset the outlet elevation of the basin with new geometry
                        !          to a default value indicating that it is an unknown value
                        !      *******************************************************************
                        LAKE_AREA(K)=0.0
                    ENDIF
                ENDDO
                L160: DO J=JLW,JHGH
                    M160: DO I=ILW,IHGH
                        K=BASIN_NUMBER(I,J)
                        IF (NEWGEOMETRY(K)) THEN
                            !      *******************************************************************
                            !       We need to look at all points on the periphery of the basin to
                            !       determine if it is an outlet, and we identify the lowest peripheral
                            !       location as the outlet
                            !      *******************************************************************
                            L170: DO  L=2,9
                                J_NEW=J+DOWNSTREAM(L,2)
                                IF (IS_Y_PERIODIC) THEN
                                    IF (J_NEW < 1) J_NEW = MY
                                    IF (J_NEW > MY) J_NEW=1
                                ELSE
                                    IF (J_NEW < 1) CYCLE
                                    IF (J_NEW > MY) CYCLE
                                ENDIF
                                I_NEW=I+DOWNSTREAM(L,1)
                                IF (IS_X_PERIODIC) THEN
                                    IF (I_NEW < 1) I_NEW = MX
                                    IF (I_NEW > MX) I_NEW=1
                                ELSE
                                    IF (I_NEW < 1) CYCLE
                                    IF (I_NEW > MX) CYCLE
                                ENDIF
                                !      *******************************************************************
                                !       Do the following if the present point is at a peripheral location
                                !       check to see if the periphery point or a contiguous point outside
                                !       the basin is lower than any identified lake outlet
                                !      *******************************************************************
                                IF ((BASIN_NUMBER(I_NEW,J_NEW) /= K)) &
                                THEN
                                    XCOMP=MAX(ELEVATION(I_NEW,J_NEW),ELEVATION(I,J))
                                    !      *******************************************************************
                                    !       If this is a potential lake outlet then there are two cases:
                                    !         1.  The basin peripery is lower than any contiguous point outside
                                    !             the basin.  If this is the case (outer is true) the outlet is
                                    !             identified with the outside point and the outlet elevation is
                                    !             equal to that of the outside point.
                                    !         2.  The basin periphery is higher than the lowest contiguous point
                                    !             outside.  The outlet is the periphery point and the lake outlet
                                    !             elevation is its elevation.  In this case we need to record where
                                    !             the flow is going to and its direction (ifrom,jfrom,kfrom)
                                    !      *******************************************************************
                                    IF (XCOMP < LAKE_OUTLET_ELEVATION(K)) THEN
                                        IF (ELEVATION(I,J) <= ELEVATION(I_NEW,J_NEW)) THEN
                                            OUTER(K)=.TRUE.
                                            I_OUTFLOW(K)=I_NEW
                                            J_OUTFLOW(K)=J_NEW
                                            LAKE_OUTLET_ELEVATION(K)=ELEVATION(I_NEW,J_NEW)
                                        ELSE
                                            OUTER(K)=.FALSE.
                                            I_OUTFLOW(K)=I
                                            J_OUTFLOW(K)=J
                                            IOUTLET=I
                                            JOUTLET=J
                                            LAKE_OUTLET_ELEVATION(K)=ELEVATION(I,J)
                                            IFROM=I_NEW
                                            JFROM=J_NEW
                                            KKFROM=L
                                            if (flow_direction(ioutlet,joutlet)/=1) FLOW_DIRECTION(IOUTLET,JOUTLET)=KKFROM
                                            XXDIST=CELL_SIZE*SQRT(FLOAT((JFROM-JOUTLET)**2 &
                                            +(IOUTLET-IFROM)**2))
                                            D8_GRADIENT(IOUTLET,JOUTLET)=(ELEVATION(IOUTLET,JOUTLET) &
                                            -ELEVATION(IFROM,JFROM))/(XXDIST)
                                        ENDIF
                                        !      *******************************************************************
                                        !       Downstream_basin is the basin into which the present basin drains
                                        !      *******************************************************************
                                        DOWNSTREAM_BASIN(K)=BASIN_NUMBER(I_NEW,J_NEW)
                                    ENDIF
                                ENDIF
                            ENDDO L170
                        ENDIF
                    ENDDO M160
                ENDDO L160
                !      *******************************************************************
                !       We have now identified the lowest outlet for the present basin
                !       so we set its outlet location.  If the basin outlet is on the
                !       periphery (not outer) then we need to redefine the drainage direction
                !       and gradient of the peripery point for flow to the outlet rather
                !       than into the present basin
                !      *******************************************************************
                !      *******************************************************************
                !       We have now done this basin - we don't need to worry about it for
                !       the time being.
                !      *******************************************************************
                DO K=1,TOTAL_BASINS
                    IF (ENCLOSED(K)) THEN
                        IF (NEWGEOMETRY(K)) THEN
                            NEWGEOMETRY(K)=.FALSE.
                            LAKE_VOLUME(K)=0.0
                            LAKE_AREA(K)=0.0
                            LOWEST_LAKE_ELEVATION(K)=1.0E+25
                            DO  J=JLW,JHGH
                                DO  I=ILW,IHGH
                                    IF (BASIN_NUMBER(I,J) == K) THEN
                                        IF ( ELEVATION(I,J) <= LAKE_OUTLET_ELEVATION(K) ) THEN
                                            LAKE_AREA(K)=LAKE_AREA(K)+CELL_SIZE*CELL_SIZE
                                            LAKE_VOLUME(K)=LAKE_VOLUME(K)+(LAKE_OUTLET_ELEVATION(K)-ELEVATION(I,J)) &
                                            *CELL_SIZE*CELL_SIZE
                                            IF ( ELEVATION(I,J) < LOWEST_LAKE_ELEVATION(K) ) LOWEST_LAKE_ELEVATION(K)=ELEVATION(I,J)
                                        ENDIF
                                    ENDIF
                                ENDDO
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
                !      ******************************************************************
                !       If we have a basin which is indicated as being enclosed but which
                !        does not actually drain into another basin, then we reset its status
                !        to not being enclosed and reset a number of indicator variables
                !      ******************************************************************
                DO K=1,TOTAL_BASINS
                    IF (ENCLOSED(K)) THEN
                        IF(DOWNSTREAM_BASIN(K) == 0) THEN
                            write(*,1511) K
                            1511     FORMAT(' NO DOWNSTREAM_BASIN, K=',I8)
                            ENCLOSED(K)=.FALSE.
                            NEEDTODO(K)=.FALSE.
                            NEWGEOMETRY(K)=.FALSE.
                            NEXTCYCLE(K)=.FALSE.
                        ENDIF
                    ENDIF
                ENDDO
                NENCLOSED=0
                NOVERFLOWS=0
                DO K=1,TOTAL_BASINS
                    BASIN_INFLUX(K)=0.0
                    IF (ENCLOSED(K)) THEN
                        IF (NEEDTODO(K)) THEN
                            DO KK=1,TOTAL_BASINS
                                IF (K /= KK) THEN
                                    IF ((DOWNSTREAM_BASIN(KK) == K).AND.OVERFLOWS(KK)) THEN
                                        BASIN_INFLUX(K)=BASIN_INFLUX(K)+BASIN_OUTFLUX(KK)
                                    ENDIF
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDIF
                ENDDO
                !      ******************************************************************
                !        Cycle through all enclosed basins.  Reset the needtodo index.
                !         then calculate the total inflow to each basin from all other
                !         basins (the variable basin_influx)
                !      ******************************************************************
                L1150: DO K=1,TOTAL_BASINS
                    IF ( ENCLOSED(K) ) THEN
                        NENCLOSED=NENCLOSED+1
                        IF (NEEDTODO(K)) THEN
                            NEEDTODO(K)=.FALSE.
                            OLDBASIN_OUTFLUX=BASIN_OUTFLUX(K)
                            NEW_BASIN_OUTFLUX(K)=0.0
                            !      ******************************************************************
                            !        Initially for each basin we assume that there is no overflow into
                            !          adjacent basins, so that overflows is false and basin_outflux is zero
                            !        The inflow into the basin (insupply) is the locally-derived drainage (local_basin_discharge)
                            !          plus the total inflow from other basins
                            !      ******************************************************************
                            OVERFLOWS(K)=.FALSE.
                            BASIN_OUTFLUX(K)=0.0
                            IOVERFLOWS=0
                            INSUPPLY=BASIN_INFLUX(K)+LOCAL_BASIN_DISCHARGE(K)
                            SUPPLYVOLUME=INSUPPLY-LAKE_AREA(K)*EVAPORATION_RATE
                            !      ******************************************************************
                            !       The available watersupply to the lake is the amount of inflow
                            !         minus the maximum possible loss from evaporation (the evaporation rate
                            !         times the surface
                            !         area of the maximum lake_surface_elevation.
                            !       Then we check to make sure that the basin has a finite storage capacity
                            !        (positive lake_volume).
                            !       If the supplyvolume is positive, then we have overflow of that magnitude
                            !
                            !      ******************************************************************
                            IF ( LAKE_VOLUME(K) > 0.0 ) THEN
                                IF ( SUPPLYVOLUME >  0.0) THEN
                                    OVERFLOWS(K)=.TRUE.
                                    NEW_BASIN_OUTFLUX(K)=SUPPLYVOLUME
                                    LAKE_SURFACE_ELEVATION(K)=LAKE_OUTLET_ELEVATION(K)
                                    NOVERFLOWS=NOVERFLOWS+1
                                    IOVERFLOWS=1
                                    !      *****************************************************************
                                    !        Then we see if the new volume of overflow equals the previous value
                                    !        if it is not, then the basin to which the flow is going has to be
                                    !        reevaluated again iteratively.  We indicate this by setting
                                    !        nextcycle to true for the basin receiving the overflow.
                                    !      ****************************************************************
                                    IF (NEW_BASIN_OUTFLUX(K) /= OLDBASIN_OUTFLUX) THEN
                                        IF (ENCLOSED(DOWNSTREAM_BASIN(K))) NEXTCYCLE(DOWNSTREAM_BASIN(K))=.TRUE.
                                    ENDIF
                                ELSE
                                    !      *********************************************************************
                                    !        If supplyvolume is negative or zero, then the basin is only partially
                                    !          filled with water and does not overflow.  We need to calculate the
                                    !          actual lake elevation (lake_surface_elevation). rootarg is the actual supply of
                                    !          water divided by the maximum amount of evaporation for an overflowing
                                    !          lake.  The actual lake level is initially assumed to be linearly
                                    !          related to rootarg, such that if rootarg is zero the lake_surface_elevation
                                    !          would equal the lowest elevation in the basin, and if it is unity
                                    !          the lake just at the verge of overflowing the lake_surface_elevation would equal
                                    !          the overflow level
                                    !      **********************************************************************
                                    ROOTARG=INSUPPLY/(LAKE_AREA(K)*EVAPORATION_RATE)
                                    IF (ROOTARG > 0.0) THEN
                                        LAKE_SURFACE_ELEVATION(K)=LOWEST_LAKE_ELEVATION(K)+  &
                                        (LAKE_OUTLET_ELEVATION(K)-LOWEST_LAKE_ELEVATION(K))* &
                                        ROOTARG
                                    ELSE
                                        LAKE_SURFACE_ELEVATION(K)=LOWEST_LAKE_ELEVATION(K)
                                    ENDIF
                                ENDIF
                            ELSE
                                !      *********************************************************************
                                !        This is the presumably rare case where the basin has no storage
                                !         capacity.  If there is a positive water input, it is just passed
                                !         downstream
                                !      *********************************************************************
                                IF (SUPPLYVOLUME > 0.0) THEN
                                    OVERFLOWS(K)=.TRUE.
                                    IOVERFLOWS=1
                                    LAKE_SURFACE_ELEVATION(K)=LAKE_OUTLET_ELEVATION(K)
                                    NOVERFLOWS=NOVERFLOWS+1
                                    NEW_BASIN_OUTFLUX(K)=BASIN_INFLUX(K)+LOCAL_BASIN_DISCHARGE(K)
                                    IF (NEW_BASIN_OUTFLUX(K) /= OLDBASIN_OUTFLUX) THEN
                                        IF (ENCLOSED(DOWNSTREAM_BASIN(K))) NEXTCYCLE(DOWNSTREAM_BASIN(K))=.TRUE.
                                    ENDIF
                                ELSE
                                    OVERFLOWS(K)=.FALSE.
                                    NEW_BASIN_OUTFLUX(K)=0.0
                                    BASIN_OUTFLUX(K)=0.0
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO L1150
                !      *******************************************************************
                !        Here we set the variable basin_outflux to equal the variable new_basin_outflux
                !         for those basins with overflow
                !      *******************************************************************
                DO K=1,TOTAL_BASINS
                    IF (OVERFLOWS(K)) THEN
                        BASIN_OUTFLUX(K)=NEW_BASIN_OUTFLUX(K)
                    ELSE
                        BASIN_OUTFLUX(K)=0.0
                    ENDIF
                ENDDO
                IF ( NENCLOSED > 0 ) THEN
                    FRACTOVERFLOWS=FLOAT(NOVERFLOWS)/FLOAT(NENCLOSED)
                ELSE
                    FRACTOVERFLOWS=0.0
                ENDIF
                !      *******************************************************************
                !       We have now completed one stage of basin agglomeration.  We need
                !       to check for the case where lake basins drain into each other
                !       if this is true then redo is set true and needtodo(k) is set true
                !       redo is false if each chain of enclosed basin outlets following downstream
                !       eventually leads to an external outlet
                !      *******************************************************************
                K=0
                KKK=0
                !      *******************************************************************
                !         Cycle through all basins
                !      *******************************************************************
                L260: DO
                    K=K+1
                    KKCOUNT=0
                    !      *******************************************************************
                    !         Exit from loop if we have looked at all basins
                    !      *******************************************************************
                    IF (K >= TOTAL_BASINS) EXIT L260 
                    !      *******************************************************************
                    !         If we find an enclosed basin we need to check that tracing successive
                    !         outlets downstream eventually leads to an external outlet, otherwise
                    !         we need an additional stage of basin agglomeration
                    !         If the newly calculated amount of flow leaving the basin does not
                    !           equal the previously calculated value (oldoverflows), and if the next basin
                    !           downstream was not overflowing, then we need to recalculate the
                    !           water balance for the downstream basin.
                    !      *******************************************************************
                    IF (OVERFLOWS(K)) THEN
                        !
                        KK=K
                        IF (OLDOVERFLOWS(K) /= BASIN_OUTFLUX(K)) THEN
                            IF (.NOT.OVERFLOWS(DOWNSTREAM_BASIN(K))) THEN
                                NEXTCYCLE(DOWNSTREAM_BASIN(K))=.TRUE.
                            ENDIF
                        ENDIF
                        !      *******************************************************************
                        !         In this loop we go downstream through successive basins receiving
                        !          drainage from the present basin.
                        !       ******************************************************************
                        L623: DO
                            !      *******************************************************************
                            !         kkk is the next basin downstream from where we are now (basin kk)
                            !      ********************************************************************
                            KKK=DOWNSTREAM_BASIN(KK)
                            !      *******************************************************************
                            !         If the next basin flows into the basin we started at (k) and if that
                            !           next basin is overflowing, we can stop our downstream march --
                            !         We have then to deal with the fact that the two basins drain into
                            !           one another.
                            !      *******************************************************************
                            IF ((DOWNSTREAM_BASIN(KKK) == K).AND.OVERFLOWS(KKK)) EXIT L623
                            !      *********************************************************************
                            !         If the next basin is not overflowing, then we can go on to the next
                            !            basin in our list.
                            !         We can also go on if the basin if we get to a basin with a 'sinkhole'
                            !           so that it cannot fill up.
                            !      **********************************************************************
                            IF (.NOT.OVERFLOWS(KKK)) CYCLE L260
                            IF (IS_Y_PERIODIC) THEN
                                IF ((I_OUTFLOW(KKK) == ILOWEST).AND.(J_OUTFLOW(KKK) == JLOWEST)) CYCLE L260
                            ELSE
                                IF (COMPLETE_RUNOFF) THEN
                                    IF ((I_OUTFLOW(KKK) == ILOWEST).AND.(J_OUTFLOW(KKK) == JLOWEST)) CYCLE L260
                                ENDIF
                            ENDIF
                            !      *********************************************************************
                            !        If we have gotten here, we need to go on to the next basin downstream
                            !      *********************************************************************
                            KK=KKK
                            KKCOUNT=KKCOUNT+1
                            IF (KKCOUNT > 2000) CYCLE L260
                        ENDDO L623
                        !      *******************************************************************
                        !       his is the case of mutually-draining basins
                        !      *******************************************************************
                        KK=K
                        KEND=KKK
                        !      *******************************************************************
                        !      Wwe set newgeometry to true for our basin because we have to combine
                        !        our basin (k) with the basin that drains into it (kend).  We reset
                        !        some of the hydrology variables for our present basin (k).
                        !      *******************************************************************
                        NEWGEOMETRY(K)=.TRUE.
                        OVERFLOWS(K)=.FALSE.
                        BASIN_OUTFLUX(K)=0.0
                        NEEDTODO(K)=.TRUE.
                        KKK=DOWNSTREAM_BASIN(KK)
                        !      *******************************************************************
                        !       If there are two mutually-draining basins, then they are arbitrarily
                        !       agglomerated by reassigning one basins points (kkk) to the other (k)
                        !       We cycle through each basin in the loop of mutually draining basins
                        !       and agglomerate then into basin (k).
                        !      *******************************************************************
                        DO  J=1,MYY
                            DO  I=1,MX
                                IF (BASIN_NUMBER(I,J) == KKK) BASIN_NUMBER(I,J)=K
                            ENDDO
                        ENDDO
                        !      *******************************************************************
                        !       If other basins drain to the reassigned drainage basin, then their
                        !       exit basin have to be reassined to the agglomerated basin
                        !      *******************************************************************
                        DO  KKKK=1,TOTAL_BASINS
                            IF (ENCLOSED(KKKK)) THEN
                                IF (KKKK == K) CYCLE
                                IF (KKKK == KKK) CYCLE
                                IF (DOWNSTREAM_BASIN(KKKK) == KKK) DOWNSTREAM_BASIN(KKKK)=K
                            ENDIF
                        ENDDO
                        !      *******************************************************************
                        !       We need to deactivate the various values of the deactivated
                        !       (agglomerated) basin (kkk)
                        !      *******************************************************************
                        ENCLOSED(KKK)=.FALSE.
                        OVERFLOWS(KKK)=.FALSE.
                        BASIN_OUTFLUX(KKK)=0.0
                        LAKE_AREA(KKK)=0.0
                        LOCAL_BASIN_DISCHARGE(KKK)=0.0
                        I_OUTFLOW(KKK)=-1
                        J_OUTFLOW(KKK)=-1
                        KK=DOWNSTREAM_BASIN(KKK)
                        DOWNSTREAM_BASIN(KKK)=0
                        LAKE_VOLUME(KKK)=0.0
                        LOWEST_LAKE_ELEVATION(KKK)=0.0
                        !      *******************************************************************
                        !       Look at additional basins if we have not examined all of them
                        !      *******************************************************************
                        IF (KKK == KEND) CYCLE L260
                        KKK=KK
                    ENDIF
                    !      *******************************************************************
                    !       Go on to the next basin
                    !      *******************************************************************
                ENDDO L260
                !      *******************************************************************
                !       If we have not completely agglomerated drainage basins so that
                !        the overflows and discharges have not changed during this iterative
                !        cycle, then we neet to indicate which enclosed basins need to be reexamined
                !        That is, if newgeometry or nextcycle is true.
                !       If any basin needs to be reexamined, then we set redo to true.
                !      *******************************************************************
                RECYCLES=.FALSE.
                REDO=.FALSE.
                DO K=1,TOTAL_BASINS
                    IF (ENCLOSED(K)) THEN
                        IF (NEWGEOMETRY(K).OR.NEXTCYCLE(K)) THEN
                            NEEDTODO(K)=.TRUE.
                            NEXTCYCLE(K)=.FALSE.
                            REDO=.TRUE.
                        ELSE
                            NEEDTODO(K)=.FALSE.
                        ENDIF
                    ELSE
                        NEEDTODO(K)=.FALSE.
                    ENDIF
                ENDDO
                !      *********************************************************************
                !       If (groan!!!) we have not yet reached a stable hydrology, we go back
                !         and do it again, at least until we reach our specified maximum number
                !         of recalculations (flowiters)
                !      *********************************************************************
                IF (REDO) THEN
                    NFLOWITERS=NFLOWITERS+1
                    IF (NFLOWITERS >= FLOWITERS) then 
                       write(*,433) nflowiters,flowiters
433                    format('too many flow iterations, nflowiters=',i6,' max_flowiters=',i6)
                      
                       EXIT L250
                     endif
                ELSE
                    WRITE(*,373) NFLOWITERS
                    WRITE(OUTHIST,373) NFLOWITERS
                    373     FORMAT(' NO MORE REDO, FLOW ITERATIONS=',I7)
                    EXIT L250
                ENDIF
            ENDDO L250
            !      *******************************************************************
            !       Now that we have calculated contributions from upstream enclosed basins
            !       we need to calculate the total lake drainage area of each basin
            !      *******************************************************************
            !      *******************************************************************
            !       We now go downstream from each enclosed basin outlet and increment the
            !       drainage area by the amount of drainage flowing through the outlet
            !      *******************************************************************
            KOVERFLOWS=0
            KLINES=0
            L180: DO K=1,TOTAL_BASINS
                IF (OVERFLOWS(K)) THEN
                    KOVERFLOWS=KOVERFLOWS+1
                    I_OLD=I_OUTFLOW(K)
                    J_OLD=J_OUTFLOW(K)
                    !      *******************************************************************
                    !       In the rare to non-existent case where the outlet of our basin is
                    !        also the outlet of the next basin downstream, then we don't need to
                    !        do anything here.
                    !     ********************************************************************
                    IF ((I_OLD == I_OUTFLOW(DOWNSTREAM_BASIN(K))).AND.  &
                    (J_OLD == J_OUTFLOW(DOWNSTREAM_BASIN(K)))) THEN
                        !WRITE(OUTHIST,28443) I_OLD,J_OLD,K,DOWNSTREAM_BASIN(K)
                        28443 FORMAT('SAME OUTLET, I_OLD=',I5,' J_OLD=',I5,' K=',I8,' TB=',I8)
                        CALL CHECK_FLOW_PATH(K)
                        CYCLE L180
                    ENDIF
                    !      *******************************************************************
                    !        This just checks to make sure that nothing terrible has gone
                    !         wrong in our calculations
                    !      *******************************************************************
                    IF ((I_OLD < 1).OR.(I_OLD > MX).OR. &
                    (J_OLD < 1).OR.(J_OLD > MY)) THEN
                        WRITE(OUTHIST,377) K,I_OLD,J_OLD
                        377             FORMAT(' ***BASIN INTEGRITY PROBLEM AT BASIN' &
                        ,I5,', I_OLD & J_OLD=', &
                        2I5)
                    ENDIF
                    !      ******************************************************************
                    !        Basically what is happening here is that we are accumulating discharges
                    !         downstream.
                    !        Within each basin we add the total cumulative discharge from upstream until
                    !          we get to a depression or sinkhole in that basin.  If the location in the
                    !          basin is below the lake_surface_elevation, then the discharge is really imaginary, but
                    !          we calculate it anyway.
                    !      ********************************************************************
                    ADDQ=BASIN_OUTFLUX(K)
                    IF(OUTER(K)) THEN
                        DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)+ADDQ
                        ROUTED_DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)
                    ELSE
                        DISCHARGE(I_OLD,J_OLD)=ADDQ
                        ROUTED_DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)
                    ENDIF
                    IF (DO_FLOW_BOUNDARIES) THEN
                        IF (FLOW_DIRECTION(I_OLD,J_OLD) == 1) THEN
                            CYCLE L180
                        ENDIF
                    ELSE
                        IF (IS_Y_PERIODIC) THEN
                            IF ((I_OLD == ILOWEST).AND. &
                            (J_OLD == JLOWEST)) CYCLE L180
                        ELSE
                            IF (J_OLD == MY) CYCLE L180
                        ENDIF
                    ENDIF
                    LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                    IF (LOCAL_DIRECTION >= 0) THEN
                        IF ((I_OLD /= ILOWEST).AND.(J_OLD /= JLOWEST)) THEN
                            IF (LOCAL_DIRECTION >= 0) THEN !            exit l190 !goto 200
                                L190: DO
                                    J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                                    I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                                    IF (IS_X_PERIODIC) THEN
                                        IF (I_NEW < 1) I_NEW = MX
                                        IF (I_NEW > MX) I_NEW=1
                                    ENDIF
                                    IF (IS_Y_PERIODIC) THEN
                                        IF (J_NEW < 1) J_NEW=MY
                                        IF (J_NEW > MY) J_NEW=1
                                    ENDIF
                                    I_OLD=I_NEW
                                    J_OLD=J_NEW
                                    IF (BASIN_NUMBER(I_OLD,J_OLD) /= DOWNSTREAM_BASIN(K)) THEN
                                        !WRITE(OUTHIST,823) K,DOWNSTREAM_BASIN(K), &
                                        !BASIN_NUMBER(I_OLD,J_OLD),I_OLD,J_OLD  &
                                        !,DISCHARGE(I_OLD,J_OLD),ADDQ
                                        823               FORMAT(' CHANGE BASIN NUMBER',I5, &
                                        ' TO BASIN'  &
                                        ,I5,/,' NOW IN BASIN',I5,' AT I=',I5,' J=',I5, &
                                        ' AREA=',G12.5,' ADDQ=',G12.5)
                                    ENDIF
                                    !      *******************************************************************
                                    !       In the rare to non-existent case where the outlet of our basin is
                                    !        also the outlet of the next basin downstream, then we don't need to
                                    !        do anything here.
                                    !     ********************************************************************
                                    IF ((I_OLD == I_OUTFLOW(DOWNSTREAM_BASIN(K)))&
                                    .AND.(J_OLD == J_OUTFLOW(DOWNSTREAM_BASIN(K)))) THEN
                                        !WRITE(OUTHIST,28442) I_OLD,J_OLD,K,DOWNSTREAM_BASIN(K)
                                        CALL CHECK_FLOW_PATH(K)
                                        28442 FORMAT(' SAME OUTLET, I_OLD=',I5,' J_OLD=',I5,' K=',I8,' TB=',I8)
                                        EXIT L190
                                    ENDIF
                                    DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)+ADDQ
                                    ROUTED_DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)
                                    LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                                    !      *******************************************************************
                                    !       Again, look for infinite loops and abort if they occur
                                    !      *******************************************************************
                                    IF(DISCHARGE(I_OLD,J_OLD) > 1.0E+20) THEN
                                        ABORT_SIMULATION=.TRUE.
                                        IABORTMIN=I_OLD-3
                                        IABORTMAX=I_OLD+3
                                        JABORTMIN=J_OLD-3
                                        JABORTMAX=J_OLD+3
                                        WRITE(OUTHIST,378) I_OLD,J_OLD,BASIN_NUMBER(I_OLD,J_OLD)
                                        378             FORMAT (' *** INFINITE CYCLE IN OUTLET COMBINE AT '&
                                        ,2I6,' ***' &
                                        ,/,' AT BASIN NUMBER: ',I5)
                                        ABORT_SIMULATION=.TRUE.
                                        RETURN
                                    ENDIF
                                    IF ((I_OLD == I_OUTFLOW(DOWNSTREAM_BASIN(K))).AND.(J_OLD == J_OUTFLOW(DOWNSTREAM_BASIN(K)))) &
                                    EXIT L190
                                    IF (LOCAL_DIRECTION < 0) EXIT L190
                                    IF (DO_FLOW_BOUNDARIES) THEN
                                        IF (LOCAL_DIRECTION == 1) THEN
                                            EXIT L190
                                        ENDIF
                                    ELSE
                                        IF (IS_Y_PERIODIC) THEN
                                            IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST))&
                                            EXIT L190 !    goto 200
                                        ELSE
                                            IF (J_OLD == MY) EXIT L190
                                        ENDIF
                                    ENDIF
                                ENDDO L190
                            ENDIF
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO L180
            !      *******************************************************************
            !       Finished with basin and lake discharge determination.  All we have
            !       left to do is to is to do some summary bookkeeping and also to
            !       more carefully calculate the lake_surface_elevation whose evaporation just
            !       balances the water inputs.  That is what is going on in the do loop
            !       labeled l8281.
            !      *******************************************************************             
        ENDIF
        ADDAREA=0
        ADDQ=0.0

        285   FORMAT(' END OF AREA CALCULATIONS')
        DO K=1,TOTAL_BASINS
            IF (ENCLOSED(K)) THEN
                IF (LOCAL_BASIN_DISCHARGE(K) == 0.0) THEN
                    IF (OVERFLOWS(K)) THEN
                        III1=1
                    ELSE
                        III1=0
                    ENDIF
                    III2=0
                    TTEMP=0.0
                    DO J=1,MY
                        DO I=1,MX
                            IF (BASIN_NUMBER(I,J) == K) THEN
                                III2=III2+1
                                TTEMP=TTEMP+DISCHARGE(I,J)
                            ENDIF
                        ENDDO
                    ENDDO
                    !WRITE(*,3765) K,III1,III2,DOWNSTREAM_BASIN(K),I_OUTFLOW(K),J_OUTFLOW(K) &
                    !,TTEMP,LAKE_AREA(K),LAKE_VOLUME(K) &
                    !,LOWEST_LAKE_ELEVATION(K),LAKE_OUTLET_ELEVATION(K),BASIN_OUTFLUX(K)
                    !WRITE(OUTHIST,3765) K,III1,III2,DOWNSTREAM_BASIN(K),I_OUTFLOW(K),J_OUTFLOW(K) &
                    !,TTEMP,LAKE_AREA(K),LAKE_VOLUME(K) &
                    !,LOWEST_LAKE_ELEVATION(K),LAKE_OUTLET_ELEVATION(K),BASIN_OUTFLUX(K)
                    3765  FORMAT(' ZERO BASIN, K=',I6,' OVERFLOWS=',I5,' NPOINTS=',I5' DOWNSTREAM_BASIN=',I6, &
                    ' IOUT=',I5,' JOUT=',I5,/, &
                    ' TOTAL Q=',G12.5,' LAKE_AREA=',G12.5,' LAKE_VOLUME=',G12.5, &
                    /,' LOWEST_LAKE_ELEVATION=',G12.5,' LAKE_OUTLET_ELEVATION=',G12.5,' BASIN_OUTFLUX=',G12.5)
                ENDIF
            ENDIF
        ENDDO
        GOTO 9281
        L8281: DO KK=1,TOTAL_BASINS
            IF (ENCLOSED(KK).AND.(.NOT.OVERFLOWS(KK))) THEN
                ELRANGE=LAKE_OUTLET_ELEVATION(KK)-LOWEST_LAKE_ELEVATION(KK)
                IF ((ELRANGE > 0.0).AND.(LAKE_VOLUME(KK) > 0.0)) THEN
                    ELOWER=LOWEST_LAKE_ELEVATION(KK)
                    BASIN_INFLUX(KK)=0.0
                    DO K=1,TOTAL_BASINS
                        IF (K /= KK) THEN
                            IF ((DOWNSTREAM_BASIN(K) == KK).AND.OVERFLOWS(K)) THEN
                                BASIN_INFLUX(KK)=BASIN_INFLUX(KK)+BASIN_OUTFLUX(K)
                            ENDIF
                        ENDIF
                    ENDDO
                    INSUPPLY=BASIN_INFLUX(KK)+LOCAL_BASIN_DISCHARGE(KK)
                    OLDAREA=0.0
                    NEEDAREA=INSUPPLY/EVAPORATION_RATE
                    LAKE_AREA(KK)=MIN(NEEDAREA,LAKE_AREA(KK))
                    IF (INSUPPLY <= 0.0) THEN
                        LAKE_SURFACE_ELEVATION(KK)=LOWEST_LAKE_ELEVATION(KK)
                    ELSE
                        L2177: DO
                            EUPPER=MIN((ELOWER+ELRANGE/10.0),LAKE_OUTLET_ELEVATION(KK))
                            LAREA=0.0
                            DO  J=JLW,JHGH
                                DO  I=JLW,JHGH
                                    IF (BASIN_NUMBER(I,J) == KK) THEN
                                        IF ( ELEVATION(I,J) <= EUPPER ) THEN
                                            LAREA=LAREA+CELL_SIZE*CELL_SIZE
                                        ENDIF
                                    ENDIF
                                ENDDO
                            ENDDO
                            IF (LAREA >= NEEDAREA) THEN
                                AAAA=(EUPPER-ELOWER)/(LAREA-OLDAREA)
                                BBBB=ELOWER-AAAA*OLDAREA
                                LAKE_SURFACE_ELEVATION(KK)=AAAA*NEEDAREA+BBBB
                                EXIT L2177
                            ELSE
                                ELOWER=EUPPER
                                OLDAREA=LAREA
                                IF (ELOWER >= LAKE_OUTLET_ELEVATION(KK)) THEN
                                    LAKE_SURFACE_ELEVATION(KK)=LAKE_OUTLET_ELEVATION(KK)
                                    EXIT L2177
                                ENDIF
                            ENDIF
                        ENDDO L2177
                    ENDIF
                ENDIF
            ENDIF
        ENDDO L8281
        9281  CONTINUE
        !WRITE(*,5668)
        5668  FORMAT(' NEW LAKELEVELS CALCULATED')
        338 FORMAT('********')
        !      ********************************************************************
        !        If discharge is not directly proportional to drainage area, then
        !        correct for non-linearity
        !      ********************************************************************
        DO J=1,MY
            DO I=1,MX
                IF (DISCHARGE(I,J) > MAXIMUMQ) MAXIMUMQ=DISCHARGE(I,J)
            ENDDO
        ENDDO
        MAXIMUMQ=MAXIMUMQ/(CELL_AREA*DISCHARGE_SCALE_FACTOR)
        IF (RESCALE_DISCHARGES)  THEN
            DISCHARGE=DISCHARGE**DISCHARGE_EXPONENT
            ROUTED_DISCHARGE=ROUTED_DISCHARGE**DISCHARGE_EXPONENT
        ENDIF
        CHANNEL_WIDTH=CHANNEL_WIDTH_CONSTANT*DISCHARGE**CHANNEL_WIDTH_EXPONENT
        !     **********************************************************************
        !      Determine sediment transport through each cell based upon elevation
        !         changes locally and upstream during last iteration.  Sediment is not
        !         routed through depressions
        !     **********************************************************************
        ERODETOADD=0.0
        ERGRAVTOADD=0.0
        !      *********************************************************************
        !       sedadd is the sediment load added to the stream flow by local bedrock
        !         erosion
        !       tsedcont sumS up the sediment flux contributions to
        !         local sediment yield  - used for debugging
        !         CHANNEL_EROSION_RATE holds the local bedrock channel erosion from the last iteration
        !      *********************************************************************
        !      *********************************************************************
        !        Sediment_yield is accumulated downstream - It is assumed that all sediment in
        !         bedrock channels moves instantaneously downstream - i.e., moves much more
        !         rapidly than the timescale for vertical elevation changes
        !        ditto for bederode
        !      *********************************************************************
        L710: DO  J=1,MYY
            M710: DO  I=1,MX
                I_OLD=I
                J_OLD=J
                IF (.NOT.IS_SEDIMENT_COVERED(I,J)) THEN
                    IF (CHANNEL_EROSION_RATE(I,J) > 0.0) THEN
                        SEDADD=CHANNEL_EROSION_RATE(I,J)*CELL_AREA
                    ELSE
                        SEDADD=0.0
                    ENDIF
                    !      *********************************************************************
                    !       If there is an incoming river from out of the matrix,
                    !         add in its sediment supply
                    !      *********************************************************************
                    IF (HAVE_INFLUENT_RIVERS.AND.IS_INFLUENT_RIVER_LOCATION(I,J)) THEN
                        DO KL=1,NUMBER_OF_INFLUENT_RIVERS
                            IF((I == I_INFLUENT_RIVER(KL)).AND.(J == J_INFLUENT_RIVER(KL)))THEN
                                SEDADD=SEDADD+INFLUENT_RIVER_SEDIMENT_LOAD(KL)
                            ENDIF
                        ENDDO
                    ENDIF
                    TSEDCONT=TSEDCONT+SEDADD
                    TOSEDCOVER=.FALSE.
                ELSE
                    CYCLE M710
                ENDIF
                SEDIMENT_YIELD(I,J)=SEDIMENT_YIELD(I,J)+SEDADD
                IF (SEDIMENT_YIELD(I,J)>ITERATION_MAXIMUM_SEDIMENT_YIELD) &
                    ITERATION_MAXIMUM_SEDIMENT_YIELD=SEDIMENT_YIELD(I,J)
                L720: DO
                    LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                    !     **********************************************************************
                    !      Branch if a drainage exit or depression
                    !     **********************************************************************
                    IF (LOCAL_DIRECTION < 2) EXIT L720
                    J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                    I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                    IF (IS_X_PERIODIC) THEN
                        IF (I_NEW < 1) I_NEW=MX
                        IF (I_NEW > MX) I_NEW=1
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (J_NEW > MY) J_NEW=1
                        IF (J_NEW < 1) J_NEW=MY
                    ENDIF
                    IF (IS_SEDIMENT_COVERED(I_NEW,J_NEW)) THEN
                        TOSEDCOVER=.TRUE.
                        !      *********************************************************************
                        !       Increments total work and transport variables
                        !      *********************************************************************
                        IF (DO_FLUVIAL_DETACHMENT) THEN
                            ERODETOADD=ERODETOADD+SEDADD*SQRT(REAL((I_OLD-I)**2+ &
                            (J_OLD-J)**2))
                            ERGRAVTOADD=ERGRAVTOADD+SEDADD*(ELEVATION(I,J)-ELEVATION(I_OLD,J_OLD))
                        ENDIF
                    ELSE
                        IF (TOSEDCOVER) EXIT L720
                    ENDIF
                    I_OLD=I_NEW
                    J_OLD=J_NEW
                    !     **********************************************************************
                    !      Otherwise add in contribution from upstream
                    !     **********************************************************************
                    SEDIMENT_YIELD(I_OLD,J_OLD)=SEDIMENT_YIELD(I_OLD,J_OLD) + SEDADD
                    IF (SEDIMENT_YIELD(I_OLD,J_OLD)>ITERATION_MAXIMUM_SEDIMENT_YIELD) &
                    ITERATION_MAXIMUM_SEDIMENT_YIELD=SEDIMENT_YIELD(I_OLD,J_OLD)
                    IF (IS_Y_PERIODIC) THEN
                        IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) EXIT L720
                    ELSE
                        IF (J_OLD == MY) EXIT L720
                    ENDIF
                ENDDO L720
            ENDDO M710
        ENDDO L710
        CALL IS_IT_SUBMERGED()
        !      *********************************************************************
        !       Print out some info for debugging purposes
        !      *********************************************************************
        TSEDCONT=TSEDCONT*TIME_INCREMENT
        !WRITE(OUTHIST,764) TSEDCONT
        IF (DO_SEDIMENT_TRANSPORT) THEN
            IF(MOD((ITERATION+1),OUTPUT_PRINT_INTERVAL) == 0) THEN
                WRITE(OUTHIST,764) TSEDCONT
                764           FORMAT(      &
                ' BED EROSION CONTRIBUTION TO SEDIMENT YIELD=',G12.5)
                WRITE(OUTHIST,436)
                436         FORMAT(/,' SEDIMENT YIELD')
                CALL SUMMARIZE_MATRIX_DATA(SEDIMENT_YIELD,THEAVG,THEMAX,THEMIN)
            ENDIF
        ENDIF
        L8554: DO KK=1,TOTAL_BASINS
            IF (ENCLOSED(KK)) THEN
                K=DOWNSTREAM_BASIN(KK)
                KKK=K
                IF (ENCLOSED(K)) THEN
                    II=I_OUTFLOW(K)
                    JJ=J_OUTFLOW(K)
                    I_OLD=II
                    J_OLD=JJ
                    LOCAL_DIRECTION=FLOW_DIRECTION(I_OLD,J_OLD)
                    IF (SUBMERGED(I_OLD,J_OLD).OR.(LOCAL_DIRECTION < 2)) CYCLE L8554
                    L8556:    DO
                        J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                        I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                        IF (IS_X_PERIODIC) THEN
                            IF (I_NEW < 1) I_NEW=MX
                            IF (I_NEW > MX) I_NEW=1
                        ENDIF
                        IF (IS_Y_PERIODIC) THEN
                            IF (J_NEW > MY) J_NEW=1
                            IF (J_NEW < 1) J_NEW=MY
                        ENDIF
                        LOCAL_DIRECTION=FLOW_DIRECTION(I_NEW,J_NEW)
                        IF (SUBMERGED(I_NEW,J_NEW).OR.(LOCAL_DIRECTION < 2) &
                        .OR.(BASIN_NUMBER(I_NEW,J_NEW) /= BASIN_NUMBER(I_OLD,J_OLD ))) &
                        THEN
                            K=BASIN_NUMBER(I_OLD,J_OLD)
                            IF (ENCLOSED(K)) THEN
                                IF  ((II >= 1).AND.(II <= MX) &
                                .AND.(I_OUTFLOW(DOWNSTREAM_BASIN(K)) >= 1) &
                                .AND.(J_OUTFLOW(K) >= 1).AND.(I_OUTFLOW(K) <= MX).AND. &
                                (J_OUTFLOW(K) <= MY).AND.(JJ >= 1).AND.(JJ <= MY)) THEN
                                    CALL DRAWLINE(II,JJ,I_OUTFLOW(K), &
                                    J_OUTFLOW(K),ROUTED_DISCHARGE(I_NEW,J_NEW))
                                ENDIF
                            ENDIF
                            EXIT L8556
                        ENDIF
                        I_OLD=I_NEW
                        J_OLD=J_NEW
                        IF (IS_Y_PERIODIC) THEN
                            IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) EXIT L8556
                            IF (J_OLD == MY) EXIT L8556
                        ENDIF
                    ENDDO L8556
                ENDIF
            ENDIF
        ENDDO L8554
        RETURN
    END !SUBROUTINE DRAINAGE_BASIN_LAKE_FLOW 
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DRAWLINE(ISTART,JSTART,IEND,JEND,FILLVALUE)
        USE ERODE_GLOBALS
        USE AREA_GLOBALS
        USE LAKE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ISTART,JSTART,IEND,JEND
        INTEGER :: INOW,JNOW,INC1, II,CRITICAL_DISTANCE,ITOTAL
        INTEGER :: IDX,JDY,IXX,JYY,IXT,JYT,IPLOT,JPLOT
        LOGICAL :: DOPLOT,DONOTHING
        REAL(4), INTENT(IN) :: FILLVALUE
        CRITICAL_DISTANCE=300
        ITOTAL=(ISTART-IEND)**2+(JSTART-JEND)**2
        IF (ITOTAL>CRITICAL_DISTANCE) RETURN
        DONOTHING=.FALSE.
        IF (DONOTHING) RETURN
        DOPLOT=.TRUE.
        IF ((ISTART == IEND).AND.(JSTART == JEND)) RETURN
        INOW=ISTART
        JNOW=JSTART
        IDX=IEND-ISTART
        JDY=JEND-JSTART
        IXX=ABS(IDX)
        JYY=ABS(JDY)
        INC1=MAX(IXX,JYY)
        INOW=ISTART
        JNOW=JSTART
        IPLOT=ISTART
        JPLOT=JSTART
        IF ((IPLOT > 0).AND.(JPLOT > 0).AND. &
        (IPLOT <= MX).AND.(JPLOT <= MY)) THEN
            IF (ROUTED_DISCHARGE(IPLOT,JPLOT) < FILLVALUE)  THEN
                ROUTED_DISCHARGE(IPLOT,JPLOT)=FILLVALUE
            ENDIF
        ENDIF
        IXT=0
        JYT=0
        DO II=1,INC1
            IXT=IXT+IXX
            JYT=JYT+JYY
            DOPLOT=.FALSE.
            IF (IXT > INC1) THEN
                DOPLOT=.TRUE.
                IXT=IXT-INC1
                IPLOT=IPLOT+SIGN(1,IDX)
            ENDIF
            IF (JYT > INC1) THEN
                DOPLOT=.TRUE.
                JYT=JYT-INC1
                JPLOT=JPLOT+SIGN(1,JDY)
            ENDIF
            IF (DOPLOT) THEN
                IF ((IPLOT > 0).AND.(JPLOT > 0).AND. &
                (IPLOT <= MX).AND.(JPLOT <= MY)) THEN
                    IF (ROUTED_DISCHARGE(IPLOT,JPLOT) < FILLVALUE) THEN
                        ROUTED_DISCHARGE(IPLOT,JPLOT)=FILLVALUE
                    ENDIF
                ENDIF
            ENDIF
        ENDDO
        RETURN
    END ! SUBROUTINE DRAWLINE
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DRAWSLINE(ISTART,JSTART,IEND,JEND,FILLVALUE)
        USE ERODE_GLOBALS
        USE AREA_GLOBALS
        USE LAKE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ISTART,JSTART,IEND,JEND
        INTEGER :: INOW,JNOW,INC1, II,CRITICAL_DISTANCE,ITOTAL
        INTEGER :: IDX,JDY,IXX,JYY,IXT,JYT,IPLOT,JPLOT
        LOGICAL :: DOPLOT
        REAL(4), INTENT(IN) :: FILLVALUE
        CRITICAL_DISTANCE=300
        ITOTAL=(ISTART-IEND)**2+(JSTART-JEND)**2
        IF (ITOTAL>CRITICAL_DISTANCE) RETURN
        RETURN
        IF ((ISTART == IEND).AND.(JSTART == JEND)) RETURN
        INOW=ISTART
        JNOW=JSTART
        IDX=IEND-ISTART
        JDY=JEND-JSTART
        IXX=ABS(IDX)
        JYY=ABS(JDY)
        INC1=MAX(IXX,JYY)
        INOW=ISTART
        JNOW=JSTART
        IPLOT=ISTART
        JPLOT=JSTART
        IF ((IPLOT > 0).AND.(JPLOT > 0).AND. &
        (IPLOT <= MX).AND.(JPLOT <= MY)) THEN
        ENDIF
        IXT=0
        JYT=0
        DO II=1,INC1
            IXT=IXT+IXX
            JYT=JYT+JYY
            DOPLOT=.FALSE.
            IF (IXT > INC1) THEN
                DOPLOT=.TRUE.
                IXT=IXT-INC1
                IPLOT=IPLOT+SIGN(1,IDX)
            ENDIF
            IF (JYT > INC1) THEN
                DOPLOT=.TRUE.
                JYT=JYT-INC1
                JPLOT=JPLOT+SIGN(1,JDY)
            ENDIF
            IF (DOPLOT) THEN
                IF ((IPLOT > 0).AND.(JPLOT > 0).AND. &
                (IPLOT <= MX).AND.(JPLOT <= MY)) THEN
                ENDIF
            ENDIF
        ENDDO
        RETURN
    END ! SUBROUTINE DRAWSLINE
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PELAGIC_DEPOSIT()
        USE ERODE_GLOBALS
        USE AREA_GLOBALS
        USE LAKE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,IE,IW,JN,JS,K
        REAL(4) :: DEPOSITDEPTH,ECREEP,WCREEP,NCREEP,SCREEP,LOCAL_DEPTH
        REAL(4) :: NECREEP,NWCREEP,SECREEP,SWCREEP,SUMPELAGIC,SUMDEPTH,NDEPTH
    !  MODIFIES: LAKE_AREA,PELAGICCHANGE,PELAGICDIFF,ELEVATION
        PELAGICCHANGE=0.0
        PELAGICEDIFF=0.0
        SUMPELAGIC=0.0
        SUMDEPTH=0.0
        NDEPTH=0.0
        LAKE_AREA=0.0
        DO J=1,MYY
            DO I=1,MX
                IF (SUBMERGED(I,J)) THEN
                    LAKE_AREA(BASIN_NUMBER(I,J))=LAKE_AREA(BASIN_NUMBER(I,J))+CELL_SIZE*CELL_SIZE
                ENDIF
            ENDDO
        ENDDO
        DO K=1,TOTAL_BASINS
            IF (PELAGIC_SEDIMENT_VOLUME(K) > 0.0) THEN
                SUMPELAGIC=SUMPELAGIC+PELAGIC_SEDIMENT_VOLUME(K)
                IF (LAKE_AREA(K) > MINIMUM_PELAGIC_LAKE_SIZE*CELL_SIZE*CELL_SIZE) THEN
 !               IF (LAKE_AREA(K) > 0.0) THEN
                    DEPOSITDEPTH=PELAGIC_SEDIMENT_VOLUME(K)/LAKE_AREA(K)
                ELSE
                    DEPOSITDEPTH=0.0
                ENDIF
                SUMDEPTH=SUMDEPTH+DEPOSITDEPTH
                NDEPTH=NDEPTH+1.0
                DO J=1,MYY
                    DO I=1,MX
                        IF ((BASIN_NUMBER(I,J) == K).AND.(SUBMERGED(I,J))) THEN
                          LOCAL_DEPTH=LAKE_SURFACE_ELEVATION(K)-ELEVATION(I,J)
                            IF (LOCAL_DEPTH>0.0) THEN
                                 PELAGICCHANGE(I,J)=MIN(DEPOSITDEPTH,LOCAL_DEPTH)
                            ENDIF
                        ENDIF
                    ENDDO
                ENDDO
            ENDIF
        ENDDO
        SUMDEPTH=SUMDEPTH/NDEPTH
        SUMPELAGIC=SUMPELAGIC/NDEPTH
        !WRITE(*,3421)SUMDEPTH,SUMPELAGIC,NDEPTH,TOTAL_BASINS
        3421 FORMAT (' AVG. PELAGIC DEPTH=',G12.5,' AVG. VOLUME=',G12.5,' N=',G12.5,' TOTAL_BASINS=',I6)
        DO J=1,MYY
            DO I=1,MX
                IF ((PELAGIC_SEDIMENT_VOLUME(BASIN_NUMBER(I,J)) > 0.0).AND.SUBMERGED(I,J)) THEN
                    IF (ELEVATION(I,J) > INITIAL_ELEVATION(I,J)) THEN
                        ECREEP=0.0
                        WCREEP=0.0
                        NCREEP=0.0
                        SCREEP=0.0
                        NWCREEP=0.0
                        SWCREEP=0.0
                        NECREEP=0.0
                        SECREEP=0.0
                        IF (IS_X_PERIODIC) THEN
                            IE=I+1
                            IF (IE > MX) IE=1
                            IW=I-1
                            IF (IW < 1) IW=MX
                        ELSE
                            IE=I+1
                            IE=MIN(IE,MX)
                            IW=I-1
                            IW=MAX(IW,1)
                        ENDIF
                        IF (IS_Y_PERIODIC) THEN
                            JS=J+1
                            IF (JS > MYY) JS=1
                            JN=J-1
                            IF (JN < 1) JN=MYY
                        ELSE
                            JS=J+1
                            JS=MIN(JS,MYY)
                            JN=J-1
                            JN=MAX(JN,1)
                        ENDIF
                        IF (ELEVATION(I,J) > ELEVATION(IE,J)) THEN
                            ECREEP=PELAGICCREEP*(ELEVATION(I,J)-ELEVATION(IE,J))
                        ENDIF
                        IF (ELEVATION(I,J) > ELEVATION(IW,J)) THEN
                            WCREEP=PELAGICCREEP*(ELEVATION(I,J)-ELEVATION(IW,J))
                        ENDIF
                        IF (ELEVATION(I,J) > ELEVATION(I,JN)) THEN
                            NCREEP=PELAGICCREEP*(ELEVATION(I,J)-ELEVATION(I,JN))
                        ENDIF
                        IF (ELEVATION(I,J) > ELEVATION(I,JS)) THEN
                            SCREEP=PELAGICCREEP*(ELEVATION(I,J)-ELEVATION(I,JS))
                        ENDIF
                        IF (ELEVATION(I,J) > ELEVATION(IE,JN)) THEN
                            NECREEP=PELAGICCREEP*(ELEVATION(I,J)-ELEVATION(IE,JN))/SQRTOFTWO
                        ENDIF
                        IF (ELEVATION(I,J) > ELEVATION(IE,JS)) THEN
                            SECREEP=PELAGICCREEP*(ELEVATION(I,J)-ELEVATION(IE,JS))/SQRTOFTWO
                        ENDIF
                        IF (ELEVATION(I,J) > ELEVATION(IW,JN)) THEN
                            NWCREEP=PELAGICCREEP*(ELEVATION(I,J)-ELEVATION(IW,JN))/SQRTOFTWO
                        ENDIF
                        IF (ELEVATION(I,J) > ELEVATION(IW,JS)) THEN
                            SWCREEP=PELAGICCREEP*(ELEVATION(I,J)-ELEVATION(IW,JS))/SQRTOFTWO
                        ENDIF
                        PELAGICEDIFF(I,J)=PELAGICEDIFF(I,J)-ECREEP*CROSS_WEIGHTING
                        PELAGICEDIFF(IE,J)=PELAGICEDIFF(IE,J)+ECREEP*CROSS_WEIGHTING
                        PELAGICEDIFF(I,J)=PELAGICEDIFF(I,J)-WCREEP*CROSS_WEIGHTING
                        PELAGICEDIFF(IW,J)=PELAGICEDIFF(IW,J)+WCREEP*CROSS_WEIGHTING
                        PELAGICEDIFF(I,J)=PELAGICEDIFF(I,J)-NCREEP*CROSS_WEIGHTING
                        PELAGICEDIFF(I,JN)=PELAGICEDIFF(I,JN)+NCREEP*CROSS_WEIGHTING
                        PELAGICEDIFF(I,J)=PELAGICEDIFF(I,J)-SCREEP*CROSS_WEIGHTING
                        PELAGICEDIFF(I,JS)=PELAGICEDIFF(I,JS)+SCREEP*CROSS_WEIGHTING
                        IF (IE /= I) THEN
                            IF (JN /= J) THEN
                                PELAGICEDIFF(I,J)=PELAGICEDIFF(I,J)-NECREEP*DIAGONAL_WEIGHTING
                                PELAGICEDIFF(IE,JN)=PELAGICEDIFF(IE,JN)+NECREEP*DIAGONAL_WEIGHTING
                            ENDIF
                        ENDIF
                        IF (IE /= I) THEN
                            IF (JS /= J) THEN
                                PELAGICEDIFF(I,J)=PELAGICEDIFF(I,J)-SECREEP*DIAGONAL_WEIGHTING
                                PELAGICEDIFF(IE,JS)=PELAGICEDIFF(IE,JS)+SECREEP*DIAGONAL_WEIGHTING
                            ENDIF
                        ENDIF
                        IF (IW /= I) THEN
                            IF (JN /= J) THEN
                                PELAGICEDIFF(I,J)=PELAGICEDIFF(I,J)-NWCREEP*DIAGONAL_WEIGHTING
                                PELAGICEDIFF(IW,JN)=PELAGICEDIFF(IW,JN)+NWCREEP*DIAGONAL_WEIGHTING
                            ENDIF
                        ENDIF
                        IF (IW /= I) THEN
                            IF (JS /= J) THEN
                                PELAGICEDIFF(I,J)=PELAGICEDIFF(I,J)-SWCREEP*DIAGONAL_WEIGHTING
                                PELAGICEDIFF(IW,JS)=PELAGICEDIFF(IW,JS)+SWCREEP*DIAGONAL_WEIGHTING
                            ENDIF
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        DO J=1,MYY
            DO I=1,MX
                ELEVATION(I,J)=ELEVATION(I,J)+PELAGICCHANGE(I,J)+PELAGICEDIFF(I,J)*TIME_INCREMENT
            ENDDO
        ENDDO
        RETURN
    END ! SUBROUTINE PELAGIC_DEPOSIT
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE CHECK_FLOW_PATH(K)
        USE ERODE_GLOBALS
        USE AREA_GLOBALS
        USE LAKE_GLOBALS
        IMPLICIT NONE
        INTEGER K,I_OLD,J_OLD,LOCAL_DIRECTION,I_NEW,J_NEW,NN
        REAL(4) :: QQQQ,QAAA
        RETURN
        IF (OVERFLOWS(K)) THEN
            I_OLD=I_OUTFLOW(K)
            J_OLD=J_OUTFLOW(K)
            LOCAL_DIRECTION=FLOW_DIRECTION(I_OLD,J_OLD)
            !WRITE(OUTHIST,6321) K,I_OLD,J_OLD,LOCAL_DIRECTION
            6321 FORMAT(/,'DEBUG FLOWROUTE, BASIN ',I5, &
            ' IOUT=',I5,' JOUT=',I5,' LOCAL_DIRECTION=',I2)
            !      *******************************************************************
            !       In the rare to non-existent case where the outlet of our basin is
            !        also the outlet of the next basin downstream, then we don't need to
            !        do anything here.
            !     ********************************************************************
            IF ((I_OLD == I_OUTFLOW(DOWNSTREAM_BASIN(K))).AND.  &
            (J_OLD == J_OUTFLOW(DOWNSTREAM_BASIN(K)))) THEN
                !WRITE(OUTHIST,28443) I_OLD,J_OLD,K,DOWNSTREAM_BASIN(K)
                28443 FORMAT(' SAME OUTLET, IOUT=',I5,' JOUT=',I5,' K=',I8,' TB=',I8,/)
                RETURN
            ENDIF
            !      *******************************************************************
            !        This just checks to make sure that nothing terrible has gone
            !         wrong in our calculations
            !      *******************************************************************
            IF ((I_OLD < 1).OR.(I_OLD > MX).OR. &
            (J_OLD < 1).OR.(J_OLD > MY)) THEN
                WRITE(OUTHIST,377) K,I_OLD,J_OLD
                377             FORMAT(' ***BASIN INTEGRITY PROBLEM AT BASIN' &
                ,I5,', I_OLD & J_OLD=', &
                2I5)
            ENDIF
            !      ******************************************************************
            !        Basically what is happening here is that we are accumulating discharges
            !         downstream.
            !        Within each basin we add the total cumulative discharge from upstream until
            !          we get to a depression or sinkhole in that basin.  if the location in the
            !          basin is below the lakelevel, then the discharge is really imaginary, but
            !          we calculate it anyway.
            !      ********************************************************************
            QAAA=BASIN_OUTFLUX(K)
            QQQQ=DISCHARGE(I_OLD,J_OLD)
            !WRITE(OUTHIST,4776) K,QQQQ,QAAA
            4776 FORMAT('K=',I6,' Q=',G12.5,' ADDQ=',G12.5)
            IF (DO_FLOW_BOUNDARIES) THEN
                IF (FLOW_DIRECTION(I_OLD,J_OLD) == 1) THEN
                    RETURN
                ENDIF
            ELSE
                IF (IS_Y_PERIODIC) THEN
                    IF ((I_OLD == ILOWEST).AND. &
                    (J_OLD == JLOWEST)) RETURN
                ELSE
                    IF (J_OLD == MY) RETURN
                ENDIF
            ENDIF
            LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
            IF (LOCAL_DIRECTION >= 0) THEN
                IF ((I_OLD /= ILOWEST).AND.(J_OLD /= JLOWEST)) THEN
                    IF (LOCAL_DIRECTION >= 0) THEN
                        L190: DO
                            J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                            I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                            IF (IS_X_PERIODIC) THEN
                                IF (I_NEW < 1) I_NEW = MX
                                IF (I_NEW > MX) I_NEW=1
                            ENDIF
                            IF (IS_Y_PERIODIC) THEN
                                IF (J_NEW < 1) J_NEW=MY
                                IF (J_NEW > MY) J_NEW=1
                            ENDIF
                            I_OLD=I_NEW
                            J_OLD=J_NEW
                            LOCAL_DIRECTION=FLOW_DIRECTION(I_OLD,J_OLD)
                            NN=BASIN_NUMBER(I_OLD,J_OLD)
                            QQQQ=DISCHARGE(I_OLD,J_OLD)
                            !WRITE(OUTHIST,4771) I_OLD,J_OLD,LOCAL_DIRECTION,BASIN_NUMBER(I_OLD,J_OLD),QQQQ
                            4771 FORMAT(' NEW POINT, I=',I5,' J=',I5,' LOCAL_DIRECTION=',I5,' BN=',I5,' Q=', &
                            G12.5)
                            IF (BASIN_NUMBER(I_OLD,J_OLD) /= DOWNSTREAM_BASIN(K)) THEN
                                !WRITE(OUTHIST,823) K,DOWNSTREAM_BASIN(K), &
                                !BASIN_NUMBER(I_OLD,J_OLD),I_OLD,J_OLD  &
                                !,DISCHARGE(I_OLD,J_OLD),QAAA
                                823 FORMAT(' CHANGE BASIN NUMBER',I5, &
                                ' TO BASIN'  &
                                ,I5,/,' NOW IN BASIN',I5,' AT I=',I5,' J=',I5, &
                                ' AREA=',G12.5,' ADDQ=',G12.5,/)
                                RETURN
                            ENDIF
                            !      *******************************************************************
                            !       In the rare to non-existent case where the outlet of our basin is
                            !        also the outlet of the next basin downstream, then we don't need to
                            !        do anything here.
                            !     ********************************************************************
                            IF ((I_OLD == I_OUTFLOW(DOWNSTREAM_BASIN(K))).AND.(J_OLD == J_OUTFLOW(DOWNSTREAM_BASIN(K)))) THEN
                                !WRITE(OUTHIST,28442) &
                                !I_OLD,J_OLD,K,DOWNSTREAM_BASIN(K),I_OUTFLOW(DOWNSTREAM_BASIN(K)),J_OUTFLOW(DOWNSTREAM_BASIN(K))
                                28442 FORMAT(' SAME OUTLETS, I_OLD=',I5,' J_OLD=',I5,' K=',I8,' TB=',I8 &
                                ,' NEXT I_OLDUT=',I5,' NEXT JOUT=',I5,/)
                                EXIT L190
                            ENDIF
                            !      *******************************************************************
                            !       Again, look for infinite loops and abort if they occur
                            !      *******************************************************************
                            IF(DISCHARGE(I_OLD,J_OLD) > 1.0E+20) THEN
                                ABORT_SIMULATION=.TRUE.
                                IABORTMIN=I_OLD-3
                                IABORTMAX=I_OLD+3
                                JABORTMIN=J_OLD-3
                                JABORTMAX=J_OLD+3
                                !WRITE(OUTHIST,378) I_OLD,J_OLD,BASIN_NUMBER(I_OLD,J_OLD)
                                378             FORMAT (' *** INFINITE CYCLE IN OUTLET COMBINE AT '&
                                ,2I6,' ***' &
                                ,/,' AT BASIN NUMBER: ',I5,/)
                                ABORT_SIMULATION=.TRUE.
                                RETURN
                            ENDIF
                            IF ((I_OLD == I_OUTFLOW(DOWNSTREAM_BASIN(K))).AND.(J_OLD == J_OUTFLOW(DOWNSTREAM_BASIN(K)))) &
                            EXIT L190 
                            IF (LOCAL_DIRECTION < 0) EXIT L190
                            IF (DO_FLOW_BOUNDARIES) THEN
                                IF (FLOW_DIRECTION(I_OLD,J_OLD) == 1) THEN
                                    EXIT L190
                                ENDIF
                            ELSE
                                IF (IS_Y_PERIODIC) THEN
                                    IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST))&
                                    EXIT L190
                                ELSE
                                    IF (J_OLD == MY) EXIT L190
                                ENDIF
                            ENDIF
                        ENDDO L190
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
        RETURN
    END ! SUBROUTINE CHECK_FLOW_PATH

