
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
    
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !  This file contains several subroutines:
    !      READ_INPUT_PARAMETERS() which reads in simulation-wide parameters such as size of
    !         simulation domain, which modules are used, boundary conditions, simulation time scale
    !         simulation duration, data output frequency, and miscellaneous parameters.  It also
    !         allocates the relevant data matrices and calls routines to read in parameters
    !         appropriate to individual modules, and reads in initial elevations
    !      INITIALIZE_VARIABLES() Sets up simulation-wide variables and structures.
    !      Several functions to generate random numbers are included, e.g. REAL (4) FUNCTION NORMAL_RANDOM_DEVIATE()
    !      BOUNDARY_CONDITIONS() is called each iteration to enforce conditions or to change parameters
    !         according to some rule.  This routine can be edited to change rules
    !      DETERMINE_ERODIBILITY() and READ_ERODIBILITY() These routines can be utilized with a direct
    !        access 3-D file of bedrock "erodibility" which determines fluvial erosion rates and weatherability
    !      MAKE_EVENT() and SETUP_EVENTS() Can be programmed to change simulation parameters and boundary
    !         conditions at specified times during the simulation
    !      FIND_OCEAN_ELEVATION() can be used to make time-varying ocean elevation as
    !         base level control.
    !      DETERMINE_CRATERING_RATE() can be used to enforce time-varying rates of new crater impacts, 
    !         usually decreasing with time
    !      ALLOCERROR() shuts down the simulation if matrix allocation fails during program initialization
    !      
    !      
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE READ_INPUT_PARAMETERS()
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        USE EOLIAN_GLOBALS
        USE LAVA_GLOBALS
        USE SEDDEBUG_GLOBALS
        USE SEDROUTE_GLOBALS
        USE AREA_GLOBALS
        USE LAKE_GLOBALS
        USE ACCRETION_GLOBALS
        USE EVENT_GLOBALS
		USE GRAVEL_MIXTURE_GLOBALS
        USE MASS_FLOW_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,IXDEBUG,IWIDTH,JWIDTH,RANDDIRUSE,STICKYUSE,PARAMS
        INTEGER :: USEWET,STATISTDO, SHADE_BORDER_INDEX, FIXED_SUNANGLE_INDEX
        INTEGER :: ISEDIMENT,ISEDROUTE,ISEDDIFFUSE,ISLOPEUSE,IREADALLUV
        INTEGER :: IVARRATEUSE,RANDTHRESHUSE,HIGHRATEUSE,WEATHER2USE
        INTEGER :: INEW_SIMULATION,DEFORMUSE
        INTEGER :: IRESINPUT,IREFLECT,IW,IE,IHORIZONTAL_LOWER_BOUNDARY,INON_ERODING_LOWER_BOUNDARY
        INTEGER :: BEDROCKHIGH,USEYPERIODIC,CRUSTUSE,USEXPERIODIC
        INTEGER :: IOTEMP3,BEDEXPLICIT,ISEDDEBUG
        INTEGER :: SMOOTHSEDUSE,NEWDIRECTIONUSE,DETACHUSE,IWATERLOWER
        INTEGER :: MX1,MY1,ISOFTCRATER,IQCONSTANT,IOCEANVAR
        INTEGER :: IDIVAVG,QQUSE,SEEPUSE,DOINRIVER,GRAVEL_PRINT
        INTEGER :: IDOCRATER,IDOLAVA,IDOEOLIAN,IDOERODE
        INTEGER :: ITAUAREA,IEVENTS,ISEEPAGEWEATHER,RADIATION_USE, GRAVEL_PRING
        INTEGER :: IDOOCEAN,ROERINGUSE, WIDTH_RELATIVE_INDEX
        INTEGER :: IUSETOTALEXPOSE,IUSENORMAL,IDOLAKES,SPATIAL_VARIATION_USE
        INTEGER :: IEXPFLOW,IDOACCRETION,EXPOSECREEPUSE,IERROR,USE_CRATER_EDGE_ABORT
        INTEGER :: USEFLOWBOUND, INVERSEUSE, TOPEXPOSEUSE, EXPOSESMOOTH, USE_REGOLITH_ABLATION
        INTEGER :: SHEAR_RATIO_USE
        INTEGER :: IDOAVALANCHE, IDO_DEPTHCREEP, DO_WRITE_BINARY_IMAGE_FILE, DO_ONLY_EOLIAN_DEPOSITION
        INTEGER :: IDOMASSFLOW,WE_DO_FLOW_VOLUME,VARIABLE_MASS_FLOW_RATE_USE
        CHARACTER (80) :: TEXTLINE

        REAL (4) :: READ_ERODIBILITY,EEEMIN
        REAL (4) :: THEMIN,THEAVG,THEMAX
        real(4) :: logsd,normvar,normmean,normsd,rrrr,LOGNORMAL_RANDOM_DEVIATE1
        REAL(4):: THEMINIMUM,ATIME
        INTEGER (4) :: IFOLD
        INTEGER :: ISHOWCALC, DO_CORRECT_BIAS,DO_CENTRAL_PEAK,EROSION_MASK_USE
        INTEGER :: IFLUXDEPEND,IMODEL_PELAGIC_DEPOSITION
        INTEGER :: ICHOSE,DO_SEDIMENT_YIELD_SCALING
        REAL (4) :: TEMP1,TEMP2,TEMP3
        INTEGER :: NRANDOM,XMX,XMY
        INTEGER :: OUTACTIVE, OUTD
        EXTERNAL READ_ERODIBILITY
        !        *******************************************************************
        !   MODIFIES:  (MOST SIMULATION GLOBAL PARAMETERS AND MATRICES)
        !   CALLS:  READ_REGOLITH_THICKNESS, SETUP_DISTANCE_WEIGHTING, FIND_MODIFICATION_RANGE
        !           SUMMARIZE_LOGICAL_MATRIX, SUMMARIZE_MATRIX_DATA
        !
        !       **********************************************************************

        !        *******************************************************************
        !        *******************************************************************
        !          r e a d  i n  t h e  s i m u l a t i o n  p a r a m e t e r s
        !        *******************************************************************
        !        *******************************************************************
        !
        !        *******************************************************************
        !        ********************** boundary conditions ************************
        !        *******************************************************************
        USE_CRATER_EVENT=.TRUE.
		PARAMS=92
        IOTEMP3=43
        DEFAULT_TIME_INCREMENT=0.001
		OPEN(PARAMS,FILE='MARSSIM_INITIAL_BOUNDARY_CONDITIONS.PRM',ACTION='READ')
        READ(PARAMS,22) TEXTLINE
        22    FORMAT(A)
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) ISEED
        write(*,822) iseed
822     format(I8)
        READ(PARAMS,*) INEW_SIMULATION
        READ(PARAMS,*) IDOERODE
        READ(PARAMS,*) IDOCRATER
        READ(PARAMS,*) IDOLAVA
        READ(PARAMS,*) IDOEOLIAN
        READ(PARAMS,*) IDOOCEAN
        READ(PARAMS,*) IDOACCRETION
        READ(PARAMS,*) IDOLAKES
        READ(PARAMS,*) IDOAVALANCHE
        READ(PARAMS,*) IDOMASSFLOW
        write(*,823) idoerode, idocrater, idolava, idoeolian, idoocean,idoaccretion,idolakes, &
        IDOAVALANCHE,IDOMASSFLOW
823     format(10(i8,' '))
        IF (INEW_SIMULATION > 0) THEN
            NEW_SIMULATION=.TRUE.
        ELSE
            NEW_SIMULATION=.FALSE.
        ENDIF
        IF (IDOCRATER > 0) THEN
            MODEL_IMPACT_CRATERING=.TRUE.
        ELSE
            MODEL_IMPACT_CRATERING=.FALSE.
        ENDIF
        IF (IDOLAVA > 0) THEN
            MODEL_LAVA_FLOWS=.TRUE.
        ELSE
            MODEL_LAVA_FLOWS=.FALSE.
        ENDIF
        IF (IDOEOLIAN > 0) THEN
            MODEL_EOLIAN_CHANGES=.TRUE.
        ELSE
            MODEL_EOLIAN_CHANGES=.FALSE.
        ENDIF
        IF (IDOERODE > 0) THEN
            FLUVIAL_AND_SLOPE_MODELING=.TRUE.
        ELSE
            FLUVIAL_AND_SLOPE_MODELING=.FALSE.
        ENDIF
        IF (IDOOCEAN > 0) THEN
            MODEL_OCEAN_LEVEL=.TRUE.
        ELSE
            MODEL_OCEAN_LEVEL=.FALSE.
        ENDIF
        IF (IDOACCRETION > 0) THEN
            MODEL_ACCRETION_AND_ABLATION=.TRUE.
        ELSE
            MODEL_ACCRETION_AND_ABLATION=.FALSE.
        ENDIF
        IF (IDOLAKES > 0) THEN
            MODEL_LAKE_EVAPORATION=.TRUE.
        ELSE
            MODEL_LAKE_EVAPORATION=.FALSE.
        ENDIF
        IF (IDOAVALANCHE > 0) THEN
            DO_AVALANCHE = .TRUE.
        ELSE
            DO_AVALANCHE = .FALSE.
        ENDIF
        IF (IDOMASSFLOW > 0) THEN
            USE_FLOW_VOLUME=.TRUE.
        ELSE
            USE_FLOW_VOLUME=.FALSE.
        ENDIF
        write(*,*)  fluvial_and_slope_modeling, model_impact_cratering, model_lava_flows,model_eolian_changes, &
            model_ocean_level, model_accretion_and_ablation, model_lake_evaporation,DO_AVALANCHE
8335        format(3(' ',i8),2(' ',g13.6))
        READ(PARAMS,*) MX
        READ(PARAMS,*) MY
        READ(PARAMS,*) MZ
        READ(PARAMS,*) INPUT_CELL_SIZE
        READ(PARAMS,*) VERTICAL_SCALING_FACTOR
        READ(PARAMS,*) CONVERT_TO_METERS
        write(*,8335) mx, my, mz
        IMMX=MX
        JMMX=MY
        LMMX=MX*MY
        RMMX=2*MAX(MX,MY)
        ALLOCATE(ELEVATION(MX,MY),INITIAL_ELEVATION(MX,MY),STAT=IERROR)
        IF (IERROR > 0) CALL ALLOCERROR()
        ALLOCATE(CUMULATIVE_ELEVATION_CHANGE(MX,MY),STAT=IERROR)
        IF (IERROR > 0) CALL ALLOCERROR()
        IF (FLUVIAL_AND_SLOPE_MODELING.OR.MODEL_EOLIAN_CHANGES.OR.MODEL_IMPACT_CRATERING.OR. &
                MODEL_ACCRETION_AND_ABLATION.OR.USE_FLOW_VOLUME) THEN
            ALLOCATE(ERODE_SLOPE(MX,MY),CHANGE_REGOLITH_STATE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
        ENDIF
        IF (FLUVIAL_AND_SLOPE_MODELING.OR.MODEL_EOLIAN_CHANGES.OR.MODEL_ACCRETION_AND_ABLATION.OR. &
                MODEL_IMPACT_CRATERING.OR.USE_FLOW_VOLUME) THEN
            ALLOCATE(D8_GRADIENT(MX,MY),CFNW(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(CUMULATIVE_CRATERING_CHANGE(MX,MY),CUMULATIVE_CRATER_EXCAVATION(MX,MY), &
                CUMULATIVE_EJECTA_DEPOSITION(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(FLOW_DIRECTION(MX,MY),IDO(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(CUMULATIVE_EOLIAN_CHANGE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
        ENDIF
        IF (FLUVIAL_AND_SLOPE_MODELING.OR.MODEL_ACCRETION_AND_ABLATION.OR.USE_FLOW_VOLUME) THEN
            ALLOCATE(CFW(MX,MY),CFN(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(CFNE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
        ENDIF
        IF (FLUVIAL_AND_SLOPE_MODELING.OR.MODEL_EOLIAN_CHANGES.OR.USE_FLOW_VOLUME.OR. &
            MODEL_ACCRETION_AND_ABLATION) THEN
            ALLOCATE(IS_SEDIMENT_COVERED(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(REGOLITH(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(IS_ROCK_SURFACE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(SEDIMENT_BASE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
        ENDIF
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            ALLOCATE(CUMULATIVE_WEATHERING(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(CUMULATIVE_SEDIMENT_DEPOSITION(MX,MY),CUMULATIVE_MASS_WASTING(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(CUMULATIVE_FLUVIAL_EROSION(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(ERODE_CHANNEL(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(RELATIVE_RESISTANCE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(CHANNEL_WIDTH(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(SEDIMENT_YIELD(MX,MY),SEDIMENT_FLUX(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(EQUILIBRIUM_GRADIENT(MX,MY),DRAINAGE_AREA(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(MAXIMUM_SEDIMENT_YIELD(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(PREVIOUS_ELEVATION(MX,MY),DISCHARGE(MX,MY),MAXIMUM_DISCHARGE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(LOCAL_DISCHARGE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            !LOCAL_DISCHARGE=1.0
            ALLOCATE(DIVERGENCE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(ERODE_REGOLITH_CHANNEL(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(DEFORMATION(MX,MY),I_OUTFLOW(LMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(J_OUTFLOW(LMMX),DOWNSTREAM_BASIN(LMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(BASIN_NUMBER(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(EROSION_DEPTH_INDEX(MX,MY),ENCLOSED(LMMX),OVERFLOWS(LMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(SUBMERGED(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(ACCELERATED_EROSION(MX,MY),DO_ACCELERATED_EROSION(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(SEDIMENT_DEPOSITED(MX,MY),IS_INFLUENT_RIVER_LOCATION(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(ERODING_LOWER_BOUNDARY(MX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(STARTING_ELEVATION(RMMX),ALLUVIAL_GRADIENT(RMMX),ACTUAL_GRADIENT(RMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR() 
            ALLOCATE(STEP_DISTANCE(RMMX),NEW_ELEVATION(RMMX),PROVISIONAL_ELEVATION(RMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(IS_BEDROCK_CHANNEL(RMMX),WATER_LEVEL(RMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(TOBEADDED(LMMX),LOCAL_BASIN_DISCHARGE(LMMX),QTOADD(LMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(NEEDTODO(LMMX),OUTER(LMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(BASIN_DRAINAGE_AREA(LMMX),LAKE_OUTLET_ELEVATION(LMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(LAKE_SURFACE_ELEVATION(LMMX),SORTING_VECTOR(LMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(NOMINAL_ERODE_SLOPE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(DONE_ONCE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(I_LOCATION(RMMX),J_LOCATION(RMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(IS_EXIT_BASIN(LMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(ROUTED_DISCHARGE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(LAKE_AREA(LMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            IF (DO_AVALANCHE) THEN
                ALLOCATE(AVALANCHE_FLUX(MX,MY),OLD_AVALANCHE_FLUX(MX,MY),STAT=IERROR)
                IF (IERROR > 0) CALL ALLOCERROR()
                ALLOCATE(CUMULATIVE_AVALANCHE_EROSION(MX,MY),STAT=IERROR)
                IF (IERROR > 0) CALL ALLOCERROR()
            ENDIF
            IF(MODEL_LAKE_EVAPORATION.OR.MODEL_IMPACT_CRATERING) THEN
                ALLOCATE(LAKE_VOLUME(LMMX),STAT=IERROR)
                IF (IERROR > 0) CALL ALLOCERROR()
                ALLOCATE(LOWEST_LAKE_ELEVATION(LMMX),NEW_BASIN_OUTFLUX(LMMX),STAT=IERROR)
                IF (IERROR > 0) CALL ALLOCERROR()
                ALLOCATE(BASIN_OUTFLUX(LMMX),OLDOVERFLOWS(LMMX),STAT=IERROR)
                IF (IERROR > 0) CALL ALLOCERROR()
                ALLOCATE(LAKEMIN(LMMX),BASINMIN(LMMX),STAT=IERROR)
                IF (IERROR > 0) CALL ALLOCERROR()
                ALLOCATE(NEXTCYCLE(LMMX),NEWGEOMETRY(LMMX),STAT=IERROR)
                IF (IERROR > 0) CALL ALLOCERROR()
                ALLOCATE(ISBORDER(LMMX),STAT=IERROR)
                IF (IERROR > 0) CALL ALLOCERROR()
                ALLOCATE(BASIN_INFLUX(LMMX),STAT=IERROR)
                IF (IERROR > 0) CALL ALLOCERROR()
            ENDIF
        ENDIF
        IF (MODEL_IMPACT_CRATERING.OR.USE_CRATER_EVENT) THEN
            ALLOCATE(CRATER_EVENT(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            CRATER_EVENT=0.0
        ENDIF
        IF (MODEL_LAVA_FLOWS) THEN
            ALLOCATE(IS_LAVA_COVERED(MX,MY),ACTIVE_LAVA_FLOW(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(LAVA_SOURCE_DIRECTION(MX,MY),ERUPTION_AGE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(ILOC(RMMX),JLOC(RMMX),LAVA_ELAPSED_TIME(RMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(LAVA_THICKNESS(RMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(LAVA_FLOW_PROBABILITY(9,RMMX),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(CUMULATIVE_LAVA_CHANGE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
        ENDIF
        IF (USE_FLOW_VOLUME) THEN
            ALLOCATE(FLOW_VOLUME_CHANGE(MX,MY),STAT = IERROR)
                       IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(FLOW_VOLUME_FLUX(MX,MY),STAT = IERROR)
                       IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(CUMULATIVE_FLOW_VOLUME_FLUX(MX,MY),STAT = IERROR)
                       IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(FLOW_VOLUME_EROSION_BASE(MX,MY),STAT = IERROR)
                       IF (IERROR > 0) CALL ALLOCERROR()
            ALLOCATE(CUMULATIVE_FLOW_VOLUME_CHANGE(MX,MY), STAT = IERROR)
                        IF (IERROR > 0) CALL ALLOCERROR()
        ENDIF
        IF (MODEL_ACCRETION_AND_ABLATION) THEN
            ALLOCATE(LAPLACE(MX,MY),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()          
        ENDIF    
        VERTICAL_SCALING=VERTICAL_SCALING_FACTOR*CONVERT_TO_METERS
        CELL_SIZE=INPUT_CELL_SIZE*CONVERT_TO_METERS
        READ(PARAMS,*) IVARRATEUSE
        IF (IVARRATEUSE > 0) THEN
            VARIABLE_EROSION_RATE = .TRUE.
            OPEN(INRATES,FILE='INRATES.DAT',ACTION='READ')
        ELSE
            VARIABLE_EROSION_RATE = .FALSE.
        ENDIF
        READ(PARAMS,*) IEVENTS
        DO_EVENTS=.FALSE.
        IF (IEVENTS > 0) THEN
            OPEN(72,FILE='EVENTS.PRM',ACTION='READ')
            READ(72,*) NUMBER_OF_EVENTS,EVENT_TYPE
            ALLOCATE(EVENT_DONE(NUMBER_OF_EVENTS),EVENT_TIMES(NUMBER_OF_EVENTS),EVENT_ITERATIONS(NUMBER_OF_EVENTS),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
            EVENT_DONE=.TRUE.
            IF (NUMBER_OF_EVENTS > 0) THEN
                DO I=1,NUMBER_OF_EVENTS
                    IF (EVENT_TYPE == 0) THEN
                        READ(72,*) EVENT_TIMES(I)
                    ELSE
                        READ(72,*) EVENT_ITERATIONS(I)
                    ENDIF
                    EVENT_DONE(I)=.FALSE.
                ENDDO
                EVENT_INDEX=1
                DO_EVENTS=.TRUE.
                CALL SETUP_EVENTS()
                IF ((.NOT.NEW_SIMULATION).AND.(EVENT_TYPE==0)) THEN
                    IF (PRESENT_TIME>0.0) THEN
                        DO I=1,NUMBER_OF_EVENTS
                            IF (EVENT_TIMES(I)<=PRESENT_TIME) THEN
                                EVENT_INDEX=EVENT_INDEX+1
                            ENDIF
                        ENDDO
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
        CHANNEL_TIMESTEP_SCALING=200.0
        READ(PARAMS,*) DEFAULT_CHANNEL_TIMESTEP
        READ(PARAMS,*) MAXIMUM_TIME_INCREMENT
        READ(PARAMS,*) MINUMUM_TIME_INCREMENT
        READ(PARAMS,*) DO_SEDIMENT_YIELD_SCALING
        USE_SEDIMENT_YIELD_SCALING=.FALSE.
        IF (DO_SEDIMENT_YIELD_SCALING>0) THEN
            USE_SEDIMENT_YIELD_SCALING=.TRUE.
        ENDIF
        READ(PARAMS,*) SEDIMENT_YIELD_TIMESTEP_SCALING
        MASS_WASTING_TIMESTEP_SCALING=100.1
        READ(PARAMS,*) ICENT
        READ(PARAMS,*) JCENT
        READ(PARAMS,*) IWIDTH
        READ(PARAMS,*) JWIDTH
        TIME_INCREMENT=MINIMUM_TIME_INCREMENT
        IWINLOW=ICENT-IWIDTH
        IWINHIGH=ICENT+IWIDTH
        JWINLOW=JCENT-JWIDTH
        JWINHIGH=JCENT+JWIDTH
        IF (IWINLOW < 1) IWINLOW=1
        IF (IWINHIGH > MX) IWINHIGH=MX
        IF (JWINLOW < 1) JWINLOW=1
        IF (JWINHIGH > MY) JWINHIGH=MY
       ALLOCATE(PLOTVALS(IWINLOW:IWINHIGH,JWINLOW:JWINHIGH),STAT=IERROR)
	   IF (IERROR > 0) CALL ALLOCERROR()
        READ(PARAMS,*) PLOT_DIAGNOSTICS

        IXDEBUG=0

        IF (IXDEBUG > 0) THEN
            DO_DEBUGGING=.TRUE.
        ELSE
            DO_DEBUGGING=.FALSE.
        ENDIF
           READ(PARAMS,*) STATISTDO
        READ(PARAMS,*)  SHADE_BORDER_INDEX
        READ(PARAMS,*) FIXED_SUNANGLE_INDEX
        READ(PARAMS,*) SUN_ANGLE_GRADIENT
        IF (STATISTDO > 0) THEN
            DO_MORPHOMETRY=.TRUE.
        ELSE
            DO_MORPHOMETRY=.FALSE.
        ENDIF
		IF (SHADE_BORDER_INDEX>0) THEN
			DO_SHADE_BORDER=.TRUE.
		ELSE
			DO_SHADE_BORDER=.FALSE.
		ENDIF
		IF (FIXED_SUNANGLE_INDEX>0) THEN
			IS_FIXED_SUNANGLE=.TRUE.
		ELSE
			IS_FIXED_SUNANGLE=.FALSE.
		ENDIF
        !        ********************************************************************
        !        inelev.dat is initial elevations
        !        ********************************************************************
        OPEN(INDATA,FILE='INELEV.DAT')
        READ(PARAMS,*) IHORIZONTAL_LOWER_BOUNDARY
        READ(PARAMS,*) BOUNDARY_LOWERING_RATE 
        READ(PARAMS,*) INON_ERODING_LOWER_BOUNDARY
        READ(PARAMS,*) USEYPERIODIC
        READ(PARAMS,*) USEXPERIODIC
        READ(PARAMS,*) USEFLOWBOUND
        IF (INON_ERODING_LOWER_BOUNDARY > 0) THEN
            NON_ERODING_LOWER_BOUNDARY=.TRUE.
        ELSE
            NON_ERODING_LOWER_BOUNDARY=.FALSE.
        ENDIF
        IF (IHORIZONTAL_LOWER_BOUNDARY > 0) THEN
            HORIZONTAL_LOWER_BOUNDARY=.TRUE.
        ELSE
            HORIZONTAL_LOWER_BOUNDARY=.FALSE.
        ENDIF
        IF (USEYPERIODIC > 0) THEN
            IS_Y_PERIODIC=.TRUE.
            MYY=MY
        ELSE
            IS_Y_PERIODIC=.FALSE.
            MYY=MY-1
        ENDIF
        IF (USEXPERIODIC > 0) THEN
            IS_X_PERIODIC=.TRUE.
        ELSE
            IS_X_PERIODIC=.FALSE.
        ENDIF
        IF (USEFLOWBOUND>0) THEN
            DO_FLOW_BOUNDARIES=.TRUE.
            IS_X_PERIODIC=.FALSE.
            IS_Y_PERIODIC=.FALSE.
            HORIZONTAL_LOWER_BOUNDARY=.FALSE.
            NON_ERODING_LOWER_BOUNDARY=.FALSE.
        ELSE
            DO_FLOW_BOUNDARIES=.FALSE.
        ENDIF
        IF (HORIZONTAL_LOWER_BOUNDARY.OR.NON_ERODING_LOWER_BOUNDARY) THEN
            REFLECTIVE_UPPER_BOUNDARY = .TRUE.
        ELSE
            REFLECTIVE_UPPER_BOUNDARY = .FALSE.
        ENDIF
        !        *******************************************************************
        !        ********************** output and recalculation timing *************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) STARTING_ITERATION
        READ(PARAMS,*) ATIME
        IF (NEW_SIMULATION) PRESENT_TIME=ATIME
        IF (NEW_SIMULATION) TOTAL_ITERATIONS = STARTING_ITERATION
        READ(PARAMS,*) MAXIMUM_ITERATION
        READ(PARAMS,*) MAXIMUM_SIMULATION_TIME
        READ(PARAMS,*) ELEVATION_PRINT_INTERVAL
        READ(PARAMS,*) OUTPUT_PRINT_INTERVAL
        READ(PARAMS,*) RECALCULATE_GRADIENT_INTERVAL
        READ(PARAMS,*) RECALCULATE_DISCHARGE_INTERVAL
        READ(PARAMS,*) WRITE_CHANGE_INTERVAL
        READ(PARAMS,*) KWRITE
        IF (KWRITE == 1)CRITICAL_SOURCE_DIVERGENCE=0.0
        IF (KWRITE == 2)CRITICAL_SOURCE_DIVERGENCE=0.1
        IF (KWRITE == 3)CRITICAL_SOURCE_DIVERGENCE=0.2
        IF (KWRITE == 4)CRITICAL_SOURCE_DIVERGENCE=0.4
        IF (KWRITE == 5)CRITICAL_SOURCE_DIVERGENCE=0.8
        IF (KWRITE == 6)CRITICAL_SOURCE_DIVERGENCE=1.6
        READ(PARAMS,*) ONEONLY
        IF (ONEONLY > 0) THEN
            WRITE_ABSOLUTE_ELEVATION = .TRUE.
        ELSE
            WRITE_ABSOLUTE_ELEVATION = .FALSE.
        ENDIF
        READ(PARAMS,*) IMAGE_OUTPUT_INTERVAL
        READ(PARAMS,*) DIVERGEINTERVAL
        READ(PARAMS,*) DO_WRITE_BINARY_IMAGE_FILE
        IF (DO_WRITE_BINARY_IMAGE_FILE>0) THEN
           WRITE_BINARY_IMAGE_FILE = .TRUE.
           OUTCRATERS=277
           OUTSUBMERGE=278
           OPEN(OUT_IMAGE_FILE, FILE='BINARY_ELEVATION.DAT',ACTION='WRITE', ACCESS='STREAM', &
               FORM='UNFORMATTED')!,RECORDTYPE='STREAM',CONVERT='LITTLE_ENDIAN')
           IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                OPEN(OUT_BINARY_DISCHARGES,FILE='BINARY_DISCHARGES.DAT',ACTION='WRITE',ACCESS='STREAM', &
                    FORM='UNFORMATTED')
            ENDIF
            IF (MODEL_IMPACT_CRATERING.OR.DO_EVENTS) THEN
                OPEN(OUTCRATERS,FILE='BINARY_CRATERS.DAT',ACTION='WRITE',ACCESS='STREAM', &
                    FORM='UNFORMATTED')
            ENDIF         
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                OPEN(OUTSUBMERGE,FILE='BINARY_SUBMERGE.DAT',ACTION='WRITE',ACCESS='STREAM', &
                    FORM='UNFORMATTED')
            ENDIF
           ALLOCATE(OUT_32_DATA(MX,MY),STAT=IERROR)
           IF (IERROR > 0) CALL ALLOCERROR()
        ELSE
            WRITE_BINARY_IMAGE_FILE=.FALSE.
        ENDIF
        IF (FLUVIAL_AND_SLOPE_MODELING.AND.(.NOT.WRITE_BINARY_IMAGE_FILE)) THEN 
            OPEN(OUTDISCHARGES,FILE='ROUTED_DISCHARGE.DAT',ACTION='WRITE')
        ENDIF
        CALL READ_FLOW_PARAMETERS
        CALL READ_BEDROCK_CHANNEL_PARAMETERS
		CALL READ_ALLUVIAL_CHANNEL_PARAMETERS
        CALL READ_WEATHERING_PARAMETERS
        CALL READ_MASS_WASTING_PARAMETERS
        !        *******************************************************************
        !        ******** random variability of simulation parameters    ***************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) RANDTHRESHUSE
        READ(PARAMS,*) RANDDISCHUSE
        READ(PARAMS,*) CRITICAL_SHEAR_VARIABILITY
        READ(PARAMS,*) DISCHARGE_COEFF_VARIATION
        READ(PARAMS,*) OMEGA_WEIGHT
        USE_RANDOM_DISCHARGE=.FALSE.
        IF (RANDTHRESHUSE > 0) THEN
            RANDOM_CRITICAL_SHEAR = .TRUE.
        ELSE
            RANDOM_CRITICAL_SHEAR = .FALSE.
            IF (RANDDISCHUSE > 0) THEN
                USE_RANDOM_DISCHARGE=.TRUE.
            ELSE
                USE_RANDOM_DISCHARGE=.FALSE.
            ENDIF
        ENDIF
        IF (OMEGA_WEIGHT < 1.0)    THEN
            OCORRECT = SQRT((EXP(1.0)-1.0)/(EXP(OMEGA_WEIGHT)-1.0))
        ELSE
            OCORRECT = 1.0
        ENDIF
        !        *******************************************************************
        !        **************** rock and surface deformation *********************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) DEFORMUSE
        READ(PARAMS,*) DEFORMSCALE
		CALL READ_GROUNDWATER_PARAMETERS
        CALL READ_EOLIAN_PARAMETERS
        CALL READ_LAVA_PARAMETERS
        CALL READ_CRATER_PARAMETERS
        !       **********************************************************************
        !
        !       Part 5:  Read in ocean parameters
        !
        !       **********************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) IOCEANVAR
        IF (MODEL_OCEAN_LEVEL.AND.(IOCEANVAR > 0)) THEN
            VARIABLE_OCEAN_ELEVATION=.TRUE.
            ALLOCATE(OCEAN_RECALCULATION_TIMES(2000),OCEAN_LEVELS(2000),STAT=IERROR)
            IF (IERROR > 0) CALL ALLOCERROR()
        ELSE
            VARIABLE_OCEAN_ELEVATION=.FALSE.
        ENDIF
        CALL READ_ACCRETION_ABLATION_PARAMETERS
		CALL READ_GRAVEL_TRANSPORT_PARAMETERS
        !        *******************************************************************
        !        **************** masking fluvial and slope erosion ****************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        USE_EROSION_MASK=.FALSE.
        READ(PARAMS,*) EROSION_MASK_USE
        !write(*,241) EROSION_MASK_USE
241 format('%%%%%%%%%%%$$$$$$$$#######  EROSION MASK %%%%%%%%%@@@@@', i6)
        IF (EROSION_MASK_USE>0) THEN
            USE_EROSION_MASK = .TRUE.
            DO_SEDIMENT_TRANSPORT=.FALSE.
            ALLOCATE(EROSION_MASK(MX,MY))
        ELSE
            USE_EROSION_MASK = .FALSE.
        ENDIF
        IF (USE_EROSION_MASK) CALL READ_EROSION_MASK
        IF (.NOT.MODEL_GROUNDWATER) SEEPAGE_WEATHERING=.FALSE.

        ! *************************************************************************
        !
        !    Read in spatial variation
        !
        ! *************************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE             
        READ(PARAMS,*) SPATIAL_VARIATION_USE
        IF (SPATIAL_VARIATION_USE>0) THEN
            DO_SPATIAL_VARIATION=.TRUE.
        ELSE
            DO_SPATIAL_VARIATION=.FALSE.
        ENDIF
        IF (DO_SPATIAL_VARIATION) THEN
            SPATIAL_MAX=-1.0E25
            SPATIAL_MIN=-SPATIAL_MAX
            OPEN(133,FILE='SPATIAL_VARIATION.DAT',ACTION='READ')
            ALLOCATE(SPATIAL_VARIATION(MX,MY),STAT=IERROR)
            IF (IERROR/=0) THEN
                WRITE(*,2046)
2046            FORMAT('ERROR ALLOCATING SPATIAL_VARIATION')
                STOP
            ENDIF
            READ(133,*) XMX,XMY
            IF ((XMX/=MX).OR.(XMY/=MY)) THEN
               WRITE(*,2045) XMX,XMY,MX,MY
2045           FORMAT('INCOMPATIBLE SPATIAL VARIATION INPUT, XMX=',I6,' XMY=',I6,' MX=',I6,' MY=',I6)
               STOP
            ENDIF
            DO I=1,MX
                DO J=1,MY
                    READ(133,*) SPATIAL_VARIATION(I,J)
                    IF (SPATIAL_VARIATION(I,J)>SPATIAL_MAX) SPATIAL_MAX=SPATIAL_VARIATION(I,J)
                    IF (SPATIAL_VARIATION(I,J)<SPATIAL_MIN) SPATIAL_MIN=SPATIAL_VARIATION(I,J)                   
                ENDDO
            ENDDO
            CLOSE(133)
        ENDIF
        ! *************************************************************************
        !
        !    Read in glacial parameters
        !
        ! *************************************************************************
        IF (USE_FLOW_VOLUME) CALL READ_MASS_FLOW_PARAMETERS()        
        !        *******************************************************************
        !        *******************************************************************
        !          w r i t e  o u t  t h e  s i m u l a t i o n  p a r a m e t e r s
        !       !!! this needs work to write out all current simulation parameters
        !        *******************************************************************
        !        *******************************************************************
        !
        !        *******************************************************************
        !        ********************** boundary conditions ************************
        !        *******************************************************************
        WRITE(OUTHIST,500)
        500 FORMAT(' %%%%%%%%%%%%%%% PLANETARY SIMULATION PROGRAM %%%%%%%%%%%%%%'/,&
        ' ******************  BOUNDARY CONDITIONS *******************')
        IF (NEW_SIMULATION) THEN
            WRITE(OUTHIST,487)
487         FORMAT (' NEW START')
         ELSE
            WRITE(OUTHIST,486)
486         FORMAT (' RESTART')
        ENDIF
            IF ((PRESENT_TIME>0).AND.DO_EVENTS) THEN
            WRITE(OUTHIST,2487)
2487 FORMAT(' EVENTS OCCURING BEFORE SIMULATION STARTING TIME ARE NOT REPEATED')            
         ENDIF
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            WRITE(OUTHIST,1487)
            1487 FORMAT(' MODELING FLUVIAL AND SLOPE PROCESSES')
        ENDIF
        IF (MODEL_IMPACT_CRATERING) THEN
            WRITE(OUTHIST,1486)
            1486 FORMAT(' MODELING IMPACT CRATERING')
        ENDIF
        IF (MODEL_LAVA_FLOWS) THEN
            WRITE(OUTHIST,1488)
            1488 FORMAT(' MODELING LAVA FLOWS')
        ENDIF
        IF (MODEL_EOLIAN_CHANGES) THEN
            WRITE(OUTHIST,1489)
            1489 FORMAT(' MODELING EOLIAN EROSION AND DEPOSITION')
        ENDIF
        IF (MODEL_OCEAN_LEVEL) THEN
            WRITE(OUTHIST,1490)
            1490 FORMAT(' MODELING VARIABLE OCEAN ELEVATION')
        ENDIF
        IF (MODEL_ACCRETION_AND_ABLATION) THEN
            WRITE(OUTHIST,1491)
            1491 FORMAT(' MODELING ACCRETION AND ABLATION')
        ENDIF
        IF (MODEL_LAKE_EVAPORATION) THEN
            WRITE(*,8476)
            WRITE(OUTHIST,8476)
            8476 FORMAT(' MODELING LAKE EVAPORATION')
        ENDIF
        IF (DO_AVALANCHE) THEN
            WRITE(*,1492)
            WRITE(OUTHIST,1492)
1492        FORMAT(' MODELING AVALANCHE EROSION')
        ENDIF
        IF (USE_FLOW_VOLUME) THEN
            WRITE(*,1493)
            WRITE(OUTHIST,1493)
1493        FORMAT(' MODELING MASS FLOWS')
        ENDIF
        WRITE(OUTHIST,501) MX, MY, MZ,CELL_SIZE, VERTICAL_SCALING,CONVERT_TO_METERS
        501 FORMAT(' MX=',I5,' MY=',I5, ' MZ=',I5,/,' HORIZONTAL CELL SIZE=',G12.5 &
            ,' VERTICAL SCALING SIZE=',G12.5,/,'CONVERSION TO METERS=',G12.5)
        IF (.NOT.VARIABLE_EROSION_RATE) THEN
            WRITE(OUTHIST,503) BOUNDARY_LOWERING_RATE
            503 FORMAT(' CONSTANT SIMULATION PARAMETERS',/,' EROSION RATE=',G12.5)
        ELSE
            WRITE(OUTHIST,2492)
2492 FORMAT(' USING TIME-VARYING PARAMETERS')
        ENDIF
        IF (DO_EVENTS) THEN
            WRITE(OUTHIST,1564) NUMBER_OF_EVENTS
1564        FORMAT('SIMULATION USES ',I6,' EPISODIC EVENTS')
        ENDIF
        WRITE(OUTHIST,504) CHANNEL_TIMESTEP_SCALING,DEFAULT_CHANNEL_TIMESTEP &
           ,MAXIMUM_TIME_INCREMENT, MINUMUM_TIME_INCREMENT
        504 FORMAT(' FRACTION OF GRADIENT TIMESTEP IN MAX EROSION RATE='&
        ,G12.5,  &
        /,' MAXIMUM CHANNEL TIMESTEP=',G12.5,/, &
        ' MAXIMUM PERMITTED TIMESTEP=',G12.5,' MINIMUM TIMESTEP=',G12.5)
        WRITE(OUTHIST,934) SEDIMENT_YIELD_TIMESTEP_SCALING,MASS_WASTING_TIMESTEP_SCALING
934     FORMAT(' SEDIMENT AND SLOPE SCALING FACTORS FOR ',&
        'MAXIMUM EROSION RATE :' &
        ,/,2G12.5)
        IF (USE_SEDIMENT_YIELD_SCALING) THEN
            WRITE(2493) SEDIMENT_YIELD_TIMESTEP_SCALING
2493        FORMAT(' SEDIMENT YIELD AFFECTS TIME INCREMENT, WITH FACTOR: ',G12.5)
        ENDIF
        IF (DO_DEBUGGING) THEN
            WRITE(OUTHIST,564) ICENT,JCENT
            564 FORMAT(' CENTER OF DEBUG WINDOW - I=',I5,' J=',I5)
        ENDIF
        IF (DO_MORPHOMETRY) THEN
            WRITE(OUTHIST,3492)
3492        FORMAT('MORPHOMETRIC PARAMETERS REPORTED')
        ENDIF
        IF (DO_SHADE_BORDER) THEN
            WRITE(OUTHIST,1565)
1565        FORMAT(' SHADED RELIEF IMAGES HAVE REPEATED PERIODIC BORDERS')
        ENDIF
        IF (IS_FIXED_SUNANGLE) THEN
            WRITE(OUTHIST,1566) SUN_ANGLE_GRADIENT
1566        FORMAT('SHADED RELIEF IMAGES HAVE SUNANGLE OF ',G12.5)
        ENDIF
        IF (NON_ERODING_LOWER_BOUNDARY) THEN
            WRITE(OUTHIST,1495)
            1495 FORMAT(' THE LOWER BOUNDARY IS NOT ACTIVELY ERODED -' &
            ' IT IS CHANGED ONLY IN THE BOUNDARY CONDITIONS ROUTINE')
        ENDIF
        IF (HORIZONTAL_LOWER_BOUNDARY) THEN
            WRITE(OUTHIST,849)
            849   FORMAT(' LOWER BOUNDARY IS FORCED TO BE LEVEL')
        ELSE
            WRITE(OUTHIST,848)
            848   FORMAT(' AN UNEVEN LOWER BOUNDARY IS ALLOWED ')
        ENDIF
       IF (DO_FLOW_BOUNDARIES) THEN
           WRITE(OUTHIST,11490)
            11490 FORMAT(' THE OUTER BOUNDARIES ARE NON-EROSIONAL BUT FLOW AND SEDIMENT MAY EXIT')
        ENDIF
        IF (IS_X_PERIODIC) THEN
            WRITE(OUTHIST,1496)
            1496 FORMAT(' THE LATERAL BOUNDARIES ARE PERIODIC')
        ELSE
            WRITE(OUTHIST,1497)
            1497 FORMAT(' THE LATERAL BOUNDARIES ARE NON-FLUX')
        ENDIF
        IF (IS_Y_PERIODIC) THEN
            WRITE(OUTHIST,1498)
            1498 FORMAT(' THE TOP AND BOTTOM BOUNDARIES ARE PERIODIC')
        ELSE
            WRITE(OUTHIST,1499)
            1499 FORMAT(' THE TOP BOUNDARY IS NON-FLUX')
        ENDIF
        !        *******************************************************************
        !        ********************** output and recalculation timing *************
        !        *******************************************************************
        WRITE(OUTHIST,514) MAXIMUM_ITERATION, MAXIMUM_SIMULATION_TIME
        514 FORMAT(' ***************  OUTPUT AND RECALCULATION TIMING ******'&
        ,/,' TOTAL ITERATIONS=',I8,' MAXIMUM SIMULATION ELAPSED TIME=',G12.5)
        WRITE(OUTHIST,1500) STARTING_ITERATION,PRESENT_TIME
        1500 FORMAT(' STARTING ITERATION NUMBER=',I6,' STARTING TIME=',G12.5)
        WRITE(OUTHIST,515)ELEVATION_PRINT_INTERVAL,OUTPUT_PRINT_INTERVAL
        515 FORMAT(' ELEVATION OUTPUT INTERVAL=',I5, &
        ' PRINT INTERVAL=',I5)
        WRITE(OUTHIST,721)WRITE_CHANGE_INTERVAL
        721 FORMAT(' DIRECTION AND GRADIENT CHANGE PRINT INTERVAL=',I5)
        WRITE(OUTHIST,516)RECALCULATE_GRADIENT_INTERVAL
        516 FORMAT(' GRADIENT REEVALUATE INTERVAL=',I5)
        WRITE(OUTHIST,517) CRITICAL_SOURCE_DIVERGENCE
        517 FORMAT(' MINIMUM DIVERGENCE FOR CHANNEL PLOTS=',G12.5)
        WRITE(OUTHIST,784) IMAGE_OUTPUT_INTERVAL,DIVERGEINTERVAL
784     FORMAT(' INTERVAL FOR PRINTING IMAGE FILE=',I5,/, &
        ' INTERVAL FOR CALCULATING DIVERGENCE=',I5)
        IF (WRITE_BINARY_IMAGE_FILE) THEN
            WRITE(OUTHIST,2784)
2784        FORMAT('CUMULATIVE BINARY FILES ARE PRODUCED')
        ENDIF
       ! WRITE(OUTHIST,1784)
1784 FORMAT(' **************** FLOW ROUTING AND DRAINAGE AREA *************')
        !        *******************************************************************
        !        ******** random variability of simulation parameters    ***************
        !        *******************************************************************
        IF (RANDOM_CRITICAL_SHEAR) THEN
            WRITE(OUTHIST,498) CRITICAL_SHEAR_VARIABILITY,OMEGA_WEIGHT
            498   FORMAT('  *************** RANDOM VARIABILITY *******************'&
            ,/,' TEMPORALLY VARIABLE FLUVIAL EROSION THRESHOLD,' &
            ,' VARIANCE=',  &
            G12.5,' OMEGA=',G12.5)
        ENDIF
        IF (USE_RANDOM_DISCHARGE) THEN
            WRITE(OUTHIST,499) DISCHARGE_COEFF_VARIATION,OMEGA_WEIGHT
            499   FORMAT('  *************** RANDOM VARIABILITY *******************'&
            ,/,' TEMPORALLY VARIABLE FLUVIAL DISCHARGE, VARIANCE=',&
            G12.5,' OMEGA=',G12.5)
        ENDIF
        !        *******************************************************************
        !        **************** rock and surface deformation *********************
        !        *******************************************************************
        IF (DEFORMUSE > 0) THEN
            DO_ROCK_DEFORMATION=.TRUE.
            WRITE(OUTHIST,776) DEFORMSCALE
            776   FORMAT(' ********* ROCK AND SURFACE DEFORMATION ****************' &
            ,/,' ROCK AND SURFACE DEFORMATION IS BEING USED WITH',/,  &
            ' DEFORMATION RATE SCALING OF ',G12.5,' WITH UPLIFT POSITIVE')
        ELSE
            DO_ROCK_DEFORMATION=.FALSE.
        ENDIF 
        ! **************************************************************************
        !    Ocean parameters
        ! **************************************************************************
        IF (VARIABLE_OCEAN_ELEVATION) THEN
            WRITE(OUTHIST,1560)
            1560 FORMAT('******************MODELING TIME-VARYING OCEAN LEVELS*************')
        ENDIF
        ! **************************************************************************
        IF (USE_EROSION_MASK) THEN
            WRITE(OUTHIST,1561)
1561        FORMAT('AN EROSION MASK IS USED TO ISOLATE AREAS SUBJECT TO EROSIONAL MODIFICATION')
        ENDIF
        !        *******************************************************************
        !         Set up initial rock resistance to unity and regolith thickness to
        !           initial_regolith_thickness
        !        *******************************************************************
            DO  J = 1, MY
                DO  I = 1, MX
                    IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                        IF (RESISTANT_SURFACE_LAYER) THEN
                            RELATIVE_RESISTANCE(I,J)=SURFACE_LAYER_RESISTANCE
                        ELSE
                            RELATIVE_RESISTANCE(I,J)=1.0
                        ENDIF
                    ENDIF
                    IF (NEW_SIMULATION.AND.FLUVIAL_AND_SLOPE_MODELING) REGOLITH(I,J)=INITIAL_REGOLITH_THICKNESS
                ENDDO
            ENDDO
        THEMINIMUM=1.0E+15
        !        *******************************************************************
        !         Set initial sediment yield to zero, is_sediment_covered to false,
        !           is_rock_surface to false, accelerated_erosion to false, and last erosion rate to zero
        !        *******************************************************************
        IF (NEW_SIMULATION) THEN
                READ(INDATA,*)MX1,MY1
            IF ((MX1 /= MX).OR.(MY1 /= MY)) THEN
                WRITE(*,1121) MX, MX1,MY,MY1
                1121      FORMAT(' INCOMPATIBLE INPUT DIMENSIONS, MX=',I5,' MX1=',I5, &
                '  MY=',I5,' MY1=',I5)
                STOP
            ENDIF
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                SEDIMENT_YIELD=0.0
                ACCELERATED_EROSION=.FALSE.
                DO_ACCELERATED_EROSION=.FALSE.
                PREVIOUS_ELEVATION=0.0
                DEFORMATION=0.0
            ENDIF
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                IS_SEDIMENT_COVERED=.FALSE.
            ENDIF
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                IS_ROCK_SURFACE=.FALSE.
            ENDIF
        ELSE
            OPEN(645,FILE='FINAL_BEDROCK.DAT',ACTION='READ')
            CALL READ_BEDROCK_LOCATIONS(645)
            CLOSE(645)
            OPEN(645,FILE='FINAL_SEDCOVER.DAT',ACTION='READ')
            CALL READ_ALLUVIAL_LOCATIONS(645)
            CLOSE(645)
            OPEN(645,FILE='FINAL_SEDBASE.DAT',ACTION='READ')
            CALL READ_SEDIMENT_BASE(645)
            CLOSE(645)
            OPEN(645,FILE='FINAL_REGOLITH.DAT',ACTION='READ')
            CALL READ_REGOLITH_THICKNESS(645)
            CLOSE(645)
            OPEN(645,FILE='FINAL_ELEVATION.DAT',ACTION='READ')
            CALL READ_ELEVATIONS(645)
            CLOSE(645)
            OPEN(645,FILE='FINAL_TIME.DAT',ACTION='READ')
            READ(645,*) TOTAL_ITERATIONS, PRESENT_TIME,EREFERENCE           
        ENDIF
         !       **********************************************************************
        !
        !       ***** More variable initialization
        !
        !       **********************************************************************
        NOISECALL=.TRUE.
        DIFFUSECALL=.TRUE.
        !        *******************************************************************
        !         Read in the initial elevations and determine the minimum elevation
        !        *******************************************************************
            IF (NEW_SIMULATION) EREFERENCE=-1.0E25
            IF (CRUSTUSE.GT.1) THEN
                OPEN(83,FILE='SURFACE_CRUST.DAT',ACTION='READ')
                READ(83,*) MX1,MY1
                IF ((MX1.NE.MX).OR.(MY1.NE.MY)) THEN
                    WRITE(*,3756)
3756                 FORMAT('SURFACE CRUST FILE HAS INCOMPATIBLE DIMENSIONS')
                     STOP
                ENDIF
            ENDIF
            DO  I = 1, MX
                DO  J = 1, MY
                    IF (NEW_SIMULATION) THEN
                        READ (INDATA,*) ELEVATION(I,J)
                        ELEVATION(I,J)=ELEVATION(I,J)*VERTICAL_SCALING
                    ENDIF
                    INITIAL_ELEVATION(I,J)=ELEVATION(I,J)
                    IF(NEW_SIMULATION) THEN
                        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                            SEDIMENT_BASE(I,J)=ELEVATION(I,J)
                        ENDIF
                    ENDIF
                    IF (FLUVIAL_AND_SLOPE_MODELING.AND.RESISTANT_SURFACE_LAYER) THEN
                        IF (CRUSTUSE.GT.1) THEN
                            READ(83,*) TEMP1
                            PREVIOUS_ELEVATION(I,J)=ELEVATION(I,J)-TEMP1
                        ELSE
                           PREVIOUS_ELEVATION(I,J)=ELEVATION(I,J)-SURFACE_LAYER_THICKNESS
                        ENDIF
                    ENDIF
                    IF (NEW_SIMULATION) THEN
                        IF (ELEVATION(I,J) > EREFERENCE) EREFERENCE=ELEVATION(I,J)
                    ENDIF
                    IF (ELEVATION(I,J) < THEMINIMUM) THEMINIMUM=ELEVATION(I,J)
                    IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                        EROSION_DEPTH_INDEX(I,J)=-1
                    ENDIF
                ENDDO
            ENDDO
            IF (CRUSTUSE.GT.1) CLOSE(83)
            WRITE(*,555) READREGOLITH
555         FORMAT('READREGOLITH STATUS=',I5)
            IF (READREGOLITH > 0) THEN
                WRITE(*,556)
556             FORMAT('CALLING READ REGOLITH')
                OPEN(IOTEMP3,FILE='INREG.DAT',ACTION='read')
                CALL READ_REGOLITH_THICKNESS(IOTEMP3)
                CLOSE(IOTEMP3)
            ENDIF
        OPEN(1488,FILE='INITIAL_ELEVATION.DAT',ACTION='WRITE')
        WRITE(1488,4487) MX,MY
4487     FORMAT(I9,' ',I9)
        DO I=1,MX
            DO J=1,MY
                WRITE(1488,4486) ELEVATION(I,J)
            ENDDO
        ENDDO
4486     FORMAT(G13.6)
        CLOSE(1488)
        WRITE(OUTHIST,221)
        221   FORMAT(' INITIAL ELEVATION')
        CALL SUMMARIZE_MATRIX_DATA(ELEVATION,THEAVG,THEMAX,THEMIN)
        IF (IREADALLUV == 1) THEN
            WRITE(OUTHIST,222)
            222       FORMAT(' INITIAL ALLUVIAL BASE')
            CALL SUMMARIZE_MATRIX_DATA(SEDIMENT_BASE,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,223)
            223       FORMAT(' INITIAL SEDIMENT COVER')
            CALL SUMMARIZE_LOGICAL_MATRIX(IS_SEDIMENT_COVERED)
        ENDIF
        !       ********************************************************************
        !        If we use temporally variable erosion rates, read in a table.
        !           The first line (number_of_parameter_changes) tells how many erosion rate intervals
        !           there are to be read.  For each rate interval the closing time for
        !           that erosion rate catagory and the erosion rate during that interval
        !           are read in, as well as the potentially time-varying simulation
        !           parameters (see the subroutine doeroderate for definition of the
        !           individual simulation_parameters values)
        !       ********************************************************************
        IF (VARIABLE_EROSION_RATE) THEN
            READ(INRATES,*) NUMBER_OF_PARAMETER_CHANGES
            WRITE(OUTHIST,738) NUMBER_OF_PARAMETER_CHANGES
            738   FORMAT(' NUMBER OF EROSION RATE CATEGORIES=',I5,/, &
            'TABLE OF LOWEST ELEVATIONS AND THEIR EROSION RATES',/)
            DO  I=1,NUMBER_OF_PARAMETER_CHANGES
                READ(INRATES,*) TIMES_FOR_PARAMETER_CHANGES(I),EROSION_RATE_VALUES(I)
                READ(INRATES,*)(SIMULATION_PARAMETERS(I,J),J=1,12)
                WRITE(OUTHIST,737) TIMES_FOR_PARAMETER_CHANGES(I),EROSION_RATE_VALUES(I)
                WRITE(OUTHIST,740) (SIMULATION_PARAMETERS(I,J),J=1,12)
                737   FORMAT (2G12.5)
                740   FORMAT(5G12.5,/,5G12.5,/,2G12.5)
            ENDDO
        ENDIF
        IF (.NOT.IS_Y_PERIODIC) THEN
            DO  I=1,MX
                IW=I-1
                IE=I+1
                IF (IS_X_PERIODIC) THEN
                    IF (IE > MX) IE=1
                    IF (IW < 1) IW=MX
                ELSE
                    IF (IE > MX) IE=MX
                    IF (IW < 1) IW=1
                ENDIF
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                    IF (HORIZONTAL_LOWER_BOUNDARY.OR.DO_FLOW_BOUNDARIES) THEN
                        ERODING_LOWER_BOUNDARY(I)=.FALSE.
                    ELSE
                        IF (NON_ERODING_LOWER_BOUNDARY) THEN
                            ERODING_LOWER_BOUNDARY(I)=.FALSE.
                        ELSE
                            IF ((ELEVATION(I,MY) <= ELEVATION(IW,MY)).AND.(ELEVATION(I,MY) <= ELEVATION(IE,MY))) &
                            THEN
                                ERODING_LOWER_BOUNDARY(I)=.FALSE.
                            ELSE
                                ERODING_LOWER_BOUNDARY(I)=.TRUE.
                            ENDIF
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO
        ENDIF
        IF (DO_SPATIAL_VARIATION) THEN
            WRITE(OUTHIST,8742)
8742        FORMAT('SPATIAL VARIATION OF DISCHARGE AND/OR WEATHERING BEING USED')
        ENDIF
		CLOSE(92)
        RETURN
    END !     SUBROUTINE READ_INPUT_PARAMETERS
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE INITIALIZE_VARIABLES()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,K,IFLAG
        REAL (4) :: XXXX,XRATE,XINC,XGRAD,RAPID_CREEP,ROCK_MASS_WASTING,SUM
        EXTERNAL RAPID_CREEP
        EXTERNAL ROCK_MASS_WASTING
        !      *******************************************************************
        !       Set up the matrix im which determines relative x,y locations of
        !       surrounding points.  The matrix id indicates local drainage
        !       direction according to the following scheme, where 1 is present
        !       location and entries are cells towards which flow is directed
        !                         relative i
        !                         -1   0   1
        !                          ---------
        !                     -1 | 7   4   8  |
        !        relative j    0 | 3   1   5  | values of the matrix flow_direction(i,j)
        !                      1 | 6   2   9  |
        !
        !      flow_direction(i,j) is given a negative sign if the d8_gradient(i,j)<=0.0 (a depression)
        !
        !      flow_direction(i,j)=1 specifies an exit point
        !
        !       The matrix downstream(k,l) takes the value of flow_direction(i,j) as the k argument
        !       and gives relative i value of the downstream location for l=1 and relative j for l=2
        !  MODIFIES: DOWNSTREAM, WEIGHTS, MP, SQRTOFTWO, FORCE_SEDIMENT_CONSERVATION, PREVIOUS_TIME_INCREMENT
        !            CUMULEXCESS, DEPOSITWORK, ERODEWORK, ERODEGRAV, DEPOSITGRAV, CRATERGRAV
        !            CRATERWORK, LAVAWORK, LAVADEPTH, EOLIANDEPTH, ELAPSEDTIME, SLOPEWORK, SLOPEGRAV
        !            LASTSEDBIAS, SEDBIAS, SEDIMENT_DEPOSITED, CUMULATIVE_ELEVATION_CHANGE,
        !            CUMULATIVE_EOLIAN_CHANGE, CUMULATIVE_LAVA_CHANGE, CUMULATIVE_CRATERING_CHANGE
        !            FLEN, TIME_INCREMENT, ITERATION, JABORTMIN, IABORTMIN, IABORTMAX, JABORTMAX
        !            CROSS_WEIGHTING, DIAGONAL_WEIGHTING, ONE_SIXTH, MAXCRIT, NCRITS, CELL_AREA
        !            PREVIOUS_DISCHARGE_COEFFICIENT, GRADMAX,FAILMAX,GRADCAT,RATECAT,DELRATE
        !            RGRADCAT, RRATECAT, RDELRATE, RFAILMAX, STICKYFACTOR, STICKYPROB, STICKYCAT, STICKYDEL
        !            PRESENT_TIME, INUMLEVELS, OCEANNEXTTIME, PRESENT_TIME
        ! CALLS: FINDOCEAN_ELEVATION

        !      *******************************************************************
        DOWNSTREAM(1,1) =  0
        DOWNSTREAM(1,2) =  0
        DOWNSTREAM(2,1) =  0
        DOWNSTREAM(2,2) =  1
        DOWNSTREAM(3,1) = -1
        DOWNSTREAM(3,2) =  0
        DOWNSTREAM(4,1) =  0
        DOWNSTREAM(4,2) = -1
        DOWNSTREAM(5,1) =  1
        DOWNSTREAM(5,2) =  0
        DOWNSTREAM(6,1) = -1
        DOWNSTREAM(6,2) =  1
        DOWNSTREAM(7,1) = -1
        DOWNSTREAM(7,2) = -1
        DOWNSTREAM(8,1) =  1
        DOWNSTREAM(8,2) = -1
        DOWNSTREAM(9,1) =  1
        DOWNSTREAM(9,2) =  1
        MP = MX * MY
        SQRTOFTWO = SQRT(2.0)
        !       weights for printing shaded-relief images
        WEIGHTS(1,1,1)=4.0/4.0
        WEIGHTS(1,1,2)=0.0/4.0
        WEIGHTS(1,1,3)=0.0/4.0
        WEIGHTS(1,1,4)=0.0/4.0
        WEIGHTS(2,1,1)=2.0/4.0
        WEIGHTS(2,1,2)=2.0/4.0
        WEIGHTS(2,1,3)=0.0/4.0
        WEIGHTS(2,1,4)=0.0/4.0
        WEIGHTS(1,2,1)=2.0/4.0
        WEIGHTS(1,2,2)=0.0/4.0
        WEIGHTS(1,2,3)=2.0/4.0
        WEIGHTS(1,2,4)=0.0/4.0
        WEIGHTS(2,2,1)=1.0/4.0
        WEIGHTS(2,2,2)=1.0/4.0
        WEIGHTS(2,2,3)=1.0/4.0
        WEIGHTS(2,2,4)=1.0/4.0
     weights_3(1,1)=1.0
     weights_3(1,3)=1.0
     weights_3(3,1)=1.0
     weights_3(3,3)=1.0
     weights_3(1,2)=4.0
     weights_3(2,1)=4.0
     weights_3(2,3)=4.0
     weights_3(3,2)=4.0
     weights_3(2,2)=-20.0
     weights_5(1,1)=0.0
     weights_5(1,2)=0.0
     weights_5(1,4)=0.0 
     weights_5(1,5)=0.0 
     weights_5(2,1)=0.0
     weights_5(2,5)=0.0 
     weights_5(4,1)=0.0 
     weights_5(4,5)=0.0 
     weights_5(5,1)=0.0
     weights_5(5,2)=0.0 
     weights_5(5,4)=0.0
     weights_5(5,5)=0.0
     weights_5(1,3)=1.0 
     weights_5(2,2)=1.0 
     weights_5(2,4)=1.0
     weights_5(3,1)=1.0 
     weights_5(3,5)=1.0 
     weights_5(4,2)=1.0 
     weights_5(4,4)=1.0
     weights_5(5,3)=1.0
     weights_5(2,3)=2.0
     weights_5(3,2)=2.0
     weights_5(4,3)=2.0 
     weights_5(3,4)=2.0 
     weights_5(3,3)=-16.0
     weights_7(1,1)=0.0 
     weights_7(1,2)=0.0 
     weights_7(1,6)=0.0 
     weights_7(1,7)=0.0
     weights_7(2,1)=0.0
     weights_7(2,7)=0.0 
     weights_7(6,1)=0.0 
     weights_7(6,7)=0.0 
     weights_7(7,1)=0.0
     weights_7(7,2)=0.0
     weights_7(7,6)=0.0
     weights_7(7,7)=0.0
     weights_7(3,3)=0.0
     weights_7(3,5)=0.0
     weights_7(5,3)=0.0
     weights_7(5,5)=0.0
     weights_7(1,3)=1.0
     weights_7(1,4)=1.0
     weights_7(1,5)=1.0
     weights_7(2,2)=1.0
     weights_7(2,6)=1.0
     weights_7(3,1)=1.0
     weights_7(3,7)=1.0
     weights_7(4,1)=1.0
     weights_7(4,7)=1.0
     weights_7(5,1)=1.0
     weights_7(5,7)=1.0
     weights_7(6,2)=1.0
     weights_7(6,6)=1.0
     weights_7(7,3)=1.0
     weights_7(7,4)=1.0
     weights_7(7,5)=1.0
     weights_7(2,3)=3.0
     weights_7(2,4)=3.0
     weights_7(2,5)=3.0
     weights_7(3,2)=3.0
     weights_7(3,6)=3.0
     weights_7(4,2)=3.0
     weights_7(4,6)=3.0
     weights_7(5,2)=3.0
     weights_7(5,6)=3.0
     weights_7(6,3)=3.0
     weights_7(6,4)=3.0
     weights_7(6,5)=3.0
     weights_7(3,4)=-7.0
     weights_7(4,3)=-7.0
     weights_7(4,5)=-7.0
     weights_7(5,4)=-7.0
     weights_7(4,4)=-24.0
     weights_9(1,1)=0.0
     weights_9(1,2)=0.0
     weights_9(1,8)=0.0
     weights_9(1,9)=0.0
     weights_9(2,1)=0.0
     weights_9(2,9)=0.0
     weights_9(3,5)=0.0
     weights_9(5,3)=0.0
     weights_9(5,7)=0.0
     weights_9(7,5)=0.0
     weights_9(8,1)=0.0
     weights_9(8,9)=0.0
     weights_9(9,1)=0.0
     weights_9(9,2)=0.0
     weights_9(9,8)=0.0
     weights_9(9,9)=0.0
     weights_9(1,4)=2.0
     weights_9(1,5)=2.0
     weights_9(1,6)=2.0
     weights_9(2,2)=2.0
     weights_9(2,8)=2.0
     weights_9(4,1)=2.0
     weights_9(4,9)=2.0
     weights_9(5,1)=2.0
     weights_9(5,9)=2.0
     weights_9(6,1)=2.0
     weights_9(6,9)=2.0
     weights_9(8,2)=2.0
     weights_9(8,8)=2.0
     weights_9(9,4)=2.0
     weights_9(9,5)=2.0
     weights_9(9,6)=2.0
     weights_9(1,3)=3.0
     weights_9(1,7)=3.0 
     weights_9(2,3)=3.0 
     weights_9(2,7)=3.0     
     weights_9(3,1)=3.0
     weights_9(3,2)=3.0 
     weights_9(3,4)=3.0 
     weights_9(3,6)=3.0 
     weights_9(3,8)=3.0
     weights_9(3,9)=3.0 
     weights_9(4,3)=2.0 
     weights_9(4,7)=2.0 
     weights_9(6,3)=2.0
     weights_9(6,7)=2.0
     weights_9(7,1)=3.0
     weights_9(7,2)=3.0 
     weights_9(7,4)=3.0 
     weights_9(7,6)=3.0 
     weights_9(7,8)=3.0
     weights_9(7,9)=3.0 
     weights_9(8,3)=3.0
     weights_9(8,7)=3.0 
     weights_9(9,3)=3.0 
     weights_9(9,7)=3.0
     weights_9(2,4)=5.0
     weights_9(2,5)=5.0
     weights_9(2,6)=5.0
     weights_9(8,4)=5.0
     weights_9(8,5)=5.0
     weights_9(8,6)=5.0
     weights_9(3,3)=5.0
     weights_9(3,7)=5.0
     weights_9(7,3)=5.0
     weights_9(7,7)=5.0
     weights_9(4,2)=5.0
     weights_9(4,8)=5.0
     weights_9(5,2)=5.0
     weights_9(5,8)=5.0
     weights_9(6,2)=5.0
     weights_9(6,8)=5.0
     weights_9(4,4)=-12.0
     weights_9(4,6)=-12.0
     weights_9(6,4)=-12.0
     weights_9(6,6)=-12.0
     weights_9(4,5)=-23.0
     weights_9(5,4)=-23.0
     weights_9(5,6)=-23.0
     weights_9(6,5)=-23.0
     weights_9(5,5)=-40.0      
        L1100: DO I=1, 2
            L1101: DO J= 1, 2
                SUM = 0.0
                L1110: DO K= 1,4
                    SUM = SUM+WEIGHTS(I,J,K)
                ENDDO L1110
                IF (SUM /= 1.0) THEN
                    WRITE(*,1120) I,J
                    1120              FORMAT(' BAD WEIGHTS AT I=',I5,' J=',I5)
                    STOP
                ENDIF
            ENDDO L1101
        ENDDO L1100
        FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF ((.NOT.IS_X_PERIODIC).OR.(.NOT.IS_Y_PERIODIC)) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (.NOT.DO_FLUVIAL_DETACHMENT) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (.NOT.DO_SEDIMENT_TRANSPORT) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (DO_ROCK_DEFORMATION) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (MODEL_EOLIAN_CHANGES) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (MODEL_LAVA_FLOWS) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (HAVE_INFLUENT_RIVERS) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (.NOT.DO_ALLUVIAL_SMOOTHING) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        RSUFFIX='.RAW'
        XSUFFIX='.PGM'
        FIRST_IMAGE=.TRUE.
        RPREFIX='BSHADE'
        GRADCUT=0.0001
        XDIRECT='RESULTS'
        X1PREFIX='UNKN'
        X2PREFIX='DQS_'
        X3PREFIX='BDQS'
        IJDEBUG=0
            IFILE1=48
            IFILE2=48
            IFILE3=48
            IFILE4=48
            XFILE1=48
            XFILE2=48
            XFILE3=48
            XFILE4=48
            RNUM1=CHAR(IFILE1)
            RNUM2=CHAR(IFILE2)
            RNUM3=CHAR(IFILE3)
            RNUM4=CHAR(IFILE4)
            PREVIOUS_TIME_INCREMENT=0.0
            CUMULEXCESS=0.0
            DEPOSITWORK=0.0
            ERODEWORK=0.0
            ERODEGRAV=0.0
            DEPOSITGRAV=0.0
            CRATERGRAV=0.0
            CRATERWORK=0.0
            LAVAWORK=0.0
            LAVADEPTH=0.0
            EOLIANDEPTH=0.0
            ELAPSEDTIME=0.0
            SLOPEWORK=0.0
            SLOPEGRAV=0.0
            FLOWWORK=0.0
            FLOWGRAV=0.0
            CUM_DEPOSIT_WORK=0.0
            CUM_ERODE_WORK=0.0
            CUM_CRATER_WORK=0.0
            CUM_LAVA_WORK=0.0
            CUM_SLOPE_WORK=0.0
            CUM_FLOW_WORK=0.0
            CUM_DEPOSIT_GRAV=0.0
            CUM_ERODE_GRAV=0.0
            CUM_CRATER_GRAV=0.0
            CUM_LAVA_GRAV=0.0
            CUM_SLOPE_GRAV=0.0
            CUM_FLOW_GRAV=0.0
            CUM_ACCRETION=0.0            
            LASTSEDBIAS=0.0
            SEDBIAS=1.0
            IF (DO_AVALANCHE) THEN
                AVALANCHE_FLUX=0.0
                OLD_AVALANCHE_FLUX=0.0
            ENDIF
            CRATERING_CALLED=.FALSE.
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                SEDIMENT_DEPOSITED=.FALSE.
            ENDIF
            L355: DO  J=1,MYY
                L356: DO  I=1,MX
                    CUMULATIVE_ELEVATION_CHANGE(I,J)=0.0
                    IF (MODEL_EOLIAN_CHANGES) CUMULATIVE_EOLIAN_CHANGE(I,J)=0.0
                    IF (MODEL_LAVA_FLOWS) CUMULATIVE_LAVA_CHANGE(I,J)=0.0
                    IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                        CUMULATIVE_CRATERING_CHANGE(I,J)=0.0
                        CUMULATIVE_WEATHERING(I,J)=0.0
                        CUMULATIVE_EJECTA_DEPOSITION(I,J)=0.0
                        CUMULATIVE_CRATER_EXCAVATION(I,J)=0.0
                        CUMULATIVE_SEDIMENT_DEPOSITION(I,J)=0.0
                        CUMULATIVE_MASS_WASTING(I,J)=0.0
                        CUMULATIVE_FLUVIAL_EROSION(I,J)=0.0
                    ENDIF
                    IF (DO_AVALANCHE) CUMULATIVE_AVALANCHE_EROSION(I,J)=0.0
                    if (USE_FLOW_VOLUME) CUMULATIVE_FLOW_VOLUME_CHANGE(I,J) = 0.0
                ENDDO L356
            ENDDO L355
        !      *******************************************************************
        !       flen(k)is the distance to the surrounding matrix point for values
        !       of the argument k, which is the appropriate value of flow_direction(i,j)
        !      *******************************************************************
        L100: DO I=1,5
            FLEN(I)=1.0
        ENDDO L100
        L110: DO I=6,9
            FLEN(I)=1.0/SQRTOFTWO
        ENDDO L110
        NMAX = MAX(MX,MY)
        LMAX = NMAX*NMAX
        TIME_INCREMENT = 1.0E-5
        IF (.NOT.FLUVIAL_AND_SLOPE_MODELING)TIME_INCREMENT=1.0
            ITERATION = 0
        !      *******************************************************************
        !       Sets range of i,j for debugging printout
        !      *******************************************************************
        JABORTMIN=NMAX+1
        IABORTMIN=NMAX+1
        IABORTMAX=0
        JABORTMAX=0
        !      *******************************************************************
        !       cross_weighting, diagonal_weighting, and one_sixth are used in divergence
        !       calculation
        !      *******************************************************************
        CROSS_WEIGHTING = 2.0/(3.0*CELL_SIZE)
        DIAGONAL_WEIGHTING = 1.0/(6.0*CELL_SIZE)
        ONE_SIXTH = 1.0/6.0
        !      *******************************************************************
        !       Maxcrit is used in writnew.f90 to determine the channel network
        !      *******************************************************************
        MAXCRIT = 5.0
        !      *******************************************************************
        !       Ncrits is also used in writnew.f90
        !      *******************************************************************
        NCRITS = 50
        CELL_AREA=CELL_SIZE*CELL_SIZE
        PREVIOUS_DISCHARGE_COEFFICIENT=DISCHARGE_COEFFICIENT
        !      ********************************************************************
        !       this sets up a table-based look-up function (rapid_creep) for use when there is
        !         a critical value for slope failure - this avoids having to take powers
        !         at each iteration and each cell direction
        !      ********************************************************************
        IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
            GRADMAX=CRITICAL_SLOPE_GRADIENT*(1.0-MAXIMUM_DIFFUSIVITY_INCREASE)**(1.0/SLOPE_GRADIENT_EXPONENT)
            FAILMAX=SLOPE_FAILURE_DIFFUSIVITY*(1.0/MAXIMUM_DIFFUSIVITY_INCREASE-1.0) 
            WRITE(OUTHIST,135) GRADMAX,FAILMAX
            135 FORMAT(' FAILURE GRADIENT=',G12.5,' MAXIMUM RATE=',G12.5)
            L130: DO I=1,101
                GRADCAT(I)=(I-1)*GRADMAX/100.0
                IF (USE_ROERING_MASS_WASTING) THEN
                    RATECAT(I)=SLOPE_FAILURE_DIFFUSIVITY*(1.0/(1.0-CRITICAL_GRADIENT_TERM &
                    *GRADCAT(I)**SLOPE_GRADIENT_EXPONENT))
                ELSE
                    RATECAT(I)=SLOPE_FAILURE_DIFFUSIVITY*(1.0/(1.0-CRITICAL_GRADIENT_TERM &
                    *GRADCAT(I)**SLOPE_GRADIENT_EXPONENT)-1.0)
                ENDIF
                WRITE(OUTHIST,150)GRADCAT(I),RATECAT(I)
            ENDDO L130
            L120: DO  I=1,100
                DELRATE(I)=(RATECAT(I+1)-RATECAT(I))/(GRADCAT(I+1)-GRADCAT(I))
            ENDDO L120
            DELRATE(101)=0.0
            WRITE(OUTHIST,160)
            160 FORMAT(' TABLE OF GRADIENTS AND RATES')
            XXXX=1.5*GRADMAX
            XINC=XXXX/75.0
            XGRAD=0.0
            L140: DO I=1,75
                XRATE=RAPID_CREEP(XGRAD)
                XGRAD=XGRAD+XINC
                WRITE(OUTHIST,150) XGRAD,XRATE
                150 FORMAT(2G12.5)
            ENDDO L140
        ENDIF
        !      ********************************************************************
        !       This sets up a table-based look-up function (rapid_mass_wasting) for use when there is
        !         a critical value for bedrock slope failure - this avoids having to take powers
        !         at each iteration and each cell direction
        !      ********************************************************************
        IF (CRITICAL_BEDROCK_GRADIENT > 0.0) THEN
            XXXX=1.0/CRITICAL_BEDROCK_GRADIENT**SLOPE_GRADIENT_EXPONENT
            RGRADMAX=CRITICAL_BEDROCK_GRADIENT*(1.0-MAXIMUM_DIFFUSIVITY_INCREASE)**(1.0/SLOPE_GRADIENT_EXPONENT)
            FAILMAX=SLOPE_FAILURE_DIFFUSIVITY*(1.0/MAXIMUM_DIFFUSIVITY_INCREASE-1.0)
            WRITE(OUTHIST,148) RGRADMAX,FAILMAX
            148 FORMAT(' ROCK FAILURE GRADIENT=',G12.5,' MAXIMUM RATE=',G12.5)
            L180: DO I=1,101
                RGRADCAT(I)=(I-1)*RGRADMAX/100.0
                IF (USE_ROERING_MASS_WASTING)THEN
                    RRATECAT(I)=SLOPE_FAILURE_DIFFUSIVITY*(1.0/(1.0-XXXX &
                    *RGRADCAT(I)**SLOPE_GRADIENT_EXPONENT))
                ELSE
                    RRATECAT(I)=SLOPE_FAILURE_DIFFUSIVITY*(1.0/(1.0-XXXX &
                    *RGRADCAT(I)**SLOPE_GRADIENT_EXPONENT)-1.0)
                ENDIF
            ENDDO L180
            L190: DO I=1,100
                RDELRATE(I)=(RRATECAT(I+1)-RRATECAT(I))/ &
                (RGRADCAT(I+1)-RGRADCAT(I))
            ENDDO L190
            RDELRATE(101)=0.0
        ENDIF
        !     **********************************************************************
        !      Determines the lookup table for probability of drainage direction
        !         change on alluvial surfaces if sticky_sediment_routing is true.  Probabilities
        !         depend upon the parameter sticky_routing_critical_value
        !     **********************************************************************
        IF  (STICKY_SEDIMENT_ROUTING) THEN
            STICKYFACTOR = - LOG(0.5E0)/(0.5*(STICKY_ROUTING_CRITICAL_VALUE-1.0))
            L201: DO I=1,11
                STICKYPROB(I)=1.0-EXP(-0.5*(I-1))
                STICKYCAT(I)=1.0+(I-1)/STICKYFACTOR
            ENDDO L201
            L202: DO I=1,10
                STICKYDEL(I)=(STICKYPROB(I+1)-STICKYPROB(I))/ &
                (STICKYCAT(I+1) -STICKYCAT(I))
            ENDDO L202
        ENDIF
        !     **********************************************************************
        !      If a variable ocean base level is used, this reads in the levels and times for recalculation
        !     **********************************************************************
        IF (VARIABLE_OCEAN_ELEVATION) THEN
            OPEN(67,FILE='OCEANLEVELS.DAT',action='read')
            I=1
            L6430: DO
                READ(67,*,IOSTAT=IFLAG) OCEAN_RECALCULATION_TIMES(I),OCEAN_LEVELS(I)
                IF (IFLAG < 0) EXIT L6430
                I=I+1
            ENDDO L6430
            INUMLEVELS=I
            OCEANNEXTTIME=0.0
            WRITE(*,7440) INUMLEVELS
            7440     FORMAT(' NO OF OCEAN_RECALCULATION_TIMES=',I6)
            CLOSE(67)
            CALL FINDOCEAN_ELEVATION(PRESENT_TIME)
        ENDIF
        RETURN
    END ! SUBROUTINE INITIALIZE_VARIABLES

    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (4) FUNCTION NORMAL_RANDOM_DEVIATE()
        IMPLICIT NONE
        !     ********************************************************************
        !      Generates a normally disributed random deviate
        !     ********************************************************************
        REAL (4) :: V1, V2, SS
        REAL (4) :: TEMP1,TEMP2
        REAL (4) :: TEMP3, RRAND
        EXTERNAL RRAND
        L10: DO
            TEMP1=RRAND()
            TEMP2=RRAND()
            V1 = 2.0*TEMP1-1.0
            V2 = 2.0*TEMP2-1.0
            SS = V1*V1 + V2*V2
            IF (SS < 1.0) EXIT L10
        ENDDO L10
        TEMP3=-2.0*LOG(SS)/SS
        IF (TEMP3 > 0.0) THEN
            SS = SQRT(TEMP3)
        ELSE
            WRITE(*,150)
            150      FORMAT(' NEGATIVE ARGUMENT TO RANDOM DEVIATE')
            STOP
        ENDIF
        NORMAL_RANDOM_DEVIATE = V1*SS
        RETURN
    END ! REAL (4) FUNCTION NORMAL_RANDOM_DEVIATE
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (4) FUNCTION LOGNORMAL_RANDOM_DEVIATE1(NORMMEAN,NORMSD)
        IMPLICIT NONE
        !     ********************************************************************
        !       Generates a lognormal_random_deviately distributed random deviate
        !     ********************************************************************
        REAL (4), INTENT(IN) :: NORMSD,NORMMEAN
        REAL (4) ::  Y,NORMAL_RANDOM_DEVIATE
        !         write(*,200) normmean,normsd
        200   FORMAT(' NORMMEAN=',G12.5,' NORMSD=',G12.5)
        Y = NORMAL_RANDOM_DEVIATE()
        LOGNORMAL_RANDOM_DEVIATE1 = EXP(NORMMEAN+NORMSD*Y)
        RETURN
    END ! REAL (4) FUNCTION LOGNORMAL_RANDOM_DEVIATE1
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (4) FUNCTION LOGNORMAL_RANDOM_DEVIATE(STDEV)
        IMPLICIT NONE
        !     ********************************************************************
        !       Generates a lognormally distributed random deviate
        !     ********************************************************************
        REAL (4), INTENT(IN) :: STDEV
        REAL (4) ::  Y,NORMAL_RANDOM_DEVIATE
        Y = NORMAL_RANDOM_DEVIATE()
        LOGNORMAL_RANDOM_DEVIATE = EXP((1.0-0.5*STDEV)*Y*STDEV)
        RETURN
    END ! REAL (4) FUNCTION LOGNORMAL_RANDOM_DEVIATE
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (4) FUNCTION EXPONENTIAL_DISTRIBUTION()
        IMPLICIT NONE
        !     ********************************************************************
        !       Generates a exponentially distributed random deviate
        !     ********************************************************************c
        REAL (4) :: STDEV
        REAL (4) :: Y
        REAL (4) :: RRAND
        EXTERNAL RRAND
        !          write(*,300)
        300   FORMAT(' EXPONENTIAL_DISTRIBUTION')
        L100: DO
            Y=RRAND()
            IF (Y /= 0.0) EXIT L100
        ENDDO L100
        STDEV=Y
        EXPONENTIAL_DISTRIBUTION=-LOG(STDEV)
        RETURN
    END ! REAL (4) FUNCTION EXPONENTIAL_DISTRIBUTION
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL(4) FUNCTION RRAND()
    REAL :: RANDNUMBER
        CALL RANDOM_NUMBER(RANDNUMBER)
        RRAND=RANDNUMBER
    RETURN
    END FUNCTION RRAND
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE BOUNDARY_CONDITIONS()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J
        REAL (4) :: UPLIFT
          !!!!!!!write(outhist,4833)
4833    format('boundary called')      
        !     **********************************************************************
        !       This subroutine applies boundary conditions at the end of the
        !         iteration.  The default is just a lowering of the
        !         lower matrix boundary.  But it could also include variable surface
        !         and subsurface uplift or depression through direct manipulation of
        !         elevation(i,j) and of the direct access resistance file or more indirectly
        !         through manipulation of the erosion_depth_index(i,j) file through, as shown below, a
        !         new matrix deformation(i,j)
        !   MODIFIES:  ELEVATION, SEDIMENT_BASE, DEFORMATION, NEWBASE
        !     **********************************************************************
        IF (.NOT.IS_Y_PERIODIC) THEN
            IF (MODEL_OCEAN_LEVEL) THEN
                DO I=1,MX
                    IF (ELEVATION(I,MY) > OCEAN_ELEVATION) THEN
                        ELEVATION(I,MY)=OCEAN_ELEVATION
                        IF (ELEVATION(I,MY) < SEDIMENT_BASE(I,MY)) SEDIMENT_BASE(I,MY)=ELEVATION(I,MY)
                        IS_ROCK_SURFACE(I,MY)=.FALSE.
                    ENDIF
                ENDDO
            ELSE
                IF (HORIZONTAL_LOWER_BOUNDARY) THEN
                    NEWBASE =TIME_INCREMENT * BOUNDARY_LOWERING_RATE
                    DO  I = 1 ,MX
                        ELEVATION(I,MY) = ELEVATION(I,MY) - NEWBASE
                    ENDDO
                ELSE
                    DO  I = 1, MX
                        IF (ERODING_LOWER_BOUNDARY(I)) THEN
                            ELEVATION(I,MY) = ELEVATION(I,MY) +TIME_INCREMENT * CFNW(I,MY-1)
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
        ENDIF
		!   **************************************************************************
		!   Force outer perimeter to remain at their initial elevation if fixed flow boundaries
        !		are specified - New: 10/14/17
		!   *************************************************************************
		IF (DO_FLOW_BOUNDARIES) THEN
		    J=1
			DO I=1,MX
			    ELEVATION(I,J)=INITIAL_ELEVATION(I,J)
			ENDDO
			J=MY
			DO I=1,MX
			    ELEVATION(I,J) = INITIAL_ELEVATION(I,J)
			ENDDO
			I=1
			DO J=1,MY
			    ELEVATION(I,J) = INITIAL_ELEVATION(I,J)
			ENDDO
			I=MX
			DO J=1,MY
			    ELEVATION(I,J)= INITIAL_ELEVATION(I,J)
			ENDDO
		ENDIF
        !     **********************************************************************
        !      Here is an example of a spatially-varying uplift rate
        !          affecting both the surface elevation and the
        !         rock structural elevation
        !         to implement this I have:
        !           declared deformation as a matrix in erode.ins
        !           declared deformscale as a global variable
        !              ncluded it in the time-dependent
        !              list of simulation parameters included in the matrix simulation_parameters
        !              and read in from file inrates.dat
        !           Write out deformation as an output matrix and read it in again if
        !              the program is restarted (much as, say, sedbase is read in)
        !           Modified the subroutine finderodibility as shown below so as to
        !             have the entry into the 3-d resistance matrix affected by the
        !             cumulative deformation
        !     **********************************************************************
        IF (DO_ROCK_DEFORMATION) THEN
            L140:   DO  J=1,MY
                IF (J <= 40) THEN
                    UPLIFT=DEFORMSCALE *TIME_INCREMENT
                ELSE
                    !   Example of linearly varying deformation
                    IF (J < 50) THEN
                        UPLIFT=(5.0-J/10.0)*DEFORMSCALE *TIME_INCREMENT
                    ELSE
                        UPLIFT=0.0
                    ENDIF
                ENDIF
                L150:   DO  I=1,MX
                    !    Example of sinusoidal deformation
                    !                  uplift=dsin(aiparam*i+biparam)*dsin(ajparam*j+bjparam)
                    !         +          *time_increment * deformscale
                    143       FORMAT(2I5,3G12.5)
                    IF (USE_3D_SLOPE_RESISTANCE) THEN
                        DEFORMATION(I,J)=DEFORMATION(I,J)+UPLIFT
                    ENDIF
                    ELEVATION(I,J)=ELEVATION(I,J)+UPLIFT
                    IF (DO_SEDIMENT_TRANSPORT) SEDIMENT_BASE(I,J)=SEDIMENT_BASE(I,J)+UPLIFT
                ENDDO L150
            ENDDO L140
        ENDIF
        RETURN
    END !  SUBROUTINE BOUNDARY_CONDITIONS
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DETERMINE_ERODIBILITY()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL (4) :: READ_ERODIBILITY
        !     ***********************************************************************
        !      Determines the erodibilty as a function of time and space.  The
        !      weatherability and the bedrock fluvial erodibilty are assumed to be
        !      proportional.  Creep rate and regolith fluvial erodibility are
        !      assumed not to vary with erodibility.
        !  MODIFIES: EROSION_DEPTH_INDEX, RELATIVE_RESISTANCE
        !  CALLS: READ_ERODIBILITY
        !     ***********************************************************************
        INTEGER :: I,J,K
        L100: DO  J=1,MY
            L101: DO  I=1,MX
                IF (IS_ROCK_SURFACE(I,J)) THEN
                    IF (DO_ROCK_DEFORMATION) THEN
                        K=((EREFERENCE-ELEVATION(I,J))+DEFORMATION(I,J))/ &
                        VERTICAL_RESISTANCE_SCALING+1
                    ELSE
                        K=(EREFERENCE-ELEVATION(I,J))/VERTICAL_RESISTANCE_SCALING+1
                    ENDIF
                ELSE
                    IF (DO_ROCK_DEFORMATION) THEN
                        K=((EREFERENCE-ELEVATION(I,J))+DEFORMATION(I,J)&
                        - REGOLITH(I,J))/VERTICAL_RESISTANCE_SCALING+1
                    ELSE
                        K=(EREFERENCE-ELEVATION(I,J)-REGOLITH(I,J))/VERTICAL_RESISTANCE_SCALING+1
                    ENDIF
                ENDIF
                IF (K > MZ) K=MZ
                IF (K < 1) K=1
                IF (K /= EROSION_DEPTH_INDEX(I,J)) THEN
                    RELATIVE_RESISTANCE(I,J)=READ_ERODIBILITY(I,J,K)
                    EROSION_DEPTH_INDEX(I,J)=K
                ENDIF
            ENDDO L101
        ENDDO L100
        RETURN
    END ! SUBROUTINE DETERMINE_ERODIBILITY
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (4) FUNCTION READ_ERODIBILITY(I,J,K)
        !     ***********************************************************************
        !      Reads a value from the erodibility "cube" input file of random deviates
        !      and scales it to a lognormal distribution with mean unity and standard
        !      deviation RESISTANCE_VARIABILITY
        !     ***********************************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER I,J,K,IREC
        REAL (4) :: RRRR
        IREC=J+MY*(I-1)+MX*MY*(K-1)
        READ(INRESIST,REC=IREC) RRRR
        IF (SCALE_3D_ROCK_ERODIBILITY) THEN
            READ_ERODIBILITY=RANDMULT*EXP(RRRR*SIGMANORM-SIGMASQ)
        ELSE
            READ_ERODIBILITY=RRRR
        ENDIF
        RETURN
    END ! REAL (4) FUNCTION READ_ERODIBILITY
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DETERMINE_EROSION_RATE()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        !     **********************************************************************
        !      Eetermines the time-varying simulation parameters during each iteration
        !         first checks to see if we have exceeded the maximum time for the
        !         current erosion rate category.  if so, then set the new simulation
        !         parameters from the matrices erosion_rate_values and simulation_parameters
        !       %%%%%% needs to be updated with newer simumaltion parameters that
        !              might be temporally varying %%%%%%%%%
        !  MODIFIES: PARAMETER_CHANGE_INDEX
        !  CALLS: PRINT_RATE_STATISTICS
        !     **********************************************************************
        IF (PRESENT_TIME > TIMES_FOR_PARAMETER_CHANGES(PARAMETER_CHANGE_INDEX)) THEN
            IF (PARAMETER_CHANGE_INDEX > NUMBER_OF_PARAMETER_CHANGES) THEN
                 PARAMETER_CHANGE_INDEX=NUMBER_OF_PARAMETER_CHANGES
            ELSE
                WRITE(OUTHIST,100) PARAMETER_CHANGE_INDEX
                100  FORMAT(' NOW USING EROSION PARAMETERS FOR PARAMETER_CHANGE_INDEX=',I5)
                CALL PRINT_RATE_STATISTICS()
                BOUNDARY_LOWERING_RATE=EROSION_RATE_VALUES(PARAMETER_CHANGE_INDEX)
                SLOPE_DIFFUSIVITY=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,1)
                BEDROCK_ERODIBILITY=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,2)
                DISCHARGE_COEFFICIENT=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,3)
                DETACHMENT_CRITICAL_SHEAR=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,4)
                TRANSPORTFACTOR=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,5)
                TRANSPORT_CRITICAL_DIM_SHEAR=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,6)
                BISTABLE_CRITICAL_SHEAR=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,7)
                LOW_EROSION_THRESHOLD=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,8)
                HIGH_EROSION_THRESHOLD=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,9)
                BISTABLE_RUNOFF_FACTOR=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,10)
                DEFORMSCALE=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,11)
                BISTABLE_BEDROCK_ERODIBILITY=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,12)
                WRITE(OUTHIST,200) BOUNDARY_LOWERING_RATE,SLOPE_DIFFUSIVITY,BEDROCK_ERODIBILITY,DISCHARGE_COEFFICIENT &
                ,DETACHMENT_CRITICAL_SHEAR,TRANSPORTFACTOR,TRANSPORT_CRITICAL_DIM_SHEAR,BISTABLE_CRITICAL_SHEAR &
                ,LOW_EROSION_THRESHOLD,HIGH_EROSION_THRESHOLD,BISTABLE_RUNOFF_FACTOR,DEFORMSCALE,BISTABLE_BEDROCK_ERODIBILITY
                200 FORMAT(' BOUNDARY LOWERING RATE=',G12.5,/, &
                ' SLOPE DIFFUSIVITY=',G12.5,/, &
                ' BEDROCK ERODIBILITY=',G12.5,/, &
                ' DISCHARGE COEFFICIENT=',G12.5,/, &
                ' DETACHMENT CRITICAL SHEAR=',G12.5,/, &
                ' TRANSPORT FACTOR=',G12.5,/, &
                ' TRANSPORT CRITICAL DIMENSIONLESS SHEAR=',G12.5,/, &
                ' BISTABLE CRITICAL SHEAR=',G12.5,/, &
                ' LOW EROSION THRESHOLD=',G12.5,/, &
                ' HIGH EROSION THRESHOLD=',G12.5,/, &
                ' BISTABLE RUNOFF FACTOR=',G12.5,/, &
                ' DEFORMSCALE=',G12.5,/, &
                ' BISTABLE BEDROCK ERODIBILITY=',G12.5)
                DISCHARGE_COEFFICIENT=WATER_DENSITY*(MANNING/CHANNEL_WIDTH_CONSTANT)**0.6*GRAVITY * &
               (9.8/GRAVITY)**0.3
                PREVIOUS_DISCHARGE_COEFFICIENT=DISCHARGE_COEFFICIENT
                PREVIOUS_CRITICAL_SHEAR=DETACHMENT_CRITICAL_SHEAR
            ENDIF
            PARAMETER_CHANGE_INDEX=PARAMETER_CHANGE_INDEX+1
        ENDIF
        RETURN
    END ! SUBROUTINE DETERMINE_EROSION_RATE

    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE MAKE_EVENT(EVENT_NUMBER)
        USE ERODE_GLOBALS
        USE EVENT_GLOBALS
        USE CRATER_GLOBALS
        USE MASS_FLOW_GLOBALS
        IMPLICIT NONE
        REAL (4) :: ZPLANE,X1,X2,Y1,Y2
        REAL (4) :: AA,BB,DD,SS,SLINE,INTLINE,ZCOMP,RRAND
        INTEGER :: I,J,JCOMP
        EXTERNAL RRAND
        INTEGER, INTENT(IN) :: EVENT_NUMBER
        LOGICAL :: IS_RANDOM
        !     **********************************************************************
        !      Does user-defined modifications of the simulation at specified times
        !        during the simulations.  those times are read in from the 'events.prm'
        !        file.  this routine is called only once for each event
        !        the times can be accessed in the matrix event_times(event_number)
        !  MODIFIES: Whatever appropriate simulation variables
        !     **********************************************************************
        WRITE(*,100) EVENT_NUMBER
        WRITE(OUTHIST,100) EVENT_NUMBER
        100   FORMAT(' EVENT NUMBER ',I5,' HAS OCCURRED')
        !     An example: changing to universal runoff
    !    COMPLETE_RUNOFF=.TRUE.
    !    DO_ALLUVIAL_REROUTING=.FALSE.
    !    RETURN
    !END
    !     **********************************************************************
    !      some other possible event assigments:
    !           change erosion threshold:
    !               detachment_critical_shear=200.0
    !           add a resistant layer:
    !               resistant_surface_layer=.true.
    !               surface_layer_thickness=10.0
    !               surface_layer_resistance=10.0
    !               do j=1,my
    !               do i=1,mx
    !                   previous_elevation(i,j)=elevation(i,j)-surface_layer_thickness
    !               enddo
    !               enddo
    !           change the discharge to a lower, but more steady value
    !               discharge_constant=1.0e-05
    !               use_random_discharge=.false.
    !               flow_fraction=1.0
    !           increase downstream discharges relative to upstream values
    !               discharge_exponent=1.0
    !           reduce creep and rate of weathering
    !               slope_diffusivity=0.002
    !               rock_weathering_rate=0.0001
    !     **********************************************************************
    !        the following commented-out section has been used to model wave-cut erosion in a simulation
    !           of erosion of a coastal-plain landscape
    !     **********************************************************************

    !           ss=0.2
    !           if (event_number == 1) then
    !              zplane=23.0
    !              x1=cell_size*(mx/10)
    !              x2=cell_size*mx
    !              y1=cell_size*my
     !             y2=cell_size*(my/20)
     !           else
    !              zplane=8.0
    !              x1=cell_size*(mx/2)
    !              x2=cell_size*mx
    !              y1=cell_size*my
    !              y2=cell_size*(my/2)
    !            endif
    !            aa=sqrt((ss**2*(y1-y2)**2)/((y1-y2)**2+(x1-x2)**2))
    !            bb=sqrt(ss**2-aa**2)
    !            dd=-aa*x1-bb*y1-zplane
    !            sline=(y1-y2)/(x1-x2)
    !            intline=y1-sline*x1
    !            !write (*,200) aa,bb,dd,sline,intline
    !            200 format('aa=',g12.5,' bb=',g12.5,' dd=',g12.5,' sline=',g12.5, &
    !            ' inter=',g12.5)
    !            do j=1,my
    !               do i=1,mx
    !                  jcomp=(sline*i*cell_size+intline)/cell_size
    !                  if (jcomp <= j) then
    !                     if (elevation(i,j) > zplane) then
    !                        elevation(i,j)=zplane+(rrand()-0.5)*0.1
    !                        is_rock_surface(i,j)=.false.
     !                       regolith(i,j)=initial_regolith_thickness
    !                        if (elevation(i,j) <= sediment_base(i,j)) then
    !                           sediment_base(i,j)=elevation(i,j)
    !                           is_sediment_covered(i,j)=.false.
    !                        endif
    !                     endif
    !                  else
    !                     zcomp=-aa*i*cell_size-bb*j*cell_size-dd
    !                     if (elevation(i,j) > zcomp) then
    !                        elevation(i,j)=zcomp+(rrand()-0.5)*0.1
    !                        is_rock_surface(i,j)=.false.
    !                        regolith(i,j)=initial_regolith_thickness
    !                        if (elevation(i,j) <= sediment_base(i,j)) then
    !                           sediment_base(i,j)=elevation(i,j)
    !                           is_sediment_covered(i,j)=.false.
    !                        endif
    !                     endif
    !                  endif
    !               enddo
    !            enddo
! This code is to make new impact craters at specific times, diameters, and x-y locations
			IF (USE_FLOW_VOLUME) THEN
			        GLEN_LAW_ALPHA=ALPHAS(EVENT_NUMBER)
                    FLOW_DIFFUSIVITY=VISCOUS(EVENT_NUMBER)
                    DELTA_HB=CRITDEPTHS(EVENT_NUMBER)
                    MASS_FLOW_EROSION_RATE=SCOURS(EVENT_NUMBER)
                    MAXIMUM_FLOW_DEPTH_EROSION=MAXDEPTHS(EVENT_NUMBER)
                    ROCK_WEATHERING_RATE=WEATHERS(EVENT_NUMBER)
			ELSE
                 IS_RANDOM=.FALSE.
                 IF (EVENT_NUMBER>NUMBER_OF_EVENTS) RETURN
                 IF (USE_REAL_CRATERS) THEN
                 WRITE(*,717) DIAMETERS(EVENT_NUMBER),MINIMUM_REAL_CRATER_DIAMETER,MAXIMUM_REAL_CRATER_DIAMETER
717              FORMAT('CRATER EVENT, DIAMETER=',G13.6,' MINREALSIZE=',G13.6,' MAXREALSIZE=',G13.6)
                 ENDIF
                 IF ((USE_REAL_CRATERS).AND.(DIAMETERS(EVENT_NUMBER)>=MINIMUM_REAL_CRATER_DIAMETER).AND. &
                     (DIAMETERS(EVENT_NUMBER)<=MAXIMUM_REAL_CRATER_DIAMETER)) THEN
                     CALL CREATE_REAL_CRATER(IS_RANDOM,DIAMETERS(EVENT_NUMBER),CRATER_X_LOCATIONS(EVENT_NUMBER) &
                         ,CRATER_Y_LOCATIONS(EVENT_NUMBER))
                 ELSE
                     CALL DO_IMPACT_CRATERING(IS_RANDOM,DIAMETERS(EVENT_NUMBER),CRATER_X_LOCATIONS(EVENT_NUMBER) &
                     ,CRATER_Y_LOCATIONS(EVENT_NUMBER))
                 ENDIF
                 CRATERING_CALLED=.TRUE.
			ENDIF
               return
               end ! SUBROUTINE MAKE_EVENT
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     SUBROUTINE FINDOCEAN_ELEVATION(TDUMMY)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL (4), INTENT(INOUT) :: TDUMMY
        REAL (4) :: EDUMMY,LASTOCEANLEVEL
        INTEGER :: I,J
        LOGICAL :: HASFOUND
    !  MODIFIES: TDUMMY, OCEAN_ELEVATION, PRESENT_TIME,MAXIMUMELEVATION, OCEANLSOPE, OCEANINT
    !            OCEANNEXTTIME
        IF (.NOT.VARIABLE_OCEAN_ELEVATION) RETURN
        LASTOCEANLEVEL=OCEAN_ELEVATION
        EDUMMY=-1.0E+25
        DO J=1,MY
            DO I=1,MX
                IF (ELEVATION(I,J) > EDUMMY) EDUMMY=ELEVATION(I,J)
            ENDDO
        ENDDO
        MAXIMUMELEVATION=EDUMMY
        L200: DO
            IF (TDUMMY >= OCEANNEXTTIME) THEN
                HASFOUND=.FALSE.
                DO I=1,INUMLEVELS-1
                    IF (TDUMMY <= OCEAN_RECALCULATION_TIMES(I+1)) THEN
                        HASFOUND=.TRUE.
                        EXIT
                    ENDIF
                ENDDO
                IF (HASFOUND) THEN
                    OCEANNEXTTIME=OCEAN_RECALCULATION_TIMES(I+1)
                    OCEANTSLOPE=(OCEAN_LEVELS(I+1)-OCEAN_LEVELS(I))/  &
                    (OCEAN_RECALCULATION_TIMES(I+1)-OCEAN_RECALCULATION_TIMES(I))
                    OCEANTINT=OCEAN_LEVELS(I)-OCEANTSLOPE*OCEAN_RECALCULATION_TIMES(I)
                ELSE
                    RETURN
                ENDIF

            ENDIF
            OCEAN_ELEVATION=OCEANTSLOPE*PRESENT_TIME+OCEANTINT
            IF (OCEAN_ELEVATION >= MAXIMUMELEVATION) THEN
                PRESENT_TIME=PRESENT_TIME+100.0
                TDUMMY=TDUMMY+100.0
                IF (PRESENT_TIME < MAXIMUM_SIMULATION_TIME) CYCLE L200
            ENDIF
            EXIT L200
        ENDDO L200
        RETURN
    END ! SUBROUTINE FINDOCEAN_ELEVATION
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LOGICAL (4) FUNCTION CHANGE_FLOW_DIRECTION(GRADIENT_RATIO)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL (4), INTENT(IN) :: GRADIENT_RATIO
        REAL (4) :: THEPROB, RRAND
        EXTERNAL RRAND
        INTEGER ICAT
        !     **********************************************************************
        !      This subroutine determines if the direction of flow on an alluvial
        !         surface remains the same or changes direction to a path with steeper
        !         gradient.  The decision is random with probability being an increasing
        !         function of the ratio of the new to the existing gradient and
        !         of the magnitude of the parameter sticky_routing_critical_value.  Uses a lookup
        !         function
        !     **********************************************************************
        IF (GRADIENT_RATIO <= 1.0) THEN
            CHANGE_FLOW_DIRECTION=.FALSE.
        ELSEIF (GRADIENT_RATIO > STICKYCAT(11)) THEN
            CHANGE_FLOW_DIRECTION=.TRUE.
        ELSE
            ICAT=INT((STICKYFACTOR*(GRADIENT_RATIO-1.0))+1.0)
            THEPROB=STICKYPROB(ICAT)+(GRADIENT_RATIO-STICKYCAT(ICAT))* &
            STICKYDEL(ICAT)
            IF (RRAND() > THEPROB) THEN
                CHANGE_FLOW_DIRECTION=.FALSE.
            ELSE
                CHANGE_FLOW_DIRECTION=.TRUE.
            ENDIF
        ENDIF
        RETURN
    END ! LOGICAL (4) FUNCTION CHANGE_FLOW_DIRECTION

 !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE SETUP_EVENTS()
       USE ERODE_GLOBALS
       USE EVENT_GLOBALS
       IMPLICIT NONE
       INTEGER :: IERROR, I
      !     **********************************************************************
       !  This routine initializes any information required in order to simulate
       !  the occurrence of events in the program, such as, for example, reading
       !  information from a file.
       !  Comment out any code except 'RETURN' if SETUP_EVENTS() is to do nothing
      !     **********************************************************************
       ! This example code reads in the diameters, x and y locations (in the simulation
       ! reference frame) of individual cratering events
	IF (USE_FLOW_VOLUME) THEN
	    ALLOCATE(VISCOUS(NUMBER_OF_EVENTS),CRITDEPTHS(NUMBER_OF_EVENTS),WEATHERS(NUMBER_OF_EVENTS) &
       ,MAXDEPTHS(NUMBER_OF_EVENTS),ALPHAS(NUMBER_OF_EVENTS),SCOURS(NUMBER_OF_EVENTS),STAT=IERROR)
       OPEN(47,FILE='BINGHAM_EVENTS.PRM',ACTION='READ')
 !      READ(47,*) EVENT_TIME_SCALE
       DO I=1,NUMBER_OF_EVENTS
           READ(47,*) ALPHAS(I),VISCOUS(I), CRITDEPTHS(I),SCOURS(I),MAXDEPTHS(I),WEATHERS(I)
     !      EVENT_TIMES(I)=EVENT_TIMES(I)*EVENT_TIME_SCALE
       ENDDO
	ELSE
       ALLOCATE(DIAMETERS(NUMBER_OF_EVENTS),CRATER_X_LOCATIONS(NUMBER_OF_EVENTS) &
       ,CRATER_Y_LOCATIONS(NUMBER_OF_EVENTS),STAT=IERROR)
       IF (IERROR > 0) CALL ALLOCERROR()
       OPEN(47,FILE='CRATER_EVENTS.PRM',ACTION='READ')
       READ(47,*) EVENT_TIME_SCALE
       DO I=1,NUMBER_OF_EVENTS
           READ(47,*) DIAMETERS(I),CRATER_X_LOCATIONS(I),CRATER_Y_LOCATIONS(I)
           EVENT_TIMES(I)=EVENT_TIMES(I)*EVENT_TIME_SCALE
       ENDDO
	ENDIF
       CLOSE(47)
       RETURN
    END ! SUBROUTINE SETUP_EVENTS
 !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DETERMINE_CRATERING_RATE(NEW_CRATER_RATE)
    USE ERODE_GLOBALS
    USE CRATER_GLOBALS
    IMPLICIT NONE
    REAL (4) :: END_RATE,END_TIME,ALPHA_TERM,NEW_CRATER_RATE
    LOGICAL SKIP
     !     **********************************************************************
       !  This routine changes the relative cratering rate from the value read
       !  in from MARSSIM.PRM to the relative value END_RATE following a negative exponential
       !  function, reaching END_RATE*IMPACT_PROBABILITY at END_TIME.
       !  this code is skipped if SKIP is set to .TRUE.
       !  The parameters for this routine may in the future be programmed to be read in
       !  from MARSSIM.PRM
      !     **********************************************************************
       SKIP=.TRUE.
       IF (SKIP) THEN
           NEW_CRATER_RATE=IMPACT_PROBABILITY
           RETURN
       ELSE
           IF (PRESENT_TIME<(0.25*MAXIMUM_SIMULATION_TIME)) THEN
               NEW_CRATER_RATE=2.67*IMPACT_PROBABILITY
           ELSE
               IF (PRESENT_TIME<(0.6*MAXIMUM_SIMULATION_TIME)) THEN
                   NEW_CRATER_RATE=0.67*IMPACT_PROBABILITY
               ELSE
                   NEW_CRATER_RATE=0.17*IMPACT_PROBABILITY
               ENDIF
           ENDIF
        
            !END_RATE=0.01
            !END_TIME=MAXIMUM_ITERATION*MINUMUM_TIME_INCREMENT
           ! ALPHA_TERM=-LOG(END_RATE)/END_TIME
            !NEW_CRATER_RATE=IMPACT_PROBABILITY*EXP(-ALPHA_TERM*PRESENT_TIME)
            RETURN
       ENDIF
       END ! SUBROUTINE DETERMINE_CRATERING_RATE
 !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE ALLOCERROR()
    USE ERODE_GLOBALS
    IMPLICIT NONE
    WRITE(OUTHIST,100)
100 FORMAT('ALLOCATION ERROR, STOPPING PROGRAM')
    WRITE(*,100)
    STOP
    END !SUBROUTINE ALLOCERROR



