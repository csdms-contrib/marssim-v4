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
    PROGRAM MARSSIM
        !USE IFPORT
        USE ERODE_GLOBALS
        USE LAKE_GLOBALS
        USE CRATER_GLOBALS
        USE EOLIAN_GLOBALS
        USE LAVA_GLOBALS
        USE ACCRETION_GLOBALS
        USE EVENT_GLOBALS
		USE GRAVEL_MIXTURE_GLOBALS
        USE GROUNDWATER_VARIABLES
        USE MASS_FLOW_GLOBALS
        IMPLICIT NONE
        INTEGER I,J,K,NCALL,SOURCEUSE,KLOW,IDOSTOP,INSTOP
        REAL (4) :: PROBVALUE
        REAL(4) :: EXPONENTIAL_DISTRIBUTION,EEEMIN,NETDIFF
        REAL(4) :: LAVAPROBTOUSE,EOPROBTOUSE,CRATERPROBTOUSE
        REAL(4) :: TOTQ,BOUNDQ, RRAND,LOCAL_PROBABILITY,ice_base_elevation,max_FLOW_VOLUME_increment
        LOGICAL IS_RANDOM
        REAL(4) :: DUMMY_D,DUMMY_X,DUMMY_Y
        LOGICAL DODETAILS
        EXTERNAL RRAND
        !       *********************************************************************
        !        This is the master program for the marssim landform evolution model.
        !        The basic algorithms, assumptions and structure of the portions of the
        !        model dealing with weathering, slope processes, fluvial erosion, fluvial transport
        !        and sedimentation, cratering, sublimation, eolian processes, and groundwater flow/sapping
        !        have been sumamrized in several publications:
        !          howard, a. d., 1994, water res. research, 30(7), 2261-85.
        !          howard, a. d., 1997, earth surf. proc. landforms, 22, 211-227.
        !          howard, a. d., 1999, in incised river channels, s. darby and a. simon, eds.,
        !                   wiley, 277-300.
        !          fagherazzi, s., howard, a. d., and wiberg, p. l.,2004, j. geophys res. 109,
        !              doi:10.1029/2003jf000091.
        !          howard, a. d., 2007, geomorphology, 91, 332-363.
        !          howard, a.d. and moore, j.m., 2008, geophys. res. lett.,35, doi:10.1028/2007gl032618
        !          barnhart, c. j., howard, a. d., and moore, j. m., 2008, j. geophys. res., 114,
        !            e01003, doi:10.1029/2008je003122
        !          howard, a.d. and tierney, h.e., 2011, geomorphology, 137, 27-40
        !          matsubara, Y. et al., 2011, j. geophys. res., 116, e04001, doi:10.1029/2010JE003739
        !          luo, w., howard, a.d., 2008, j. geophys. res., 113,E05002 doi:10.1029/2007je002981
        !          howard, a.d. et al., 2012, icarus, 220, 268-276
        !          matsubara, y. et al., 2018, j. geophys. res., 123, joi:10.1029/2018JE005572
        !          howard, a.d., et al., 2016, icarus, 270, 100-113
        !          moore, j.m. et al., 2017, icarus, 287, 320-333
        !
        !        References to publications related to the more specialized planetary processes are
        !        given in individual subroutines.
        !
        !        Portability issues:  The program is almost entirely written in standard fortran 90.
        !          Most data file output and input is from ascii files and should be identical on all architectures.
        !          However, the program
        !            also writes out binary image files (Photoshop 'raw' files), mostly shaded relief images.
        !          These routines, "imagewrite", "shadewrite" , "write_lava+ages",
        !            and "colorwrite" might have to be modified for other compilers or architectures.
        !          These use the non-standard
        !            output control character '$', which does not terminate a record.
        !          Why is FORTRAN so backward in terms of output format control???
        !          For standard compilers such as gfortran the next best solution is to change the format statement from:
        !                 130 format(a1,$)
        !          to:
        !                 130 format(a1)
        !          and the corresponding write statements from:
        !                 write(5,130) achar
        !          to:
        !                 write(5,130,advance='no')
        !          The disadvantage of this is that a record terminator (CR & LF) or (LF) is appended to the end
        !            of the file and you need to deal with this when importing into, say,
        !            Photoshop by telling it to ignore the last 1 or 2 bytes.
        !          All image-writing routines that use the $ formatting. are in "read_and_write_data_files.f90".
        !          A replacement routine using the above approach with advance='no' is provided as "alternate_read_and_write.f90".
        !          Probably the best solution is to link in a standard image library, such as TIFF and to rewrite
        !            the routines to directly write formatted image files.
        !          But I'm happy with my PGI and INTEL compilers and raw images.
        !        This version of 04/2020 has several changes from the previous version, including
        !           * The cratering module can use a database of martian fresh crater instead of a parameterized
        !               shape file. This makes craters more visually appealing, but also indroduces important
        !               aspects of real craters, including rim and ejecta variability, including interior rim slumps,
        !               central peaks and pits. The morphology is based upon the MOLA shot data for the individual
        !               crater.  A disadvantage is that the DEM for the craters must be created for the simulation
        !               resolution.  The MOLA database and the DEM generation routine is included with the program
        !               distribution.  Due to limitations of the MOLA data for small craters and the rarity of large
        !               fresh craters, the range of diameters for which the database is useful is about 15-100 km.
        !               More documentation is found in the cratering module
        !           * A new mass flow module is included.  It models deeper mass flows using either Bingham or Glen's Law
        !               rheology.  The flow occurs within the surficial regolith layer.  Regolith can be generated by
        !               accretion or by rock weathering, or as a uniform initial layer.  The model assumes the shallow
        !               flow approximation and continuous deformation. It thus does not handle flows with appreciable
        !               momentum, such as avalanches, or spatially varying rheological or density variations.  It can,
        !               account for a specified rigid surface layer or a critical yeild stress in Bingham flow.  For
        !               Bingham flow either a high yeild stress, shallow regolith, or low surface gradients can cause
        !               flow to cease.  The mass flow can optionally erode the underlying bedrock, converting it into
        !               regolith. The "regolith" might be conceptualized as glacial ice or deformable particulates as in
        !               earth flows.
        !           * Fluvial bedrock erosion can utilize the Sklar & Dietrich bed abrasion model instead of stream
        !               power.
        !           * In the previous version the governing parameters were in a single parameter file. In the present
        !               version the are several parameter files governing specific process suites. Each process module
        !               now reads the parameter file governing the processes, so that the variables and their use
        !               within the program are more directly compared.
        !           * There is the option to output pereodic binary output of state variable matrices in addition to
        !               ascii files
        !           ^ The program calculates and outputs extimates of rhe relative "work" done by various process,
        !               both in terms of lateral and vertical mass transfer
        !           * At the close of the simulation, summary matrices are printed out showing the net changes produced
        !               by the simulation, including overall net elevation change, elevation changes produced by
        !               particular processes, net weathering, etc.
        !           * In addition to periodic output of state variable matrices, matrices of the final state variable
        !               matrices are also produced, which are useful for restarting a continuing simulation or
        !               summarizing the final evolution state.
        !           * The shaded relief output files have been changed from a RAW file of bytes to a PBM file, which
        !               can be opened directly without needing to specify image dimensions in Photoshop or Gimp
        !        ********************************************************************
        !        define input - output file unit numbers
        !   MODIFIES:  (MOST INPUT-OUTPUT FILE NUMBERS), SUBMERGED, ITERATION, PREVIOUS_ELEVATION, ELAPSEDTIME,
        !              WRITE_OUT_CHANNELS, PRESENT_TIME, EVENT_DONE, EVENT_INDEX, DO_EVENTS
        !              SUBGRADCHANGE, ABSGRADCHANGE, NUMIDCHANGE
        !   CALLS:  READ_INPUT_PARAMETERS, SETUP_FLUVIAL_SLOPE_EROSION, WRITE_SHADED_RELIEF_IMAGE,
        !           WRITE_COLOR_SHADED_RELIEF_IMAGE, WRITE_IMAGE, DO_FLUVIAL_AND_SLOPE, DEPOSIT_ICE,
        !           REPORT_MAX, FINDOCEAN_ELEVATION, WRITE_DEBUG_INFO, CHANNEL_PROPERTIES,WRITE_FIRST_DATA_MATRIX
        !           PRINT_MORPHOMETRY, CALCULATE_TOPO_DIVERGENCE, SUMMARIZE_CHANNELS, MAKE_EVENT
        !           WRITE_SECOND_DATA_MATRIX, WRITE_SOURDED_DISCHARGE, WRITE_SUBMERGED_LOCATIONS, WRITE_REPORT
        !           PRINT_SIMULATION_INFORMATION, DO_LAVA_FLOWS, WRITE_LAVA_INFO, WRITE_LAVA_AGES, DO_EOLOIAN_CHANGE
        !           DO_IMPACT_CRATERING, DO_EXPOSURE_DEPENDENT_CREEP, DO_ACCRETION_ABLATION, FINALIZE_FLUVIAL_SLOPE_EROSION
        !           WRITE_ALLUVIAL_LOCATIONS, WRITE_BEDROCK_LOCATIONS, WRITE_EROSION_DEPTH_INDEX, WRITE_DEFORMATION
        !           WRITE_SEDIMENT_BASE, WRITE_LAKE_INFO, WRITE_REGOLITH_THICKNESS, WRITE_ROCK_RESISTANCE
        !           PRINT_SIMULATION_INFORMATION, WRITE_ACCELERATED_EROSION_STATE, WRITE_GROUNDWATER_ELEVATION
        !           WRITE_GROUNDWATER_FLOW
        !
        !        ********************************************************************
        PARAMETER_CHANGE_INDEX=1
        INDATA=15
        INRESIST=13
        OUTDATA=19
        OUTHIST=17
        OUTCHAN=18
        OUTRECORD=11
        OUTSUMMARY=12
        OUTSOURCE=14
        OUTREPORT=30
        OUTCRATER=73
        OUTSEDDEBUG=74   
        OUTLAKE=75
        INRATES=32
        OUTIMGDAT=34
        OUTROCK=35
        OUTRESIST=36
        OUTGRAD=26
        OUTALLUV=25
        OUTEROSION_DEPTH_INDEX=45
        OUTELEV=47
        IOTEMP1=21
        IOTEMP2=22
        OUTBASE=23
		BINGFLUX=276
        OUTWORK=277
        OUTDEFORM=48
        OUTHIGH=49
        OUTREGOLITH=51
        OUTDISCHARGES=477
        OUT_BINARY_DISCHARGES=476
        OUTAVALANCHE=87
        OUT_IMAGE_FILE=88
        EEEMIN=0.0
        DODETAILS=.FALSE.
        IS_RANDOM=.TRUE.
        DUMMY_D=0.0
        DUMMY_X=0.0
        DUMMY_Y=0.0
        TOTAL_RANDOM_CRATERS=0
        INSTOP=31
        ITERATION=0
        FRAC2312 = 23.0E0/12.0E0
        FRAC43 = 4.0E0/3.0E0
        FRAC512 = 5.0E0/12.0E0
        FRAC23 = 2.0E0/3.0E0
        FRAC112 = 1.0E0/12.0E0
        !        ********************************************************************
        !        Simulation parameters are read from several files - see the individual
        !            source files.
        !        ********************************************************************
        !      ```````````
        OPEN (82,FILE='DEBUG.PRN',ACTION='WRITE')
        !      ``````````

        !        ********************************************************************
        !        basin.lst is descriptive file of simulation parameters and
        !        basin statistics during simulation
        !        ********************************************************************
        OPEN(OUTHIST,FILE='BASIN.LST',ACTION='WRITE')
        !        ********************************************************************
        !        channel.dat is list of drainage network loations and flow directions
        !        ********************************************************************
        OPEN(OUTCHAN,FILE='CHANNEL.DAT',ACTION='WRITE')
        !        ********************************************************************
        !        record.dat is historical record of progression of simulation
        !        towards steady state
        !        ********************************************************************
        OPEN(OUTRECORD,FILE='RECORD.DAT',ACTION='WRITE')
        !        ********************************************************************
        !        report.dat is historical record of relief & erosion rate
        !        ********************************************************************
        OPEN(OUTREPORT,FILE='REPORT.PRN',ACTION='WRITE')
        !        ********************************************************************
        !        summary.dat is like basin.lst without text
        !        ********************************************************************
        OPEN(OUTSUMMARY,FILE='SUMMARY.DAT',ACTION='WRITE')
        !        ********************************************************************
        !        source.dat is output file of channel source locations
        !        ********************************************************************
        OPEN(OUTSOURCE,FILE='SOURCE.DAT',ACTION='WRITE')
        !       **********************************************************************
        !        crater.dat is information about modeled crater events
        !        statistics.prm is information about the progress of the simulation
        !       **********************************************************************
        OPEN(OUTCRATER,FILE='CRATER.DAT',ACTION='WRITE')
        OPEN(78,FILE='STATISTICS.PRN',ACTION='WRITE')
        OPEN(OUTWORK,FILE='WORK_STATISTICS.DAT',ACTION='WRITE')
        !        ********************************************************************
        !        Read in initial parameters & elevations and summarize
        !        ********************************************************************
        CALL READ_INPUT_PARAMETERS()
        !        ********************************************************************
        !        Allocate and initialize simulation data
        !        ********************************************************************
        CALL SETUP_FLUVIAL_SLOPE_EROSION()                                                                      
        IF (FLUVIAL_AND_SLOPE_MODELING) SUBMERGED=.FALSE.
        SUBMERGED=.FALSE.
        !       **********************************************************************
        !        Write the raw image file for initial conditions
        !       **********************************************************************
            CALL WRITE_IMAGE()
            IF (USE_SOLAR_EROSION) THEN
                CALL WRITE_SHADED_RELIEF_IMAGE()
            ELSE
                CALL WRITE_SHADED_RELIEF_IMAGE()
            ENDIF
        ELOWEST=-1.0E+25
            KLOW=1
        ELAPSED_TIME=0.0 
        if (is_gravel_mixture) CALL write_output_gravel
        IF ((WRITE_INITIAL_EXPOSURE>0).AND.MODEL_EOLIAN_CHANGES) CALL WRITE_EXPOSURE
        !       **********************************************************************
        !        THIS IS THE MASTER ITERATION CYCLE
        !       **********************************************************************
        L200: DO  K=KLOW,MAXIMUM_ITERATION
            ITERATION=ITERATION+1
            TOTAL_ITERATIONS=TOTAL_ITERATIONS+1
            CALL GRADIENT_AND_FLOW_DIRECTION()
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                IF (BISTABLE_FLUVIAL_EROSION) THEN
                ELSE
                    IF(RESISTANT_SURFACE_LAYER) THEN
                    ELSE
                        PREVIOUS_ELEVATION=ELEVATION
                    ENDIF
                ENDIF
                IF ((ITERATION==1).AND.USE_SLOPE_DEPENDENT_RUNOFF) THEN
                    CALL SLOPE_RUNOFF_CHARACTERISTICS()
                ENDIF
                !       **********************************************************************
                !       Call Fluvial and slope erosion
                !       **********************************************************************
                CALL DO_FLUVIAL_AND_SLOPE()
                !       **********************************************************************
                !       Call sublimation and condensation routine
                !       **********************************************************************

                IF (DODETAILS) THEN
                    CALL REPORT_MAX()
                ENDIF
                CALL FINDOCEAN_ELEVATION(PRESENT_TIME)
                !       **********************************************************************
                !       Stop the simulation if we have exceeded the maximum simulation time
                !       **********************************************************************
                IF (PRESENT_TIME >= MAXIMUM_SIMULATION_TIME) EXIT L200
                !       **********************************************************************
                !       Write various reports at specified times
                !       **********************************************************************
                IF (DO_DEBUGGING) THEN
                    IF ((MAXIMUM_ITERATION-ITERATION) < 4) CALL WRITE_DEBUG_INFO()
                ENDIF
                IF (MODEL_ACCRETION_AND_ABLATION.AND.USE_SOLAR_EROSION) CALL DEPOSIT_ICE() 
            ENDIF
         
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                ELAPSEDTIME=ELAPSEDTIME+TIME_INCREMENT
            ELSE
                ELAPSEDTIME=ELAPSEDTIME+1.0
            ENDIF
            IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) == 0) THEN
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                    CALL CHANNEL_PROPERTIES()
                    IF ((.NOT.WRITE_BINARY_IMAGE_FILE).AND.USE_DISCHARGE) CALL WRITE_DISCHARGES()
                ENDIF
                CALL WRITE_ELEVATION_MATRIX()
            ENDIF
            IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) == 0) THEN
                IF (DO_MORPHOMETRY) THEN
                    write (*, 7012)
7012 format('print_morphometry')
                    CALL PRINT_MORPHOMETRY()
                    write(*,7013)
7013 format('topo_divergence')
                    CALL CALCULATE_TOPO_DIVERGENCE()
                    WRITE_OUT_CHANNELS=.TRUE.
                    write(*,7014)
7014 format('summarize_channels')
                    CALL SUMMARIZE_CHANNELS()
                ENDIF
            ENDIF
            !       ********************************************************************
            !        If it is time to write out a raw image of elevations, do so
            !       ***********************************************************************
            IF (MOD(ITERATION,IMAGE_OUTPUT_INTERVAL) == 0) THEN
                CALL WRITE_IMAGE()
                IF (USE_SOLAR_EROSION) THEN
                    CALL WRITE_SHADED_RELIEF_IMAGE()
                ELSE
                    CALL WRITE_SHADED_RELIEF_IMAGE()
                    if (FLUVIAL_AND_SLOPE_MODELING.AND.USE_DISCHARGE.AND. &
                    (.NOT.WRITE_BINARY_IMAGE_FILE)) &
                        CALL WRITE_LAKE_INFO(OUTLAKE)
                ENDIF
                IF(WRITE_BINARY_IMAGE_FILE) THEN
                    CALL OUTPUT_BINARY_DATA()
                ELSE
                    IF (FLUVIAL_AND_SLOPE_MODELING.AND.USE_DISCHARGE) CALL WRITE_ROUTED_DISCHARGE()
                    IF (FLUVIAL_AND_SLOPE_MODELING.AND.USE_DISCHARGE) CALL WRITE_SUBMERGED_LOCATIONS()
                    IF (IS_GRAVEL_MIXTURE) CALL WRITE_SEDIMENT_FLUX()
                ENDIF
            ENDIF
            !        ********************************************************************
            !        If it is time to print erosion rate statistics do so
            !        ********************************************************************
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                IF(MOD(ITERATION,WRITE_CHANGE_INTERVAL) == 0) THEN
                    CALL WRITE_REPORT()
                ENDIF
                IF (MOD(ITERATION,OUTPUT_PRINT_INTERVAL) == 0) THEN
                    IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) /= 0) THEN
                        CALL PRINT_SIMULATION_INFORMATION()
                    ENDIF
                ENDIF
                IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) == 0) THEN
                    CALL WRITE_ALLUVIAL_LOCATIONS(OUTALLUV)
                    CALL WRITE_BEDROCK_LOCATIONS(OUTROCK)
                    CALL WRITE_EROSION_DEPTH_INDEX(OUTEROSION_DEPTH_INDEX)
                    CALL WRITE_DEFORMATION(OUTDEFORM)
                    CALL WRITE_SEDIMENT_BASE()
                    IF (USE_DISCHARGE) CALL WRITE_LAKE_INFO(OUTLAKE)
                    CALL WRITE_REGOLITH_THICKNESS()
                    CALL WRITE_ROCK_RESISTANCE()
                    CALL PRINT_SIMULATION_INFORMATION()
                    CALL WRITE_ACCELERATED_EROSION_STATE(OUTHIGH)
                    IF (MODEL_GROUNDWATER) THEN
                        CALL WRITE_GROUNDWATER_ELEVATION(IOTEMP1)
                        CALL WRITE_GROUNDWATER_FLOW(IOTEMP1)
                    ENDIF
                    IF (DO_AVALANCHE) THEN
                        CALL WRITE_AVALANCHE_FLUX()
                    ENDIF
                    IF(USE_FLOW_VOLUME) THEN
                        IF (THIS_IS_BINGHAM_FLOW) THEN    !ALAN OCT2019                                                           ! added by (orkan) March 31-2014
                            CALL WRITE_MASS_FLUX(1) !ALAN OCT2019
                        ELSE
                            CALL WRITE_MASS_FLUX(2)
                        ENDIF
                    ENDIF !ALAN OCT2019
                     CALL WRITE_FINAL_STATE()
                ENDIF
            ENDIF
             !       **********************************************************************
             !       Call lava flow procedure
             !       **********************************************************************
            IF (MODEL_LAVA_FLOWS) THEN
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                    LAVAPROBTOUSE=LAVA_EVENT_PROBABILITY*TIME_INCREMENT
                ELSE
                    LAVAPROBTOUSE=LAVA_EVENT_PROBABILITY
                ENDIF
                IF (LAVAPROBTOUSE < 1.0) THEN
                    PROBVALUE=RRAND()
                    IF (PROBVALUE <= LAVAPROBTOUSE) THEN
                        PROBVALUE=RRAND()
                        SOURCEUSE=INT(PROBVALUE*NUMBER_OF_LAVA_SOURCES+1.0)
                        IF (SOURCEUSE < 1) SOURCEUSE=1
                        IF (SOURCEUSE > NUMBER_OF_LAVA_SOURCES) SOURCEUSE=NUMBER_OF_LAVA_SOURCES
                        WRITE(*,7011) SOURCEUSE
                        7011           FORMAT(' LAVAFLOW FROM SOURCE ',I3)
                        CALL DO_LAVA_FLOWS(SOURCEUSE)
                        IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) == 0) THEN
                            CALL WRITE_LAVA_INFO()
                            CALL WRITE_LAVA_AGES()
                        ENDIF

                    ENDIF
                ELSE
                    NCALL=INT(EXPONENTIAL_DISTRIBUTION()*LAVAPROBTOUSE+0.5)
                    IF (NCALL >= 1) THEN
                        DO  I=1,NCALL
                            PROBVALUE=RRAND()
                            SOURCEUSE=INT(PROBVALUE*NUMBER_OF_LAVA_SOURCES+1.0)
                            IF (SOURCEUSE < 1) SOURCEUSE=1
                            IF (SOURCEUSE > NUMBER_OF_LAVA_SOURCES) SOURCEUSE=NUMBER_OF_LAVA_SOURCES
                            WRITE(*,7011)
                            CALL DO_LAVA_FLOWS(SOURCEUSE)
                        ENDDO
                    ENDIF
                ENDIF
                IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) == 0) THEN
                    CALL WRITE_LAVA_INFO()
                    CALL WRITE_LAVA_AGES()
                ENDIF
            ENDIF
            !       **********************************************************************
            !       Call eolian deposition and erosion
            !       **********************************************************************
            IF (MODEL_EOLIAN_CHANGES) THEN
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                    EOPROBTOUSE=EOLIAN_EVENT_PROBABILITY*TIME_INCREMENT
                ELSE
                    EOPROBTOUSE=EOLIAN_EVENT_PROBABILITY
                ENDIF
                IF (EOPROBTOUSE < 1.0) THEN
                    PROBVALUE=RRAND()
                    IF (PROBVALUE <= EOPROBTOUSE) THEN
                        WRITE(*,6012)EOPROBTOUSE, PROBVALUE
                        6012                FORMAT(' EOLIAN, PROB=',g12.5,' PR. ACT=',g12.5)
                        CALL DO_EOLIAN_CHANGE()
                        IF (.NOT.FLUVIAL_AND_SLOPE_MODELING) THEN
                            PRESENT_TIME=PRESENT_TIME+EOLIAN_TIME_INCREMENT
                        ENDIF
                    ENDIF
                ELSE
                    NCALL=INT(EXPONENTIAL_DISTRIBUTION()*EOPROBTOUSE+0.5)
                    IF (NCALL>5) NCALL=5
                    IF (NCALL >= 1) THEN
                        DO  I=1,NCALL
                            WRITE(*,6013)NCALL,I
                            6013 FORMAT(' NO. CALLS=',I5,' N=',i5)
                            CALL DO_EOLIAN_CHANGE()
                            IF (.NOT.FLUVIAL_AND_SLOPE_MODELING) THEN
                                PRESENT_TIME=PRESENT_TIME+EOLIAN_TIME_INCREMENT
                            ENDIF
                        ENDDO
                    ENDIF
                ENDIF
            ENDIF
            !       **********************************************************************
            !       Call impact cratering
            !       **********************************************************************
            IF (MODEL_IMPACT_CRATERING) THEN
                CALL DETERMINE_CRATERING_RATE(LOCAL_PROBABILITY)
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                    CRATERPROBTOUSE=LOCAL_PROBABILITY*TIME_INCREMENT
                ELSE
                    CRATERPROBTOUSE=LOCAL_PROBABILITY
                ENDIF
                IF (CRATERPROBTOUSE < 1.0) THEN
                    PROBVALUE=RRAND()
                    IF (PROBVALUE <= CRATERPROBTOUSE) THEN
                        IS_RANDOM=.TRUE.
                        CALL DO_IMPACT_CRATERING(IS_RANDOM,DUMMY_D,DUMMY_X,DUMMY_Y)
                        CRATERING_CALLED=.TRUE.
                    ENDIF
                ELSE
                    NCALL=INT(EXPONENTIAL_DISTRIBUTION()*CRATERPROBTOUSE+0.5)
                    IF (NCALL >= 1) THEN
                        DO  I=1,NCALL
                            IS_RANDOM=.TRUE.
                            CALL DO_IMPACT_CRATERING(IS_RANDOM,DUMMY_D,DUMMY_X,DUMMY_Y)
                        ENDDO
                    ENDIF
                ENDIF
            ENDIF
 !
!       **********************************************************************
!
!   Glacial Subroutine Added here - 02/26/2014 (orkan)
!
!       **********************************************************************
            IF (USE_FLOW_VOLUME) THEN
                IF (THIS_IS_BINGHAM_FLOW) THEN
                    WRITE(*,3666)
                    3666     FORMAT(' CALL BINGHAM GLACIAL FLOW')
                    CALL DO_BINGHAM_MASS_FLOW()
                    DO J=1,MY
                        DO I=1,MX
                            IF (REGOLITH(I,J) <= 0.0) THEN
                                REGOLITH(I,J) = 0.0
                                ENDIF
                        ENDDO
                    ENDDO
                ELSE  ! if this is not Bingham flow
                    WRITE(*,3667)
                    3667     FORMAT(' CALLing Glen Law FLOW')
                    CALL DO_GLEN_LAW_FLOW()
                        DO J=1,MY
                        DO I=1,MX
                        IF (REGOLITH(I,J) <= 0.0) THEN
                            REGOLITH(I,J) = 0.0
                            ENDIF
                        ENDDO
                        ENDDO
                ENDIF
            max_FLOW_VOLUME_INCREMENT=-1.0E25
            if (use_variable_FLOW_VOLUME_rate) then
                time_increment=TARGET_MASS_FLOW_CHANGE_RATE/max_FLOW_VOLUME_CHANGE
                if (time_increment>maximum_time_increment) time_increment=maximum_time_increment
                if (time_increment<minumum_time_increment) time_increment=minumum_time_increment
                LAST_TIME_INCREMENT=TIME_INCREMENT
            endif
                    DO J=1,MY
                    DO I=1,MX
                    FLOW_VOLUME_INCREMENT = TIME_INCREMENT*FLOW_VOLUME_CHANGE(I,J)
                    if (abs(FLOW_VOLUME_INCREMENT)>max_FLOW_VOLUME_INCREMENT) max_FLOW_VOLUME_INCREMENT=abs(FLOW_VOLUME_INCREMENT)
                    ELEVATION(I,J) = ELEVATION(I,J) + FLOW_VOLUME_INCREMENT
                    REGOLITH(I,J) = REGOLITH(I,J) + FLOW_VOLUME_INCREMENT
                    CUMULATIVE_FLOW_VOLUME_CHANGE(i,j)=CUMULATIVE_FLOW_VOLUME_CHANGE(I,J)+FLOW_VOLUME_INCREMENT
                    IF (REGOLITH(I,J) <= 0.0E0) THEN
                        REGOLITH(I,J) = 0.0E0
                        IS_SEDIMENT_COVERED(I,J) = .FALSE.
                        IS_ROCK_SURFACE(I,J) = .TRUE.
                        ENDIF
                    ENDDO
                    ENDDO
                    FLOWWORK=FLOWWORK+FLOWTOADD*TIME_INCREMENT
                    FLOWGRAV=FLOWGRAV+FLGRAVTOADD+TIME_INCREMENT
                    write(*,7344) max_FLOW_VOLUME_INCREMENT
7344                format('max glacial increment=',g12.5)
                ENDIF  ! end of mass flow if statement
        ! *************************************************************
        !   End of Glacial Subroutine - 02/26/2014  (orkan)
        ! *************************************************************
            !       **********************************************************************
            !       Call accretion and ablation
            !       **********************************************************************
            IF (MODEL_ACCRETION_AND_ABLATION) THEN
                !WRITE(*,3221)
                3221     FORMAT(' CALL ACCRETION')
                IF (EXPOSURE_DEPENDENT_CREEP) THEN
                    CALL DO_EXPOSURE_DEPENDENT_CREEP()
                ELSE
                    CALL DO_ACCRETION_ABLATION()
 
                ENDIF
            ENDIF
            !       **********************************************************************
            !       End accretion and ablation
            !       **********************************************************************      
            NETDIFF=0.0
            DO J=1,MY
                DO I=1,MX
                    NETDIFF=NETDIFF+ELEVATION(I,J)-INITIAL_ELEVATION(I,J)
                ENDDO
            ENDDO
            NETDIFF=NETDIFF/(MX*MY)
            WRITE(OUTHIST,798) ITERATION,PRESENT_TIME,MAXIMUM_ELEVATION_CHANGE, &
            ITERATION_MAXIMUM_SEDIMENT_YIELD, NETDIFF
            WRITE(*,798) ITERATION,PRESENT_TIME,MAXIMUM_ELEVATION_CHANGE, &
            ITERATION_MAXIMUM_SEDIMENT_YIELD, NETDIFF
            798           FORMAT(' I=',I9,' T=',G12.5,' MC=',G12.5  &
            ,' MX SED=',G12.5,' NETDIF=',F12.5)
            !       **********************************************************************
            !       Call specified events at specified times
            !       **********************************************************************
            IF (DO_EVENTS) THEN
                IF (EVENT_TYPE == 1) THEN
                    IF(K >= EVENT_ITERATIONS(EVENT_INDEX)) THEN
                        IF (.NOT.EVENT_DONE(EVENT_INDEX)) THEN
                            CALL MAKE_EVENT(EVENT_INDEX)
                            EVENT_DONE(EVENT_INDEX)=.TRUE.
                            EVENT_INDEX=EVENT_INDEX+1
                            IF ((EVENT_INDEX > 10).OR. &
                            (EVENT_INDEX > NUMBER_OF_EVENTS))  &
                            DO_EVENTS=.FALSE.
                        ENDIF
                    ENDIF
                ELSE
                    DO WHILE(PRESENT_TIME >= EVENT_TIMES(EVENT_INDEX))
                        IF (.NOT.EVENT_DONE(EVENT_INDEX)) THEN
                            CALL MAKE_EVENT(EVENT_INDEX)
                            EVENT_DONE(EVENT_INDEX)=.TRUE.
                            EVENT_INDEX=EVENT_INDEX+1
                            IF (EVENT_INDEX > NUMBER_OF_EVENTS) THEN
                                DO_EVENTS=.FALSE.
                                EXIT
                            ENDIF
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
            IF (PRESENT_TIME >= MAXIMUM_SIMULATION_TIME) EXIT L200
            IF (.NOT.FLUVIAL_AND_SLOPE_MODELING) THEN
                IF(MOD(ITERATION,10) == 0) THEN
                    OPEN(INSTOP,FILE='ERODE.STOP')
                    READ(INSTOP,*) IDOSTOP
                    CLOSE(INSTOP)
                    IF (IDOSTOP > 0) THEN
                        CALL WRITE_IMAGE()
                        IF (USE_SOLAR_EROSION) THEN
                            CALL WRITE_SHADED_RELIEF_IMAGE()
                        ELSE
                            CALL WRITE_SHADED_RELIEF_IMAGE()
                        ENDIF
                        IF(WRITE_BINARY_IMAGE_FILE) CALL OUTPUT_BINARY_DATA()
                        IF (MODEL_LAVA_FLOWS) THEN
                            CALL WRITE_LAVA_INFO()
                            CALL WRITE_LAVA_AGES()
                        ENDIF
                      IF(USE_FLOW_VOLUME) THEN
                        IF (THIS_IS_BINGHAM_FLOW) THEN    !ALAN OCT2019                                                           ! added by (orkan) March 31-2014
                            CALL WRITE_MASS_FLUX(1) !ALAN OCT2019
                        ELSE
                            CALL WRITE_MASS_FLUX(2)
                        ENDIF
                       ENDIF
                        CALL WRITE_NET_CHANGE_MATRICES()
                        CLOSE(OUTHIST)
                        CALL WRITE_ELEVATION_MATRIX()
                        CALL WRITE_FINAL_STATE()
                        STOP
                    ENDIF
                ENDIF
            ENDIF
        ENDDO L200
        !       **********************************************************************
        !       END OF MASTER ITERACTION CYCLING
        !       DO END-OF-SIMULATION BOOKKEEPING
        !       **********************************************************************
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            IF (HORIZONTAL_LOWER_BOUNDARY.AND.(.NOT.IS_Y_PERIODIC)) THEN
                BOUNDQ=0.0
                TOTQ=0.0
                DO I=1,MX
                    BOUNDQ=BOUNDQ+DISCHARGE(I,MY-1)
                    DO J=1,MY-1
                        TOTQ=TOTQ+CONVERT_DISCHARGE*CELL_AREA
                    ENDDO
                ENDDO
                WRITE(OUTHIST,4398) TOTQ,BOUNDQ,CONVERT_DISCHARGE,CELL_AREA
                WRITE(*,4398) TOTQ,BOUNDQ,CONVERT_DISCHARGE,CELL_AREA
                4398  FORMAT('TOTQ=',G12.5,' BOUNDQ=',G12.5, 'CONVERT_DISCHARGE=',G12.5,'CELL_AREA=',G12.5)
            ENDIF
            IF(MOD(ITERATION,WRITE_CHANGE_INTERVAL) /= 0) THEN
                IF (CHANGECOUNT > 0.0) THEN
                    SUMGRADCHANGE=SUMGRADCHANGE/  &
                    (CHANGECOUNT*(MX-1)*(MY-1))
                    ABSGRADCHANGE=ABSGRADCHANGE/  &
                    (CHANGECOUNT*(MX-1)*(MY-1))
                    NUMBIDCHANGE=NUMBIDCHANGE/(CHANGECOUNT* &
                    (MX-1)*(MY-1))
                    WRITE(OUTRECORD,580) PRESENT_TIME,SUMGRADCHANGE, &
                    ABSGRADCHANGE, &
                    NUMBIDCHANGE
                ENDIF
            ENDIF
            580             FORMAT(4(' ',E15.7))
            IF(MOD(ITERATION,WRITE_CHANGE_INTERVAL) /= 0) THEN
                CALL WRITE_REPORT()
            ENDIF
            IF (MOD(ITERATION,OUTPUT_PRINT_INTERVAL) /= 0) THEN
                IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) /= 0) THEN
                    CALL CHANNEL_PROPERTIES()
                    CALL WRITE_ELEVATION_MATRIX()
                     IF (.NOT.WRITE_BINARY_IMAGE_FILE) CALL WRITE_DISCHARGES()
                    CALL WRITE_ALLUVIAL_LOCATIONS(OUTALLUV)
                    CALL WRITE_BEDROCK_LOCATIONS(OUTROCK)
                    CALL WRITE_EROSION_DEPTH_INDEX(OUTEROSION_DEPTH_INDEX)
                    CALL WRITE_DEFORMATION(OUTDEFORM)
                    CALL WRITE_SEDIMENT_BASE()
                    CALL WRITE_REGOLITH_THICKNESS()
                    CALL WRITE_ROCK_RESISTANCE()
                    CALL PRINT_SIMULATION_INFORMATION()
                    IF(USE_FLOW_VOLUME) THEN
                        IF (THIS_IS_BINGHAM_FLOW) THEN    !ALAN OCT2019                                                           ! added by (orkan) March 31-2014
                            CALL WRITE_MASS_FLUX(1) !ALAN OCT2019
                        ELSE
                            CALL WRITE_MASS_FLUX(2)
                        ENDIF
                    ENDIF
                    CALL WRITE_ACCELERATED_EROSION_STATE(OUTHIGH)
                    IF (MODEL_GROUNDWATER) THEN
                        CALL WRITE_GROUNDWATER_ELEVATION(IOTEMP1)
                        CALL WRITE_GROUNDWATER_FLOW(IOTEMP1)
                    ENDIF
                    IF (DO_AVALANCHE) THEN
                        CALL WRITE_AVALANCHE_FLUX()
                    ENDIF
                    IF (DO_MORPHOMETRY) THEN
                                            write (*, 8012)
8012 format('print_morphometry')
                    CALL PRINT_MORPHOMETRY()
                    write(*,8013)
8013 format('topo_divergence')
                    CALL CALCULATE_TOPO_DIVERGENCE()
                    WRITE_OUT_CHANNELS=.TRUE.
                    write(*,8014)
8014 format('summarize_channels')
                        CALL SUMMARIZE_CHANNELS()
                    ENDIF
                ENDIF
                IF (MOD(ITERATION,IMAGE_OUTPUT_INTERVAL) /= 0) THEN
                    CALL WRITE_IMAGE()
                    IF (USE_SOLAR_EROSION) THEN
                        CALL WRITE_SHADED_RELIEF_IMAGE()
                    ELSE
                        CALL WRITE_SHADED_RELIEF_IMAGE()
                        IF (.NOT.WRITE_BINARY_IMAGE_FILE) CALL WRITE_LAKE_INFO(OUTLAKE)
                    ENDIF
                    IF(WRITE_BINARY_IMAGE_FILE) THEN
                        CALL OUTPUT_BINARY_DATA()
                    ELSE
                        IF (FLUVIAL_AND_SLOPE_MODELING.AND.USE_DISCHARGE) CALL WRITE_ROUTED_DISCHARGE()
                        IF (FLUVIAL_AND_SLOPE_MODELING.AND.USE_DISCHARGE) CALL WRITE_SUBMERGED_LOCATIONS()
                    ENDIF
                ENDIF
                IF (MODEL_LAVA_FLOWS) THEN
                    IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) /= 0) THEN
                       CALL WRITE_LAVA_INFO()
                       CALL WRITE_LAVA_AGES()
                    ENDIF
                ENDIF
            ENDIF
        ELSE
            CALL WRITE_FINAL_STATE()
            CALL WRITE_NET_CHANGE_MATRICES()
        ENDIF
        WRITE(OUTHIST,32289) ELAPSED_TIME
        WRITE(*,32289) ELAPSED_TIME
32289   FORMAT('ELAPSED TIME=',G13.6)
        IF (PLOT_DIAGNOSTICS.GT.0) CALL WRITE_PLOTS()
        CALL FINALIZE_FLUVIAL_SLOPE_EROSION()
        CLOSE(OUTHIST)
        STOP
    END ! PROGRAM MARSSIM
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE SETUP_FLUVIAL_SLOPE_EROSION
        USE ERODE_GLOBALS
        USE GROUNDWATER_VARIABLES
        IMPLICIT NONE
        REAL(4) :: LOGSD,NORMVAR,NORMMEAN,NORMSD,RRRR,LOGNORMAL_RANDOM_DEVIATE1,NORMAL_RANDOM_DEVIATE
        INTEGER :: INSTOP, I,J
        LOGICAL :: TEMPORARY_TRUE
        !      ********************************************************************
        !   MODIFIES: CFW, ABORT_SIMULATION,TEMPQCONSTANT,LASTQCONSTANT, DISCHARGE_SCALE_FACTOR
        !             EVAPORATION_SCALE, EVAPORATION_RATE
        !
        !   CALLS: INITIALIZE_VARIABLES, CALCULATE_DIVERGENCE, EXPONENTIAL_HYDR_COND_GRNDWTR
        !          FIND_GROUNDWATER_FLUX, CONSTANT_HYDR_COND_GRNDWTR,DRAINAGE_BASIN_LAKE_FLOW
        !          DRAINAGE_BASIN_AREA_FLOW, DETERMINE_EROSION_RATE
        !
        !      ********************************************************************
        !      Set initial parameter values
        !      ********************************************************************
        CALL INITIALIZE_VARIABLES()
        !      ********************************************************************
        !      Abort_simulation is true if run is to be aborted
        !      ********************************************************************
        ABORT_SIMULATION = .FALSE.
        INSTOP=31
        !     **********************************************************************
        !      Determine topographic divergence
        !     **********************************************************************
        IF (FLUVIAL_AND_SLOPE_MODELING) CALL CALCULATE_DIVERGENCE()
        !      ********************************************************************
        !      Determine flow directions and gradients, and reset sed transport
        !         matrix
        !      ********************************************************************
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            CFW=0.0
        ENDIF
        !      ********************************************************************
        !      Calculate groundwater flow
        !      ********************************************************************
        IF (MODEL_GROUNDWATER) THEN
            FIRST_GROUNDWATER_CALL=.TRUE.
            IF (EXPONENTIAL_PERMEABILITY_DECAY) THEN
                CALL EXPONENTIAL_HYDR_COND_GRNDWTR()
            ELSE
                CALL CONSTANT_HYDR_COND_GRNDWTR()
            ENDIF
            CALL FIND_GROUNDWATER_FLUX()
            CALL WRITE_GROUNDWATER_ELEVATION(IOTEMP1)
            CALL WRITE_GROUNDWATER_FLOW(IOTEMP1)
            FIRST_GROUNDWATER_CALL=.FALSE.
        ENDIF
        !      ********************************************************************
        !      Determine discharge scaling
        !      ********************************************************************
        IF (USE_RANDOM_DISCHARGE) THEN
            LOGSD=DISCHARGE_CONSTANT*DISCHARGE_COEFF_VARIATION*OCORRECT
            NORMVAR=LOG((LOGSD**2+DISCHARGE_CONSTANT**2)/DISCHARGE_CONSTANT**2)
            NORMMEAN=LOG(DISCHARGE_CONSTANT**2/SQRT(DISCHARGE_CONSTANT**2+LOGSD**2))
            NORMSD=SQRT(NORMVAR)
            WRITE(*,8211) NORMMEAN,NORMSD
            8211   FORMAT(' NORMMEAN=',G12.5,' NORMSD=',G12.5)
            RRRR=LOGNORMAL_RANDOM_DEVIATE1(NORMMEAN,NORMSD)
            TEMPQCONSTANT=LASTQCONSTANT*(1.0-OMEGA_WEIGHT)  &
            +OMEGA_WEIGHT*RRRR
            LASTQCONSTANT=TEMPQCONSTANT
        ELSE
            TEMPQCONSTANT=DISCHARGE_CONSTANT
            LASTQCONSTANT=DISCHARGE_CONSTANT
        ENDIF
        !      ********************************************************************
        !      Determine discharges and lake levels
        !      ********************************************************************
        DISCHARGE_SCALE_FACTOR=(TEMPQCONSTANT)**(1.0/DISCHARGE_EXPONENT)
        STEADY_DISCHARGE_FACTOR=DISCHARGE_CONSTANT/TEMPQCONSTANT
        CALL GRADIENT_AND_FLOW_DIRECTION()
        IF (FLUVIAL_AND_SLOPE_MODELING.AND.USE_DISCHARGE) THEN
            IF (DO_SEDIMENT_TRANSPORT.OR.DO_FLUVIAL_DETACHMENT) THEN
                IF (MODEL_LAKE_EVAPORATION.AND.USE_DISCHARGE) THEN
                   ! EVAPORATION_SCALE=(EVAPORATION_MEAN+EVAPORATION_STANDARD_DEVIATION*NORMAL_RANDOM_DEVIATE())
                    EVAPORATION_SCALE=LOGNORMAL_RANDOM_DEVIATE1(LOG_EVAP_MEAN,LOG_EVAP_SD)
                    IF (EVAPORATION_SCALE < 0.0) EVAPORATION_SCALE=0.0
                    EVAPORATION_RATE=DISCHARGE_SCALE_FACTOR*EVAPORATION_SCALE
                    WRITE(*,8298) EVAPORATION_SCALE,EVAPORATION_RATE,DISCHARGE_SCALE_FACTOR
                    8298 FORMAT('EVAPORATION_SCALE=',G12.5,' EVAPORATION_RATE=',G12.5,' DSCHARGE_SCALE_FACTOR=',G12.5)
                   WRITE(OUTHIST,8298) EVAPORATION_SCALE,EVAPORATION_RATE,DISCHARGE_SCALE_FACTOR
                    CALL DRAINAGE_BASIN_LAKE_FLOW()
                ELSE
                    CALL DRAINAGE_BASIN_AREA_FLOW()
                ENDIF
            ENDIF
        ENDIF
        !      ********************************************************************
        !      Stop simulation if something has gone wrong in subroutine drbasins
        !      ********************************************************************
        IF (ABORT_SIMULATION) THEN
            WRITE(OUTHIST,1777)
            WRITE(*,1777)
            1777    FORMAT(' ABORTING DUE TO DRBASINS ERROR')
            STOP
        ENDIF
        !     **********************************************************************
        !      Find present value of erosion rate and simulation parameters if they
        !         are time-varying
        !     **********************************************************************
        IF (VARIABLE_EROSION_RATE) CALL DETERMINE_EROSION_RATE()
        RETURN
    END ! SUBROUTINE SETUP_FLUVIAL_SLOPE_EROSION
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DO_FLUVIAL_AND_SLOPE()
        USE GRAVEL_MIXTURE_GLOBALS
        USE ERODE_GLOBALS
        USE GROUNDWATER_VARIABLES
		USE MASS_FLOW_GLOBALS
        IMPLICIT NONE
        INTEGER I,J,ILU,IUM
        REAL(4) :: LOGSD,NORMVAR,NORMMEAN,NORMSD,RRRR,LOGNORMAL_RANDOM_DEVIATE1,NORMAL_RANDOM_DEVIATE
        REAL(4) :: THEAVG,THEMIN,THEMAX,FRACTION_ACTIVE
        INTEGER INSTOP,IDOSTOP
        !      ********************************************************************
        !      This routine is called each iteration to do fluvial and slope erosion
        !  MODIFIES:  ABORT_SIMULATION,TEMPQCONSTANT,LASTQCONSTANT,DISCHARGE_SCALE_FACTOR
        !             EVAPORATION_RATE, EVAPORATION_SCALE,, PREVIOUS_TIME_INCREMENT,CHANGECOUNT
        !
        !  CALLS:  DETERMINE_ERODIBILITY, DETERMINE_EROSION_RATE, EXPONENTIAL_HYDR_COND_GRNDWTR
        !          CONSTANT_HYDR_COND_GRNDWTR, DRAINAGE_BASIN_LAKE_FLOW
        !          DRAINAGE_BASIN_AREA_FLOW, CALCULATE_DIVERGENCE, IS_IT_SUBMERGED
        !          PELAGIC_DEPOSIT, DO_THE_EROSION, PRINT_MORPHOMETRY, SUMMARIZE_CHANNELS
        !          CALCULATE_TOPO_CONVERGENCE, WRITE_ELEVATION_MATRIX, WRITE_ALLUVIAL_LOCATIONS
        !          WRITE_BEDROCK_LOCATIONS, WRITE_LAKE_INFO, WRITE_EROSION_DEPTH_INDEX
        !          WRITE_DEFORMATION, WRITE_SEDIMENT_BASE, WRITE_REGOLITH_THICKNESS
        !          WRITE_ROCK_RESISTANCE, PRINT_SIMULATION_INFORMATION, WRITE_ACCELERATED_EROSION_STATE
        !          WRITE_SHADED_RELIEF_IMAGE, WRITE_IMAGE, WRITE_GROUNDWATER_ELEVATION,
        !          WRITE_GROUNDWATER_FLOW, WRITE_SECOND_DATA_MATRIX

        !      ********************************************************************
        !      Abort_simulation is true if run is to be aborted
        !      ********************************************************************
        ABORT_SIMULATION = .FALSE.
        INSTOP=31
        !     **********************************************************************
        !      Find matrix of erodibilities if spatio-temporal variability of bedrock
        !             resistance is used
        !     **********************************************************************
        IF (USE_3D_SLOPE_RESISTANCE) THEN
            CALL DETERMINE_ERODIBILITY()
        ENDIF
        !     **********************************************************************
        !      Determine amount of bedrock weathering for both regolith-mantled and
        !             bedrock slopes
        !     **********************************************************************
        ICLAST=ICENT
        JCLAST=JCENT
        !     **********************************************************************
        !      Find present value of erosion rate and simulation parameters if they
        !         are time-varying
        !     **********************************************************************
        IF (VARIABLE_EROSION_RATE) CALL DETERMINE_EROSION_RATE()
        !      ********************************************************************
        !      Calculate groundwater flow
        !      ********************************************************************
        IF (MODEL_GROUNDWATER) THEN
            IF (MOD(ITERATION,SEEPAGE_ITERATION_INTERVAL) == 0) THEN
                IF (EXPONENTIAL_PERMEABILITY_DECAY) THEN
                    CALL EXPONENTIAL_HYDR_COND_GRNDWTR()
                ELSE
                    CALL CONSTANT_HYDR_COND_GRNDWTR()
                ENDIF
            ENDIF
        ENDIF
        !      ********************************************************************
        !      Determine discharge scaling
        !      ********************************************************************
        IF (MOD(ITERATION,RECALCULATE_GRADIENT_INTERVAL) == 0) THEN
            IF (USE_RANDOM_DISCHARGE) THEN
                LOGSD=DISCHARGE_CONSTANT*DISCHARGE_COEFF_VARIATION*OCORRECT
                NORMVAR=LOG((LOGSD**2+DISCHARGE_CONSTANT**2)/DISCHARGE_CONSTANT**2)
                NORMMEAN=LOG(DISCHARGE_CONSTANT**2/SQRT(DISCHARGE_CONSTANT**2+LOGSD**2))
                NORMSD=SQRT(NORMVAR)
                RRRR=LOGNORMAL_RANDOM_DEVIATE1(NORMMEAN,NORMSD)
                TEMPQCONSTANT=LASTQCONSTANT*(1.0-OMEGA_WEIGHT)  &
                +OMEGA_WEIGHT*RRRR
                LASTQCONSTANT=TEMPQCONSTANT
            ELSE
                TEMPQCONSTANT=DISCHARGE_CONSTANT
                LASTQCONSTANT=DISCHARGE_CONSTANT
            ENDIF
            DISCHARGE_SCALE_FACTOR=(TEMPQCONSTANT)**(1.0/DISCHARGE_EXPONENT)
            STEADY_DISCHARGE_FACTOR=DISCHARGE_CONSTANT/TEMPQCONSTANT
            4337      FORMAT('QS=',G12.5)
        !      ********************************************************************
        !      Determine discharges and lake levels
        !      ********************************************************************
            IF (FLUVIAL_AND_SLOPE_MODELING.AND.USE_DISCHARGE) THEN
                IF (DO_SEDIMENT_TRANSPORT.OR.DO_FLUVIAL_DETACHMENT) THEN
                  IF(MOD(ITERATION,RECALCULATE_DISCHARGE_INTERVAL) == 0) THEN
                    IF (MODEL_LAKE_EVAPORATION) THEN
                        IF (MOD(ITERATION,NCALCEVAP) == 0) THEN
                            IF(MOD(ITERATION,NCHANGEEVAP) == 0) THEN
                                EVAPORATION_SCALE=LOGNORMAL_RANDOM_DEVIATE1(LOG_EVAP_MEAN,LOG_EVAP_SD)
                                IF (EVAPORATION_SCALE < 0.0) EVAPORATION_SCALE=0.0
                                WRITE(*,8298) EVAPORATION_SCALE
                                WRITE(OUTHIST,8298) EVAPORATION_SCALE
                                8298 FORMAT('EVAPORATION_RATE=',G12.5)
                            ENDIF
                            EVAPORATION_RATE=DISCHARGE_SCALE_FACTOR*EVAPORATION_SCALE
                            CALL DRAINAGE_BASIN_LAKE_FLOW()
                        ENDIF
                    ELSE
                        CALL DRAINAGE_BASIN_AREA_FLOW()
                    ENDIF
                  ENDIF
                ENDIF
            ENDIF
        ENDIF
        IF (MOD(ITERATION,DIVERGEINTERVAL) == 0) THEN
            if (fluvial_and_slope_modeling) then
                 if ((integer_divergence_range==5).or.(integer_divergence_range==7).or. &
                     (integer_divergence_range==9)) then 
                     call new_calculate_divergence(INTEGER_DIVERGENCE_RANGE)
                     !write(*,12345) Integer_divergence_range
12345                format ('div, n=',I5)
                 else 
                     CALL CALCULATE_DIVERGENCE()
                     !write(*,23456)
23456                format('div3')
                 endif
          endif
        ENDIF
        !      ********************************************************************
        !      Do the actual erosion
        !      ********************************************************************
        IF (MODEL_PELAGIC_DEPOSITION) PELAGIC_SEDIMENT_VOLUME=0.0
        If (IS_GRAVEL_MIXTURE) then
           !
            time_increment=minumum_time_increment

            call DRAINAGE_BASIN_AREA_FLOW()
            aterm_avg=0.0
            avg_change=0.0
            a1_avg=0.0
            a2_avg=0.0
            agterm_avg=0.0
            gd1_avg=0.0
            gd2_avg=0.0
            Fnew_avg=0.0
            F_diff=0.0
            net1_avg=0.0
            net2_avg=0.0
            el_avg=0.0
            sf_avg=0.0
            ew_avg=0.0
            DO GRAVEL_ITERATION=1, STEP_BY_MARSSIM_ITERATION
                CALL GRAVEL_MIXTURE_TRANSPORT
            END DO
            PRESENT_TIME = PRESENT_TIME +TIME_INCREMENT

            If (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) == 0) then
                CALL write_output_gravel
            END IF

        ELSE
            CALL DO_THE_EROSION()
            IF (FLUVIAL_AND_SLOPE_MODELING.AND.DO_FLUVIAL_DETACHMENT) THEN
                IF (DETACHMENT_CRITICAL_SHEAR>0.0) THEN
                    FRACTION_ACTIVE=FLOAT(NUMBER_ABOVE_THRESHOLD)/ &
                        (FLOAT(NUMBER_BELOW_THRESHOLD)+FLOAT(NUMBER_ABOVE_THRESHOLD))
                    WRITE (*,734) FRACTION_ACTIVE
734                 FORMAT('FRACTION OF BEDROCK CHANNELS ABOVE THRESHOLD=',G13.6)
                ENDIF
            ENDIF
        END IF
        IF (FLUVIAL_AND_SLOPE_MODELING.AND.USE_DISCHARGE) THEN
                IF (DO_SEDIMENT_TRANSPORT.OR.DO_FLUVIAL_DETACHMENT) THEN
                   IF (MODEL_PELAGIC_DEPOSITION) THEN
                      CALL IS_IT_SUBMERGED()
                      CALL PELAGIC_DEPOSIT()
                   ENDIF
                ENDIF
        ENDIF
        PREVIOUS_TIME_INCREMENT=TIME_INCREMENT
        IF (.NOT.ABORT_SIMULATION) THEN
            !      ********************************************************************
            !      If it is time to calculate basin statistics do so
            !      ********************************************************************
            IF (MOD(ITERATION,OUTPUT_PRINT_INTERVAL) == 0) THEN
                IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) /= 0) THEN
                    IF (DO_MORPHOMETRY) THEN
                        write(*,7743)
7743 format('print morphometry')
                        CALL PRINT_MORPHOMETRY()
                        write (*,7744)
7744 format('topo_divergence')
                        CALL CALCULATE_TOPO_DIVERGENCE()
                        WRITE_OUT_CHANNELS=.FALSE.
                        write(*,7745)
7745 format('summarize channels')
                        CALL SUMMARIZE_CHANNELS()
                    ENDIF
                ENDIF
            ENDIF
            !      ********************************************************************
            !      Initialize simulation summary variables
            !      ********************************************************************
            IF (ITERATION < 2) THEN
                NUMBIDCHANGE=0.0
                SUMGRADCHANGE=0.0
                ABSGRADCHANGE=0.0
                CHANGECOUNT=0.0
            ENDIF
            !      ********************************************************************
            !      If it is time to summarize basin changes then do so
            !      ********************************************************************
            IF(MOD(ITERATION,WRITE_CHANGE_INTERVAL) == 0) THEN
                IF (CHANGECOUNT > 0.0) THEN
                    SUMGRADCHANGE=SUMGRADCHANGE/  &
                    (CHANGECOUNT*(MX-1)*(MY-1))
                    ABSGRADCHANGE=ABSGRADCHANGE/  &
                    (CHANGECOUNT*(MX-1)*(MY-1))
                    NUMBIDCHANGE=NUMBIDCHANGE/(CHANGECOUNT*(MX-1)*(MY-1))
                    WRITE(OUTRECORD,580) PRESENT_TIME,SUMGRADCHANGE, &
                    ABSGRADCHANGE, &
                    NUMBIDCHANGE
                    580             FORMAT(4(' ',E15.7))
                ENDIF
                NUMBIDCHANGE=0.0
                SUMGRADCHANGE=0.0
                ABSGRADCHANGE=0.0
                CHANGECOUNT=0.0
            ENDIF
            !      ********************************************************************
            !      Changecount is number of iterations over which changes are
            !        summarized
            !      ********************************************************************
            CHANGECOUNT=CHANGECOUNT+1.0
            !      ********************************************************************
            !      Check to see if operator intervention is requesting a stop
            !      ********************************************************************
            IF(MOD(ITERATION,10) == 0) THEN
                OPEN(INSTOP,FILE='ERODE.STOP')
                READ(INSTOP,*) IDOSTOP
                CLOSE(INSTOP)
                IF (IDOSTOP > 0) THEN
                    IF(MOD(ITERATION,WRITE_CHANGE_INTERVAL) /= 0) THEN
                        IF (CHANGECOUNT > 0.0) THEN
                            SUMGRADCHANGE=SUMGRADCHANGE/  &
                            (CHANGECOUNT*(MX-1)*(MY-1))
                            ABSGRADCHANGE=ABSGRADCHANGE/  &
                            (CHANGECOUNT*(MX-1)*(MY-1))
                            NUMBIDCHANGE=NUMBIDCHANGE/(CHANGECOUNT* &
                            (MX-1)*(MY-1))
                            WRITE(OUTRECORD,580) PRESENT_TIME,SUMGRADCHANGE, &
                            ABSGRADCHANGE,  &
                            NUMBIDCHANGE
                        ENDIF
                    ENDIF
                    IF(MOD(ITERATION,WRITE_CHANGE_INTERVAL) /= 0) THEN
                        CALL WRITE_REPORT()
                    ENDIF
                    IF (MOD(ITERATION,OUTPUT_PRINT_INTERVAL) /= 0) THEN
                        IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) /= 0) THEN
                            IF (FLUVIAL_AND_SLOPE_MODELING.AND.USE_DISCHARGE) CALL CHANNEL_PROPERTIES()
                            CALL WRITE_ELEVATION_MATRIX()
                            IF (USE_DISCHARGE) CALL WRITE_ALLUVIAL_LOCATIONS(OUTALLUV)
                            CALL WRITE_BEDROCK_LOCATIONS(OUTROCK)
                            CALL WRITE_EROSION_DEPTH_INDEX(OUTEROSION_DEPTH_INDEX)
                            CALL WRITE_DEFORMATION(OUTDEFORM)
                            CALL WRITE_SEDIMENT_BASE()
                            CALL WRITE_REGOLITH_THICKNESS()
                            CALL WRITE_ROCK_RESISTANCE()
                            CALL PRINT_SIMULATION_INFORMATION()
                            CALL WRITE_ACCELERATED_EROSION_STATE(OUTHIGH)
		                    IF(USE_FLOW_VOLUME) THEN
                                IF (THIS_IS_BINGHAM_FLOW) THEN    !ALAN OCT2019                                                           ! added by (orkan) March 31-2014
                                    CALL WRITE_MASS_FLUX(1) !ALAN OCT2019
                                ELSE
                                    CALL WRITE_MASS_FLUX(2)
                                ENDIF
                            ENDIF
                            IF (MODEL_GROUNDWATER) THEN
                                CALL WRITE_GROUNDWATER_ELEVATION(IOTEMP1)
                                CALL WRITE_GROUNDWATER_FLOW(IOTEMP1)
                            ENDIF
                            IF (DO_AVALANCHE) THEN
                                CALL WRITE_AVALANCHE_FLUX()
                            ENDIF
                            IF (DO_MORPHOMETRY) THEN
                        write(*,8743)
8743 format('print morphometry')
                        CALL PRINT_MORPHOMETRY()
                        write (*,8744)
8744 format('topo_divergence')
                        CALL CALCULATE_TOPO_DIVERGENCE()
                        WRITE_OUT_CHANNELS=.TRUE.
                        write(*,8745)
8745 format('summarize channels')
                                CALL SUMMARIZE_CHANNELS()
                            ENDIF
                            IF (MOD(ITERATION,IMAGE_OUTPUT_INTERVAL) /= 0) THEN
                                CALL WRITE_IMAGE()
                                CALL WRITE_SHADED_RELIEF_IMAGE()
                                IF (.NOT.WRITE_BINARY_IMAGE_FILE) CALL WRITE_LAKE_INFO(OUTLAKE)
                                IF(WRITE_BINARY_IMAGE_FILE) CALL OUTPUT_BINARY_DATA()
                            ENDIF
                        ENDIF
                    ENDIF
                    CLOSE(OUTHIST)
                    CLOSE(OUTCHAN)
                    CLOSE(OUTRECORD)
                    CLOSE(OUTSUMMARY)
                    CLOSE(OUTSOURCE)
                    CLOSE(OUTREPORT)
                    CLOSE(IOTEMP1)
                    CLOSE(IOTEMP2)
                    CLOSE(OUTCRATER)
                    IF (WRITE_BINARY_IMAGE_FILE) THEN
                        CLOSE(OUT_IMAGE_FILE)
                        IF (FLUVIAL_AND_SLOPE_MODELING.AND.COMPLETE_RUNOFF) CLOSE(OUTSUBMERGE)
                        IF (MODEL_IMPACT_CRATERING.OR.DO_EVENTS) CLOSE(OUTCRATERS)
                        IF (FLUVIAL_AND_SLOPE_MODELING.AND.(.NOT.WRITE_BINARY_IMAGE_FILE)) CLOSE(OUTDISCHARGES)
                    ENDIF
                    CALL WRITE_NET_CHANGE_MATRICES()
                    CALL WRITE_FINAL_STATE()
                    STOP
                ENDIF
            ENDIF

            !      ********************************************************************
            !      Recycle for another iteration if everything is ok and we are
            !           not done
            !      ********************************************************************
            IF ((MAXGRADIENT < 10000.0*CELL_SIZE) &
            .AND.(.NOT.ABORT_SIMULATION)) RETURN
        ENDIF
        IF (ABORT_SIMULATION) THEN
            !      ********************************************************************
            !      Something is wrong - print out summary information and abort
            !      ********************************************************************
            WRITE(OUTHIST,501)
            501     FORMAT('*** ABORTED DUE TO SINKHOLES ***')
            IF (FLUVIAL_AND_SLOPE_MODELING.AND.USE_DISCHARGE) CALL CHANNEL_PROPERTIES()
            CALL WRITE_ELEVATION_MATRIX()
            IF (USE_DISCHARGE) CALL WRITE_ALLUVIAL_LOCATIONS(OUTALLUV)
            CALL WRITE_BEDROCK_LOCATIONS(OUTROCK)
            CALL WRITE_EROSION_DEPTH_INDEX(OUTEROSION_DEPTH_INDEX)
            CALL WRITE_DEFORMATION(OUTDEFORM)
            CALL WRITE_SEDIMENT_BASE()
            CALL WRITE_REGOLITH_THICKNESS()
            CALL WRITE_ROCK_RESISTANCE()
            CALL PRINT_SIMULATION_INFORMATION()
            CALL WRITE_ACCELERATED_EROSION_STATE(OUTHIGH)
            CALL WRITE_SHADED_RELIEF_IMAGE()
            IF ((.NOT.WRITE_BINARY_IMAGE_FILE).AND.USE_DISCHARGE) CALL WRITE_LAKE_INFO(OUTLAKE)
            IF(WRITE_BINARY_IMAGE_FILE) CALL OUTPUT_BINARY_DATA()
            CALL WRITE_IMAGE()
            IF (MODEL_GROUNDWATER) THEN
                CALL WRITE_GROUNDWATER_ELEVATION(IOTEMP1)
                CALL WRITE_GROUNDWATER_FLOW(IOTEMP1)
            ENDIF
            IF (DO_AVALANCHE) THEN
                CALL WRITE_AVALANCHE_FLUX()
            ENDIF
            IF (DO_MORPHOMETRY) THEN
                CALL PRINT_MORPHOMETRY()
                WRITE_OUT_CHANNELS=.TRUE.
                CALL CALCULATE_TOPO_DIVERGENCE()
                CALL SUMMARIZE_CHANNELS()
            ENDIF
            ILU=-1
            IUM=MX+1
            IABORTMIN = IABORTMIN-1
            IF (IABORTMIN < 1) THEN
                IUM=MX-(1-IABORTMIN)
                IABORTMIN=1
            ENDIF
            IABORTMAX = IABORTMAX+1
            IF (IABORTMAX  >  MX) THEN
                ILU=(IABORTMAX-MX)
                IABORTMAX=MX
            ENDIF
            JABORTMIN = JABORTMIN-1
            IF (JABORTMIN < 1) JABORTMIN=1
            JABORTMAX = JABORTMAX+1
            IF (JABORTMAX  >  MY) JABORTMAX =MY
            WRITE(OUTHIST,702)
            702     FORMAT('*****  AREAS ****')
            !             write(outhist,703)
            703     FORMAT('I=     ')
            WRITE(OUTHIST,704)(I,I=IUM,MX), &
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            704     FORMAT(' I=    ',20I5,//)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,705)J,(DRAINAGE_AREA(I,J),I=IUM,MX),  &
                (DRAINAGE_AREA(I,J),I=IABORTMIN,IABORTMAX),(DRAINAGE_AREA(I,J),I=1,ILU)
                705         FORMAT(20I5)
            ENDDO
            WRITE(OUTHIST,706)
            706     FORMAT(//)
            WRITE(OUTHIST,712)
            712     FORMAT('*****  BASINS ****')
            !             write(outhist,713)
            713     FORMAT('I=     ')
            WRITE(OUTHIST,714)(I,I=IUM,MX),&
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            714     FORMAT(' I=    ',20I5,//)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,715)J,(BASIN_NUMBER(I,J),I=IUM,MX), &
                (BASIN_NUMBER(I,J),I=IABORTMIN,IABORTMAX), &
                (BASIN_NUMBER(I,J),I=1,ILU)
                715         FORMAT(20I5)
            ENDDO
            WRITE(OUTHIST,716)
            716     FORMAT(//)
            WRITE(OUTHIST,502)
            502     FORMAT('*****  DIRECTIONS ****')
            503     FORMAT('I=     ')
            WRITE(OUTHIST,504)(I,I=IUM,MX), &
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            504     FORMAT(' I=    ',20I5,//)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,505)J,(FLOW_DIRECTION(I,J),I=IUM,MX), &
                (FLOW_DIRECTION(I,J),I=IABORTMIN,IABORTMAX),(FLOW_DIRECTION(I,J),I=1,ILU)
                505         FORMAT(20I5)
            ENDDO
            WRITE(OUTHIST,506)
            506     FORMAT(//)
            WRITE(OUTHIST,507)
            507     FORMAT ('*****  ELEVATIONS ****')
            WRITE(OUTHIST,508)(I,I=IUM,MX), &
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            508     FORMAT(' I= ',10I8,//)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,509)J,(ELEVATION(I,J),I=IUM,MX), &
                (ELEVATION(I,J),I=IABORTMIN,IABORTMAX),(ELEVATION(I,J),I=1,ILU)
                509         FORMAT(I4,6(G12.4),/,3(G12.4))
            ENDDO
            WRITE(OUTHIST,506)
            WRITE(OUTHIST,599)
            599     FORMAT ('*****  GRADIENTS ****')
            WRITE(OUTHIST,503)
            WRITE(OUTHIST,508)(I,I=IUM,MX),    &
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,509)J,(D8_GRADIENT(I,J),I=IUM,MX), &
                (D8_GRADIENT(I,J),I=IABORTMIN,IABORTMAX),(D8_GRADIENT(I,J),I=1,ILU)
            ENDDO
            WRITE(OUTHIST,506)
            WRITE(OUTHIST,587)
            587     FORMAT ('*****  RESISTANCE ****')
            WRITE(OUTHIST,588)(I,I=IUM,MX), &
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            588     FORMAT(' I= ',10I8,//)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,589)J,(RELATIVE_RESISTANCE(I,J),I=IUM,MX), &
                (RELATIVE_RESISTANCE(I,J),I=IABORTMIN,IABORTMAX),&
                (RELATIVE_RESISTANCE(I,J),I=1,ILU)
                589         FORMAT(I4,6(G12.4),/,3(G12.4))
            ENDDO
            WRITE(OUTHIST,506)
            WRITE(OUTHIST,597)
            597     FORMAT ('*****  REGOLITH ****')
            WRITE(OUTHIST,528)(I,I=IUM,MX), &
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            528     FORMAT(' I= ',10I8,//)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,529)J,(REGOLITH(I,J),I=IUM,MX), &
                (REGOLITH(I,J),I=IABORTMIN,IABORTMAX),&
                (REGOLITH(I,J),I=1,ILU)
                529         FORMAT(I4,6(G12.4),/,3(G12.4))
            ENDDO
            WRITE(OUTHIST,506)
                  CALL WRITE_NET_CHANGE_MATRICES()
        ENDIF
        IF (MAXGRADIENT > 10000.0) THEN
            WRITE(OUTHIST,510)
            510     FORMAT('***ABORTED DUE TO HIGH GRADIENTS***')
            IF (FLUVIAL_AND_SLOPE_MODELING.AND.USE_DISCHARGE) CALL CHANNEL_PROPERTIES()
            CALL WRITE_ELEVATION_MATRIX()
            IF (USE_DISCHARGE) CALL WRITE_ALLUVIAL_LOCATIONS(OUTALLUV)
            CALL WRITE_BEDROCK_LOCATIONS(OUTROCK)
            CALL WRITE_EROSION_DEPTH_INDEX(OUTEROSION_DEPTH_INDEX)
            CALL WRITE_DEFORMATION(OUTDEFORM)
            CALL WRITE_SEDIMENT_BASE()
            CALL WRITE_REGOLITH_THICKNESS()
            CALL WRITE_ROCK_RESISTANCE()
            CALL PRINT_SIMULATION_INFORMATION()
            CALL WRITE_ACCELERATED_EROSION_STATE(OUTHIGH)
            CALL WRITE_SHADED_RELIEF_IMAGE()
            IF (.NOT.WRITE_BINARY_IMAGE_FILE) CALL WRITE_LAKE_INFO(OUTLAKE)
            IF (WRITE_BINARY_IMAGE_FILE) CALL OUTPUT_BINARY_DATA()
            CALL WRITE_IMAGE()
            IF (MODEL_GROUNDWATER) THEN
                CALL WRITE_GROUNDWATER_ELEVATION(IOTEMP1)
                CALL WRITE_GROUNDWATER_FLOW(IOTEMP1)
            ENDIF
            IF (DO_AVALANCHE) THEN
                    CALL WRITE_AVALANCHE_FLUX()
            ENDIF
            IF (DO_MORPHOMETRY) THEN
                CALL PRINT_MORPHOMETRY()
                WRITE_OUT_CHANNELS=.TRUE.
                CALL CALCULATE_TOPO_DIVERGENCE()
                CALL SUMMARIZE_CHANNELS()
            ENDIF
            WRITE(OUTHIST,433)
            433     FORMAT(/,' ELEVATION')
            CALL SUMMARIZE_MATRIX_DATA(ELEVATION,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,434)
            434     FORMAT(/,' ERODECHANNEL')
            CALL SUMMARIZE_MATRIX_DATA(ERODE_CHANNEL,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,435)
            435     FORMAT(/,' SEDIMENT DIVERGENCE')
            CALL SUMMARIZE_MATRIX_DATA(CFW,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,436)
            436     FORMAT(/,' SEDIMENT YIELD')
            CALL SUMMARIZE_MATRIX_DATA(SEDIMENT_YIELD,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,437)
            437     FORMAT(/,' ERODE_SLOPE')
            CALL SUMMARIZE_MATRIX_DATA(ERODE_SLOPE,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,438)
            438     FORMAT(/,' GRADIENT')
            CALL SUMMARIZE_MATRIX_DATA(D8_GRADIENT,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,439)
            439     FORMAT(' SEDBASE')
            CALL SUMMARIZE_MATRIX_DATA(SEDIMENT_BASE,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,741)
            741     FORMAT(/,' BEDROCK')
            CALL SUMMARIZE_LOGICAL_MATRIX(IS_ROCK_SURFACE)
            WRITE(OUTHIST,742)
            742     FORMAT(/,' REGOLITH')
            CALL SUMMARIZE_MATRIX_DATA(REGOLITH,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,743)
            743     FORMAT(/,' RESISTANCE')
            CALL SUMMARIZE_MATRIX_DATA(RELATIVE_RESISTANCE,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,441)
            441     FORMAT(/,' IS_SEDIMENT_COVERED')
            CALL SUMMARIZE_LOGICAL_MATRIX(IS_SEDIMENT_COVERED)
445         FORMAT (' TOTAL DELTA VOLUME CHANGE=',G12.5)
                  CALL WRITE_NET_CHANGE_MATRICES()
        ENDIF
        STOP
    END ! SUBROUTINE DO_FLUVIAL_AND_SLOPE
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE FINALIZE_FLUVIAL_SLOPE_EROSION()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER I,J,IPY,IPX
        REAL(4) :: EEEMIN
        LOGICAL WRITEDETAIL
        !     ****************************************************************
        !      Close down shop
        !   CALLS: GRAD_DISCH_WRITE
        !     **********************************************************************
        WRITEDETAIL=.FALSE.
        WRITE(OUTHIST,511) ITERATION
        511 FORMAT(' FINAL ITERATION=',I9)
        IF (FLUVIAL_AND_SLOPE_MODELING.AND.WRITEDETAIL.AND.USE_DISCHARGE) THEN
            OPEN(77,FILE='STATE.DAT')
            IF (IS_X_PERIODIC) THEN
                IPX=1
            ELSE
                IPX=0
            ENDIF
            IF (IS_Y_PERIODIC) THEN
                IPY=1
            ELSE
                IPY=0
            ENDIF
            WRITE(77,998) MX,MY,IPX,IPY
            WRITE(77,996) CELL_SIZE
            DO I=1,MX
                DO J=1,MY
                    WRITE(77,997) FLOW_DIRECTION(I,J),ELEVATION(I,J),DRAINAGE_AREA(I,J), &
                    DISCHARGE(I,J),D8_GRADIENT(I,J),BASIN_NUMBER(I,J)
                ENDDO
            ENDDO
            CLOSE(77)
            998   FORMAT(I6,' ',I6,' ',I6,' ',I6)
            996   FORMAT(G13.6)
            997   FORMAT(I6,4(' ',G13.6),' ',I7)
            CALL GRAD_DISCH_WRITE()
        ENDIF
        IF (MODEL_IMPACT_CRATERING) THEN
               WRITE(OUTHIST,8334) TOTAL_RANDOM_CRATERS
8334           FORMAT('TOTAL RANDOM CRATERS SIMULATED=',I7)
        ENDIF
        CLOSE(OUTCHAN)
        CLOSE(OUTRECORD)
        CLOSE(OUTSUMMARY)
        CLOSE(OUTSOURCE)
        CLOSE(OUTREPORT)
        CLOSE(IOTEMP1)
        CLOSE(IOTEMP2)
        CLOSE(OUTCRATER)
        CLOSE(OUTHIST)
        CLOSE(OUTWORK)
         IF (.NOT.WRITE_BINARY_IMAGE_FILE) CLOSE(OUTDISCHARGES)
        IF (WRITE_BINARY_IMAGE_FILE) THEN
            CLOSE(OUT_IMAGE_FILE)
            IF (FLUVIAL_AND_SLOPE_MODELING.AND.USE_DISCHARGE.AND. &
            (.NOT.WRITE_BINARY_IMAGE_FILE)) CLOSE(OUTDISCHARGES)
            IF (FLUVIAL_AND_SLOPE_MODELING.AND.COMPLETE_RUNOFF.AND.USE_DISCHARGE) CLOSE(OUTSUBMERGE)
            IF (MODEL_IMPACT_CRATERING.OR.DO_EVENTS) CLOSE(OUTSUBMERGE)
        ENDIF
        CALL WRITE_NET_CHANGE_MATRICES()
        CALL WRITE_FINAL_STATE()
        999   EEEMIN=1.0E+36
        DO  J=1,MY
            DO  I=1,MX
                IF (ELEVATION(I,J) < EEEMIN) EEEMIN=ELEVATION(I,J)
            ENDDO
        ENDDO
        RETURN
    END ! SUBROUTINE FINALIZE_FLUVIAL_SLOPE_EROSION
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE REPORT_MAX()
        USE ERODE_GLOBALS
        INTEGER SUBDET,SEDDET
        REAL(4) :: DIF1
        IF (SUBMERGED(IDET,JDET)) THEN
            SUBDET=1
        ELSE
            SUBDET=0
        ENDIF
        IF (IS_SEDIMENT_COVERED(IDET,JDET)) THEN
            SEDDET=1
        ELSE
            SEDDET=0
        ENDIF
        DIF1=ELEVATION(IDET,JDET)-PREVIOUS_ELEVATION(IDET,JDET)
        WRITE(OUTHIST,3992) IDET,JDET,SUBDET,SEDDET, FLOW_DIRECTION(IDET,JDET) &
        ,IDO(IDET,JDET),DISCHARGE(IDET,JDET),D8_GRADIENT(IDET,JDET),MAXIMUM_DISCHARGE(IDET,JDET) &
        ,MAXIMUM_ELEVATION_CHANGE,ELEVATION(IDET,JDET),ERODE_CHANNEL(IDET,JDET),ERODE_SLOPE(IDET,JDET) &
        ,CFN(IDET,JDET),CFNE(IDET,JDET),CFW(IDET,JDET),CFNW(IDET,JDET) &
        ,SEDIMENT_YIELD(IDET,JDET),MAXIMUM_SEDIMENT_YIELD(IDET,JDET) &
        ,PREVIOUS_ELEVATION(IDET,JDET),ERODE_REGOLITH_CHANNEL(IDET,JDET),SEDIMENT_FLUX(IDET,JDET) &
        ,EQUILIBRIUM_GRADIENT(IDET,JDET),DIF1
        3992  FORMAT('I=',I4,' J=',I4,' SUB=',I1,' SED=',I1,' ID=',I2 &
        ,' IDO=',I2,' QV=',G12.5,' S=',G12.5,'QM=',G12.5,/, &
        'MC=',G12.5,' ELEVATION=',G12.5,' ECH=',G12.5,' ESL=',G12.5, &
        ' CFN=',G12.5,' CFNE=',G12.5,' CFW=',G12.5,' CFNW=',G12.5 &
        ,/,'SY=',G12.5,' SYM=',G12.5,' PREVIOUS_ELEVATION=',G12.5,' ERG=',G12.5, &
        ' SQ=',G12.5,' SE=',G12.5,' EDIFF=',G12.5)
        RETURN
    END ! SUBROUTINE REPORT_MAX



