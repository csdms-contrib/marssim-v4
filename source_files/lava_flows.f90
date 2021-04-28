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
	SUBROUTINE READ_LAVA_PARAMETERS
	USE ERODE_GLOBALS
	USE LAVA_GLOBALS
	IMPLICIT NONE
	CHARACTER (80) :: TEXTLINE
	INTEGER :: LAVA_PARAMS, I 
	    LAVA_PARAMS=93
		OPEN(LAVA_PARAMS,FILE='LAVA_FLOW_PARAMETERS.PRM',ACTION='READ')
        !       ********************************************************************
        !
        !       Part 3:  Read in lavaflow parameters
        !       
        !       ********************************************************************
        !        *******************************************************************
        !         number_of_lava_sources is the number of lava sources in the simulation domain
        !           i_lava and j_lava hold the i and j coordinates of these sources
        !        *******************************************************************
        READ(LAVA_PARAMS,22) TEXTLINE
22      FORMAT(A)
        WRITE(*,22) TEXTLINE
        READ(LAVA_PARAMS,*) NUMBER_OF_LAVA_SOURCES
        IF (MODEL_LAVA_FLOWS) THEN
            OPEN(222,FILE="LAVA_SOURCES.TXT",ACTION='READ')
            DO I=1,NUMBER_OF_LAVA_SOURCES
                READ(222,*) I_LAVA(I),J_LAVA(I)
            ENDDO
            CLOSE(222)
        ENDIF
        !        *******************************************************************
        !         lava_event_probability is the probability of a lava flow event during a single
        !               iteration
        !        *******************************************************************
        READ(LAVA_PARAMS,*) LAVA_EVENT_PROBABILITY
        IF (FLUVIAL_AND_SLOPE_MODELING) LAVA_EVENT_PROBABILITY=LAVA_EVENT_PROBABILITY/DEFAULT_TIME_INCREMENT
        !        *******************************************************************
        !         minimum_lava_flow_slope is the minimum gradient for lava flow at the flow source
        !        *******************************************************************
        READ(LAVA_PARAMS,*) MINIMUM_LAVA_FLOW_SLOPE
        !        *******************************************************************
        !         lava_flow_thickness is the assumed thickness of individual lava flow deposits
        !          -- assumed to be spatially uniforn and temporally constant
        !        *******************************************************************
        READ(LAVA_PARAMS,*) LAVA_FLOW_THICKNESS
        !        *******************************************************************
        !         minimum_flow_thickness is the minimum thickness of lava flow that can flow into
        !           an adjoining cell
        !        *******************************************************************
        READ(LAVA_PARAMS,*) MINIMUM_FLOW_THICKNESS
        !        *******************************************************************
        !         new_segment_interval is the number if iterations between reassessment of
        !           which matrix points are potential sources for new flow segments
        !        *******************************************************************
        READ(LAVA_PARAMS,*) NEW_SEGMENT_INTERVAL
        !        *******************************************************************
        !         source_change_interval is the number of iterations between changeover
        !           between different flow sources
        !        *******************************************************************
        READ(LAVA_PARAMS,*) SOURCE_CHANGE_INTERVAL
        !        *******************************************************************
        !         eruption_stop_probability is the probability, per iteration, that the existing flow
        !           will solidify and stop being active_lava_flow -- presumably due to
        !           interruption of the flow source.  If the flow becomes inlactive
        !           a new flow starts from the source
        !        *******************************************************************
        READ(LAVA_PARAMS,*) ERUPTION_STOP_PROBABILITY
        !        *******************************************************************
        !         this is the lower limit of probability for a cell to be a source for
        !           a new flow segment.  if the probability drops below this value the
        !           cell is no longer considered a possible flow source
        !        *******************************************************************
        READ(LAVA_PARAMS,*) NO_FLOW_PROBABILITY
        !        *******************************************************************
        !         lava_gradient_weight is a parameter that determines how much the gradient
        !           between the edge of a flow and a neighboring point determines
        !           the probability of flow in that direction.
        !        *******************************************************************
        READ(LAVA_PARAMS,*) LAVA_GRADIENT_WEIGHT
        !        *******************************************************************
        !         lava_duration_weight is a parameter that determines how rapidly a new cell
        !           diminishes in probability that it can be the source of a flow
        !           to a neighboring cell -- Presumably represents effects of cooling
        !           of emplaced lava
        !        *******************************************************************
        READ(LAVA_PARAMS,*) LAVA_DURATION_WEIGHT
		CLOSE (LAVA_PARAMS)
        !        *******************************************************************
        IBACK(2)=4
        IBACK(3)=5
        IBACK(4)=2
        IBACK(5)=3
        IBACK(6)=8
        IBACK(7)=9
        IBACK(8)=6
        IBACK(9)=7	
		        !       ************************************************************************
        !        Write out lavaflow parameters
        !       ************************************************************************
        IF (MODEL_LAVA_FLOWS) THEN
            WRITE(OUTHIST,3411)
            3411 FORMAT('**********************LAVA FLOW PARAMETERS*********************')
            WRITE(OUTHIST,3412) NUMBER_OF_LAVA_SOURCES
            3412   FORMAT(' NUMBER_OF_LAVA_SOURCES=',I5)
            DO  I=1,NUMBER_OF_LAVA_SOURCES
                WRITE(OUTHIST,3414) I,I_LAVA(I),J_LAVA(I)
                3414   FORMAT(' I=',I5,' I_LAVA=',I5,' J_LAVA=',I5)
            ENDDO
            WRITE(OUTHIST,3430) MINIMUM_LAVA_FLOW_SLOPE
            3430   FORMAT(' MINIMUM_LAVA_FLOW_SLOPE=',G12.5)
            WRITE(OUTHIST,3440) LAVA_FLOW_THICKNESS
            3440    FORMAT(' LAVA_FLOW_THICKNESS=',G12.5)
            WRITE(OUTHIST,3445) MINIMUM_FLOW_THICKNESS
            3445      FORMAT(' MINIMUM_FLOW_THICKNESS=',G12.5)
            WRITE(OUTHIST,3455) NEW_SEGMENT_INTERVAL
            3455   FORMAT(' NEW_SEGMENT_INTERVAL=',I8)
            WRITE(OUTHIST,3456) SOURCE_CHANGE_INTERVAL
            3456   FORMAT(' SOURCE_CHANGE_INTERVAL=',I8)
            WRITE(OUTHIST,3480) ERUPTION_STOP_PROBABILITY
            3480    FORMAT(' PROGCLOG=',G12.5)
            WRITE(OUTHIST,3485) NO_FLOW_PROBABILITY
            3485    FORMAT(' NO_FLOW_PROBABILITY=',G12.5)
            WRITE(OUTHIST,3501) LAVA_GRADIENT_WEIGHT
            3501    FORMAT(' LAVA_GRADIENT_WEIGHT=',G12.5)
            WRITE(OUTHIST,3510) LAVA_DURATION_WEIGHT
            3510    FORMAT(' LAVA_DURATION_WEIGHT=',G12.5)
        ENDIF
	RETURN
	END !SUBROUTINE READ_LAVA_PARAMETERS
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
    SUBROUTINE DO_LAVA_FLOWS(NNNN)
        !      *******************************************************************
        !       this subroutine models lava flows originating from defined vents.  the activity of the
        !         vents can be stochastic.  lava flows downhill but the probability of flow into a neighboring
        !         cell depends upon the topographic gradient, the lava flow thickness, and, inversely, to the
        !         length of time since the possible source cell has been actively flowing.
        !       the model is outlined briefly in howard, a. d., abstract 1112, lunar and planetary science
        !         conference xxx, 1999.
        !        the model is based several
        !      assumptions: 1) individual flows are thin compared to
        !      their length and breadth, and have a spatially and
        !      temporally constant thickness; 2) emplacement of
        !      individual flows occurs over an extended period of
        !      time; 3) most new emplacement of lava occurs at the
        !      leading edge of the flow - cooler, earlier deposited
        !      parts of the flow cool they are less likely to be a source
        !      for new lava; 4) lava is erupted from relatively
        !      few, localized, long-lived source vents;
        !      5) flow events cease either due to
        !      interruption at the feeding source or by ponding at its
        !      terminus; 6) erosion by lava flows is negligible.
        !         one to several cells serve as source vents. lava emplacement
        !      occurs during each iteration at the edges of the flow.
        !      a new cell at the boundary of the flow is chosen during
        !      each iteration to receive lava deposition of thickness
        !      h. each cell, i, of the active flow is examined as
        !      a potential source for the flow extension. at each
        !      source cell, i, the surrounding eight cells, j, are examined
        !      as potential receiving site. the probability
        !      function, pij, is given by:
        !      pij=exp[-c*ai]*(1-exp[-z*sij]),
        !      where c and z are parameters, ai is the elapsed time
        !      since cell i was inundated with lava, measured in iterations,
        !      and sij is the gradient from cell i to j. thus
        !      the probability of a cell on the flow serving as a source
        !      diminishes with time since its emplacement. the
        !      probability is zero if the neighboring cell is part of the
        !      current flow or if the gradient to the cell is uphill. a
        !      particular source cell and neighboring deposition cell
        !      is chosen randomly in proportion to the probabilities
        !      pij. the selected cell is flooded with lava to depth h
        !      and becomes a new potential source cell. this process
        !      of cell selection and flooding occurs until there is either
        !      no positive probability pij for any active part of
        !      the existing flow, or a flow interruption event occurs.
        !      a finite probability, px, exists that the flow interruption
        !      will occur after each iteration. if an interruption
        !      occurs, a new flow is started from the vent.!
        !      *******************************************************************
        USE ERODE_GLOBALS
        USE LAVA_GLOBALS
        IMPLICIT NONE
        REAL (4) :: RANDTEMP
        INTEGER :: INX,JNX,LLL,KKK,I,J,K,ITERFLOWSTART,NNNN
        REAL(4) :: MINEL,PROBCOMP
        REAL (4) :: NLACTIVE1,NLACTIVE2
        REAL(4) :: TVOLUME,DTEMP, RRAND
        EXTERNAL RRAND
        !      *******************************************************************
        !       istart and jstart are the i and j coordinates of the current lava
        !         flow vent source
        !  MODIFIES:  (ALL LAVA_GLOBALS), IS_SEDIMENT_COVERED, SEDIMENT_BASE, IS_ROCK_SURFACE
        !              ELEVATION, REGOLITH
        !  CALLS: FIND_ACTIVE_LAVA_SITES, FIND_NEXT_LAVA_SITE, FIND_LAVA_START_PLACE
        !      *******************************************************************
        ISTART=I_LAVA(NNNN)
        JSTART=J_LAVA(NNNN)
        TVOLUMESUM=0.0
        TVOLUMESQ=0.0
        TLACTIVESUM=0.0
        TLACTIVESQ=0.0
        TACTSAMP=0.0
        TVOLUMESAMP=0.0
        WRITE (OUTHIST,3478) ISTART,JSTART
        3478  FORMAT (' VENT AT I=',I3,' J=',I3)
        !      *******************************************************************
        !       Sets up initial values of simulation parameters:
        !         lava(,) is true if lava has ever been deposited on that site
        !         active_lava_flow(,) is true if a given cell is a potential source for lava
        !           to a neighboring point
        !         lava_source_direction(,) is the direction pointer to the cell that contributed
        !           lava to the current cell
        !         eruption_age(,) measures the time since new lava was last deposited in a
        !           cell from a neighboring point
        !           -- it doesn't include the possibility that a cell may have
        !           flowing lava below a cooled surface, presumably as a lava tunnel
        !      *******************************************************************
        IS_LAVA_COVERED=.FALSE.
        ACTIVE_LAVA_FLOW=.FALSE.
        LAVA_SOURCE_DIRECTION=0
        ERUPTION_AGE=0
        !      *******************************************************************
        !        Information about active_lava_flow source cells is stored in several vectors
        !          iloc() and jloc() give the i and j coordinates of each active_lava_flow
        !            cell
        !          lava_elapsed_time() is the number of iterations since an active_lava_flow flow entered
        !            the cell
        !        This sets up initial values
        !      *******************************************************************
        LAVA_THICKNESS=0.0
        LAVA_ELAPSED_TIME=0.0
        !      *******************************************************************
        !        number_of_active_lava_cells is the number of active_lava_flow cells -- not really the total
        !          flow length -- initially zero of course
        !      *******************************************************************
        NUMBER_OF_ACTIVE_LAVA_CELLS=0
        !      *******************************************************************
        !       find_lava_start_place determines the initial flow direction from the source
        !         cell at location istart, jstart.  Location inx, jnx is the location
        !         of the neighboring cell to receive the lava inundation.  lll is the
        !         flow direction to the cell receiving new lava
        !      *******************************************************************
        CALL FIND_LAVA_START_PLACE(ISTART,JSTART,LLL,INX,JNX)
        WRITE(OUTHIST,510)LLL,INX,JNX
        510      FORMAT('FIND_LAVA_START_PLACE, LLL=',I5,' INX=',I5,' JNX=',I5)
        !      *******************************************************************
        !        The following code assures that the source cell elevation is always
        !          higher than the cell receiving the lava
        !      *******************************************************************
        IF (LLL <= 5) THEN
            MINEL=ELEVATION(INX,JNX)+MINIMUM_LAVA_FLOW_SLOPE*CELL_SIZE
        ELSE
            MINEL=ELEVATION(INX,JNX)+MINIMUM_LAVA_FLOW_SLOPE*CELL_SIZE/0.7071068
        ENDIF
        MINEL=MINEL+LAVA_FLOW_THICKNESS
        IF (ELEVATION(ISTART,JSTART) < MINEL) THEN
            ELEVATION(ISTART,JSTART)=MINEL
            DTEMP=MINEL-ELEVATION(ISTART,JSTART)
            CUMULATIVE_LAVA_CHANGE(ISTART,JSTART)=CUMULATIVE_LAVA_CHANGE(ISTART,JSTART) &
            +DTEMP
        ENDIF
        !      *******************************************************************
        !       Increment the elevation of the cell receiving the lava by lava_flow_thickness
        !      *******************************************************************
        ELEVATION(INX,JNX)=ELEVATION(INX,JNX)+LAVA_FLOW_THICKNESS
        CUMULATIVE_LAVA_CHANGE(INX,JNX)=CUMULATIVE_LAVA_CHANGE(INX,JNX)  &
        + LAVA_FLOW_THICKNESS
        LAVAWORK=LAVAWORK+LAVA_FLOW_THICKNESS*    &
        SQRT(REAL((INX-ISTART)**2+(JNX-JSTART)**2))
        LAVAGRAV=LAVAGRAV+LAVA_FLOW_THICKNESS

        !      *******************************************************************
        !       Both the source cell and the cell receiving the lava are indicated
        !        to be lava covered and active_lava_flow.  The total elapsed time since inundation
        !        is zero. The number of active_lava_flow cells (number_of_active_lava_cells) is incremented by one.
        !        The i and j locations of the newly inundated cell are recorded, as
        !        is the direction back to the source cell.  The total time the inundated
        !        cell has been active_lava_flow is also set to zero.
        !      *******************************************************************
        IS_LAVA_COVERED(ISTART,JSTART)=.TRUE.
        ACTIVE_LAVA_FLOW(ISTART,JSTART)=.TRUE.
        !      *******************************************************************
        !         interaction with erosion program
        !      *******************************************************************
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            IS_SEDIMENT_COVERED(ISTART,JSTART)=.FALSE.
            SEDIMENT_BASE(INX,JNX)=ELEVATION(INX,JNX)
            SEDIMENT_BASE(ISTART,JSTART)=ELEVATION(ISTART,JSTART)
            IS_ROCK_SURFACE(ISTART,JSTART)=.TRUE.
            REGOLITH(ISTART,JSTART)=-ROCK_WEATHERING_RATE
            IS_SEDIMENT_COVERED(INX,JNX)=.FALSE.
            IS_ROCK_SURFACE(INX,JNX)=.TRUE.
            REGOLITH(INX,JNX)=-ROCK_WEATHERING_RATE
        ENDIF
        ERUPTION_AGE(ISTART,JSTART)=0
        IS_LAVA_COVERED(INX,JNX)=.TRUE.
        ACTIVE_LAVA_FLOW(INX,JNX)=.TRUE.
        NUMBER_OF_ACTIVE_LAVA_CELLS=NUMBER_OF_ACTIVE_LAVA_CELLS+1
        ILOC(NUMBER_OF_ACTIVE_LAVA_CELLS)=INX
        JLOC(NUMBER_OF_ACTIVE_LAVA_CELLS)=JNX
        LAVA_THICKNESS(NUMBER_OF_ACTIVE_LAVA_CELLS)=LAVA_FLOW_THICKNESS
        LAVA_SOURCE_DIRECTION(INX,JNX)=IBACK(LLL)
        LAVA_ELAPSED_TIME(NUMBER_OF_ACTIVE_LAVA_CELLS)=0
        ITERFLOWSTART=1
        !      *******************************************************************
        !       This is the main loop  -- Inundation of a single cell with lava
        !         occurs during each iteration
        !      *******************************************************************
        L100: DO LAVAITER=1,SOURCE_CHANGE_INTERVAL
            IF (MOD(LAVAITER,100) == 0) THEN
                WRITE(OUTHIST,410) LAVAITER,NUMBER_OF_ACTIVE_LAVA_CELLS
                410           FORMAT(' I=',I8,' FL=',I5)
            ENDIF
            !      *******************************************************************
            !        We next determine if the existing flow will become inlactive (clogged)
            !         by drawing a random number and comparing it with the probability
            !         parameter, eruption_stop_probability
            !      *******************************************************************
            RANDTEMP=RRAND()
            PROBCOMP=RANDTEMP
            IF (PROBCOMP > ERUPTION_STOP_PROBABILITY) THEN
                !      *******************************************************************
                !        This is the case where there is no clogging of the active_lava_flow flow.
                !         We therefore look at all active_lava_flow cells to see which one will serve
                !         as the source for inundating a neighboring cell and which cell will
                !         be inundated.  inx and jnx are the i and j coordinates of the cell
                !         to be inundated.  kkk is the index of the source cell in the active_lava_flow
                !         cell vectors, and lll is the flow direction from the source cell to
                !         the inundated cell
                !      *******************************************************************
                CALL FIND_NEXT_LAVA_SITE(LLL,KKK,INX,JNX)
                !      *******************************************************************
                !        If kkk<=0 then there is no cell selected -- in this case we will
                !          start a new flow from the source.  this also occurs if we get
                !          too many active_lava_flow cells
                !      *******************************************************************
                IF ((KKK > 0).AND.(NUMBER_OF_ACTIVE_LAVA_CELLS < RMMX))THEN
                    !      *******************************************************************
                    !        This is the case where we do have lava inundation of a cell.
                    !        values of several variables are set as indicated above, and the
                    !        elevations are incremented to reflect the lava deaccess
                    !      *******************************************************************
                    EOLD=ELEVATION(INX,JNX)
                    ELEVATION(INX,JNX)=AMIN1((ELEVATION(INX,JNX)+LAVA_FLOW_THICKNESS) &
                    ,ELEVATION(ILOC(KKK),JLOC(KKK)))
                    DTEMP=ELEVATION(INX,JNX)-EOLD
                    CUMULATIVE_LAVA_CHANGE(INX,JNX)=CUMULATIVE_LAVA_CHANGE(INX,JNX) &
                    +DTEMP
                    LAVAWORK=LAVAWORK+(ELEVATION(INX,JNX)-EOLD)*    &
                    SQRT(REAL((INX-ISTART)**2+(JNX-JSTART)**2))
                    IS_LAVA_COVERED(INX,JNX)=.TRUE.
                    !      *******************************************************************
                    !         Interaction with erosion program
                    !      *******************************************************************
                    IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                        IS_SEDIMENT_COVERED(INX,JNX)=.FALSE.
                        SEDIMENT_BASE(INX,JNX)=ELEVATION(INX,JNX)

                        IS_ROCK_SURFACE(INX,JNX)=.TRUE.
                        REGOLITH(INX,JNX)=-ROCK_WEATHERING_RATE
                    ENDIF
                    !      *******************************************************************
                    NUMBER_OF_ACTIVE_LAVA_CELLS=NUMBER_OF_ACTIVE_LAVA_CELLS+1
                    ILOC(NUMBER_OF_ACTIVE_LAVA_CELLS)=INX
                    JLOC(NUMBER_OF_ACTIVE_LAVA_CELLS)=JNX
                    LAVA_THICKNESS(NUMBER_OF_ACTIVE_LAVA_CELLS)=ELEVATION(INX,JNX)-EOLD
                    LAVA_SOURCE_DIRECTION(INX,JNX)=IBACK(LLL)
                    ACTIVE_LAVA_FLOW(INX,JNX)=.TRUE.
                    ERUPTION_AGE(INX,JNX)=0
                    LAVA_ELAPSED_TIME(NUMBER_OF_ACTIVE_LAVA_CELLS)=0
                    !      *******************************************************************
                    !        The elapsed time for all cells except for the newly inundated
                    !         cell is incremented
                    !      *******************************************************************
                    DO  K=1,NUMBER_OF_ACTIVE_LAVA_CELLS-1
                        LAVA_ELAPSED_TIME(K)=LAVA_ELAPSED_TIME(K)+1
                    ENDDO
                ELSE
                    !      *******************************************************************
                    !        This is the case where the existing flow is abandoned and a new
                    !          flow is initiated.  the pressure code is not used.
                    !      *******************************************************************
                    WRITE(OUTHIST,1420) KKK
                    1420  FORMAT(' NEW FLOW STARTED, KKK=',I5)
                    !      *******************************************************************
                    !        All active_lava_flow cells are dropped, a new direction of flow from the
                    !         source is selected, and the first lava deaccess for the new
                    !         flow occurs -- all parallel to code for the start of the program
                    !      *******************************************************************
                    CALL FIND_LAVA_START_PLACE(ISTART,JSTART,LLL,INX,JNX)
                    IF (LLL <= 5) THEN
                        MINEL=ELEVATION(INX,JNX)+MINIMUM_LAVA_FLOW_SLOPE*CELL_SIZE
                    ELSE
                        MINEL=ELEVATION(INX,JNX)+MINIMUM_LAVA_FLOW_SLOPE*CELL_SIZE/0.7071068
                    ENDIF
                    ACTIVE_LAVA_FLOW=.FALSE.
                    ACTIVE_LAVA_FLOW(ISTART,JSTART)=.TRUE.
                    ERUPTION_AGE(ISTART,JSTART)=0
                    MINEL=MINEL+LAVA_FLOW_THICKNESS
                    IF (ELEVATION(ISTART,JSTART) < MINEL) THEN
                        ELEVATION(ISTART,JSTART)=MINEL
                        DTEMP=MINEL-ELEVATION(ISTART,JSTART)
                        CUMULATIVE_LAVA_CHANGE(ISTART,JSTART)=CUMULATIVE_LAVA_CHANGE(ISTART,JSTART) &
                        +DTEMP
                    ENDIF
                    ELEVATION(INX,JNX)=ELEVATION(INX,JNX)+LAVA_FLOW_THICKNESS
                    CUMULATIVE_LAVA_CHANGE(INX,JNX)=CUMULATIVE_LAVA_CHANGE(INX,JNX) &
                    + LAVA_FLOW_THICKNESS
                    LAVAWORK=LAVAWORK+LAVA_FLOW_THICKNESS* &
                    SQRT(REAL((INX-ISTART)**2+(JNX-JSTART)**2))
                    IS_LAVA_COVERED(INX,JNX)=.TRUE.
                    !      *******************************************************************
                    !         Interaction with erosion program
                    !      *******************************************************************
                    IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                        IS_SEDIMENT_COVERED(INX,JNX)=.FALSE.
                        IS_ROCK_SURFACE(INX,JNX)=.TRUE.
                        REGOLITH(INX,JNX)=-ROCK_WEATHERING_RATE
                        SEDIMENT_BASE(INX,JNX)=ELEVATION(INX,JNX)
                        SEDIMENT_BASE(ISTART,JSTART)=ELEVATION(ISTART,JSTART)
                        IS_SEDIMENT_COVERED(ISTART,JSTART)=.FALSE.
                        IS_ROCK_SURFACE(ISTART,JSTART)=.TRUE.
                        REGOLITH(ISTART,JSTART)=-ROCK_WEATHERING_RATE
                    ENDIF
                    !      *******************************************************************
                    NUMBER_OF_ACTIVE_LAVA_CELLS=1
                    ILOC(NUMBER_OF_ACTIVE_LAVA_CELLS)=INX
                    JLOC(NUMBER_OF_ACTIVE_LAVA_CELLS)=JNX
                    LAVA_THICKNESS(NUMBER_OF_ACTIVE_LAVA_CELLS)=LAVA_FLOW_THICKNESS
                    ACTIVE_LAVA_FLOW(INX,JNX)=.TRUE.
                    ERUPTION_AGE(INX,JNX)=0
                    LAVA_SOURCE_DIRECTION(INX,JNX)=IBACK(LLL)
                    LAVA_ELAPSED_TIME(NUMBER_OF_ACTIVE_LAVA_CELLS)=0
                    TVOLUME=DFLOAT(LAVAITER-ITERFLOWSTART)
                    TVOLUMESUM=TVOLUMESUM+TVOLUME
                    TVOLUMESQ=TVOLUMESQ+TVOLUME*TVOLUME
                    TVOLUMESAMP=TVOLUMESAMP+1.0
                    ITERFLOWSTART=LAVAITER+1
                ENDIF
                !      *******************************************************************
                !        This is the case where flow gets blocked and a new flow is started:
                !        All active_lava_flow cells are dropped, a new direction of flow from the
                !         source is selected, and the first lava deaccess for the new
                !         flow occurs -- all parallel to code for the start of the program
                !      *******************************************************************
            ELSE
                WRITE(OUTHIST,450)
                450            FORMAT('NEW_LAVA_FLOW')
                CALL FIND_LAVA_START_PLACE(ISTART,JSTART,LLL,INX,JNX)
                IF (LLL <= 5) THEN
                    MINEL=ELEVATION(INX,JNX)+MINIMUM_LAVA_FLOW_SLOPE*CELL_SIZE
                ELSE
                    MINEL=ELEVATION(INX,JNX)+MINIMUM_LAVA_FLOW_SLOPE*CELL_SIZE/0.7071068
                ENDIF
                ACTIVE_LAVA_FLOW=.FALSE.
                MINEL=MINEL+LAVA_FLOW_THICKNESS
                IF (ELEVATION(ISTART,JSTART) < MINEL) THEN
                    ELEVATION(ISTART,JSTART)=MINEL
                    DTEMP=MINEL-ELEVATION(ISTART,JSTART)
                    CUMULATIVE_LAVA_CHANGE(ISTART,JSTART)=CUMULATIVE_LAVA_CHANGE(ISTART,JSTART)  &
                    +DTEMP
                ENDIF
                ELEVATION(INX,JNX)=ELEVATION(INX,JNX)+LAVA_FLOW_THICKNESS
                CUMULATIVE_LAVA_CHANGE(INX,JNX)=CUMULATIVE_LAVA_CHANGE(INX,JNX) &
                + LAVA_FLOW_THICKNESS
                LAVAWORK=LAVAWORK+LAVA_FLOW_THICKNESS*   &
                SQRT(REAL((INX-ISTART)**2+(JNX-JSTART)**2))
                IS_LAVA_COVERED(INX,JNX)=.TRUE.
                !      *******************************************************************
                !         Interaction with erosion program
                !      *******************************************************************
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                    IS_SEDIMENT_COVERED(INX,JNX)=.FALSE.
                    SEDIMENT_BASE(INX,JNX)=ELEVATION(INX,JNX)
                    IS_ROCK_SURFACE(INX,JNX)=.TRUE.
                    REGOLITH(INX,JNX)=-ROCK_WEATHERING_RATE
                    IS_SEDIMENT_COVERED(ISTART,JSTART)=.FALSE.
                    SEDIMENT_BASE(ISTART,JSTART)=ELEVATION(ISTART,JSTART)
                    IS_ROCK_SURFACE(ISTART,JSTART)=.TRUE.
                    REGOLITH(ISTART,JSTART)=-ROCK_WEATHERING_RATE
                ENDIF
                !      *******************************************************************
                ACTIVE_LAVA_FLOW(ISTART,JSTART)=.TRUE.
                ACTIVE_LAVA_FLOW(INX,JNX)=.TRUE.
                ERUPTION_AGE(ISTART,JSTART)=0
                ERUPTION_AGE(INX,JNX)=0
                NUMBER_OF_ACTIVE_LAVA_CELLS=1
                ILOC(NUMBER_OF_ACTIVE_LAVA_CELLS)=INX
                JLOC(NUMBER_OF_ACTIVE_LAVA_CELLS)=JNX
                LAVA_THICKNESS(NUMBER_OF_ACTIVE_LAVA_CELLS)=LAVA_FLOW_THICKNESS
                LAVA_SOURCE_DIRECTION(INX,JNX)=IBACK(LLL)
                LAVA_ELAPSED_TIME(NUMBER_OF_ACTIVE_LAVA_CELLS)=0
                TVOLUME=DFLOAT(LAVAITER-ITERFLOWSTART)
                TVOLUMESUM=TVOLUMESUM+TVOLUME
                TVOLUMESQ=TVOLUMESQ+TVOLUME*TVOLUME
                TVOLUMESAMP=TVOLUMESAMP+1.0
                ITERFLOWSTART=LAVAITER+1
            ENDIF
            !      *******************************************************************
            !       At intervals of new_segment_interval iterations the list of active_lava_flow cells
            !         is reevaluated.  also, the length of time since last deaccess
            !         on each cell is updated.
            !         Also statistics are collected about the number of active_lava_flow cells
            !         before and after reevaluation
            !      *******************************************************************
            IF (MOD(LAVAITER,NEW_SEGMENT_INTERVAL) == 0) THEN
                NLACTIVE1=0.0
                DO I=1,MX
                    DO J=1,MY
                        IF (ACTIVE_LAVA_FLOW(I,J)) THEN
                            NLACTIVE1=NLACTIVE1+1.0
                        ENDIF
                    ENDDO
                ENDDO
                CALL FIND_ACTIVE_LAVA_SITES()
                NLACTIVE2=0.0
                DO I=1,MX
                    DO J=1,MY
                        IF (ACTIVE_LAVA_FLOW(I,J)) THEN
                            NLACTIVE2=NLACTIVE2+1.0
                        ENDIF
                    ENDDO
                ENDDO
                NLACTIVE1=(NLACTIVE1+NLACTIVE2)/2.0
                TLACTIVESUM=TLACTIVESUM+NLACTIVE1
                TLACTIVESQ=TLACTIVESQ+NLACTIVE1*NLACTIVE1
                TACTSAMP=TACTSAMP+1.0
                DO  J=1,MY
                    DO  I=1,MX
                        ERUPTION_AGE(I,J)=ERUPTION_AGE(I,J)+1
                    ENDDO
                ENDDO
            ENDIF
        ENDDO L100
        998   IF (TACTSAMP > 1.0) THEN
            TLACTIVESQ=SQRT((TLACTIVESQ-TLACTIVESUM*TLACTIVESUM/TACTSAMP)  &
            /(TACTSAMP-1.0))
            TLACTIVESUM=TLACTIVESUM/TACTSAMP
            WRITE (OUTHIST,444) TLACTIVESUM,TLACTIVESQ,TACTSAMP
            444   FORMAT(' MEAN NO. OF LACTIVE POINTS=',G12.5,' S.D.=',G12.5,  &
            ' N=',G12.5)
        ENDIF
        IF (TVOLUMESAMP > 1.0) THEN
            TVOLUMESQ=SQRT((TVOLUMESQ-TVOLUMESUM*TVOLUMESUM/TVOLUMESAMP) &
            /(TVOLUMESAMP-1.0))
            TVOLUMESUM=TVOLUMESUM/TVOLUMESAMP
            WRITE (OUTHIST,443) TVOLUMESUM,TVOLUMESQ,TVOLUMESAMP
            443   FORMAT(' MEAN AREA OF FLOWS (IN CELLS)=',G12.5,' S.D.=',G12.5, &
            ' N=',G12.5)
        ENDIF
        RETURN
    END ! SUBROUTINE DO_LAVA_FLOWS(NNNN)
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE FIND_ACTIVE_LAVA_SITES()
        !      *******************************************************************
        !        This subroutine adjusts the list of active_lava_flow cells, dropping those
        !         which have a low probability of fathering a new flow, but keeping
        !         those which are feeding active_lava_flow cells
        !      *******************************************************************
        USE ERODE_GLOBALS
        USE LAVA_GLOBALS
        IMPLICIT NONE
        INTEGER :: I_OLD,J_OLD,IIO,JJO,I_NEW,J_NEW,KGET,KPUT,NEWFLOW,LLL,K,L
        REAL(4) :: PSUM
        !      *******************************************************************
        !        Initially reset all cells to inlactive
        !      *******************************************************************
        ACTIVE_LAVA_FLOW=.FALSE.
        !      *******************************************************************
        !        However, the flow vent remains active_lava_flow
        !      *******************************************************************
        ACTIVE_LAVA_FLOW(ISTART,JSTART)=.TRUE.
        !      *******************************************************************
        !        Here we cycle through all previously active_lava_flow cells and look for ones
        !         that still have a reasonable probability of fathering a new flow
        !      *******************************************************************
        L110: DO K=1,NUMBER_OF_ACTIVE_LAVA_CELLS
            I_OLD=ILOC(K)
            J_OLD=JLOC(K)
            PSUM=0.0
            !      *******************************************************************
            !        For each formerly active_lava_flow cell, total the probabilities of fathering
            !         a new flow over all flow directions
            !      *******************************************************************
            DO L=2,9
                PSUM=PSUM+LAVA_FLOW_PROBABILITY(L,K)
            ENDDO
            !      *******************************************************************
            !        If that probability exceeds the parameter no_flow_probability, then that cell
            !         is retained as an active_lava_flow cell
            !      *******************************************************************
            IF (PSUM > NO_FLOW_PROBABILITY) THEN
                ACTIVE_LAVA_FLOW(I_OLD,J_OLD)=.TRUE.
                IIO=I_OLD
                JJO=J_OLD
                !      *******************************************************************
                !       In addition, we need to keep the train of source cells from the active_lava_flow
                !        cell up to the source also active_lava_flow.  This code goes "upstream" along
                !        the flow path, making each feeding cell also active_lava_flow
                !        -- we can stop if the cell is
                !        already active_lava_flow because we have already identified all upstream cells from
                !        that location from some other still-active_lava_flow cell
                !      *******************************************************************
                L160: DO
                    I_NEW=IIO + DOWNSTREAM(LAVA_SOURCE_DIRECTION(IIO,JJO),1)
                    J_NEW=JJO + DOWNSTREAM(LAVA_SOURCE_DIRECTION(IIO,JJO),2)
                    IF (IS_X_PERIODIC) THEN
                        IF (I_NEW < 1) I_NEW=MX
                        IF (I_NEW > MX) I_NEW=1
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (J_NEW < 1) J_NEW=MY
                        IF (J_NEW > MY) J_NEW=1
                    ENDIF
                    IF (ACTIVE_LAVA_FLOW(I_NEW,J_NEW)) THEN
                        EXIT L160
                    ELSE
                        ACTIVE_LAVA_FLOW(I_NEW,J_NEW)=.TRUE.
                    ENDIF
                    IIO=I_NEW
                    JJO=J_NEW
                ENDDO L160
            ENDIF
        ENDDO L110
        !      *******************************************************************
        !       This code condenses the vectors containing variables about active_lava_flow
        !         cells by overwriting the entries of inlactive cells -- also the
        !         total number of active_lava_flow cells (number_of_active_lava_cells) is updated.  Initially
        !         number_of_active_lava_cells is the previous value of the number of active_lava_flow cells.
        !         kput and newflow are dummy variables
        !      *******************************************************************
        KPUT=1
        NEWFLOW=1
        DO KGET=1,NUMBER_OF_ACTIVE_LAVA_CELLS
            IF (ACTIVE_LAVA_FLOW(ILOC(KGET),JLOC(KGET))) THEN
                ILOC(KPUT)=ILOC(KGET)
                JLOC(KPUT)=JLOC(KGET)
                LAVA_ELAPSED_TIME(KPUT)=LAVA_ELAPSED_TIME(KGET)
                LAVA_THICKNESS(KPUT)=LAVA_THICKNESS(KGET)
                DO  LLL=2,9
                    LAVA_FLOW_PROBABILITY(LLL,KPUT)=LAVA_FLOW_PROBABILITY(LLL,KGET)
                ENDDO
                KPUT=KPUT+1
                NEWFLOW=NEWFLOW+1
            ENDIF
        ENDDO
        NUMBER_OF_ACTIVE_LAVA_CELLS=NEWFLOW-1
        RETURN
    END ! SUBROUTINE FIND_ACTIVE_LAVA_SITES()
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE FIND_LAVA_START_PLACE(IST,JST,LLL,I_NEW,J_NEW)
        !      *******************************************************************
        !       This subroutine finds the direction and cell coordinates of the
        !        first cell adjacent to the source cell to receive lava flow.
        !       The procedure determines the initial flow direction from the source
        !         cell at location ist, jst.  Location i_new, j_new is the location
        !         of the neighboring cell to receive the lava inundation.  lll is the
        !         flow direction to the cell receiving new lava
        !      *******************************************************************
        USE ERODE_GLOBALS
        USE LAVA_GLOBALS
        IMPLICIT NONE
        REAL(4) :: SPROB(9),EMIN,LOCALGRAD,PROBSUM,PROBCOMPARE, RRAND
        REAL (4) :: RANDTEMP
        EXTERNAL RRAND
        INTEGER :: I_NEW,J_NEW,L,III,JJJ,LLL,IST,JST
        !      *******************************************************************
        !       This loops through each cell surrounding the source cell and
        !         calculates probabilites sprob() for each flow direction, l.
        !      *******************************************************************
        L110: DO L=2,9
            SPROB(L)=0.0
            J_NEW=JST+DOWNSTREAM(L,2)
            I_NEW=IST+DOWNSTREAM(L,1)
            IF (IS_X_PERIODIC) THEN
                IF (I_NEW < 1) I_NEW=MX
                IF (I_NEW > MX) I_NEW=1
            ENDIF
            IF (IS_Y_PERIODIC) THEN
                IF (J_NEW < 1) J_NEW=MY
                IF (J_NEW > MY) J_NEW=1
            ENDIF
            !      *******************************************************************
            !       This is the gradient to the adjacent cell
            !      *******************************************************************
            LOCALGRAD=(ELEVATION(IST,JST)-ELEVATION(I_NEW,J_NEW))/CELL_SIZE
            IF (LOCALGRAD > 0.0) THEN
                IF (L > 5) THEN
                    LOCALGRAD=LOCALGRAD*0.7071068
                ENDIF
                !      *******************************************************************
                !        The probability increases as the gradient increases, with 1 as the
                !          maximum probability
                !      *******************************************************************
                SPROB(L)=(1.0-EXP(-LAVA_GRADIENT_WEIGHT*LOCALGRAD))
            ELSE
                !      *******************************************************************
                !        But the probability is zero to any cell lying above the source
                !      *******************************************************************
                SPROB(L)=0.0
            ENDIF
        ENDDO L110
        !      *******************************************************************
        !        Here is where we select that cell that is to receive the flow,
        !         with the probability of selection being proportional to the
        !         probability, sprob(), of flow to that direction
        !      *******************************************************************
        PROBSUM=0.0
        DO  L=2,9
            PROBSUM=PROBSUM+SPROB(L)
        ENDDO
        IF (PROBSUM > 0.0) THEN
            RANDTEMP=RRAND()
            PROBCOMPARE=PROBSUM*RANDTEMP
            PROBSUM=0.0
            L130: DO L=2,9
                PROBSUM=PROBSUM+SPROB(L)
                IF (PROBSUM >= PROBCOMPARE) THEN
                    J_NEW=JST+DOWNSTREAM(L,2)
                    I_NEW=IST+DOWNSTREAM(L,1)
                    IF (IS_X_PERIODIC) THEN
                        IF (I_NEW < 1) I_NEW=MX
                        IF (I_NEW > MX) I_NEW=1
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (J_NEW < 1) J_NEW=MY
                        IF (J_NEW > MY) J_NEW=1
                    ENDIF
                    LLL=L
                    EXIT L130
                ENDIF
            ENDDO L130
        ELSE
            !      *******************************************************************
            !       If all surrounding cells are higher than the source cell, then we
            !        pick the lowest of the surrounding cells-- presumably the lava fills
            !        the depression until it spills out
            !      *******************************************************************
            EMIN=1.0E+25
            DO L=2,9
                J_NEW=JSTART+DOWNSTREAM(L,2)
                I_NEW=ISTART+DOWNSTREAM(L,1)
                IF (IS_X_PERIODIC) THEN
                    IF (I_NEW < 1) I_NEW=MX
                    IF (I_NEW > MX) I_NEW=1
                ENDIF
                IF (IS_Y_PERIODIC) THEN
                    IF (J_NEW < 1) J_NEW=MY
                    IF (J_NEW > MY) J_NEW=1
                ENDIF
                IF (ELEVATION(I_NEW,J_NEW) < EMIN) THEN
                    EMIN=ELEVATION(I_NEW,J_NEW)
                    III=I_NEW
                    JJJ=J_NEW
                    LLL=L
                ENDIF
            ENDDO
            I_NEW=III
            J_NEW=JJJ
        ENDIF
        RETURN
    END ! SUBROUTINE FIND_LAVA_START_PLACE(IST,JST,LLL,I_NEW,J_NEW)

    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE FIND_NEXT_LAVA_SITE(LLL,KKK,INX,JNX)
        !      *******************************************************************
        !       This subroutine finds the next cell to receive lava flows.
        !         We therefore look at all active_lava_flow cells to see which one will serve
        !         as the source for inundating a neighboring cell and which cell will
        !         be inundated.  inx and jnx are the i and j coordinates of the cell
        !         to be inundated.  kkk is the index of the source cell in the active_lava_flow
        !         cell vectors, and lll is the flow direction from the source cell to
        !         the inundated cell
        !      *******************************************************************
        USE ERODE_GLOBALS
        USE LAVA_GLOBALS
        IMPLICIT NONE
        INTEGER :: KKK,K,L,I_OLD,J_OLD,I_NEW,J_NEW,INX,JNX,LLL
        REAL(4) :: LOCALGRAD,PROBSUM,PROBCOMPARE, RRAND
        REAL (4) :: RANDTEMP
        EXTERNAL RRAND
        !      *******************************************************************
        !        We cycle through all active_lava_flow cells
        !      *******************************************************************
        L100: DO K=1,NUMBER_OF_ACTIVE_LAVA_CELLS
            I_OLD=ILOC(K)
            J_OLD=JLOC(K)
            L110: DO L=2,9
                LAVA_FLOW_PROBABILITY(L,K)=-1.0
                J_NEW=J_OLD+DOWNSTREAM(L,2)
                I_NEW=I_OLD+DOWNSTREAM(L,1)
                IF (IS_X_PERIODIC) THEN
                    IF (I_NEW < 1) I_NEW=MX
                    IF (I_NEW > MX) I_NEW=1
                ELSE
                    IF (I_NEW < 1) CYCLE
                    IF (I_NEW > MX) CYCLE
                ENDIF
                IF (IS_Y_PERIODIC) THEN
                    IF (J_NEW < 1) J_NEW=MY
                    IF (J_NEW > MY) J_NEW=1
                ELSE
                    IF (J_NEW < 1) CYCLE
                    IF (J_NEW > MY) CYCLE
                ENDIF
                !      *******************************************************************
                !        We don't flood cells that are already active_lava_flow
                !      *******************************************************************
                111   FORMAT(' I_NEW=',I6,' J_NEW=',I6,' L=',I5,' K=',I5)
                IF (ACTIVE_LAVA_FLOW(I_NEW,J_NEW)) LAVA_FLOW_PROBABILITY(L,K)=0.0
                !      *******************************************************************
                !       The probability, lava_flow_probability(,) of flooding a neighboring cell in the l direction
                !        from the k'th active_lava_flow cell increases as the gradient to the adjacent cell
                !        increases but decreases as the length of time since the cell was
                !        last flooded with lava increases
                !      *******************************************************************
                IF (LAVA_FLOW_PROBABILITY(L,K) /= 0.0) THEN
                    LOCALGRAD=(ELEVATION(I_OLD,J_OLD)-ELEVATION(I_NEW,J_NEW))/CELL_SIZE
                    IF ((LOCALGRAD > 0.0).AND.(LAVA_THICKNESS(K) >= MINIMUM_FLOW_THICKNESS)) &
                    THEN
                        IF (L > 5) THEN
                            LOCALGRAD=LOCALGRAD*0.7071068
                        ENDIF
                        LAVA_FLOW_PROBABILITY(L,K)=(1.0-EXP(-LAVA_GRADIENT_WEIGHT*LOCALGRAD))* &
                        EXP(-LAVA_DURATION_WEIGHT*LAVA_ELAPSED_TIME(K))
                    ELSE
                        LAVA_FLOW_PROBABILITY(L,K)=0.0
                    ENDIF
                ENDIF
            ENDDO L110
        ENDDO L100
        !      *******************************************************************
        !       Here we determine the summ of all probabilites for all active_lava_flow cells
        !         and possible flood directions
        !      *******************************************************************
        PROBSUM=0.0
        DO  K=1,NUMBER_OF_ACTIVE_LAVA_CELLS
            DO  L=2,9
                PROBSUM=PROBSUM+LAVA_FLOW_PROBABILITY(L,K)
            ENDDO
        ENDDO
        !      *******************************************************************
        !       We need to check whether there is at least one suitable cell to
        !        receive lava deaccess
        !      *******************************************************************
        IF (PROBSUM > 0.0) THEN
            !      *******************************************************************
            !       If there is one or more suitable active_lava_flow cells and flow directions,
            !        we select one with probability proportional to the lava_flow_probability(,) variable
            !      *******************************************************************
            RANDTEMP=RRAND()
            PROBCOMPARE=PROBSUM*RANDTEMP
            PROBSUM=0.0
            L130: DO  K=1,NUMBER_OF_ACTIVE_LAVA_CELLS
                M130: DO  L=2,9
                    PROBSUM=PROBSUM+LAVA_FLOW_PROBABILITY(L,K)
                    IF (PROBSUM >= PROBCOMPARE) THEN
                        KKK=K
                        LLL=L
                        I_OLD=ILOC(KKK)
                        J_OLD=JLOC(KKK)
                        JNX=J_OLD+DOWNSTREAM(L,2)
                        INX=I_OLD+DOWNSTREAM(L,1)
                        IF (IS_X_PERIODIC) THEN
                            IF (INX < 1) INX=MX
                            IF (INX > MX) INX=1
                        ENDIF
                        IF (IS_Y_PERIODIC) THEN
                            IF (JNX < 1) JNX=MY
                            IF (JNX > MY) JNX=1
                        ENDIF
                        EXIT L130
                    ENDIF
                ENDDO M130
            ENDDO L130
            !      *******************************************************************
            !        If there is no suitable cell to receive the lava, then we indicate
            !        that by making the source cell, kkk, negative
            !      *******************************************************************
        ELSE
            KKK=-1
            INX=-1
            JNX=-1
        ENDIF
        RETURN
    END ! SUBROUTINE FIND_NEXT_LAVA_SITE(LLL,KKK,INX,JNX)
