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

    !      ********************************************
    !          Include file for basin erosion program
    !              defines global variables
    !      ********************************************
    MODULE ERODE_GLOBALS
        IMPLICIT NONE
        SAVE
        !         immx and jmmx are the maximum x and y size of the simulation
        !         lmmx should equal immx*jmmx
        INTEGER(4) :: IMMX,JMMX,LMMX,RMMX
        REAL(4) ::  CUT_RATIO,BIAS_PARAMETER
        REAL(4) ::  SECONDS_PER_YEAR,SEDSETTLE,WATER_DENSITY,SKLAR_FACTOR
        PARAMETER (CUT_RATIO=0.001)
        PARAMETER (BIAS_PARAMETER=0.2)
        PARAMETER (SECONDS_PER_YEAR=31.536E+6,SEDSETTLE=0.66667)
        INTEGER(4) :: INDATA,INRESIST,OUTDATA,OUTHIST,OUTCHAN,OUTSUBMERGE,OUTCRATERS
        INTEGER(4) :: OUTHIST_GLACIAL_FLUX_BINGHAM                      ! added by (orkan) 03/30/2014
        INTEGER(4) :: OUTICE                                            ! added by (orkan) 04/10/2014,
        INTEGER(4) :: OUTHIST_MATERIAL_VOLUMES       																									 																							  																										  
        INTEGER(4) :: DOWNSTREAM(9,2),OUTRECORD,OUTSUMMARY,WRITE_CHANGE_INTERVAL
        INTEGER(4) :: MX,MY,MZ,OUTSOURCE,OUTREPORT,OUTDEFORM,OUTDISCHARGES,OUT_BINARY_DISCHARGES
        INTEGER(4) :: IABORTMAX,IABORTMIN,JABORTMAX,JABORTMIN
        INTEGER(4) :: MP,EVENT_INDEX,NUMBER_OF_EVENTS
        INTEGER(4) :: IDUMMY,ICLAST,JCLAST,NUMBER_ABOVE_THRESHOLD,NUMBER_BELOW_THRESHOLD
        INTEGER(4) :: MAXIMUM_ITERATION,ELEVATION_PRINT_INTERVAL,OUTPUT_PRINT_INTERVAL,RECALCULATE_GRADIENT_INTERVAL
        INTEGER(4) :: CRITICAL_GRADIENT_USE,ITERATION,WRITETYPE(13)
        INTEGER(4) :: ONEONLY,STARTING_ITERATION,TOTAL_ITERATIONS
        INTEGER(4) :: LAYERUSE,VARIABLE_ROCK_RESISTANCE_USE,CHANNELVARUSE
        INTEGER(4) :: NCRITS,NMAX,LMAX,NEWERODE,KWRITE,RANDDISCHUSE
        INTEGER(4) :: TOTAL_BASINS,EVENT_TYPE
        INTEGER(4) :: IWINLOW,IWINHIGH,JWINLOW,JWINHIGH,OUTSEDDEBUG , INTEGER_DIVERGENCE_RANGE
        INTEGER(4) :: ICENT,JCENT,INRATES,PARAMETER_CHANGE_INDEX,NUMBER_OF_PARAMETER_CHANGES
        INTEGER(4) :: OUTIMGDAT,IMAGE_OUTPUT_INTERVAL,OUTROCK,OUTHIGH,OUTWORK
        INTEGER(4) :: DIVERGEINTERVAL,MYY,OUTREGOLITH,OUTLAKE,USE_DIVERGENCE_DEPENDENT_RUNOFF
        INTEGER(4) :: ILOWEST,JLOWEST,INUMLEVELS,IDET,JDET
        INTEGER(2) ::  IFILE1,IFILE2,IFILE3,XFILE1,XFILE2,XFILE3,IFILE4,XFILE4
        INTEGER(4) :: OUTRESIST,OUTGRAD,OUTALLUV,OUTEROSION_DEPTH_INDEX,OUTCONT
        INTEGER(4) :: OUTELEV,IOTEMP1,IOTEMP2,OUTBASE,OUTCRATER,OUTAVALANCHE
        INTEGER(4) :: IDEBUG(500),JDEBUG(500),IJDEBUG,NCALCEVAP
        INTEGER(4) :: I_INFLUENT_RIVER(10),J_INFLUENT_RIVER(10),NUMBER_OF_INFLUENT_RIVERS,BINGFLUX
        INTEGER(4) ::  ISEED,NCHANGEEVAP,PLOT_DIAGNOSTICS, MAX_OUT_DISTANCE
        INTEGER(4) :: ICASE(20), TOTAL_RANDOM_CRATERS,OUT_IMAGE_FILE,WRITE_INITIAL_EXPOSURE
        INTEGER(4) :: RECALCULATE_DISCHARGE_INTERVAL,READREGOLITH
        REAL(4) ::  DISCHARGE_CONSTANT,DISCHARGE_EXPONENT, TRANSPORTFACTOR,MANNING,CONVERT_TO_METERS 
        REAL(4) ::  SEDIMENT_POROSITY,SEDIMENT_SPECIFIC_GRAVITY,GRAVITY,INPUT_CELL_SIZE,VERTICAL_SCALING_FACTOR
        REAL(4) ::  TIMES_FOR_PARAMETER_CHANGES(50),EROSION_RATE_VALUES(50),ELOWEST, HIGH_DIVERGENCE_RUNOFF,LOW_DIVERGENCE_RUNOFF
        REAL(4) ::  DEFAULT_CHANNEL_TIMESTEP,FLOW_FRACTION,GRAIN_SIZE,REGOLITH_ERODIBILITY_FACTOR
        REAL(4) ::  SQRTOFTWO,XSEED,ISEED1,SLOPE_FAILURE_DIFFUSIVITY,GRADMAX,FAILMAX
        REAL(4) ::  RATECAT(102),GRADCAT(102),DELRATE(102),CELL_AREA
        REAL(4) ::  BOUNDARY_LOWERING_RATE,BEDROCK_DISCHARGE_EXPONENT,BEDROCK_GRADIENT_EXPONENT
        REAL(4) :: SLOPE_DIFFUSIVITY,BEDROCK_ERODIBILITY
        REAL(4) ::  DETACHMENT_CRITICAL_SHEAR,CRITICAL_SLOPE_GRADIENT,SLOPE_GRADIENT_EXPONENT,TIME_INCREMENT,PRESENT_TIME
        REAL(4) ::  CHANNEL_TIMESTEP_SCALING,MAXIMUM_ELEVATION_CHANGE,GRADAVERAGE,MAXGRADIENT,TMULT
        REAL(4) ::  NEWBASE,MAXIMUM_DIFFUSIVITY_INCREASE,CRITICAL_SOURCE_DIVERGENCE,CRITICAL_GRADIENT_TERM
        REAL(4) ::  MAXIMUM_CHANNEL_TIMESTEP,NUMBIDCHANGE,SUMGRADCHANGE,ABSGRADCHANGE
        REAL(4) ::  CHANGECOUNT,DELTA_FORESET_GRADIENT,DISCHARGE_COEFFICIENT,BISTABLE_RUNOFF_FACTOR
        REAL(4) ::  RESISTANCE_VARIABILITY,SIMULATION_PARAMETERS(50,20)
        REAL(4) ::  SURFACE_LAYER_RESISTANCE,LAYERMAX,LAYERMIN,SEDIMENT_DISCHARGE_FACTOR
        REAL(4) ::  CROSS_WEIGHTING,DIAGONAL_WEIGHTING,MAXCRIT,FLEN(9),ONE_SIXTH
        REAL(4) ::  AREAFACTOR,HYDRAULIC_CONDUCTIVITY,EVAPORATION_MEAN,EVAPORATION_STANDARD_DEVIATION,EVAPORATION_SCALE
        REAL(4) ::  CHANNEL_WIDTH_CONSTANT,CHANNEL_WIDTH_EXPONENT,DISCHARGE_COEFF_VARIATION,MINUMUM_TIME_INCREMENT
        REAL(4) ::  SEDIMENT_1_EXPONENT,SEDIMENT_2_EXPONENT,SEDIMENT_GRADIENT_EXPONENT,SEDIMENT_CONSTANT
        REAL(4) ::  TRANSPORT_CRITICAL_DIM_SHEAR,SEDIMENT_TRANSPORT_EXPONENT,CRITICAL_SHEAR_VARIABILITY
        REAL(4) ::  TEMPQCONSTANT,LASTQCONSTANT,MINIMUM_BEDROCK_GRADIENT,LOG_EVAP_MEAN,LOG_EVAP_SD
        REAL(4) ::  OCEANNEXTTIME,EVAPORATION_RATE,ELAPSED_TIME,DIVERGENCE_SCALE_PARAMETER
        REAL(4) ::  OCEANTSLOPE,OCEANTINT,MAXIMUM_SIMULATION_TIME
        REAL(4) ::  SEDIMENT_YIELD_TIMESTEP_SCALING,MASS_WASTING_TIMESTEP_SCALING,TIMECHANGE
        REAL(4) ::  ROUTETHRSHLD,OCEAN_ELEVATION,ALVCREEPFAC,DEPRESSFAC
        REAL(4) ::  STICKYFACTOR,STICKY_ROUTING_CRITICAL_VALUE,STICKYPROB(11),STICKYCAT(11)
        REAL(4) ::  EROSION_RATE_CHANGE_LAG,BISTABLE_CRITICAL_SHEAR,DISCHFACT,HIGH_EROSION_THRESHOLD
        REAL(4) ::  ROCK_WEATHERING_RATE,CRITICAL_BEDROCK_GRADIENT,RGRADMAX,STICKYDEL(11)
        REAL(4) ::  RDELRATE(102),RGRADCAT(102),RRATECAT(102),REGOLITH_CRITICAL_SHEAR_FACTOR
        REAL(4) ::  WEATHER_DIVERGENCE,WEATHER_DECAY_RATE,OMEGA_WEIGHT,MINIMUM_PELAGIC_LAKE_SIZE
        REAL(4) ::  OCORRECT,PREVIOUS_DISCHARGE_COEFFICIENT,PREVIOUS_CRITICAL_SHEAR
        REAL(4) ::  EREFERENCE,VERTICAL_RESISTANCE_SCALING,RANDMULT,SIGMANORM,SIGMASQ
        REAL(4) ::  WEATHERING_DECAY_2,WEATHERING_TERM_2,RUNOFF_CONSTANT_1,RUNOFF_CONSTANT_2
        REAL(4) ::  RUNOFF_CONSTANT_3,DIVERGENCE_FOR_MEAN_RUNOFF
        REAL(4) ::  DEFORMSCALE,BISTABLE_BEDROCK_ERODIBILITY,CELL_SIZE,WEIGHTS(2,2,4)
        REAL(4) ::  SURFACE_LAYER_THICKNESS,GRADCUT,CUMULEXCESS,LAST_TIME_INCREMENT
        REAL(4) ::  EOLIAN_EVENT_PROBABILITY,LAVA_EVENT_PROBABILITY,IMPACT_PROBABILITY,MAXIMUM_TIME_INCREMENT
        REAL(4) ::  RAINDEPTH,RAINSTD,LOW_EROSION_THRESHOLD,MINIMUM_TIME_INCREMENT
        REAL(4) ::  DEPOSITWORK,ERODEWORK,CRATERWORK,LAVAWORK,LAVADEPTH
        REAL(4) ::  EOLIANDEPTH,ELAPSEDTIME,ERODETOADD,SLOPETOADD
        REAL(4) ::  INITIAL_REGOLITH_THICKNESS,SLOPEWORK,SLGRAVTOADD
        REAL(4) ::  PREVIOUS_TIME_INCREMENT,BEDLOAD_FRACTION,EFFECTIVE_DISCHARGE_RATIO
        REAL(4) ::  INFLUENT_RIVER_DISCHARGE(10),INFLUENT_RIVER_SEDIMENT_LOAD(10)
        REAL(4) ::  DEPOSITGRAV,ERODEGRAV,CRATERGRAV,ERGRAVTOADD
        REAL(4) ::  SLOPEGRAV,VERTICAL_SCALING,ROCK_TENSILE_STRENGTH,GROUNDWATER_SCALE_FACTOR
        REAL(4) ::  SEDBIAS,LASTSEDBIAS,SKLAR_MULT,EPOWER,X_RATIO
        REAL(4) ::  VEGETATION_UPLAND_RESISTANCE,VEGETATION_CHANNEL_RESISTANCE
        REAL(4) ::  VEGETATION_AREA_MINIMUM,VEGETATION_AREA_MAXIMUM,VEGETATION_FACTOR_SLOPE
        REAL(4) ::  MEDIAN_SLOPE,SLOPE_RUNOFF_SCALE_FACTOR
        REAL(4) ::  VEGETATION_FACTOR_INTERCEPT,MAXIMUMELEVATION
        REAL(4) ::  SEEPAGE_WEATHERING_SCALING,SEEPAGE_WEATHERING_EXPONENT,DISCHARGE_SCALE_FACTOR
        REAL(4) ::  WASHLOAD_FRACTION,PELAGICCREEP,WEATHER_MULT,FLOWWORK,FLOWGRAV,FLOWTOADD,FLGRAVTOADD
        REAL(4) ::  AVGRUNOFF,AVGSEEP,AVGSEEPFRACT,NSEEPAVG,NFRACTAVG
        REAL(4) ::  DEPOSITFRACTION,STEADY_DISCHARGE_FACTOR, SEDIMENT_RUNOFF_FACTOR
        REAL(4) ::  EJECTA_FRACTION_RETAINED,VOLUME_CHANGE_COEFFICIENT,CONVERT_DISCHARGE
        REAL(4) ::  AVALANCHE_RATE_CONSTANT,AVALANCHE_SLOPE_EXPONENT,AVALANCHE_FLUX_EXPONENT,AVALANCHE_CRITICAL_VALUE
        REAL(4) ::  CREEP_RATE_HALF_DEPTH, CREEP_DEPTH_CONSTANT, SUN_ANGLE_GRADIENT,ITERATION_MAXIMUM_SEDIMENT_YIELD
        REAL(4) ::  RUNOFF_SCALE_BEDROCK, RUNOFF_SCALE_DEEP_REGOLITH, RUNOFF_DEPTH_DECAY_RATE
        REAL(4) ::  MINIMUM_ROUTING_GRADIENT 
        REAL(4) ::  CUM_DEPOSIT_WORK,CUM_ERODE_WORK,CUM_CRATER_WORK,CUM_LAVA_WORK,CUM_SLOPE_WORK,CUM_FLOW_WORK
        REAL(4) ::  CUM_DEPOSIT_GRAV,CUM_ERODE_GRAV,CUM_CRATER_GRAV,CUM_LAVA_GRAV,CUM_SLOPE_GRAV,CUM_FLOW_GRAV
        REAL(4) ::  CUM_ACCRETION,ACCRETE_WORK,LAVAGRAV
        REAL(4) ::  FRAC2312, FRAC43, FRAC512, FRAC23, FRAC112   ! added by (orkan) April 4
        REAL(4) ::  WEIGHTS_3(3,3),WEIGHTS_5(5,5),WEIGHTS_7(7,7),WEIGHTS_9(9,9)
        REAL(4) ::  THE_DIAG_WEIGHT(100),THE_CROSS_WEIGHT(100)
        LOGICAL(4) ::  DONE,WRITE_ABSOLUTE_ELEVATION,USE_CRITICAL_SLOPE_CRADIENT,ABORT_SIMULATION,NO_LAKES_EXIST
        LOGICAL(4) ::  USELAYER,USE_3D_SLOPE_RESISTANCE,USE_SLOPE_DEPENDENT_RUNOFF
        LOGICAL(4) ::  WRITE_OUT_CHANNELS,USEBEDRILL,LDUMMY,DO_DEBUGGING
        LOGICAL(4) ::  DO_SEDIMENT_TRANSPORT,DO_SEDIMENT_DIFFUSION,DO_SEDIMENT_ROUTING,COMPLETE_RUNOFF,DO_MODEL_SLOPES
        LOGICAL(4) ::  DO_MORPHOMETRY,USE_AN_OCEAN,DO_DEMON_FLOW_ROUTING,USE_RANDOM_FLOWROUTING,STICKY_SEDIMENT_ROUTING
        LOGICAL(4) ::  VARIABLE_EROSION_RATE,RANDOM_CRITICAL_SHEAR,USE_RANDOM_DISCHARGE,BISTABLE_FLUVIAL_EROSION
        LOGICAL(4) ::  TWO_TERM_WEATHERING,DIVERGENCE_DEPENDENT_RUNOFF,NEW_SIMULATION
        LOGICAL(4) ::  DO_ROCK_DEFORMATION,SCALE_3D_ROCK_ERODIBILITY,NO_FLUX_LOWER_BOUNDARY,HORIZONTAL_LOWER_BOUNDARY
        LOGICAL(4) ::  NON_ERODING_LOWER_BOUNDARY,USE_BISTABLE_BEDROCK,IS_Y_PERIODIC,RESISTANT_SURFACE_LAYER
        LOGICAL(4) ::  EXPLICIT_CHANNEL_BED_STATE,DO_ALLUVIAL_SMOOTHING,DO_ALLUVIAL_REROUTING
        LOGICAL(4) ::  MODEL_LAVA_FLOWS,MODEL_IMPACT_CRATERING,MODEL_EOLIAN_CHANGES,FLUVIAL_AND_SLOPE_MODELING,DO_FLUVIAL_DETACHMENT
        LOGICAL(4) ::  ROUTE_REGOLITH_OVER_ROCK,WRITE_SEDIMENT_DIAGNOSTICS,SEDSWITCH
        LOGICAL(4) ::  VARIABLE_OCEAN_ELEVATION,HAVE_INFLUENT_RIVERS, USE_SEDIMENT_DEPENDENT_RUNOFF
        LOGICAL(4) ::  MODEL_OCEAN_LEVEL,DISCHARGE_WEATHERING , DO_SHADE_BORDER, IS_FIXED_SUNANGLE
        LOGICAL(4) ::  USE_ROERING_MASS_WASTING,USE_SKLAR_BED_ABRASION,USE_WHIPPLE_BED_ABRASION
        LOGICAL(4) ::  MODEL_ACCRETION_AND_ABLATION,USE_TOTAL_EXPOSURE,USE_DIVERGENCE_EXPOSURE,USE_DISCHARGE
        LOGICAL(4) :: IS_X_PERIODIC,FORCE_SEDIMENT_CONSERVATION,RESCALE_DISCHARGES
        LOGICAL(4) :: PERMEABILITY_RESCALING,DEFAULT_EOLIAN_PROCESS,SEEPAGE_WEATHERING,MODEL_LAKE_EVAPORATION
        LOGICAL(4) :: MODEL_PELAGIC_DEPOSITION,DO_FLOW_BOUNDARIES,CRATERING_CALLED
        LOGICAL(4) :: USE_SEDIMENT_YIELD_SCALING, USE_REGOLITH_DEPTH_DEPENDENT_RUNOFF
        LOGICAL(4) :: USE_SKLAR_SATURATION,USE_EROSION_MASK,MODEL_GROUNDWATER,DO_SPATIAL_VARIATION, FIRST_IMAGE
        LOGICAL(4) :: REFLECTIVE_UPPER_BOUNDARY,USE_SPATIAL_RUNOFF,USE_SPATIAL_WEATHERING
        !      parameters for groundwater flow and sapping
          LOGICAL(4) ::  VARIABLE_VEGETATION_RESISTANCE,DO_EVENTS
        LOGICAL(4) ::  IS_GRAVEL_MIXTURE, DO_REGOLITH_ABLATION, USE_SHEAR_RATIO
        LOGICAL(4) ::  DO_AVALANCHE,DEPTH_DEPENDENT_CREEP
        LOGICAL(4) :: WRITE_BINARY_IMAGE_FILE,USE_FLOW_VOLUME,THIS_IS_BINGHAM_FLOW
 
        REAL(4) :: ALLUVIUM_SMOOTHING_FACTOR
        REAL(4) :: DEFAULT_TIME_INCREMENT
           !      end parameters for groundwater flow and sapping
        CHARACTER (LEN=6) :: RPREFIX
        CHARACTER (LEN=4) :: X1PREFIX, X2PREFIX , X3PREFIX, RSUFFIX,XSUFFIX
        CHARACTER (LEN=30):: X1FILENAME ,X2FILENAME ,X3FILENAME
        CHARACTER (LEN=1) :: XNUM1 ,XNUM2, XNUM3, RNUM3,RNUM2, RNUM1,RNUM4,XNUM4
        CHARACTER (LEN=7) :: XDIRECT
        CHARACTER (LEN=25) :: RFILENAME,XFILENAME
        REAL(4),ALLOCATABLE, DIMENSION(:) :: BASIN_DRAINAGE_AREA,LAKE_OUTLET_ELEVATION
        REAL(4),ALLOCATABLE, DIMENSION(:) :: LAKE_SURFACE_ELEVATION,SORTING_VECTOR
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: ELEVATION,D8_GRADIENT,ERODE_CHANNEL,ERODE_SLOPE
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: RELATIVE_RESISTANCE,CFW,CFN,CFNE,CFNW
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: CHANNEL_WIDTH,SEDIMENT_BASE,SEDIMENT_YIELD
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: SEDIMENT_FLUX,EQUILIBRIUM_GRADIENT,DRAINAGE_AREA
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: CUMULATIVE_ELEVATION_CHANGE, NOMINAL_ERODE_SLOPE
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: CUMULATIVE_EOLIAN_CHANGE,CUMULATIVE_LAVA_CHANGE
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: CUMULATIVE_CRATERING_CHANGE,MAXIMUM_SEDIMENT_YIELD
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: PREVIOUS_ELEVATION,DISCHARGE,MAXIMUM_DISCHARGE,INITIAL_ELEVATION
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: DIVERGENCE,ERODE_REGOLITH_CHANNEL
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: REGOLITH,DEFORMATION, LOCAL_DISCHARGE
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: FILTERED_GROUNDWATER_FLUX,TRANSMISSIVITY_TERM
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: WATER_ELEVATION,GROUNDWATER_FLUX,SPATIAL_VALUE
        REAL(4),ALLOCATABLE, DIMENSION(:) :: OCEAN_RECALCULATION_TIMES, OCEAN_LEVELS
        REAL(4),ALLOCATABLE, DIMENSION(:) :: PELAGIC_SEDIMENT_VOLUME
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: PELAGICCHANGE,PELAGICEDIFF
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: PLOTVALS
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: CUMULATIVE_WEATHERING,LAPLACE
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: CUMULATIVE_EJECTA_DEPOSITION, CUMULATIVE_CRATER_EXCAVATION
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: CUMULATIVE_SEDIMENT_DEPOSITION, CUMULATIVE_MASS_WASTING
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: CUMULATIVE_FLUVIAL_EROSION,CUMULATIVE_AVALANCHE_EROSION
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: CUMULATIVE_FLOW_VOLUME_CHANGE
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: AVALANCHE_FLUX, OLD_AVALANCHE_FLUX
        REAL(4),ALLOCATABLE, DIMENSION(:,:) :: OUT_32_DATA
        REAL(4), ALLOCATABLE,DIMENSION(:,:) :: REFERENCE_ELEVATION
        REAL(4), ALLOCATABLE,DIMENSION(:,:) :: SPATIAL_VARIATION
        INTEGER(4),ALLOCATABLE, DIMENSION(:) :: I_OUTFLOW,J_OUTFLOW,DOWNSTREAM_BASIN
        INTEGER(4),ALLOCATABLE, DIMENSION(:,:) :: BASIN_NUMBER,FLOW_DIRECTION,IDO,EROSION_DEPTH_INDEX
        INTEGER(4),ALLOCATABLE, DIMENSION(:,:) :: CHANGE_REGOLITH_STATE
        LOGICAL(4),ALLOCATABLE, DIMENSION(:) :: ENCLOSED,OVERFLOWS
        LOGICAL(4),ALLOCATABLE, DIMENSION(:,:) :: IS_SEDIMENT_COVERED,SUBMERGED,EROSION_MASK
        LOGICAL(4),ALLOCATABLE, DIMENSION(:,:) :: ACCELERATED_EROSION,DO_ACCELERATED_EROSION,IS_ROCK_SURFACE
        LOGICAL(4),ALLOCATABLE, DIMENSION(:,:) :: SEDIMENT_DEPOSITED,IS_INFLUENT_RIVER_LOCATION
        LOGICAL(4),ALLOCATABLE, DIMENSION(:,:) :: FIXED_HEAD,CRATER_RESET
        LOGICAL(4),ALLOCATABLE, DIMENSION(:) :: ERODING_LOWER_BOUNDARY,IS_EXIT_BASIN
        LOGICAL(4),ALLOCATABLE, DIMENSION(:,:) :: IS_ICE    
    END MODULE ERODE_GLOBALS

        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MODULE OCEAN_GLOBALS
        IMPLICIT NONE
        SAVE
        REAL(4),ALLOCATABLE,DIMENSION(:,:) :: ELL2,FLUX_X,FLUX_Y
        LOGICAL(4),ALLOCATABLE,DIMENSION(:,:) :: TURB
        REAL(4), PARAMETER :: D=1.0
        REAL(4), PARAMETER :: DL=5.0
    END MODULE OCEAN_GLOBALS
           !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MODULE EOLIAN_GLOBALS
        IMPLICIT NONE
        SAVE
        REAL(4) :: ZERO_PERCENT_EXPOSURE,EXPOSURE_50_PERCENT,EXPOSURE_90_PERCENT
        REAL(4) :: MINIMUM_EOLIAN_DEPOSIT_RATE,MAXIMUM_EOLIAN_DEPOSIT_RATE
        REAL(4) :: EXPOSURE_10_PERCENT,DEPOSITRATE
        REAL(4) :: EOLIAN_CONSTANT_1,EOLIAN_CONSTANT_2,EOLIAN_CONSTANT_3
        REAL(4) :: DISTANCE_DECAY_FACTOR,WEIGHTING_DECAY_FACTOR,RATE0
        INTEGER(4) :: NITER,NPRINT,MAXIMUM_WEIGHT_DISTANCE,WEIGHTING_CALCULATION_DISTANCE,ICNNT,JCNNT
        REAL(4) :: EOLIAN_TIME_INCREMENT,WEIGHTX(350),WEIGHTD(350)
        REAL(4) :: MAXWEATHERRATE
        REAL(4) :: DISTANCE_WEIGHTING(51,51),DISTANCE_VALUE(51,51)
        LOGICAL :: ONLY_EOLIAN_DEPOSITION
    END MODULE EOLIAN_GLOBALS
            !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MODULE AREA_GLOBALS
        IMPLICIT NONE
        SAVE
        REAL(4), ALLOCATABLE, DIMENSION(:) :: TOBEADDED,LOCAL_BASIN_DISCHARGE,QTOADD
        LOGICAL(4), ALLOCATABLE, DIMENSION(:) :: NEEDTODO,OUTER
        REAL(4) :: RUNOFF_AVAL,RUNOFF_BVAL,SPATIAL_MIN,SPATIAL_MAX
        REAL(4) :: WEATHERING_AVAL,WEATHERING_BVAL
    END MODULE AREA_GLOBALS
	!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
    MODULE LAKE_GLOBALS
        IMPLICIT NONE
        SAVE
        REAL(4), ALLOCATABLE, DIMENSION(:,:) :: ROUTED_DISCHARGE
        REAL(4), ALLOCATABLE, DIMENSION(:) :: LAKE_AREA,LAKE_VOLUME
        REAL(4), ALLOCATABLE, DIMENSION(:) :: LOWEST_LAKE_ELEVATION,NEW_BASIN_OUTFLUX,BASIN_INFLUX
        REAL(4), ALLOCATABLE, DIMENSION(:) :: BASIN_OUTFLUX,OLDOVERFLOWS
        REAL(4), ALLOCATABLE, DIMENSION(:) :: LAKEMIN,BASINMIN
        LOGICAL(1), ALLOCATABLE,DIMENSION(:) :: NEXTCYCLE,NEWGEOMETRY
        LOGICAL(1), ALLOCATABLE,DIMENSION(:) :: ISBORDER
    END MODULE LAKE_GLOBALS	
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

MODULE GRAVEL_MIXTURE_GLOBALS
IMPLICIT NONE
SAVE
! Choose Parker (1) ou Wilcock Crowe (2)

Integer ::  EQUATION_SELECTOR,GD
INTEGER :: GRAIN_ARRAY_SIZE, STRAIN_CURVE_SIZE, STEP_BY_MARSSIM_ITERATION

Character(Len=100) elev_WRITE, Dsg_WRITE, SED_FLUX_WRITE
INTEGER ::  OUTSEDFLUX, OUTDSGS, OUTD90, OUTSAND, GRAVEL_ITERATION
REAL(4) ::  abrade_term,aterm_avg,a1term,a2term,a1_avg,a2_avg,agterm_avg,avg_change,gd1_avg,gd2_avg,net1_avg,net2_avg,delay_weight
REAL(4) :: el_avg,sf_avg,ew_avg
    REAL(4), ALLOCATABLE, DIMENSION(:) :: Fnew_avg,F_diff
REAL(4), ALLOCATABLE, DIMENSION(:) :: STRAIN_CURVE_PO, STRAIN_CURVE_OO, STRAIN_CURVE_SO


REAL(4) :: ROUGHNESS_FACTOR, ACTIVE_LAYER_FACTOR, MANNING_COEF, UPWIND_COEF, AGGRADATION_COEF
REAL(4) :: INTERMITENCY, SEDIMENT_CONSTANTE,SLOPE_MUD_FRACTION, ABRASION_COEFFICIENT

REAL(4) :: time, phisgo, sas, omega
REAL(4) :: ym1, y, x, xm1, yp1, xp1, yp2, xp2



REAL(4), ALLOCATABLE, DIMENSION(:) :: psi, ds, plf, Fs
REAL(4), ALLOCATABLE, DIMENSION(:) :: di, plff, FfI, Ffs
!ds and psi calculated from di, plf(I, J)/FI(I, J)/Fs(I, J): are the proportion of grain of size between di(I, J) and di(i+1)



REAL(4), ALLOCATABLE, DIMENSION(:,:) :: ACTIVE_LAYER_THICKNESS, ACTIVE_LAYER_THICKNESSold, BED_ELEVATION_CHANGE_RATE
REAL(4), ALLOCATABLE, DIMENSION(:,:) :: SEDIMENT_FLUX_SLOPES, SEDIMENT_FLUX_RIVER
REAL(4), ALLOCATABLE, DIMENSION(:,:) :: dsgs, D90_SIZE, D50_SIZE, fracsand, SHIELD_NUMBER, taus50, H
!ACTIVE_LAYER_THICKNESS: Thickness of the active layer

REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: F, pl
!F: store the grain distribution of the active layer for each interval
LOGICAL :: PRINT_GRAVEL_DETAILS,WIDTH_RELATIVE

 END MODULE GRAVEL_MIXTURE_GLOBALS
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	MODULE GROUNDWATER_VARIABLES
	IMPLICIT NONE
	SAVE
	    LOGICAL(4) ::  SHOW_GROUNDWATER_CALCULATIONS,FIRST_GROUNDWATER_CALL
        LOGICAL(4) ::  USE_GROUNDWATER_FLOW,SEEPAGE_AVERAGING,USEIMPLICIT
        LOGICAL(4) ::  EXPONENTIAL_PERMEABILITY_DECAY
        REAL(4) ::  GROUNDWATER_DEPTH_SCALE,MAXIMUM_GROUNDWATER_ERROR,GROUNDWATER_RECHARGE_RATE
        REAL(4) ::  GROUNDWATER_FLOW_FRACTION,GROUNDWATER_RELAXATION_CONST
        REAL(4) ::  INITIAL_GROUNDWATER_DEPTH
        INTEGER(4) :: MAXIMUM_GROUNDWATER_ITERATIONS ,SEEPAGE_ITERATION_INTERVAL
        REAL(4) :: DARCIES,VISCOSITY,METRIC_PERMEABILITY,YEARLY_RECHARGE
    END !GROUNDWATER_VARIABLES 
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MODULE CRATER_GLOBALS
        IMPLICIT NONE
        SAVE
        INTEGER(4) :: KSIZE
        INTEGER(4) :: TOPFILE,SHADEFILE,BASEFILE
        INTEGER(4) :: NIMPACTS,ITERNO,ITOTHITS,DIFFINTERVAL
        INTEGER(4) :: NOISEINTERVAL,NUMBER_REAL_CRATERS
        REAL(4) :: CRATER_FREQUENCY_EXPONENT,ALINV,SMALLEST_MODELED_CRATER,LARGEST_MODELED_CRATER,CCON,BCON
        INTEGER(4) :: IMIN,IMAX,JMIN,JMAX,RCATS(10),OCATS(10)
        INTEGER(4) :: COUNT2,COUNT10
        INTEGER(4) :: IMINR,IMAXR,JMINR,JMAXR
        REAL(4) :: DIAMETER,XMIN,YMIN,XMAX,YMAX,XCENTER,YCENTER
        REAL(4) :: CRATER_DEPTH,RIM_HEIGHT,INTERIOR_SHAPE_EXPONENT,EXTERIOR_SHAPE_EXPONENT
        REAL(4) :: RADIUS,MAXRANGE,ELEVBASE,SMALLEST_POSSIBLE_CRATER
        REAL(4) :: EJECTA_THICKNESS_VARIABILITY,NOISESD,INHERITANCE_PARAMETER,FREQUENCY_CUTOFF_SCALING
        REAL(4) :: DCATS(10),LARGE_CRATER_DEPTH_SCALE,LARGE_CRATER_DEPTH_EXPONENT
        REAL(4) :: LARGE_CRATER_RIM_SCALE,LARGE_CRATER_RIM_EXPONENT,INITGRAD
        REAL(4) :: SMALL_CRATER_DEPTH_SCALE,SMALL_CRATER_DEPTH_EXPONENT,SMALL_CRATER_RIM_SCALE,SMALL_CRATER_RIM_EXPONENT
        REAL(4) :: SMALL_CRATER_SHAPE_SCALE,SMALL_CRATER_SHAPE_EXPONENT,LARGE_CRATER_SHAPE_SCALE,LARGE_CRATER_SHAPE_EXPONENT
        REAL*4  :: PEAK_HEIGHT_MULT,PEAK_HEIGHT_EXPONENT,PEAK_DIAMETER_MULT,PEAK_DIAMETER_EXPONENT
        REAL*4  :: CENTRAL_PEAK_HEIGHT,CENTRAL_PEAK_DIAMETER,CENTRAL_PEAK_VOLUME,INHERIT_EXPONENT		
        REAL(4) :: SIZES(1000),XLOCS(1000),YLOCS(1000),TRANSITION_DIAMETER,MINIMUM_HARD_DIAMETER
        REAL(4) :: MAXIMUM_RIM_GRADIENT,MINIMUM_REAL_CRATER_DIAMETER,MAXIMUM_REAL_CRATER_DIAMETER,RADIUS_MAX_INHERIT
        REAL(4) :: RADIUS_MAX_USE,COSINE_POWER
        REAL(4), ALLOCATABLE, DIMENSION(:,:) :: CRATER_EVENT
		REAL(4), ALLOCATABLE, DIMENSION(:,:) :: BWEIGHT
        LOGICAL(4) :: RANDOM_EJECTA_THICKNESS,MICRONOISE,NOISECALL,DIFFUSECALL,CORRECT_BIAS,MAKE_CENTRAL_PEAK
        LOGICAL(4) :: DO_EJECTA_WRAPAROUND, IS_REGOLITH_CRATER,CRATER_EDGE_ABORT,USE_REAL_CRATERS,USE_CRATER_EVENT
 
        REAL(4) ::  CRATERDIAM(400)
        CHARACTER(160) :: CRATER_DATABASE_LOCATION,CRATERFILENAMES(400)
    END MODULE CRATER_GLOBALS
    
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
		
    MODULE EVENT_GLOBALS
        IMPLICIT NONE
        SAVE
        LOGICAL (4), ALLOCATABLE,DIMENSION(:) :: EVENT_DONE
        REAL(4),ALLOCATABLE, DIMENSION(:) :: EVENT_TIMES
        INTEGER, ALLOCATABLE, DIMENSION(:) :: EVENT_ITERATIONS
        ! ANY VARIABLES NECESSARY TO MODEL DISCRETE EVENTS CAN BE ADDED HERE
        REAL(4) :: EVENT_TIME_SCALE
        REAL(4), ALLOCATABLE, DIMENSION (:) :: DIAMETERS, CRATER_X_LOCATIONS,CRATER_Y_LOCATIONS
		REAL(4), ALLOCATABLE, DIMENSION (:) :: VISCOUS, CRITDEPTHS,MAXDEPTHS,ALPHAS,SCOURS,WEATHERS
    END MODULE EVENT_GLOBALS
         !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        MODULE LAVA_GLOBALS
        IMPLICIT NONE
        SAVE
        !      ********************************************
        !          Include file for lava flow program
        !              defines global variables
        !      ********************************************
        INTEGER(4) ::  NUMBER_OF_ACTIVE_LAVA_CELLS,ISTART,JSTART
        INTEGER(4) :: LAVAMAXITER,LAVAITER
        INTEGER(4) :: NEW_SEGMENT_INTERVAL,IBACK(9),SOURCECOUNT
        INTEGER(4) :: NUMBER_OF_LAVA_SOURCES,I_LAVA(100),J_LAVA(100)
        INTEGER(4) :: SOURCE_CHANGE_INTERVAL
        INTEGER(4) :: ICONTINUE
        LOGICAL(1),ALLOCATABLE, DIMENSION(:,:) :: IS_LAVA_COVERED,ACTIVE_LAVA_FLOW
        INTEGER(4),ALLOCATABLE, DIMENSION(:,:) :: LAVA_SOURCE_DIRECTION,ERUPTION_AGE
        REAL(4), ALLOCATABLE, DIMENSION(:) :: LAVA_THICKNESS
        REAL(4), ALLOCATABLE, DIMENSION(:,:) :: LAVA_FLOW_PROBABILITY
        INTEGER(4), ALLOCATABLE, DIMENSION(:) :: LAVA_ELAPSED_TIME,ILOC,JLOC
        REAL(4) :: LAVA_DURATION_WEIGHT,LAVA_GRADIENT_WEIGHT,EOLD
        REAL(4) :: MINIMUM_LAVA_FLOW_SLOPE,ERUPTION_STOP_PROBABILITY,LAVA_FLOW_THICKNESS
        REAL(4) :: NO_FLOW_PROBABILITY,TLACTIVESUM,TLACTIVESQ,TACTSAMP
        REAL(4) :: TVOLUMESUM,TVOLUMESQ,TVOLUMESAMP,MINIMUM_FLOW_THICKNESS
    END MODULE LAVA_GLOBALS	
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MODULE SEDROUTE_GLOBALS
        IMPLICIT NONE
        SAVE
        REAL(4) :: DEFAULT_DRAINAGE_AREA,DEFAULT_SEDIMENT_YIELD,SEDIMENT_CARRYOVER
        INTEGER(4) :: ISTART,JSTART,IEND,JEND,INOW,JNOW,IDX
        INTEGER(4) :: K,IABCD,IXT,JYT,IIII,JDY,IXX,JYY
        INTEGER(4) :: INCRX,INCRY,EE,INC1,INC2,IDMAX,ICOUNT
        LOGICAL(4) :: IS_DEPRESSION,IS_THERE,MUST_SEARCH,DO_BACKTRACK
        LOGICAL(4), ALLOCATABLE,DIMENSION(:,:) :: DONE_ONCE
        INTEGER(4), ALLOCATABLE,DIMENSION(:) :: I_LOCATION,J_LOCATION
    END MODULE SEDROUTE_GLOBALS	
	
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MODULE SEDDEBUG_GLOBALS
        IMPLICIT NONE
        SAVE
        REAL(4) :: UPSTREAM_ELEVATION,ACTUAL_START_GRADIENT,ALLUVIAL_START_GRADIENT,SEDIMENT_VOLUME,PROVISIONAL_START_GRADIENT
        REAL(4) :: PROVISIONAL_END_GRADIENT,EL1,EL2,START_STEP_DISTANCE,MAXCHG,END_ELEVATION
        INTEGER(4) :: KKK,I_NEW,J_NEW,I_OLD,J_OLD,IMAX,JMAX
        REAL(4), ALLOCATABLE, DIMENSION(:) :: STARTING_ELEVATION,ALLUVIAL_GRADIENT,ACTUAL_GRADIENT
        REAL(4), ALLOCATABLE, DIMENSION(:) :: STEP_DISTANCE,NEW_ELEVATION,PROVISIONAL_ELEVATION,WATER_LEVEL
        LOGICAL(4) , ALLOCATABLE, DIMENSION(:) :: IS_BEDROCK_CHANNEL
    END MODULE SEDDEBUG_GLOBALS
   	        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MODULE ACCRETION_GLOBALS
        IMPLICIT NONE
        SAVE
        LOGICAL :: USE_INVERSE_EXPOSURE, EXPOSURE_DEPENDENT_CREEP, USE_EXPOSURE_SMOOTHING,USE_SOLAR_EROSION,USE_TOP_EXPOSURE
        LOGICAL :: DO_NORMAL_ACCRETION, DO_ACCRETE_REGOLITH
        REAL(4) :: RAD_CONST, RAD_DUST_FACTOR,RAD_THRESH_CONVEXITY,RAD_DEPOSIT_RATE,ACCRETION_RATE,TOTALRADERODE
        REAL(4), ALLOCATABLE, DIMENSION(:,:) :: EXPOSE_RATE,SMOOTH_EXPOSE
    END MODULE ACCRETION_GLOBALS 
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MODULE MASS_FLOW_GLOBALS
        IMPLICIT NONE
        SAVE
        REAL(4) :: GLEN_LAW_ALPHA, FLOW_DIFFUSIVITY, DELTA_HB, AVG_FLOW_VOLUME_FLUX, FRAC_BINGHAM_ACTIVE    ! modified by (orkan) March 25
        REAL(4) :: FLOW_VOLUME_INCREMENT, AVG_FLOW_VOLUME_CHANGE, MAXIMUM_FLOW_DEPTH_EROSION,bingham_maximum_thickness
        REAL(4) :: MASS_FLOW_EROSION_RATE,NUM_FLUX !  MASS_FLOW_EROSION_RATE ADDED BY ORKAN 12/19/2014
        REAL(4) :: FLUX_EXPONENT, DEPTH_EXPONENT, SINE_EXPONENT, MASS_FLOW_CRITICAL_VALUE,CRITICAL_FLOW_THICKNESSS
        REAL(4) :: STERM_AVERAGE,STERM_MAX,STERM_MIN,STERM_NUMBER,FLOW_VOLUME_DENSITY,CRITICAL_SCOUR_THICKNESS
        real(4) :: TARGET_MASS_FLOW_CHANGE_RATE,  MAX_FLOW_VOLUME_FLUX,MAX_FLOW_VOLUME_CHANGE
        LOGICAL :: scheme1,scheme2,scheme3,use_variable_FLOW_VOLUME_rate
        integer :: flow_scheme
        REAL(4), ALLOCATABLE, DIMENSION(:,:) :: FLOW_VOLUME_CHANGE, FLOW_VOLUME_EROSION_BASE 
        REAL(4), ALLOCATABLE, DIMENSION(:,:) :: FLOW_VOLUME_FLUX,CUMULATIVE_FLOW_VOLUME_FLUX                  ! modified by (orkan) 01/05/2015
    END MODULE MASS_FLOW_GLOBALS





