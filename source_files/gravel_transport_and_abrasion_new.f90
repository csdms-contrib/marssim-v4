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
    !This program is a subroutine for the MARSSIM program. It model the transport of a sediment mixture.
!This program is adapted from the Gary Parker model AgDegNormGravMixPW
!It allows the modelisation of agradated or degradated landscapes but don't acuaratly deals with mixed phenomenon. In order to do it the grain size of the agradated material must be stored in a matrix.
!For more detail about the formulation of the problem refer to "surface-based bedload transport relation for gravel-river", Gary parker (1990)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE READ_GRAVEL_TRANSPORT_PARAMETERS
  USE ERODE_GLOBALS
  USE GRAVEL_MIXTURE_GLOBALS
  IMPLICIT NONE
  INTEGER :: PARAMS_MIXTURE, GRAVEL_PRINT,WIDTH_RELATIVE_INDEX,DO_GRAVEL_MIXTURE,IERROR
  INTEGER :: PARAMS_TEXT
  CHARACTER(80) :: TEXTLINE
        PARAMS_MIXTURE=93
        PARAMS_TEXT=94
!        *******************************************************************
!        *******************************************************************
!                        Gravel mixture transport parameters
!        *******************************************************************
!        *******************************************************************

    OPEN(PARAMS_MIXTURE,FILE='GRAVEL_MIXTURE.PRM',ACTION='READ')

    READ(PARAMS_MIXTURE,22) TEXTLINE
22  FORMAT(A)    
    WRITE(*,22) TEXTLINE
    READ(PARAMS_MIXTURE,*) DO_GRAVEL_MIXTURE
    IF (DO_GRAVEL_MIXTURE>0) THEN
        OPEN(PARAMS_TEXT,FILE='GRAVEL_VALUES.TXT',ACTION='READ')
    ENDIF
    READ(PARAMS_MIXTURE,*) GRAVEL_PRINT
    READ(PARAMS_MIXTURE,*) WIDTH_RELATIVE_INDEX
    
    IF (GRAVEL_PRINT>0) then
        PRINT_GRAVEL_DETAILS=.true.
    else
        PRINT_GRAVEL_DETAILS=.false. 
    endif
    IF (WIDTH_RELATIVE_INDEX>0) then
        WIDTH_RELATIVE=.true.
    else
        WIDTH_RELATIVE=.false. 
    endif
    IF (DO_GRAVEL_MIXTURE>0) then
        DO_FLUVIAL_DETACHMENT=.false. !may be temporary restriction
        DO_SEDIMENT_TRANSPORT=.false. ! may be more integrated in fluvial transport in future
        DO_SEDIMENT_ROUTING=.false.
        IS_GRAVEL_MIXTURE=.true.
        ! 1 for Parker formula, 2 for Wilcock-Crowe
        READ(PARAMS_MIXTURE,*) EQUATION_SELECTOR
        ! number of terms in Parker formula for grain size specific transport
        READ(PARAMS_MIXTURE,*) STRAIN_CURVE_SIZE
        ALLOCATE(STRAIN_CURVE_PO(STRAIN_CURVE_SIZE), STRAIN_CURVE_OO(STRAIN_CURVE_SIZE), &
		     STRAIN_CURVE_SO(STRAIN_CURVE_SIZE),STAT=IERROR)
		IF (IERROR > 0) CALL ALLOCERROR()
        ! read the parameters for the "strain curves"
        DO GD=1,STRAIN_CURVE_SIZE
            READ(PARAMS_TEXT,*) STRAIN_CURVE_PO(GD), STRAIN_CURVE_OO(GD), STRAIN_CURVE_SO(GD)
        END DO
        ! read how many grain size cumulative values will be read in
        READ(PARAMS_MIXTURE,*) GRAIN_ARRAY_SIZE
        ! These variables are for the grain size bins between the cumulative values (and thus one less than GRAIN_ARRAY_SIZE)
! ds(GD) are the geometric mean grain sizes, in mm
! psi(GD) are the corresponding phi (logarithmic) grain sizes
! plf(GD) is the initial proportions of the upstream feed
! F(i,j,GD) are the initial surface layer gravel size proportions for the active layer for each of the spatial I,J cells
! Fs(GD) is the subsurface grain size proportions, assumed to be areally constant - and is assumed to also represent sediment input to 
!   the channel from surrounding slopes
        ALLOCATE(psi(GRAIN_ARRAY_SIZE-1), ds(GRAIN_ARRAY_SIZE-1), plf(GRAIN_ARRAY_SIZE-1), Fs(GRAIN_ARRAY_SIZE-1),STAT=IERROR)
		IF (IERROR > 0) CALL ALLOCERROR()
        ALLOCATE(di(GRAIN_ARRAY_SIZE), plff(GRAIN_ARRAY_SIZE), FfI(GRAIN_ARRAY_SIZE), Ffs(GRAIN_ARRAY_SIZE),STAT=IERROR)
		IF (IERROR > 0) CALL ALLOCERROR()
! read gravel grain size boundaries Di and cumulative size percentages for inital gravel inflow from upstream, plff, 
!    surface layer, FfI, and subsurface layer, Ffs       
        DO GD=1,GRAIN_ARRAY_SIZE
            READ(PARAMS_TEXT,*) di(GD), plff(GD), FfI(GD), Ffs(GD)
        END DO
! These are parameters in the Parker model        
        READ(PARAMS_MIXTURE,*) ROUGHNESS_FACTOR
        READ(PARAMS_MIXTURE,*) ACTIVE_LAYER_FACTOR
        READ(PARAMS_MIXTURE,*) MANNING_COEF
        READ(PARAMS_MIXTURE,*) UPWIND_COEF
        READ(PARAMS_MIXTURE,*) AGGRADATION_COEF
! intermittency is the fraction of the year with sediment transport
! Slope_mud_fraction is what percentage of the substrate and surrounding slopes are transported in suspension and thus not accounted
!   for in the gravel budget
! The abrasion_coefficent determines the comminution rate of the gravel during transport
        READ(PARAMS_MIXTURE,*) INTERMITENCY
        READ(PARAMS_MIXTURE,*) SLOPE_MUD_FRACTION
        READ(PARAMS_MIXTURE,*) ABRASION_COEFFICIENT
! How many sub-iterations of gravel transport are done per main program iteration
        READ(PARAMS_MIXTURE,*) STEP_BY_MARSSIM_ITERATION
		CLOSE(PARAMS_MIXTURE)
        CLOSE(PARAMS_TEXT)
        ! Convert the abrasion coefficient into millimeters      
        ABRASION_COEFFICIENT=ABRASION_COEFFICIENT/1000.0
        IF (DO_GRAVEL_MIXTURE>0) THEN
            WRITE(OUTHIST,1001)
1001        FORMAT('***********************MIXED-SIZE GRAVEL TRANSPORT AND EROSION IS MODELED***************')
            WRITE(OUTHIST,1002) 
1002        FORMAT('GRAVEL TRANSPORT AS IMPLEMENTED CAN ONLY BE USED FOR VERY SMALL GRID SIZE AND SMALL TIME STEPS',/, &
            'This program is adapted from the Gary Parker model AgDegNormGravMixPW',/, &
            'AS PRESENTLY IMPLEMENTED IT DOES NOT PERMIT MODELING OF MARSSIM FLUVIAL EROSION OR SEDIMENT ROUTING')
            IF (PRINT_GRAVEL_DETAILS) THEN
                WRITE(OUTHIST,1003)
1003            FORMAT('DETAILED PRINTOUT IS SELECTED')
            ENDIF
            IF (WIDTH_RELATIVE) THEN
                WRITE(OUTHIST,1004)
1004            FORMAT(' GRAVEL CHANNEL WIDTH IS SMALLER THAN GRID SIZE')
            ENDIF
            IF (EQUATION_SELECTOR==1) THEN
                WRITE(OUTHIST,1005)
1005            FORMAT('PARKER GRAVEL TRANSPORT RELATIONSHIP IS USED')
            ELSE
                WRITE(OUTHIST,1006)
1006            FORMAT('CROWE-WILCOCK GRAVEL TRANSPORT RELATIONSHIP IS USED')
            ENDIF
            WRITE(OUTHIST,1007) STRAIN_CURVE_SIZE, GRAIN_ARRAY_SIZE
1007        FORMAT('SIZE OF STRAIN CURVE ARRAY=',I6,' NUMBER OF GRAIN SIZES=',I6)
            WRITE(OUTHIST,1008) ROUGHNESS_FACTOR,ACTIVE_LAYER_FACTOR,MANNING_COEF,UPWIND_COEF,AGGRADATION_COEF, &
                INTERMITENCY,SLOPE_MUD_FRACTION,STEP_BY_MARSSIM_ITERATION
                1008 FORMAT('ROUGHNESS_FACTOR=',G12.5,/, &
                     'ACTIVE_LAYER_FACTOR=',G12.5,/, &
                     'MANNING_COEFFICIENT=',G12.5,/, &
                     'UPWIND_COEFFICIENT=',I6,/, &
                     'AGGRADATION_COEFFICIENT=',G12.5,/, &
                     'INTERMITTENCY=',G12.5,/, &
                     'SEDIMENT_CONSTANT=',G12.5,/, &
                     'MUD FRACTION=',G12.5,/, &
                     'NUMBER OF STEPS PER MASTER ITERATION=',I6)                     
  
 
            write(outhist,7392) abrasion_coefficient
7392        format('GRAVEL_ABRASION_COEFFICIENT='g13.6)
        ENDIF
! Variables are explained in the gravel transport module        
        ALLOCATE(ACTIVE_LAYER_THICKNESS(MX,MY), ACTIVE_LAYER_THICKNESSold(MX,MY), BED_ELEVATION_CHANGE_RATE(MX,MY))
        ALLOCATE(SEDIMENT_FLUX_SLOPES(MX,MY), SEDIMENT_FLUX_RIVER(MX,MY))
        ALLOCATE(dsgs(MX,MY), D90_SIZE(MX,MY), D50_SIZE(MX,MY), fracsand(MX,MY), SHIELD_NUMBER(MX,MY), taus50(MX,MY), H(MX,MY))
        ALLOCATE(F(MX,MY,GRAIN_ARRAY_SIZE-1), pl(MX,MY,GRAIN_ARRAY_SIZE-1))
        ALLOCATE(Fnew_avg(GRAIN_ARRAY_SIZE-1),F_diff(GRAIN_ARRAY_SIZE-1))
! Initialize
        BED_ELEVATION_CHANGE_RATE(:,:)=0.
        SEDIMENT_FLUX_RIVER(:,:)=0.
        SEDIMENT_FLUX_SLOPES(:,:)=0.
        dsgs(:,:)=0.
        taus50(:,:)=0.
        H(:,:)=0.
        D90_SIZE(:,:)=0.
        D50_SIZE(:,:)=0.
        fracsand(:,:)=0.
        SHIELD_NUMBER(:,:)=0.
        F(:,:,:)=0.
        pl(:,:,:)=0.
        ACTIVE_LAYER_THICKNESS(:,:)=0.
        ACTIVE_LAYER_THICKNESSold(:,:)=0.
        psi(:)=0.
        ds(:)=0.0
        plf(:)=0.
        Fs(:)=0.0
        OUTSEDFLUX = 81
        OUTDSGS = 82
        OUTD90 = 83
        OUTSAND = 84
        
        CALL GRAVEL_MIXTURE_INITIALIZE
    ELSE
        IS_GRAVEL_MIXTURE=.false.
    END IF  
  RETURN
  END !SUBROUTINE READ_GRAVEL_TRANSPORT_PARAMETERS
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GRAVEL_MIXTURE_INITIALIZE
    USE ERODE_GLOBALS
    USE GRAVEL_MIXTURE_GLOBALS
    Implicit none
!    INTEGER :: GD
! Gravel grain size boundaries Di and cumulative size percentages for inital gravel inflow from upstream, plff,
!    surface layer, FfI, and subsurface layer, Ffs
! Here these are converted to grain size proportions for each of the GRAIN_ARRAY_SIZE-1 phi intervals
! ds(GD) are the geometric mean grain sizes, in mm
! psi(GD) are the corresponding phi (logarithmic) grain sizes
! plf(GD) is the initial proportions of the upstream feed
! F(i,j,GD) are the initial surface layer gravel size proportions for the active layer for each of the spatial I,J cells
! Fs(GD) is the subsurface grain size proportions, assumed to be areally constant - and is assumed to also represent sediment input to
!   the channel from surrounding slopes

    Do GD=1, GRAIN_ARRAY_SIZE-1
        ds(GD)=(di(GD)*di(GD+1))**0.5
        psi(GD)=LOG(ds(GD))

        plf(GD)=(plff(GD)-plff(GD+1))/100.
        F(:,:,GD)=(FfI(GD)-FfI(GD+1))/100.
        Fs(GD)=(Ffs(GD)-Ffs(GD+1))/100.
    END DO
    !find sand fraction (ie : the fraction of witch grain size <=2mm), it is basically what is left over after accounting for the
    ! gravel size proportions
    fracsand=0.
    Do GD=1, GRAIN_ARRAY_SIZE-1
        If (di(GD)>2 .and. di(GD+1)<=2) Then
            fracsand=(FfI(GD+1)+(FfI(GD)-FfI(GD+1))/(LOG(di(GD))-LOG(di(GD+1)))*(LOG(2.)-LOG(di(GD+1))))/100.
        END IF
    END DO
    FfI(1)=100.
    dsgs(1,1)=0.
    ! here we calculate the initial geometric median grain size, dsgs, the D90 and median (D50) grain sizes, which are calculated
    !   for each i,j location.
    Do GD=1, GRAIN_ARRAY_SIZE-1
        dsgs(:,:)=dsgs(1,1)+F(:,:,GD)*psi(GD)
        If (FfI(GD)>=90 .and. FfI(GD+1)<90) Then
            D90_SIZE(:,:)=Exp(LOG(di(GD+1))+(LOG(di(GD))-LOG(di(GD+1)))/(FfI(GD)-FfI(GD+1))*(90-FfI(GD+1)))
        End If

        If (FfI(GD)>=50 .and. FfI(GD+1)<50) Then
            D50_SIZE(:,:)=Exp(LOG(di(GD+1))+(LOG(di(GD))-LOG(di(GD+1)))/(FfI(GD)-FfI(GD+1))*(50-FfI(GD+1)))
        End If
    END DO
    !here the formulation of gary parker was a bit different (ie : dsgs(:,:)=Exp(dsgs(1,1)*LOG(2))
    ! this reexpresses the median grain sizes in mm.
    dsgs(:,:)=Exp(dsgs(1,1))
END ! SUBROUTINE GRAVEL_MIXTURE_INITIALIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GRAVEL_MIXTURE_TRANSPORT   ! This is the main program for the calculation.
    USE ERODE_GLOBALS
    USE LAKE_GLOBALS
    USE CRATER_GLOBALS
    USE EOLIAN_GLOBALS
    USE SEDDEBUG_GLOBALS
    USE AREA_GLOBALS
    USE SEDROUTE_GLOBALS
    USE LAVA_GLOBALS
    USE ACCRETION_GLOBALS
    USE GRAVEL_MIXTURE_GLOBALS
    Implicit none
        !!!!!!!!write(outhist,4833)
4833    format('gravel called')    
    CALL Find_Shields_Stresses
    CALL Find_Load
    CALL Find_New_ELEVATION

END SUBROUTINE ! SUBROUTINE GRAVEL_MIXTURE_TRANSPORT!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Find_Shields_Stresses()
! Slopes, Shield numbers and depths are computed here
    USE ERODE_GLOBALS
    USE GRAVEL_MIXTURE_GLOBALS
    Implicit none
    INTEGER :: I, J, IADD, JADD
    REAL(4), DIMENSION(MX,MY) :: WATER_FLUX, SLOPE, ks, count
! calculate discharge per unit channel width, WATER_FLUX
    Do I=1, MX
        Do J=1, MY
            WATER_FLUX(I,J)=DISCHARGE(I,J)/CHANNEL_WIDTH(I,J)
        END Do
    END DO

    SLOPE(:,:)=0.
    count(:,:)=0.

    !channel gradient (SLOPE) calculation
    Do I=1, MX
        Do J=1, MY
            If (FLOW_DIRECTION(I,J)>1) then
                IADD = I+DOWNSTREAM(FLOW_DIRECTION(I,J),1)
                JADD = J+DOWNSTREAM(FLOW_DIRECTION(I,J),2)
                SLOPE(I,J)=(ELEVATION(I,J)-ELEVATION(IADD,JADD))/&
                    (sqrt((DOWNSTREAM(FLOW_DIRECTION(I,J),1)*CELL_SIZE)**2.+&
                    (DOWNSTREAM(FLOW_DIRECTION(I,J),2)*CELL_SIZE)**2.))

            !If the node is a depresion slope is zero
            Else IF (FLOW_DIRECTION(I,J)/=1) then
                SLOPE(I,J)=0.
            ELSE

            !If the node is an exit node the slope is a mean of the slopes of incoming nodes
                IADD = I+DOWNSTREAM(FLOW_DIRECTION(I,J),1)
                JADD = J+DOWNSTREAM(FLOW_DIRECTION(I,J),2)
                SLOPE(IADD,JADD)=SLOPE(I,J)+SLOPE(IADD,JADD)
                count(IADD,JADD)=count(IADD,JADD)+1.
            END IF
        END DO
    END DO

    Do I=1, MX
        Do J=1, MY
            If (FLOW_DIRECTION(IADD,JADD)==1) then
                SLOPE(IADD,JADD)=SLOPE(IADD,JADD)/count(IADD,JADD)
            END IF

            !shield stress calculation, and calculate roughness height, ks, for each location
            ks(I, J)=ROUGHNESS_FACTOR*D90_SIZE(I, J)/1000.
            SHIELD_NUMBER(I, J)=((ks(I, J)**(1./3.))*(WATER_FLUX(I, J)**2.) &
            /(MANNING_COEF**2.)/GRAVITY)**(3./10.)

            !tau50 and H aren't used in the program and come from the original Gary program
            !taus50(I, J)=SHIELD_NUMBER(I, J)

            If (SLOPE(I,J) <= 0) Then
                SHIELD_NUMBER(I, J)=0
                !taus50(I, J)=0
                !H(I, J)=0
            Else
                SHIELD_NUMBER(I, J)=SHIELD_NUMBER(I, J)*(SLOPE(I,J) ** (7./10.))*1000. /&
                            (SEDIMENT_SPECIFIC_GRAVITY*dsgs(I, J))
                !taus50(I, J)=taus50(I, J)*(SLOPE(I,J)**(7./10.))*1000./&
                !            (SEDIMENT_SPECIFIC_GRAVITY*D50_SIZE(I, J))
               ! H(I, J)=((ks(I, J)**(1./3.))*(WATER_FLUX(I, J)**2.)/(MANNING_COEF**2.)/&
               !             GRAVITY/SLOPE(I, J))**(3./10.)

            END IF
        END DO
    END DO

END SUBROUTINE ! SUBROUTINE Find_Shields_Stresses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Find_Load()
! This SUBROUTINEroutine computes the magnitude and size distribution of the load.
    USE ERODE_GLOBALS
    USE GRAVEL_MIXTURE_GLOBALS
    Implicit none
    INTEGER :: I, J
    REAL(4), DIMENSION(GRAIN_ARRAY_SIZE-1) :: b, plt
    REAL(4) :: Gsum, arg, SHEAR_VELOCITY, taussrg, GG, GGwc
! This routine calculates the sediment_flux_river (transport volume per unit time per unit channel width)
!   using either the Parker surface-based transport relationship or the Wilcock-Crowe relationship
! It also calculates the fractional distribution, pl(I,J,GD), of the sediment load in transport for each
!   grain size, GD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Parker bedload relation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    If (EQUATION_SELECTOR==1) Then
        Do I=1, MX
            Do J=1, MY
                sas=0
                Do GD=1, GRAIN_ARRAY_SIZE-1
                    sas=sas+((psi(GD)-LOG(dsgs(I,J))/LOG(2.))**2.)*F(I, J, GD)
                END DO
                sas=sas**0.5

                phisgo=SHIELD_NUMBER(I, J)/0.0386


                If (phisgo<0.1) Then
                    SEDIMENT_FLUX_RIVER(I, J)=0
                    pl(I, J, :)=1/(GRAIN_ARRAY_SIZE-1)
                Else
                    Gsum=0
                    Call find_omega2

                    Do GD=1, GRAIN_ARRAY_SIZE-1
                        arg=omega*phisgo*(ds(GD)/dsgs(I,J))**(-0.0951)
                        plt(GD)=F(I, J, GD)*GG(arg)
                        Gsum=Gsum+plt(GD)
                    END DO

                    pl(I,J,:)=plt(:)/Gsum
                    SHEAR_VELOCITY=(SEDIMENT_SPECIFIC_GRAVITY*GRAVITY*(dsgs(I, J)/1000.)*SHIELD_NUMBER(I, J))**0.5
                    SEDIMENT_FLUX_RIVER(I,J)=Gsum*0.00218*(SHEAR_VELOCITY**3.)/SEDIMENT_SPECIFIC_GRAVITY/GRAVITY

                END IF
            END DO
        END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wilcock-Crowe bedload relation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Else
        Do I=1, MX
            Do J=1, MY
                taussrg=0.021+0.015*Exp(-20.*fracsand(I, J))
                phisgo=SHIELD_NUMBER(I, J)/taussrg
                If (phisgo < 0.01) Then
                    SEDIMENT_FLUX_RIVER(I, J)=0
                    Do GD=1, GRAIN_ARRAY_SIZE-1
                        pl(I, J, GD)=1/(GRAIN_ARRAY_SIZE-1)
                    END DO
                Else
                    Gsum=0
                    Do GD=1, GRAIN_ARRAY_SIZE-1
                        b(GD)=0.67/(1+Exp(1.5-ds(GD)/dsgs(I, J)))
                        arg=phisgo*(ds(GD)/dsgs(I, J))**(-b(GD))
                        plt(GD)=F(I, J, GD)*GGwc(arg)
                        Gsum=Gsum+plt(GD)
                    END DO
                    pl(I, J,:)=plt(:)/Gsum

                    SHEAR_VELOCITY=(SEDIMENT_SPECIFIC_GRAVITY*GRAVITY*(dsgs(I, J)/1000.)*SHIELD_NUMBER(I, J))**0.5
                    SEDIMENT_FLUX_RIVER(I, J)=Gsum*(SHEAR_VELOCITY**3.)/SEDIMENT_SPECIFIC_GRAVITY/GRAVITY
                END IF
            END DO
        END DO
    END IF



    return
END SUBROUTINE ! SUBROUTINE Find_Load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Function GG(xarg)
! Computes the function G in the bedload formula of Parker
    USE ERODE_GLOBALS
    USE GRAVEL_MIXTURE_GLOBALS
    Implicit none
    REAL(4), intent(in) :: xarg
    REAL(4) :: GG
    If (xarg < 1) Then
        GG=xarg**14.2
    Else If (xarg <= 1.59) Then
        GG=Exp(14.2*(xarg-1)-9.28*(xarg-1)**2.)
    Else
        GG=5474.*(1.-0.853/xarg)**4.5
    END IF
    return
END Function ! Function GG


Function GGwc(xarg)
! Computes the function GGwc in the bedload formula of Wilcock-Crowe
    Implicit none
    REAL(4), intent(in) :: xarg
    REAL(4) :: GGwc
    If (xarg < 1.35) Then
        GGwc=0.002*xarg**7.5
    Else
        GGwc=14.*(1-0.894/(xarg**0.5))**4.5
    END IF
    Return
END Function ! Function GGwc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE find_omega1
! This SUBROUTINEroutine is part of the calculation of the bedload rate and size distribution
! using the Parker relation.
    USE ERODE_GLOBALS
    USE GRAVEL_MIXTURE_GLOBALS
    Implicit none
    Real :: omegaoak, sigoak
    REAL :: a1, a2, a3, a4
    INTEGER :: kk


    If (phisgo > STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-2) .and. phisgo <= STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-1)) Then
        omegaoak=STRAIN_CURVE_OO(STRAIN_CURVE_SIZE-2)+&
                (STRAIN_CURVE_OO(STRAIN_CURVE_SIZE-1)-STRAIN_CURVE_OO(STRAIN_CURVE_SIZE-2))/&
                (STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-1)-STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-2))*&
                (phisgo-STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-2))

        sigoak=STRAIN_CURVE_SO(STRAIN_CURVE_SIZE-2)+&
            (STRAIN_CURVE_SO(STRAIN_CURVE_SIZE-1)-STRAIN_CURVE_SO(STRAIN_CURVE_SIZE-2))/&
            (STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-1)-STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-2))*&
            (phisgo-STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-2))


    Else If (phisgo > STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-1) .and. phisgo <= STRAIN_CURVE_PO(STRAIN_CURVE_SIZE)) Then
        omegaoak=STRAIN_CURVE_OO(STRAIN_CURVE_SIZE-1)+&
            (STRAIN_CURVE_OO(STRAIN_CURVE_SIZE)-STRAIN_CURVE_OO(STRAIN_CURVE_SIZE-1))/&
            (STRAIN_CURVE_PO(STRAIN_CURVE_SIZE)-STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-1))*&
            (phisgo-STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-1))

        sigoak=STRAIN_CURVE_SO(STRAIN_CURVE_SIZE-1)+&
            (STRAIN_CURVE_SO(STRAIN_CURVE_SIZE)-STRAIN_CURVE_SO(STRAIN_CURVE_SIZE-1))/&
            (STRAIN_CURVE_PO(STRAIN_CURVE_SIZE)-STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-1))*&
            (phisgo-STRAIN_CURVE_PO(STRAIN_CURVE_SIZE-1))


    Else If (phisgo > STRAIN_CURVE_PO(STRAIN_CURVE_SIZE)) Then
        omegaoak=1.501
        sigoak=0.453

    Else If (phisgo <= STRAIN_CURVE_PO(1)) Then
        omegaoak=STRAIN_CURVE_OO(1)
        sigoak=STRAIN_CURVE_SO(1)

    Else
        L155: Do kk=1,STRAIN_CURVE_SIZE-2
            If (phisgo > STRAIN_CURVE_PO(kk) .And. phisgo <= STRAIN_CURVE_PO(kk+1)) Then
                xm1=STRAIN_CURVE_PO(kk-1)
                x=STRAIN_CURVE_PO(kk)
                xp1=STRAIN_CURVE_PO(kk+1)
                xp2=STRAIN_CURVE_PO(kk+2)
                ym1=STRAIN_CURVE_OO(kk-1)
                y=STRAIN_CURVE_OO(kk)
                yp1=STRAIN_CURVE_OO(kk+1)
                yp2=STRAIN_CURVE_OO(kk+2)

                !cubic_interpolator(omegaoak)
                a1=ym1
                a2=(y-a1)/(x-xm1)
                a3=(yp1-a1-a2*(xp1-xm1))/(xp1-xm1)/(xp1-x)
                a4=yp2-a1-a2*(xp2-xm1)-a3*(xp2-xm1)*(xp2-x)
                a4=a4/(xp2-xm1)/(xp2-x)/(xp2-xp1)

                omegaoak=a1+a2*(phisgo-xm1)+a3 *(phisgo-xm1)*(phisgo-x)+a4*(phisgo-xm1)*(phisgo-x)*(phisgo-xp1)


                ym1=STRAIN_CURVE_SO(kk-1)
                y=STRAIN_CURVE_SO(kk)
                yp1=STRAIN_CURVE_SO(kk+1)
                yp2=STRAIN_CURVE_SO(kk+2)

                !cubic_interpolator(sigoak)
                a1=ym1
                a2=(y-a1)/(x-xm1)
                a3=(yp1-a1-a2*(xp1-xm1))/(xp1-xm1)/(xp1-x)
                a4=yp2-a1-a2*(xp2-xm1)-a3*(xp2-xm1)*(xp2-x)
                a4=a4/(xp2-xm1)/(xp2-x)/(xp2-xp1)
                sigoak=a1+a2*(phisgo-xm1)+a3 *(phisgo-xm1)*(phisgo-x)+a4*(phisgo-xm1)*(phisgo-x)*(phisgo-xp1)

                exit L155
            END IF
        End Do L155
    END IF
    omega=1.+sas/sigoak*(omegaoak-1.)
    return
END SUBROUTINE ! SUBROUTINE find_omega1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FIND_OMEGA2
! This SUBROUTINEroutine is part of the calculation of the bedload rate and size distribution
! using the Parker relation.
    USE ERODE_GLOBALS
    USE GRAVEL_MIXTURE_GLOBALS
    Implicit none
    logical foundit
    REAL(4) :: omegaoak, sigoak
    REAL(4) :: a1, a2, a3, a4
    INTEGER :: kk
    If (phisgo > strain_curve_po(strain_curve_size-1)) Then
            omegaoak = strain_curve_oo(strain_curve_size-1) + (strain_curve_oo(strain_curve_size) - &
            strain_curve_oo(strain_curve_size-1)) / (strain_curve_po(strain_curve_size) - &
            strain_curve_po(strain_curve_size-1)) &
                * (phisgo - strain_curve_po(strain_curve_size-1))
            sigoak = strain_curve_so(strain_curve_size-1) + (strain_curve_so(strain_curve_size) - &
            strain_curve_so(strain_curve_size-1)) / (strain_curve_po(strain_curve_size) - &
            strain_curve_po(strain_curve_size-1)) &
                     * (phisgo - strain_curve_po(strain_curve_size-1))
        Else
            If (phisgo <= 0.7639) Then
                omegaoak = 1.011
                sigoak = 0.8157
            Else
                foundit = .False.
                kk = 1
                L300: Do
                    kk = kk + 1
                    If ((phisgo > strain_curve_po(kk)).and.(phisgo <= strain_curve_po(kk + 1))) Then
                        foundit = .True.
                        exit L300
                    endif
                    if (kk>(strain_curve_size-2)) then
                        kk=strain_curve_size-2
                        exit L300
                    endif
                    !cycle L300
                enddo L300
                !Lstrain_curve_oop Until foundit
                xm1 = strain_curve_po(kk - 1)
                x = strain_curve_po(kk)
                xp1 = strain_curve_po(kk + 1)
                xp2 = strain_curve_po(kk + 2)
                ym1 = strain_curve_oo(kk - 1)
                y = strain_curve_oo(kk)
                yp1 = strain_curve_oo(kk + 1)
                yp2 = strain_curve_oo(kk + 2)
                !cubic_interpolator(omegaoak)
                a1=ym1
                a2=(y-a1)/(x-xm1)
                a3=(yp1-a1-a2*(xp1-xm1))/(xp1-xm1)/(xp1-x)
                a4=yp2-a1-a2*(xp2-xm1)-a3*(xp2-xm1)*(xp2-x)
                a4=a4/(xp2-xm1)/(xp2-x)/(xp2-xp1)
                omegaoak=a1+a2*(phisgo-xm1)+a3 *(phisgo-xm1)*(phisgo-x)+a4*(phisgo-xm1)*(phisgo-x)*(phisgo-xp1)
!               cubic_interpolator(sigoak)
                ym1 = strain_curve_so(kk - 1)
                y = strain_curve_so(kk)
                yp1 = strain_curve_so(kk + 1)
                yp2 = strain_curve_so(kk + 2)
                a1=ym1
                a2=(y-a1)/(x-xm1)
                a3=(yp1-a1-a2*(xp1-xm1))/(xp1-xm1)/(xp1-x)
                a4=yp2-a1-a2*(xp2-xm1)-a3*(xp2-xm1)*(xp2-x)
                a4=a4/(xp2-xm1)/(xp2-x)/(xp2-xp1)
                sigoak=a1+a2*(phisgo-xm1)+a3 *(phisgo-xm1)*(phisgo-x)+a4*(phisgo-xm1)*(phisgo-x)*(phisgo-xp1)
            End If
        End If
    omega=1.+sas/sigoak*(omegaoak-1.)
    return
END SUBROUTINE ! SUBROUTINE FIND_OMEGA2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Find_New_ELEVATION()
! This SUBROUTINEroutine implements the Exner equation of sediment continuity to find the
! new bed and grain size distribution one ITERATION step later.
    USE ERODE_GLOBALS
    USE GRAVEL_MIXTURE_GLOBALS
    Implicit none
    INTEGER ::  I,J, IADD, JADD, KK
    REAL(4) :: Sum, netchange,xsum, divide_term
    REAL(4), DIMENSION(MX,MY) :: SEDIMENT_FLUX_IN_TOT, SEDIMENT_FLUX_OUT_TOT
    REAL(4), DIMENSION(MX,MY,GRAIN_ARRAY_SIZE-1) :: qGD1dev_IN_TOT, qGD1dev_OUT_TOT
    REAL(4), DIMENSION(GRAIN_ARRAY_SIZE-1) :: qGD1dev, qGD2dev,aterm,Fprime
    REAL(4), DIMENSION(MX,MY) :: ELEVATIONnew, INFLUX_TOT, Effective_Width
    REAL(4), DIMENSION(MX,MY, GRAIN_ARRAY_SIZE-1) :: Fnew, FIexc
    REAL(4), DIMENSION(GRAIN_ARRAY_SIZE) :: Fft


    REAL(4) :: b_elavg,b_sfi,b_sfo,b_ecr,b_alt,b_ew,b_f(grain_array_size-1),b_d90,b_dsgs

    ! SEDIMENT_FLUX_IN_TOT(i,j) and SEDIMENT_FLUX_OUT_TOT(i,j) bookeep total sediment transport per unit time entering and leaving
    !   a simulation cell (and thus utilizes SEDIMENT_RIVER_FLUX(i,j)*CHANNEL_WIDTH(I,J))

    SEDIMENT_FLUX_IN_TOT(:,:)=0.
    SEDIMENT_FLUX_OUT_TOT(:,:)=0.

! In this formulation we assume an effective channel width, effective_width(i,j), which depends upon the fraction of material in
!    the subsurface that is assumed to be in the size range carried as suspended load, which is assumed to be transported
!    outside the boundaries of the simulation domain without redeposition.
! The effective width varies from the width the channel (channel_width) to the width of the simulation cell, cell_size, as the
!    parameter slope_mud_fraction goes from 1.0 to 0.0.  This is a bit of a kludge, because even whtn the slope_mud_fraction is
!    1.0 it is assumed that there is sufficient delivery of gravel from the subsurface and surroundings to maintain a gravel layer
!    when the channel bed lowers.  When the slope_mud_fraction is 0.0 the program assumes that the channel must transport (and abrade)
!    sediment from the entire cell.  Gravel transport rates, however, are determined by the channel discharge and channel width in Find_Load().
! The effective_width also accounts for sediment porosity.
!
! As presently formulated the gravel transport module only IMPLICITLY accounts for delivery of sediment from slopes and the slopes are not
!    explicitly modeled.  In particular, in the MARSSIM.PRM file the FLUVIAL_AND_SLOPE switch should be 1, but DO_FLUVIAL_DETACHMENT should be 0
!    and DO_SEDIMENT_TRANSPORT should be 0.  It also assumes a fixed, horizontal lower (Y) boundary, so that the last line under BOUNDARY CONDITIONS
!    in MARSSIM.PRM should be 1 1 0 0 0
! Parameters for the gravel transport are in the GRAVEL_MIXTURE.PRM file. The first number in that file (DO_GRAVEL_MIXTURE) determines whether
!    this routine will be used (1 == Yes and 0 == No).
!
! The program also uses the Parker (1991) model of "Selective sorting and abrasion of river gravel. I: Theory", J. Hydr. Eng., 131-149.  WIth formulas
!    as presented in Cui,Y. (2007), Water Resourc. Res, 43, W10436, doi:10.1029/2006WR005330.



    BED_ELEVATION_CHANGE_RATE(:,:)=0.
    ACTIVE_LAYER_THICKNESS(:,:)=ACTIVE_LAYER_FACTOR*D90_SIZE(:,:)/1000.0

    !for each node the incoming sediment flux and the outflux is calculated
    Do I=1, MX
        Do J=1, MY
             IF (WIDTH_RELATIVE) THEN
                effective_width(i,j)=(slope_mud_fraction*(amin1(0.0,channel_width(i,j)-cell_size))+ &
                cell_size)*(1.0-sediment_porosity)
             ELSE
                effective_width(i,j)=(1.0-slope_mud_fraction)*cell_size*(1.0-sediment_porosity)
             ENDIF
             ew_avg=ew_avg+effective_width(i,j)
            If (FLOW_DIRECTION(I,J)/=1) then
                !if the node is a depresion there is no out flux of sediment, everything is deposited
                IF (FLOW_DIRECTION(I,J)<=0) then
                    SEDIMENT_FLUX_OUT_TOT(I,J)=0.
                Else

                    IADD = I+DOWNSTREAM(FLOW_DIRECTION(I,J),1)
                    JADD = J+DOWNSTREAM(FLOW_DIRECTION(I,J),2)

                    SEDIMENT_FLUX_OUT_TOT(I,J)=(SEDIMENT_FLUX_RIVER(I, J))*CHANNEL_WIDTH(I,J)
                    SEDIMENT_FLUX_IN_TOT(IADD,JADD)=SEDIMENT_FLUX_IN_TOT(IADD,JADD)+&
                        SEDIMENT_FLUX_OUT_TOT(I,J)

                End If
            End If
        ENd DO
    ENd DO

    ! Here the mass balance is converted to a rate of bed elevation change BED_ELEVATION_CHANGE_RATE
    ! Mass lost by abrasion (assumed to become silt carried in suspension) is calculated as the abrade_term and is proportional
    ! to the SEDIMENT_FLUX_RIVER(i,j) divided by the effective width(i,j) and the abrasion contributes to overall elevation change
    ! Note that a positive elevation change rate corresponds to depositon.

    Do I=1, MX
        Do J=1, MY

            IF (FLOW_DIRECTION(I,J)/=1) then
                netchange=(SEDIMENT_FLUX_OUT_TOT(I,J)-SEDIMENT_FLUX_IN_TOT(I,J))/(cell_size*effective_width(i,j))
                net1_avg=net1_avg+netchange
                abrade_term=6.0*ABRASION_COEFFICIENT*SEDIMENT_FLUX_RIVER(I,J)*CHANNEL_WIDTH(i,j)/effective_width(i,j)
                netchange=netchange+abrade_term
                net2_avg=net2_avg+netchange
                aterm_avg=aterm_avg+abrade_term

                    BED_ELEVATION_CHANGE_RATE(I,J)=netchange

            !if the node is an output node we asume that no change in elevation occur
            else
                netchange=0.0
                BED_ELEVATION_CHANGE_RATE(I,J)=0.
            END IF
! FIexc is the size distribution of sediment either added to the active layer from the subsurface [Fs] (or implicitly laterally from slopes)
!   if the channel is eroding (B_E_C_R>0) or a mixture of the sediment in transport [pl] and the active layer [F] if deposition is occurring (B_E_C_R<0)
            If (BED_ELEVATION_CHANGE_RATE(I,J) > 0) Then
                FIexc(I,J,:)=Fs(:)
            Else
                FIexc(I,J,:)=AGGRADATION_COEF*F(I,J,:)+(1-AGGRADATION_COEF)*pl(I,J,:)
            END IF
            avg_change=avg_change+BED_ELEVATION_CHANGE_RATE(I,J)
        ENd DO
    ENd DO


    !This part aims to calculates the new grain size distribution
    ! the scheme is the same than for the bed elevation, calculating an output and an input distribution

    ! gGD1dev_in_tot(i,j,GD) and gGD1dev_out_tot(i,j,GD) are the grain-size-specific input and output flux of sediment from the cell
    INFLUX_TOT(:,:)=0.
    qGD1dev_IN_TOT(:,:,:)=0.
    qGD1dev_OUT_TOT(:,:,:)=0.

    Do I=1, MX
        Do J=1, MY
             el_avg=el_avg+elevation(i,j)
             sf_avg=sf_avg+sediment_flux_river(i,j)

            If (FLOW_DIRECTION(I,J)/=1) Then
                If (FLOW_DIRECTION(I,J)<=0) then
                   qGD1dev_OUT_TOT(I,J,:)=0.
                else
                  IADD = I+DOWNSTREAM(FLOW_DIRECTION(I,J),1)
                  JADD = J+DOWNSTREAM(FLOW_DIRECTION(I,J),2)

                qGD1dev_OUT_TOT(I,J,:)=SEDIMENT_FLUX_RIVER(I,J)*CHANNEL_WIDTH(I,J)*pl(I,J,:)

                qGD1dev_IN_TOT(IADD,JADD,:)=qGD1dev_IN_TOT(IADD,JADD,:)+qGD1dev_OUT_TOT(I,J,:)

                End If
            End If

            !If the node is an output node we assume that The new grain size distribution will be the mean of all the input sediment flux ponderate by the value of the flux
            If (FLOW_DIRECTION(IADD,JADD)==1) Then
                qGD1dev_IN_TOT(IADD,JADD,:)=qGD1dev_IN_TOT(IADD,JADD,:)+F(I,J,:)* &
                SEDIMENT_FLUX_RIVER(I, J)*channel_width(i,j)

            ENd IF
        END DO
    END DO
! Fprime(GD) is the areal fraction of the given grain size on the channel bed -- this is used for abrasion calculations
! qGD1dev(GD) is the grain-size-specific downstream flux divergence for the cell
! qGD2dev(GD) is grain-size-specific rate of addition or removal of sediment to/from the active layer
! Here we also calculate the net elevation change, taking into account the time increment (specified in years), the elevation change
! measured in meters per second, how many sub-iterations of transport are done per year, and the fraction of the year with flow events

    Do I=1, MX
        Do J=1, MY
            xsum=0.0
            Fprime=0.
            aterm=0.
            do GD=1,GRAIN_ARRAY_SIZE-1
                xsum=xsum+F(I,J,GD)/ds(GD)**0.5
            enddo
            do GD=1,GRAIN_ARRAY_SIZE-1
                Fprime(GD)=F(I,J,GD)/ds(GD)**0.5/xsum
            enddo
            qGD2dev(:)=FIexc(I,J,:)*BED_ELEVATION_CHANGE_RATE(I,J)
            do GD=1,GRAIN_ARRAY_SIZE-1
                  gd2_avg=gd2_avg+qGD2dev(GD)
            enddo

            If (FLOW_DIRECTION(I,J)/=1) then
                qGD1dev(:)=(qGD1dev_OUT_TOT(I,J,:)-qGD1dev_IN_TOT(I,J,:))/ (cell_size*effective_width(i,j))
            else
                qGD1dev(:)=0.0
            END IF
             do GD=1,GRAIN_ARRAY_SIZE-1
                  gd1_avg=gd1_avg+qGD1dev(GD)
            enddo
            ELEVATIONnew(I, J)=ELEVATION(I,J)+TIME_INCREMENT*SECONDS_PER_YEAR/&
                STEP_BY_MARSSIM_ITERATION*(-BED_ELEVATION_CHANGE_RATE(I,J)*INTERMITENCY)
! here we calculate grain-size-specific abrasion as presented in Parker(1991) and Cui(2007)
! It involves two terms, both weighted by the ABRASION_COEFFICIENT and the sediment flux rate
! One term is for abrasion abrasion of the load and the bed and loss to silt, and the other accounts for the
!   resulting grain-size fining
                      !
            do GD=1,GRAIN_ARRAY_SIZE-1
                a1term=3.0*(pl(I,J,GD)+Fprime(GD))
                a1_avg=a1_avg+a1term
                if (GD==1) then
                    a2term=(pl(I,J,GD))+(Fprime(GD))/(log(2.)*(PSI(2)-PSI(1)))
                else
                a2term=(pl(I,J,GD)-pl(I,J,GD-1))+(Fprime(GD)-Fprime(GD-1))/(log(2.)*(PSI(GD)-PSI(GD-1)))
                endif
                a2_avg=a2_avg+a2term
                aterm(GD)=SEDIMENT_FLUX_RIVER(I,J)*ABRASION_COEFFICIENT*CHANNEL_WIDTH(I,J)* &
                (a1term+a2term)/effective_width(i,j)
                agterm_avg=agterm_avg+aterm(GD)
            enddo
! Here we calculate the new grain size distribution of the gravel in the channel bed

        If (FLOW_DIRECTION(I,J)/=1) then
            Fnew(I,J,:)=F(I,J,:)+TIME_INCREMENT*SECONDS_PER_YEAR/&
                STEP_BY_MARSSIM_ITERATION*INTERMITENCY*(-qGD1dev(:)+qGD2dev(:)-aterm(:))/&
                    ACTIVE_LAYER_THICKNESS(I, J)

        Else

                Fnew(I,J,:)=F(I,J,:)
           !
        END IF
! This adds the second term into the grain size change
            If (ITERATION>0) Then
                if (flow_direction(i,j)/=1) then
                Fnew(I,J,:)=Fnew(I,J,:)+(FIexc(I,J,:)-F(I,J,:))/ACTIVE_LAYER_THICKNESS(I, J)*&
                        (ACTIVE_LAYER_THICKNESS(I, J)-ACTIVE_LAYER_THICKNESSold(I, J))
                Fnew_avg(:)=Fnew_avg(:)+Fnew(I,J,:)
                F_diff(:)=F_diff(:)+Fnew(I,J,:)-F(I,J,:)
                endif
            END IF
        END DO
    END DO
! Just some diagnistics
    IF (GRAVEL_ITERATION==STEP_BY_MARSSIM_ITERATION) THEN
    a1_avg=a1_avg/((GRAIN_ARRAY_SIZE-1)*MX*STEP_BY_MARSSIM_ITERATION*(MY-1))
    a2_avg=a2_avg/((GRAIN_ARRAY_SIZE-1)*MX*STEP_BY_MARSSIM_ITERATION*(MY-1))
    agterm_avg=agterm_avg/((GRAIN_ARRAY_SIZE-1)*MX*STEP_BY_MARSSIM_ITERATION*(MY-1))
    gd1_avg=gd1_avg/((GRAIN_ARRAY_SIZE-1)*MX*STEP_BY_MARSSIM_ITERATION*(MY-1))
    gd2_avg=gd2_avg/((GRAIN_ARRAY_SIZE-1)*MX*STEP_BY_MARSSIM_ITERATION*(MY-1))
    Fnew_avg(:)=Fnew_avg(:)/(MX*STEP_BY_MARSSIM_ITERATION*(MY-1))
    F_diff(:)=F_diff(:)/(MX*STEP_BY_MARSSIM_ITERATION*(MY-1))
    aterm_avg=aterm_avg/(MX*STEP_BY_MARSSIM_ITERATION*(MY-1))
    net1_avg=net1_avg/(MX*STEP_BY_MARSSIM_ITERATION*(MY-1))
    net2_avg=net2_avg/(MX*STEP_BY_MARSSIM_ITERATION*(MY-1))
    avg_change=avg_change/(MX*STEP_BY_MARSSIM_ITERATION*(MY-1))
    el_avg=el_avg/(MX*STEP_BY_MARSSIM_ITERATION*(MY))
    sf_avg=sf_avg/(MX*STEP_BY_MARSSIM_ITERATION*(MY))
    ew_avg=ew_avg/(MX*STEP_BY_MARSSIM_ITERATION*(MY))
    if (PRINT_GRAVEL_DETAILS.OR.(MOD(ITERATION,IMAGE_OUTPUT_INTERVAL) == 0)) THEN
    write(outhist,433) aterm_avg,a1_avg,a2_avg,gd1_avg,gd2_avg,agterm_avg,avg_change, &
    net1_avg,net2_avg,el_avg,sf_avg,ew_avg
433 format('aterm_avg=',g13.6,' a1_avg=',g13.6,' a2_avg=',g13.6,' gd1_avg=',g13.6, &
' gd2_avg=',g13.6,/,'  agterm_avg=',g13.6,' avg_change=',g13.6,/, &
           'net1_avg=',g13.6,' net2_avg=',g13.6,/,'el_avg=',g13.6,' sf_avg=',g13.6,' ew_avg',g13.6)
    write(outhist,434)
434 format('grain size, fractional percentage, change magnitude')
    DO GD=1,GRAIN_ARRAY_SIZE-1
        write(outhist,435) ds(GD),Fnew_avg(GD),F_Diff(GD)
435     format(3(' ',G13.6))
    ENDDO
    ENDIF
    ENDIF
! Actually change the elevations
    ELEVATION(:,:)=ELEVATIONnew(:,:)
! Keep a record of the active layer thickness to be able to calculate the rate of change of active layer thickness in the next iteration
    ACTIVE_LAYER_THICKNESSold(:,:)=ACTIVE_LAYER_THICKNESS(:,:)
! Here the grain size distribution is changed based upon the above calculations.  Also dsgs,D90, and D50 are calculated.  The active layer
!   grain size distribution is renormalized and the sand fraction is calculated
    Do I=1, MX
        Do J=1, MY
            Sum=0

            F(I, J,:)=Fnew(I,J,:)
            Do GD=1, GRAIN_ARRAY_SIZE-1
                If (F(I, J, GD) < 0) Then
                    F(I, J, GD)=0
                END IF
                Sum=Sum+F(I, J, GD)
            END DO
            F(I, J,:)=F(I, J,:)/Sum

            dsgs(I, J)=0
            Fft(1)=100
            Do GD=1, GRAIN_ARRAY_SIZE-1
                dsgs(I, J)=dsgs(I, J)+F(I, J, GD)*psi(GD)
                Fft(GD+1)=Fft(GD)-F(I, J, GD)*100.
                If (Fft(GD)>=90. .and. Fft(GD+1)<90.) Then
                    D90_SIZE(I, J)=Exp(LOG(di(GD+1))+(LOG(di(GD))-LOG(di(GD+1)))/(Fft(GD)-Fft(GD+1))*(90.-Fft(GD+1)))
                End If
                If (Fft(GD)>=50. .and. Fft(GD+1)<50.) Then
                    D50_SIZE(I, J)=Exp(LOG(di(GD+1))+(LOG(di(GD))-LOG(di(GD+1)))/(Fft(GD)-Fft(GD+1))*(50.-Fft(GD+1)))
                End If
            END DO

            dsgs(I, J)=Exp(dsgs(I, J))

            fracsand(I, J)=0.
            Do GD=1, GRAIN_ARRAY_SIZE-1
                If (di(GD)>2. .and. di(GD+1)<=2.) Then
                    fracsand(I, J)=(Fft(GD+1)+(Fft(GD)-Fft(GD+1))/(LOG(di(GD))-LOG(di(GD+1))) &
                    *(LOG(2.)-LOG(di(GD+1))))/100.
                END IF
            END DO
        END DO
    END DO
    ! more diagnostics
     IF (GRAVEL_ITERATION==STEP_BY_MARSSIM_ITERATION) THEN
    b_elavg=0.0
    b_sfi=0.0
    b_sfo=0.0
    b_ecr=0.0
    b_alt=0.0
    b_ew=0.0
    b_d90=0.0
    b_dsgs=0.0
    b_f=0.0
    do i=1,mx
        b_elavg=b_elavg+elevation(i,my)
        b_sfi=b_sfi+sediment_flux_in_tot(i,my)
        b_sfo=b_sfo+sediment_flux_out_tot(i,my)
        b_ecr=b_ecr+bed_elevation_change_rate(i,my)
        b_alt=b_alt+active_layer_thickness(i,my)
        b_ew=b_ew+effective_width(i,my)
        b_d90=b_d90+d90_size(i,my)
        b_dsgs=b_dsgs+dsgs(i,my)
        do gd=1,GRAIN_ARRAY_SIZE-1
            b_f(gd)=b_f(gd)+f(i,my,gd)
        enddo
    enddo
    b_elavg=b_elavg/mx
    b_sfi=b_sfi/mx
    b_sfo=b_sfo/mx
    b_ecr=b_ecr/mx
    b_alt=b_alt/mx
    b_ew=b_ew/mx
    b_d90=b_d90/mx
    b_dsgs=b_dsgs/mx
    b_f(:)=b_f(:)/mx
    write(outhist,7131)
7131 format(/,'MY')
    write(outhist,7132) b_elavg,b_sfi,b_sfo,b_ecr,b_alt,b_ew,b_d90,b_dsgs
7132 format('b_elavg=',g13.6,' b_sfi=',g13.6,' b_sfo=',g13.6,' b_ecr=',g13.6,/,' b_alt=',g13.6, &
' b_ew=',g13.6,' b_d90=',g13.6,' b_dsgs=',g13.6)
    write(outhist,7133)
7133 format('grain size frequencies')
     do gd=1,grain_array_size-1
         write(outhist,7134) b_f(GD)
7134 format('     ',g13.6)
     enddo
    b_elavg=0.0
    b_sfi=0.0
    b_sfo=0.0
    b_ecr=0.0
    b_alt=0.0
    b_ew=0.0
    b_d90=0.0
    b_dsgs=0.0
    b_f=0.0
    do i=1,mx
        b_elavg=b_elavg+elevation(i,my-1)
        b_sfi=b_sfi+sediment_flux_in_tot(i,my-1)
        b_sfo=b_sfo+sediment_flux_out_tot(i,my-1)
        b_ecr=b_ecr+bed_elevation_change_rate(i,my-1)
        b_alt=b_alt+active_layer_thickness(i,my-1)
        b_ew=b_ew+effective_width(i,my-1)
        b_d90=b_d90+d90_size(i,my-1)
        b_dsgs=b_dsgs+dsgs(i,my-1)
        do gd=1,GRAIN_ARRAY_SIZE-1
            b_f(gd)=b_f(gd)+f(i,my-1,gd)
        enddo
    enddo
    b_elavg=b_elavg/mx
    b_sfi=b_sfi/mx
    b_sfo=b_sfo/mx
    b_ecr=b_ecr/mx
    b_alt=b_alt/mx
    b_ew=b_ew/mx
    b_d90=b_d90/mx
    b_dsgs=b_dsgs/mx
    b_f(:)=b_f(:)/mx
    write(outhist,7135)
7135 format(/,'MY-1')
    write(outhist,7132) b_elavg,b_sfi,b_sfo,b_ecr,b_alt,b_ew,b_d90,b_dsgs
!7132 format('b_elavg=',g13.6,' b_sfi=',g13.6,' b_sfo=',g13.6,' b_ecr=',g13.6,/,' b_alt=',g13.6,' b_ew=',g13.6,' b_d90=',g13.6,' b_dsgs=',g13.6)
    write(outhist,7133)
!7133 format('grain size frequencies')
     do gd=1,grain_array_size-1
         write(outhist,7134) b_f(GD)
!7134 format('     ',g13.6)
     enddo
    endif
    return
END SUBROUTINE ! SUBROUTINE Find_New_ELEVATION
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE write_output_gravel
! report on simulation progress
    USE ERODE_GLOBALS
    USE GRAVEL_MIXTURE_GLOBALS
    Implicit none
    INTEGER :: I, J
    REAL(4) :: outarray(mx,my)
    Character(Len=10) :: time_name
        OPEN(OUTSEDFLUX,FILE='SEDFLUX.DAT',action='write', position='append')
        OPEN(OUTDSGS,file='DSGS.DAT',action='write', position='append')
        open(OUTD90,file='D90.DAT',action='write', position='append')
        OPEN(OUTSAND,file='SANDFRACTION.DAT',action='write', position='append')

    do j=1,my
        do i=1,mx
            outarray(i,j)=sediment_flux_river(i,j)*channel_width(i,j)
        enddo
    enddo
    CALL REAL_WRITE(outarray,ITERATION,OUTSEDFLUX)
    CALL REAL_WRITE(D50_SIZE,ITERATION,OUTDSGS)
    CALL REAL_WRITE(D90_SIZE,ITERATION,OUTD90)
    CALL REAL_WRITE(FRACSAND,ITERATION,OUTSAND)
    close(OUTSEDFLUX)
    close(OUTDSGS)
    close(OUTD90)
    close(OUTSAND)
    RETURN
   write(time_name,'(I5.5)') ITERATION
  elev_WRITE=Trim(adjustl('Results/ELEVATION/step'))//Trim(adjustl(time_name))//Trim(adjustl('.csv'))
 Dsg_WRITE=Trim(adjustl('Results/DSGS/step'))//Trim(adjustl(time_name))//Trim(adjustl('.csv'))
    SED_FLUX_WRITE=Trim(adjustl('Results/SEDIMENT_FLUX/step'))//Trim(adjustl(time_name))//Trim(adjustl('.csv'))


    Open (UNIT=9, file=elev_WRITE, form='formatted')
    Open (UNIT=10, file=Dsg_WRITE, form='formatted')
    if (ITERATION/=0) then
        Open (UNIT=11, file=SED_FLUX_WRITE, form='formatted')
    End if

    Do I=1,MX
        Do J=1,MY
            Write(9,*) I, J, ELEVATION(I, J)
            Write(10,*) I, J, D50_SIZE(I,J)
            if (ITERATION/=0) then
                Write(11,*) I, J, SEDIMENT_FLUX_RIVER(I, J)
            End if
        End Do
        Write(9,*)
        Write(10,*)
        if (ITERATION/=0) then
            Write(11,*) I, J, SEDIMENT_FLUX_RIVER(I, J)
        End if
    End Do

    close (unit=9)
    close (unit=10)
    if (ITERATION/=0) then
        close (unit=11)
    End if


end SUBROUTINE ! SUBROUTINE write_output_gravel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE REAL_WRITE(INARRAY,INUMB,IUNIT)
USE ERODE_GLOBALS
USE GRAVEL_MIXTURE_GLOBALS
INTEGER I,J,IUNIT,INUMB
REAL(4) :: INARRAY(MX,MY)
WRITE(IUNIT,100) MX,MY,ITERATION
100 FORMAT(3(I9,' '))
DO I=1,MX
    DO J=1,MY
        WRITE(IUNIT,200) INARRAY(I,J)
200 FORMAT(G13.6)
    ENDDO
ENDDO
RETURN
END SUBROUTINE REAL_WRITE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SUBROUTINE WRITE_SEDIMENT_FLUX()
        USE ERODE_GLOBALS
        USE GRAVEL_MIXTURE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J, ITEMP
        LOGICAL :: CONNECTBASIN
        INTEGER(4) :: IFILE
        INTEGER (1) :: BVAL
        CHARACTER :: ACHAR
        REAL(4) :: MAXELEV,MINELEV,AAAA,BBBB
        REAL(4) :: ZMAX,ZMIN,QLOCAL
        IFILE=53
        CONNECTBASIN=.TRUE.
        RFILENAME='SEDFLX'//RNUM3//RNUM2//RNUM1//RSUFFIX
        OPEN(IFILE,FILE=RFILENAME)
        MAXELEV=-1.0E+25
        MINELEV=-MAXELEV
        DO J=1,MY
            DO I=1,MX
                IF (SEDIMENT_FLUX_RIVER(I,J) > 0.0) THEN
                    QLOCAL=LOG10(SEDIMENT_FLUX_RIVER(I,J)*CHANNEL_WIDTH(I,J))
                ENDIF
                IF (QLOCAL > MAXELEV) MAXELEV=QLOCAL
                IF (QLOCAL < MINELEV) MINELEV=QLOCAL
            ENDDO
        ENDDO
        ZMAX=MAXELEV
        ZMIN=MINELEV
        WRITE(OUTHIST,200) ZMAX,ZMIN
        WRITE(*,200) ZMAX,ZMIN
        200   FORMAT('QS ZMAX=',G12.5,' ZMIN=',G12.5)
        WRITE(OUTHIST,200) ZMAX,ZMIN
        AAAA=256.0/(MAXELEV-MINELEV)
        BBBB=-AAAA*MINELEV
        DO J=1,MY
            DO I=1,MX
                IF (SEDIMENT_FLUX_RIVER(I,J) > 0.0) THEN
                    QLOCAL=LOG10(SEDIMENT_FLUX_RIVER(I,J)*CHANNEL_WIDTH(I,J))
                ELSE
                    QLOCAL=LOG10(ZMIN)
                ENDIF
                ITEMP=INT(AAAA*QLOCAL+BBBB)
                IF (ITEMP < 0) ITEMP=0
                IF (ITEMP > 255) ITEMP=255
                BVAL=ITEMP
                ACHAR=CHAR(BVAL)
                WRITE(IFILE,130) ACHAR
                130   FORMAT(A1,$)
            ENDDO
        ENDDO
        CLOSE(53)
        RETURN
    END ! SUBROUTINE WRITE_SEDIMENT_FLUX
