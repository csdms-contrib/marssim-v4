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
    SUBROUTINE CHANNEL_PROPERTIES
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,K,L,NSTEP,MAXLENGTH,NLENGTH
        INTEGER :: I_OLD,J_OLD,I_NEW,J_NEW,LOCAL_DIRECTION,ISUB,JSTART,JEND,JSTEP
        INTEGER :: ILONGEST,JLONGEST,IDIST
        REAL(4) :: A1TERM,A2TERM,XAREA,QSSS,TRANSSTAGE
        REAL(4) :: BEDLOAD_FLUX,SEDSATFACTOR,ROUTEDIV
        REAL(4) :: QE,QW,QN,QS,QNE,QNW,QSE,QSW,QCOMP
        REAL(4) :: LOCAL_INVERSE_RESISTANCE,LOCAL_CRITICAL_SHEAR,LOCAL_BEDROCK_ERODIBILITY
        REAL(4) :: LOCAL_REGOLITH_CRITICAL_SHEAR,LOCAL_REGOLITH_ERODIBILITY,NOMINAL_CRITICAL_SHEAR
        REAL(4) :: TEMP,AFACTOR,GFACTOR,SHEAR_STRESS_FACTOR
        INTEGER :: IE,IW,JN,JS
        INTEGER :: KSED,KROCK,OUTCHPROP
        !   ******************************************************************************
        !    This subroutine picks regularly-spaced locations within the simulation matrix
        !       and follows the path downstream to the channel end in a depression or simulation boundary.
        !       a number of local properties are written for each location along the stream path.
        !   ******************************************************************************
        !   CALLS: LOCAL_VALUES,SEDIMENT_TRANSPORT_FLUX,
        NSTEP=64 !   picks every 64th point - can change for higher or lower resolution
        OUTCHPROP=88
        RETURN
        IF (IS_Y_PERIODIC) THEN
            JSTART=1
            JEND=MY
            JSTEP=NSTEP
        ELSE
            JSTART=1
            JEND=1
            JSTEP=1
        ENDIF
        OPEN(OUTCHPROP,FILE='CHPROP.PRN',POSITION='APPEND')
        IF (COMPLETE_RUNOFF.AND.(.NOT.IS_Y_PERIODIC)) THEN
            J=JSTART
            DO I=1,MX,NSTEP
                MAXLENGTH=0
                DO K=1,NSTEP
                    DO L=1,NSTEP
                        I_OLD=MIN(MX,I+K-1)
                        J_OLD=MIN(MY,J+L-1)
                        NLENGTH=1
                        L1110: DO
                            LOCAL_DIRECTION=FLOW_DIRECTION(I_OLD,J_OLD)
                            IF (LOCAL_DIRECTION < 2) EXIT
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
                            IF (FLOW_DIRECTION(I_NEW,J_NEW) < 2) EXIT
                            I_OLD=I_NEW
                            J_OLD=J_NEW
                            NLENGTH=NLENGTH+1
                        ENDDO L1110
                        IF (NLENGTH > MAXLENGTH) THEN
                            ILONGEST=I+K-1
                            JLONGEST=J+L-1
                            MAXLENGTH=NLENGTH
                        ENDIF
                    ENDDO
                ENDDO
                I_OLD=ILONGEST
                J_OLD=JLONGEST
                IDIST=0
                L1210: DO
                    XAREA=DISCHARGE(I_OLD,J_OLD)
                    A1TERM=(XAREA)**SEDIMENT_1_EXPONENT
                    A2TERM=(XAREA)**SEDIMENT_2_EXPONENT
                    CALL LOCAL_VALUES(I_OLD,J_OLD,LOCAL_INVERSE_RESISTANCE,NOMINAL_CRITICAL_SHEAR, &
                    LOCAL_CRITICAL_SHEAR, LOCAL_BEDROCK_ERODIBILITY, LOCAL_REGOLITH_CRITICAL_SHEAR, &
                    LOCAL_REGOLITH_ERODIBILITY)
                    IF (DISCHARGE(I_OLD,J_OLD) > 0.0) THEN
                        AFACTOR = DISCHARGE(I_OLD,J_OLD)**BEDROCK_DISCHARGE_EXPONENT
                    ELSE
                        AFACTOR=0.0
                    ENDIF
                    GFACTOR = D8_GRADIENT(I_OLD,J_OLD)**BEDROCK_GRADIENT_EXPONENT
                    TEMP = AFACTOR*GFACTOR
                    SHEAR_STRESS_FACTOR=DISCHARGE_COEFFICIENT*TEMP
                    CALL SEDIMENT_TRANSPORT_FLUX(QSSS,TRANSSTAGE,D8_GRADIENT(I_OLD,J_OLD),A1TERM,A2TERM)
                    QSSS=QSSS/CHANNEL_WIDTH(I_OLD,J_OLD)
                    BEDLOAD_FLUX=AMAX1(0.0,SEDIMENT_YIELD(I_OLD,J_OLD)/CHANNEL_WIDTH(I_OLD,J_OLD))
                    IF (QSSS > 0.0) THEN
                        SEDSATFACTOR=AMAX1(0.0,(1.0-(BEDLOAD_FLUX/QSSS)))
                    ELSE
                        SEDSATFACTOR=0.0
                    ENDIF
                    BEDLOAD_FLUX=BEDLOAD_FLUX*BEDLOAD_FRACTION
                    IF (IS_SEDIMENT_COVERED(I_OLD,J_OLD)) THEN
                        KSED=1
                    ELSE
                        KSED=0
                    ENDIF
                    IF (IS_ROCK_SURFACE(I_OLD,J_OLD)) THEN
                        KROCK=1
                    ELSE
                        KROCK=0
                    ENDIF
                    IF (SUBMERGED(I_OLD,J_OLD)) THEN
                        ISUB=1
                    ELSE
                        ISUB=0
                    ENDIF
                    IF (CFN(I_OLD,J_OLD) < 1.0E+24) THEN
                        ROUTEDIV=CFN(I_OLD,J_OLD)
                    ELSE
                        ROUTEDIV=0.0
                    ENDIF
                    LOCAL_DIRECTION=FLOW_DIRECTION(I_OLD,J_OLD)
                    WRITE(OUTCHPROP,300) ILONGEST,JLONGEST,I_OLD,J_OLD,LOCAL_DIRECTION,IDIST &
                    ,KSED,KROCK &
                    ,ELEVATION(I_OLD,J_OLD),D8_GRADIENT(I_OLD,J_OLD),EQUILIBRIUM_GRADIENT(I_OLD,J_OLD) &
                    ,ERODE_CHANNEL(I_OLD,J_OLD),ERODE_SLOPE(I_OLD,J_OLD),DISCHARGE(I_OLD,J_OLD),BEDLOAD_FLUX,QSSS &
                    ,TRANSSTAGE,SEDSATFACTOR  &
                    ,ISUB,IDO(I_OLD,J_OLD),ROUTEDIV,SEDIMENT_BASE(I_OLD,J_OLD),ERODE_REGOLITH_CHANNEL(I_OLD,J_OLD) &
                    ,CFNE(I_OLD,J_OLD),CFNW(I_OLD,J_OLD),SHEAR_STRESS_FACTOR,LOCAL_CRITICAL_SHEAR &
                    ,SEDIMENT_BASE(I_OLD,J_OLD),CFN(I_OLD,J_OLD),LAKE_OUTLET_ELEVATION(BASIN_NUMBER(I_OLD,J_OLD)) &
                    ,I_OUTFLOW(BASIN_NUMBER(I_OLD,J_OLD)),J_OUTFLOW(BASIN_NUMBER(I_OLD,J_OLD))
                    IDIST=IDIST+1
                    QCOMP=0.0
                    IE=I_OLD+1
                    IF (IS_X_PERIODIC) THEN
                        IF (IE > MX) IE=1
                    ENDIF
                    IW=I_OLD-1
                    IF (IS_X_PERIODIC) THEN
                        IF (IW < 1) IW=MX
                    ENDIF
                    JN=J_OLD-1
                    JS=J_OLD+1
                    IF (JN < 1) THEN
                    ELSE
                        IF (JS >= MY) THEN
                            EXIT L1210
                        ENDIF
                        QN=DISCHARGE(I_OLD,JN)
                        IF (QN > QCOMP)THEN
                            I_NEW=I_OLD
                            J_NEW=JN
                            QCOMP=QN
                        ENDIF
                        QNW=DISCHARGE(IW,JN)
                        IF (QNW > QCOMP)THEN
                            I_NEW=IW
                            J_NEW=JN
                            QCOMP=QNW
                        ENDIF
                        QNE=DISCHARGE(IE,JN)
                        IF (QNE > QCOMP)THEN
                            I_NEW=IE
                            J_NEW=JN
                            QCOMP=QNE
                        ENDIF
                    ENDIF
                    QE=DISCHARGE(IE,J_OLD)
                    IF (QE > QCOMP) THEN
                        I_NEW=IE
                        J_NEW=J_OLD
                        QCOMP=QE
                    ENDIF
                    QW=DISCHARGE(IW,J_OLD)
                    IF (QW > QCOMP) THEN
                        I_NEW=IW
                        J_NEW=J_OLD
                        QCOMP=QW
                    ENDIF
                    QS=DISCHARGE(I_OLD,JS)
                    IF (QS > QCOMP)THEN
                        I_NEW=I_OLD
                        J_NEW=JS
                        QCOMP=QS
                    ENDIF
                    QSW=DISCHARGE(IW,JS)
                    IF (QSW > QCOMP)THEN
                        I_NEW=IW
                        J_NEW=JS
                        QCOMP=QSW
                    ENDIF
                    QSE=DISCHARGE(IE,JS)
                    IF (QSE > QCOMP)THEN
                        I_NEW=IE
                        J_NEW=JS
                        QCOMP=QSE
                    ENDIF
                    IF (LOCAL_DIRECTION < 2) THEN
                        I_NEW=I_OUTFLOW(BASIN_NUMBER(I_OLD,J_OLD))
                        J_NEW=J_OUTFLOW(BASIN_NUMBER(I_OLD,J_OLD))
                        I_OLD=I_NEW
                        J_OLD=J_NEW
                        IF (J_OLD >= (MY-1)) EXIT L1210
                    ELSE
                        IF (DISCHARGE(I_OLD,J_OLD) >= QCOMP) THEN
                            IF (JN < 1) JN=1
                            IF (JS > MY) JS=MY
                            WRITE(*,1222) I_OLD,J_OLD,DISCHARGE(I_OLD,J_OLD),QCOMP &
                            ,QN,QNE,QE,QSE,QS,QSW,QW,QNW, &
                            FLOW_DIRECTION(I_OLD,J_OLD),FLOW_DIRECTION(I_OLD,JN),FLOW_DIRECTION(IE,JN),FLOW_DIRECTION(IE,J_OLD), &
                            FLOW_DIRECTION(IE,JS),FLOW_DIRECTION(I_OLD,JS),FLOW_DIRECTION(IW,JS), &
                            FLOW_DIRECTION(IW,J_OLD),FLOW_DIRECTION(IW,JN)
                            WRITE(OUTHIST,1222) I_OLD,J_OLD,DISCHARGE(I_OLD,J_OLD),QCOMP  &
                            ,QN,QNE,QE,QSE,QS,QSW,QW,QNW, &
                            FLOW_DIRECTION(I_OLD,J_OLD),FLOW_DIRECTION(I_OLD,JN),FLOW_DIRECTION(IE,JN),FLOW_DIRECTION(IE,J_OLD),  &
                            FLOW_DIRECTION(IE,JS),FLOW_DIRECTION(I_OLD,JS),FLOW_DIRECTION(IW,JS), &
                            FLOW_DIRECTION(IW,J_OLD),FLOW_DIRECTION(IW,JN)
                            1222  FORMAT('BUMMER, I_OLD=',I5,' J_OLD=',I5,' Q=',G12.5,' QCOMP=',G12.5, &
                            /,' Q CLOCK FROM N=',8(G12.5,' ')   &
                            ,/,'ID=',I5,' IDS CLOCK FROM N=',8(I5,' '))
                            EXIT L1210 !  goto 1200
                        ENDIF
                        I_OLD=I_NEW
                        J_OLD=J_NEW
                    ENDIF
                    IF (IDIST <= (2*MY)) EXIT L1210
                ENDDO L1210
            ENDDO
        ELSE
            DO I=1,MX,NSTEP
                DO J=JSTART,JEND,JSTEP
                    MAXLENGTH=0
                    DO K=1,NSTEP
                        DO L=1,NSTEP
                            I_OLD=I+K-1
                            J_OLD=J+L-1
                            NLENGTH=1
                            L110: DO
                                LOCAL_DIRECTION=FLOW_DIRECTION(I_OLD,J_OLD)
                                IF (LOCAL_DIRECTION < 2) EXIT L110
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
                                IF (FLOW_DIRECTION(I_NEW,J_NEW) < 2) EXIT L110
                                I_OLD=I_NEW
                                J_OLD=J_NEW
                                NLENGTH=NLENGTH+1
                            ENDDO L110
                            IF (NLENGTH > MAXLENGTH) THEN
                                ILONGEST=I+K-1
                                JLONGEST=J+L-1
                                MAXLENGTH=NLENGTH
                            ENDIF
                        ENDDO
                    ENDDO
                    I_OLD=ILONGEST
                    J_OLD=JLONGEST
                    IDIST=0
                    L210: DO
                        XAREA=DISCHARGE(I_OLD,J_OLD)
                        A1TERM=(XAREA)**SEDIMENT_1_EXPONENT
                        A2TERM=(XAREA)**SEDIMENT_2_EXPONENT
                        CALL LOCAL_VALUES(I_OLD,J_OLD,LOCAL_INVERSE_RESISTANCE,NOMINAL_CRITICAL_SHEAR, &
                        LOCAL_CRITICAL_SHEAR, LOCAL_BEDROCK_ERODIBILITY, LOCAL_REGOLITH_CRITICAL_SHEAR,  &
                        LOCAL_REGOLITH_ERODIBILITY)
                        IF (DISCHARGE(I_OLD,J_OLD) > 0.0) THEN
                            AFACTOR = DISCHARGE(I_OLD,J_OLD)**BEDROCK_DISCHARGE_EXPONENT
                        ELSE
                            AFACTOR=0.0
                        ENDIF
                        GFACTOR = D8_GRADIENT(I_OLD,J_OLD)**BEDROCK_GRADIENT_EXPONENT
                        TEMP = AFACTOR*GFACTOR
                        SHEAR_STRESS_FACTOR=DISCHARGE_COEFFICIENT*TEMP
                        CALL SEDIMENT_TRANSPORT_FLUX(QSSS,TRANSSTAGE,D8_GRADIENT(I_OLD,J_OLD),A1TERM,A2TERM)
                        QSSS=QSSS/CHANNEL_WIDTH(I_OLD,J_OLD)
                        BEDLOAD_FLUX=AMAX1(0.0,SEDIMENT_YIELD(I_OLD,J_OLD)/CHANNEL_WIDTH(I_OLD,J_OLD))
                        IF (QSSS > 0.0) THEN
                            SEDSATFACTOR=AMAX1(0.0,(1.0-(BEDLOAD_FLUX/QSSS)))
                        ELSE
                            SEDSATFACTOR=0.0
                        ENDIF
                        BEDLOAD_FLUX=BEDLOAD_FLUX*BEDLOAD_FRACTION
                        IF (IS_SEDIMENT_COVERED(I_OLD,J_OLD)) THEN
                            KSED=1
                        ELSE
                            KSED=0
                        ENDIF
                        IF (IS_ROCK_SURFACE(I_OLD,J_OLD)) THEN
                            KROCK=1
                        ELSE
                            KROCK=0
                        ENDIF
                        IF (SUBMERGED(I_OLD,J_OLD)) THEN
                            ISUB=1
                        ELSE
                            ISUB=0
                        ENDIF
                        IF (CFN(I_OLD,J_OLD) < 1.0E+24) THEN
                            ROUTEDIV=CFN(I_OLD,J_OLD)
                        ELSE
                            ROUTEDIV=0.0
                        ENDIF
                        WRITE(OUTCHPROP,300) ILONGEST,JLONGEST,I_OLD,J_OLD,LOCAL_DIRECTION,IDIST &
                        ,KSED,KROCK  &
                        ,ELEVATION(I_OLD,J_OLD),D8_GRADIENT(I_OLD,J_OLD),EQUILIBRIUM_GRADIENT(I_OLD,J_OLD) &
                        ,ERODE_CHANNEL(I_OLD,J_OLD),ERODE_SLOPE(I_OLD,J_OLD),DISCHARGE(I_OLD,J_OLD),BEDLOAD_FLUX,QSSS &
                        ,TRANSSTAGE,SEDSATFACTOR  &
                        ,ISUB,IDO(I_OLD,J_OLD),ROUTEDIV,SEDIMENT_BASE(I_OLD,J_OLD),ERODE_REGOLITH_CHANNEL(I_OLD,J_OLD) &
                        ,CFNE(I_OLD,J_OLD),CFNW(I_OLD,J_OLD),SHEAR_STRESS_FACTOR,LOCAL_CRITICAL_SHEAR &
                        ,SEDIMENT_BASE(I_OLD,J_OLD),CFN(I_OLD,J_OLD),LAKE_OUTLET_ELEVATION(BASIN_NUMBER(I_OLD,J_OLD)) &
                        ,I_OUTFLOW(BASIN_NUMBER(I_OLD,J_OLD)),J_OUTFLOW(BASIN_NUMBER(I_OLD,J_OLD))
                        300         FORMAT(8(I6,','),10(G12.5,',') &
                        ,2(I6,','),10(G12.5,','),I5,',',I5)
                        LOCAL_DIRECTION=FLOW_DIRECTION(I_OLD,J_OLD)
                        IDIST=IDIST+1
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
                        IF (FLOW_DIRECTION(I_NEW,J_NEW) < 2) EXIT L210
                        I_OLD=I_NEW
                        J_OLD=J_NEW
                    ENDDO L210
                ENDDO
            ENDDO
        ENDIF
        CLOSE(OUTCHPROP)
        RETURN
    END !  SUBROUTINE CHANNEL_PROPERTIES
