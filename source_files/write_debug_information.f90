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
    SUBROUTINE PRINT_DEBUGGING_DATA()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,IMID,JMID
        !     **********************************************************************
        !      writes out debugging information.  modify as appropriate
        !     **********************************************************************
        !          write(outhist,500) indata
        500 FORMAT(' INDATA ',I8)
        !          write(outhist,501) params
        501 FORMAT(' PARAMS ',I8)
        !          write(outhist,502) inresist
        502 FORMAT(' INRESIST ',I8)
        !          write(outhist,503) outdata
        503 FORMAT(' OUTDATA ',I8)
        !          write(outhist,504) outhist
        504 FORMAT(' OUTHIST ',I8)
        !          write(outhist,505) outchan
        505 FORMAT(' OUTCHAN',I8)
        IF (ITERATION < 2) THEN
            !          write(outhist,508) mx
            508 FORMAT(' MX ',I8)
            !          write(outhist,509) my
            509 FORMAT(' MY ',I8)
            !         write(outhist,515) iabortmax
            515 FORMAT(' IABORTMAX ',I8)
            !         write(outhist,516) iabortmin
            516 FORMAT(' IABORTMIN ',I8)
            !         write(outhist,517) jabortmax
            517 FORMAT(' JABORTMAX ',I8)
            !         write(outhist,518) jabortmin
            518 FORMAT(' JABORTMIN',I8)
            WRITE(OUTHIST,527) MP
            527 FORMAT(' MP ',I8)
            !         write(outhist,537) idummy
            537 FORMAT(' IDUMMY',I8)
            !          write(outhist,538) maximum_iteration
            538 FORMAT(' MAXIMUM_ITERATION ',I9)
            !          write(outhist,539) elevation_print_interval
            539 FORMAT(' ELEVATION_PRINT_INTERVAL ',I8)
            !          write(outhist,540) output_print_interval
            540 FORMAT(' OUTPUT_PRINT_INTERVAL ',I8)
            !          write(outhist,541) recalculate_gradient_interval
            541 FORMAT(' RECALCULATE_GRADIENT_INTERVAL',I8)
            !          write(outhist,542) critical_gradient_use
            542 FORMAT(' CRITICAL_GRADIENT_USE ',I8)
        ENDIF
        WRITE(OUTHIST,543) ITERATION
        543 FORMAT(' ITERATION ',I9)
        IF (ITERATION < 3) THEN
            !          write(outhist,544) (writetype(i),i=1,10)
            544 FORMAT(' WRITETYPE',10I2)
            !          write(outhist,545) oneonly
            545 FORMAT(' ONEONLY ',I8)
            !          write(outhist,546) minplotarea
            546 FORMAT(' MINPLOTAREA ',I8)
            !          write(outhist,547) starting_iteration
            547 FORMAT(' STARTING_ITERATION ',I9)
            !          write(outhist,548) total_iterations
            548 FORMAT(' TOTAL_ITERATIONS',I9)
            !          write(outhist,549) powertouse
            549 FORMAT(' POWERTOUSE ',I8)
            !          write(outhist,550) layeruse
            550 FORMAT(' LAYERUSE ',I8)
            !          write(outhist,551) VARIABLE_ROCK_RESISTANCE_USE
            551 FORMAT(' VARIABLE_ROCK_RESISTANCE_USE ',I8)
            !          write(outhist,552) channelvaruse
            552 FORMAT(' CHANNELVARUSE',I8)
            !          write(outhist,553) ncrits
            553 FORMAT(' NCRITS ',I8)
            !          write(outhist,554) nmax
            554 FORMAT(' NMAX ',I8)
            !          write(outhist,555) lmax
            555 FORMAT(' LMAX',I8)
            !         write(outhist,562) sqrtoftwo
            562 FORMAT(' SQRTOFTWO ',E12.5)
            !         write(outhist,563) seed
            563 FORMAT(' SEED ',E12.5)
            !         write(outhist,564) iseed1
            564 FORMAT(' ISEED1',E12.5)
            !          write(outhist,565) boundary_lowering_rate
            565 FORMAT(' BOUNDARY_LOWERING_RATE ',E12.5)
            !          write(outhist,566) bedrock_discharge_exponent
            566 FORMAT(' BEDROCK_DISCHARGE_EXPONENT ',E12.5)
            !          write(outhist,567) bedrock_gradient_exponent
            567 FORMAT(' BEDROCK_GRADIENT_EXPONENT ',E12.5)
            !          write(outhist,568) slope_diffusivity
            568 FORMAT(' SLOPE_DIFFUSIVITY ',E12.5)
            !          write(outhist,569) bedrock_erodibility
            569 FORMAT(' BEDROCK_ERODIBILITY',E12.5)
            !          write(outhist,570) detachment_critical_shear
            570 FORMAT(' DETACHMENT_CRITICAL_SHEAR ',E12.5)
            !          write(outhist,571) critical_slope_gradient
            571 FORMAT(' CRITICAL_SLOPE_GRADIENT ',E12.5)
            !          write(outhist,572) slope_gradient_exponent
            572 FORMAT(' SLOPE_GRADIENT_EXPONENT ',E12.5)
        ENDIF
        WRITE(OUTHIST,573)TIME_INCREMENT
        573 FORMAT('TIME_INCREMENT ',E12.5)
        WRITE(OUTHIST,574) PRESENT_TIME
        574 FORMAT(' PRESENT_TIME',E12.5)
        WRITE(OUTHIST,575) CHANNEL_TIMESTEP_SCALING
        575 FORMAT(' CHANNEL_TIMESTEP_SCALING ',E12.5)
        WRITE(OUTHIST,576) MAXIMUM_ELEVATION_CHANGE
        576 FORMAT(' MAXIMUM_ELEVATION_CHANGE ',E12.5)
        WRITE(OUTHIST,577) GRADAVERAGE
        577 FORMAT(' GRADAVERAGE ',E12.5)
        WRITE(OUTHIST,578) MAXGRADIENT
        578 FORMAT(' MAXGRADIENT ',E12.5)
        WRITE(OUTHIST,579) TMULT
        579 FORMAT(' TMULT',E12.5)
        WRITE(OUTHIST,580) NEWBASE
        580 FORMAT(' NEWBASE ',E12.5)
        WRITE(OUTHIST,581) MAXIMUM_DIFFUSIVITY_INCREASE
        581 FORMAT(' MAXIMUM_DIFFUSIVITY_INCREASE ',E12.5)
        IF(ITERATION < 2) THEN
            !          write(outhist,582) critical_source_divergence
            582 FORMAT(' CRITICAL_SOURCE_DIVERGENCE ',E12.5)
            !          write(outhist,583) critical_gradient_term
            583 FORMAT(' CRITICAL_GRADIENT_TERM ',E12.5)
            !          write(outhist,584) regolith_erodibility_factor
            584 FORMAT(' REGOLITH_ERODIBILITY_FACTOR',E12.5)
        ENDIF
        WRITE(OUTHIST,585) MAXIMUM_CHANNEL_TIMESTEP
        585 FORMAT(' MAXIMUM_CHANNEL_TIMESTEP ',E12.5)
        IF (ITERATION < 2) THEN
            !          write(outhist,588) RESISTANCE_VARIABILITY
            588 FORMAT(' RESISTANCE_VARIABILITY ',E12.5)
            !          write(outhist,589) channelvarscale
            589 FORMAT(' CHANNELVARSCALE',E12.5)
            !          write(outhist,590) surface_layer_resistance
            590 FORMAT(' SURFACE_LAYER_RESISTANCE ',E12.5)
            !          write(outhist,591) layermax
            591 FORMAT(' LAYERMAX ',E12.5)
            !          write(outhist,592) layermin
            592 FORMAT(' LAYERMIN',E12.5)
            !          write(outhist,593) cross_weighting
            593 FORMAT(' CROSS_WEIGHTING ',E12.5)
            !          write(outhist,594) diagonal_weighting
            594 FORMAT(' DIAGONAL_WEIGHTING ',E12.5)
            !          write(outhist,595) maxcrit
            595 FORMAT(' MAXCRIT ',E12.5)
            !         write(outhist,596) one_sixth
            596 FORMAT(' ONE_SIXTH',E12.5)
        ENDIF
        IF (DONE ) THEN
            WRITE(OUTHIST,597)
            597 FORMAT(' DONE IS TRUE')
        ELSE
            WRITE(OUTHIST,598)
            598 FORMAT(' DONE  IS FALSE')
        ENDIF
        IF (ITERATION < 2) THEN
            IF (WRITE_ABSOLUTE_ELEVATION ) THEN
                !          write(outhist,599)
                599 FORMAT(' WRITE_ABSOLUTE_ELEVATION IS TRUE')
            ELSE
                !          write(outhist,600)
                600 FORMAT(' WRITE_ABSOLUTE_ELEVATION  IS FALSE')
            ENDIF
            IF (USE_CRITICAL_SLOPE_CRADIENT ) THEN
                !          write(outhist,601)
                601 FORMAT(' USE_CRITICAL_SLOPE_CRADIENT IS TRUE')
            ELSE
                !          write(outhist,602)
                602 FORMAT(' USE_CRITICAL_SLOPE_CRADIENT  IS FALSE')
            ENDIF
        ENDIF
        IF (ABORT_SIMULATION) THEN
            WRITE(OUTHIST,603)
            603 FORMAT(' ABORT_SIMULATION IS TRUE')
        ELSE
            WRITE(OUTHIST,604)
            604 FORMAT(' ABORT_SIMULATION IS FALSE')
        ENDIF
        IF (ITERATION < 2) THEN
            IF (USELAYER ) THEN
                !          write(outhist,605)
                605 FORMAT(' USELAYER IS TRUE')
            ELSE
                !          write(outhist,606)
                606 FORMAT(' USELAYER  IS FALSE')
            ENDIF
            IF (USE_3D_SLOPE_RESISTANCE ) THEN
                !          write(outhist,607)
                607 FORMAT(' USE_3D_SLOPE_RESISTANCE IS TRUE')
            ELSE
                !          write(outhist,608)
                608 FORMAT(' USE_3D_SLOPE_RESISTANCE  IS FALSE')
            ENDIF
        ENDIF
        IMID=MX/2
        JMID=MY/2
        WRITE(OUTHIST,611)
        611 FORMAT(/,' ID  MATRIX',/)
        DO  J= JMID-3,JMID+3
            WRITE(OUTHIST,612)(FLOW_DIRECTION(I,J),I=IMID-3,IMID+3)
            612 FORMAT(8I6)
        ENDDO
        WRITE(OUTHIST,613)
        613 FORMAT(/,' DRAINAGE AREA  MATRIX',/)
        DO J= JMID-3,JMID+3
            WRITE(OUTHIST,614)(DRAINAGE_AREA (I,J),I=IMID-3,IMID+3)
            614 FORMAT(8I6)
        ENDDO
        WRITE(OUTHIST,615)
        615 FORMAT(/,' EL  MATRIX',/)
        DO J= JMID-3,JMID+3
            WRITE(OUTHIST,616)(ELEVATION(I,J),I=IMID-3,IMID+3)
            616 FORMAT(7(G11.4,' '))
        ENDDO
        WRITE(OUTHIST,617)
        617 FORMAT(/,' D8_GRADIENT  MATRIX',/)
        DO  J= JMID-3,JMID+3
            WRITE(OUTHIST,618)(D8_GRADIENT(I,J),I=IMID-3,IMID+3)
            618 FORMAT(7(G11.4,' '))
        ENDDO
        IF(ITERATION > 0) THEN
            WRITE(OUTHIST,619)
            619 FORMAT(/,' ECHANNEL  MATRIX',/)
            DO J= JMID-3,JMID+3
                WRITE(OUTHIST,620)(ERODE_CHANNEL(I,J),I=IMID-3,IMID+3)
                620 FORMAT(7(G11.4,' '))
            ENDDO
            WRITE(OUTHIST,621)
            621 FORMAT(/,' ERODE_SLOPE  MATRIX',/)
            DO  J= JMID-3,JMID+3
                WRITE(OUTHIST,622)(ERODE_SLOPE (I,J),I=IMID-3,IMID+3)
                622 FORMAT(7(G11.4,' '))
            ENDDO
            IF (USELAYER) THEN
                WRITE(OUTHIST,623)
                623 FORMAT(/,' RELATIVE_RESISTANCE  MATRIX',/)
                DO J= JMID-3,JMID+3
                    WRITE(OUTHIST,624)(RELATIVE_RESISTANCE (I,J),I=IMID-3,IMID+3)
                    624 FORMAT(7(G11.4,' '))
                ENDDO
            ENDIF
            WRITE(OUTHIST,625)
            625 FORMAT(/,' CFW MATRIX',/)
            DO J= JMID-3,JMID+3
                WRITE(OUTHIST,626)(CFW(I,J),I=IMID-3,IMID+3)
                626 FORMAT(7(G11.4,' '))
            ENDDO
            WRITE(OUTHIST,725)
            725 FORMAT(/,' CFN MATRIX',/)
            DO J= JMID-3,JMID+3
                WRITE(OUTHIST,726)(CFN(I,J),I=IMID-3,IMID+3)
                726 FORMAT(7(G11.4,' '))
            ENDDO
            WRITE(OUTHIST,727)
            727 FORMAT(/,' CFNE MATRIX',/)
            DO J= JMID-3,JMID+3
                WRITE(OUTHIST,728)(CFNE(I,J),I=IMID-3,IMID+3)
                728 FORMAT(7(G11.4,' '))
            ENDDO
            WRITE(OUTHIST,729)
            729 FORMAT(/,' CFNW  MATRIX',/)
            DO  J= JMID-3,JMID+3
                WRITE(OUTHIST,730)(CFNW(I,J),I=IMID-3,IMID+3)
                730 FORMAT(7(G11.4,' '))
            ENDDO
        ENDIF
        WRITE (OUTHIST,701) TOTAL_BASINS
        701 FORMAT (' THERE ARE ',I5,' BASINS')
        DO I=1,TOTAL_BASINS
            IF (ENCLOSED(I)) THEN
                WRITE (OUTHIST,702)I,I_OUTFLOW(I),J_OUTFLOW(I),DOWNSTREAM_BASIN(I),BASIN_DRAINAGE_AREA(I) &
                ,LAKE_OUTLET_ELEVATION(I)
                702 FORMAT ('CLOSED BASIN ',I5,'; OUTLET AT ',2I5,' DRAINS TO ',I5,/,&
                ' AREA OF ',I7,' OUTLET ELEVATION ',G15.7)
            ENDIF
        ENDDO
        RETURN
    END ! SUBROUTINE PRINT_DEBUGGING_DATA
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PRINT_INTEGER_MATRIX_DATA(IMAT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J
        INTEGER :: IMAT(IMMX,JMMX)
        DO  J=JWINLOW,JWINHIGH
            WRITE (OUTHIST,110) (IMAT(I,J),I=IWINLOW,IWINHIGH)
            110 FORMAT(20I5)
        ENDDO
        RETURN
    END ! SUBROUTINE PRINT_INTEGER_MATRIX_DATA
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PRINT_BASIN_INFORMATION()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: N
        DO  N=1,TOTAL_BASINS
            IF ((I_OUTFLOW(N) <= IWINHIGH).AND.(I_OUTFLOW(N) >= IWINLOW).AND.&
            (J_OUTFLOW(N) >= JWINLOW).AND.(J_OUTFLOW(N) <= JWINHIGH)) THEN
                WRITE(OUTHIST,100) N,I_OUTFLOW(N),J_OUTFLOW(N)
                100   FORMAT(' OUTLET OF BASIN ',I5,' AT I=',I5,' J=',I5)
            ENDIF
        ENDDO
        RETURN
    END !  SUBROUTINE PRINT_BASIN_INFORMATION
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PRINT_LOGICAL_MATRIX_DATA(LMAT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J
        LOGICAL (4) :: LMAT(IMMX,JMMX)
        INTEGER :: IROW(IMMX)
        L100: DO  J=JWINLOW,JWINHIGH
            L50: DO I=IWINLOW,IWINHIGH
                IF (LMAT(I,J)) THEN
                    IROW(I)=1
                ELSE
                    IROW(I)=0
                ENDIF
            ENDDO L50
            WRITE (OUTHIST,110) (IROW(I),I=IWINLOW,IWINHIGH)
            110 FORMAT (20I4)
        ENDDO L100
        RETURN
    END ! SUBROUTINE PRINT_LOGICAL_MATRIX_DATA
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PRINT_REAL_MATRIX_DATA(RMAT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: IIII,I,J,IBASE,IRANGE
        REAL(4) :: RMAT(IMMX,JMMX)
        IRANGE=IWINHIGH-IWINLOW
        L100: DO J=JWINLOW,JWINHIGH
            IIII=IWINLOW+MIN(9,IRANGE)
            WRITE (OUTHIST,110) (RMAT(I,J),I=IWINLOW,IIII)
            110 FORMAT(10(G10.3,' '))
            IF (IIII >= IWINHIGH) CYCLE
            130 IBASE=IIII+1
            IRANGE=IWINHIGH-IBASE
            IIII=IBASE+MIN(9,IRANGE)
            WRITE (OUTHIST,120) (RMAT(I,J),I=IBASE,IIII)
            120 FORMAT('     ',10(G10.3,' '))
            IF (IIII < IWINHIGH) CYCLE
        ENDDO L100
        RETURN
    END !  SUBROUTINE PRINT_REAL_MATRIX_DATA
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE SUMMARIZE_MATRIX_DATA(AMATRIX,AVG,MAXVAL,MINVAL)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL(4) :: AMATRIX(IMMX,JMMX)
        REAL(4) :: SUMX,SUMXX,SUMN,MAXVAL,MINVAL,AAAA,AVG,STD,TEMP
        INTEGER I,J
        SUMX=0.0
        SUMXX=0.0
        SUMN=0.0
        MAXVAL=-1.0E25
        MINVAL=1.0E25
        L100: DO I=1,MX
            M100: DO J=1,MY
                AAAA=AMATRIX(I,J)
                SUMX=SUMX+AAAA
                SUMXX=SUMXX+AAAA*AAAA
                SUMN=SUMN+1.0
                IF (AAAA > MAXVAL) MAXVAL=AAAA
                IF (AAAA < MINVAL) MINVAL=AAAA
            ENDDO M100
        ENDDO L100
        IF (SUMN > 1.0) THEN
            AVG=SUMX/SUMN
            TEMP=SUMXX-SUMN*AVG*AVG
            IF (TEMP > 0.0) THEN
                STD=SQRT((TEMP)/(SUMN-1.0))
            ELSE
                STD=0.0
            ENDIF
        ELSE
            AVG=0.0
            STD=0.0
        ENDIF
        WRITE(OUTHIST,110) SUMX,AVG,STD,MINVAL,MAXVAL
        110 FORMAT(' SUM=',G12.4,' AVG=',G12.4,' STD=',G12.4,' MIN=',G12.4, &
        ' MAX=',G12.4)
        RETURN
    END ! SUBROUTINE SUMMARIZE_MATRIX_DATA(
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE SUMMARIZE_REGOLITH_DATA(AMATRIX)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL(4) :: AMATRIX(IMMX,JMMX)
        REAL(4) :: SUMX,SUMXX,SUMN,MAXVAL,MINVAL,AAAA,AVG,STD,TEMP
        INTEGER :: I,J
        WRITE(OUTHIST,510)
        510   FORMAT(' REGOLITH SLOPES')
        SUMX=0.0
        SUMXX=0.0
        SUMN=0.0
        MAXVAL=-1.0E25
        MINVAL=1.0E25
        DO  I=1,MX
            DO  J=1,MY
                IF ( IS_ROCK_SURFACE(I,J) ) CYCLE
                AAAA=AMATRIX(I,J)
                SUMX=SUMX+AAAA
                SUMXX=SUMXX+AAAA*AAAA
                SUMN=SUMN+1.0
                IF (AAAA > MAXVAL) MAXVAL=AAAA
                IF (AAAA < MINVAL) MINVAL=AAAA
            ENDDO
        ENDDO
        IF (SUMN > 1.0) THEN
            AVG=SUMX/SUMN
            TEMP=SUMXX-SUMN*AVG*AVG
            IF (TEMP > 0.0) THEN
                STD=SQRT((TEMP)/(SUMN-1.0))
            ELSE
                STD=0.0
            ENDIF
        ELSE
            AVG=0.0
            STD=0.0
        ENDIF
        WRITE(OUTHIST,110) SUMX,AVG,STD,MINVAL,MAXVAL,SUMN
        110 FORMAT(' SUM=',G12.4,' AVG=',G12.4,' STD=',G12.4,' MIN=',G12.4,&
        ' MAX=',G12.4,' N=',G15.6)
        WRITE(OUTHIST,520)
        520   FORMAT(' BEDROCK SLOPES')
        SUMX=0.0
        SUMXX=0.0
        SUMN=0.0
        MAXVAL=-1.0E25
        MINVAL=1.0E25
        DO  I=1,MX
            DO  J=1,MY
                IF (.NOT. IS_ROCK_SURFACE(I,J) ) CYCLE
                AAAA=AMATRIX(I,J)
                SUMX=SUMX+AAAA
                SUMXX=SUMXX+AAAA*AAAA
                SUMN=SUMN+1.0
                IF (AAAA > MAXVAL) MAXVAL=AAAA
                IF (AAAA < MINVAL) MINVAL=AAAA
            ENDDO
        ENDDO
        IF (SUMN > 1.0) THEN
            AVG=SUMX/SUMN
            TEMP=SUMXX-SUMN*AVG*AVG
            IF (TEMP > 0.0) THEN
                STD=SQRT((TEMP)/(SUMN-1.0))
            ELSE
                STD=0.0
            ENDIF
        ELSE
            AVG=0.0
            STD=0.0
        ENDIF
        WRITE(OUTHIST,110) SUMX,AVG,STD,MINVAL,MAXVAL,SUMN

        RETURN
    END ! SUBROUTINE SUMMARIZE_REGOLITH_DATA
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE SUMMARIZE_LOGICAL_MATRIX(AMATRIX)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        LOGICAL :: AMATRIX(IMMX,JMMX)
        REAL(4) :: SUMX,SUMN,PROP
        INTEGER :: I,J
        SUMX=0.0
        SUMN=0.0
        DO  I=1,MX
            DO  J=1,MY
                IF (AMATRIX(I,J)) SUMX=SUMX+1.0
                SUMN=SUMN+1.0
            ENDDO
        ENDDO
        IF (SUMN > 0.0) THEN
            PROP=SUMX/SUMN
        ELSE
            PROP=0.0
        ENDIF
        WRITE(OUTHIST,110) SUMN,PROP
        110 FORMAT(' NO. POINTS=',G12.4,' PROPORTION TRUE=',G12.4)
        RETURN
    END !  SUBROUTINE SUMMARIZE_LOGICAL_MATRIX
