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
    SUBROUTINE MOMENTS(ARRAY,VALUE)
        IMPLICIT NONE
        REAL(4), INTENT(INOUT) :: ARRAY(5)
        REAL(4), INTENT(IN) ::VALUE
        !      *********************************************************************
        !       calculates and increments the sum of the 1st four powers of the
        !         passed value in the passed array as well as the total number of
        !         observations
        !      *********************************************************************
        ARRAY(1)=ARRAY(1)+VALUE
        ARRAY(2)=ARRAY(2)+VALUE*VALUE
        ARRAY(3)=ARRAY(3)+VALUE*VALUE*VALUE
        ARRAY(4)=ARRAY(4)+VALUE*VALUE*VALUE*VALUE
        ARRAY(5)=ARRAY(5)+1.0
        RETURN
    END ! SUBROUTINE MOMENTS
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE RESET_MOMENTS(ARRAY)
        !      *********************************************************************
        !       Initializes (to zero) a passed array of powers
        !      *********************************************************************
        IMPLICIT NONE
        REAL(4), INTENT(INOUT) :: ARRAY(5)
        INTEGER I
        DO  I=1,5
            ARRAY(I)=0.0
        ENDDO
        RETURN
    END ! SUBROUTINE RESET_MOMENTS
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE CALCULATE_MOMENTS(ARRAY,AVG,VARI,SKEW,KURT,TOTAL)
        !      *********************************************************************
        !       Calculates average,variance,skewness,kurtosis and total observations
        !         from the passed array of sum of powers of the variable
        !      *********************************************************************
        IMPLICIT NONE
        REAL(4), INTENT(IN) :: ARRAY(5)
        REAL(4), INTENT(OUT) :: AVG,VARI,SKEW,KURT,TOTAL
        TOTAL=ARRAY(5)
        IF (TOTAL > 0.0) THEN
            AVG=ARRAY(1)/ARRAY(5)
        ELSE
            AVG=0.0
        ENDIF
        IF (TOTAL > 1.0) THEN
            VARI=(ARRAY(2)-TOTAL*AVG*AVG)/(TOTAL-1.0)
            SKEW=(ARRAY(3)-3.0*AVG*ARRAY(2)  &
            + 2.0*TOTAL*AVG*AVG*AVG)/(TOTAL-1.0)
            KURT=(ARRAY(4)-4.0*AVG*ARRAY(3)+6.0*ARRAY(2)*AVG*AVG &
            -3.0*TOTAL*AVG*AVG*AVG*AVG)/(TOTAL-1.0)
        ELSE
            VARI=0.0
            SKEW=0.0
            KURT=0.0
        ENDIF
        RETURN
    END ! SUBROUTINE CALCULATE_MOMENTS
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PRINT_MOMENTS(AVG,VARI,SKEW,KURT,TOTAL)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        !      *********************************************************************
        !       Prints the moments calculated by calcmoments()
        !      *********************************************************************
        REAL(4), INTENT(IN) :: AVG,VARI,SKEW,KURT,TOTAL
        WRITE(OUTHIST,500)AVG,VARI,SKEW,KURT,TOTAL
        500 FORMAT(' AVG.=',G12.5,'VAR.=',G12.5,' SKEW=',G12.5,/, &
        '   KURT.=',G12.5,' TOT. OBS.=',G12.5)
        RETURN
    END ! SUBROUTINE PRINT_MOMENTS
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE FIND_TOPOGRAPHIC_EXTREMA(I,J,XSUMMIT,XSINK,XSADDLE)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        LOGICAL, INTENT(OUT) :: XSUMMIT,XSINK,XSADDLE
        INTEGER, INTENT(IN) :: I,J
        INTEGER :: IW,IE,JN,JS,IIII,ICOUNT
        REAL(4) :: ELOCAL
        ELOCAL=ELEVATION(I,J)
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
        XSUMMIT=.TRUE.
        IF (ELEVATION(I,JN) >= ELOCAL) THEN
            XSUMMIT=.FALSE.
            GOTO 100
        ENDIF
        IF (ELEVATION(I,JS) >= ELOCAL) THEN
            XSUMMIT=.FALSE.
            GOTO 100
        ENDIF
        IF (ELEVATION(IW,JN) >= ELOCAL) THEN
            XSUMMIT=.FALSE.
            GOTO 100
        ENDIF
        IF (ELEVATION(IW,JS) >= ELOCAL) THEN
            XSUMMIT=.FALSE.
            GOTO 100
        ENDIF
        IF (ELEVATION(IE,JN) >= ELOCAL) THEN
            XSUMMIT=.FALSE.
            GOTO 100
        ENDIF
        IF (ELEVATION(IE,JS) >= ELOCAL) THEN
            XSUMMIT=.FALSE.
            GOTO 100
        ENDIF
        IF (ELEVATION(IW,J) >= ELOCAL) THEN
            XSUMMIT=.FALSE.
            GOTO 100
        ENDIF
        IF (ELEVATION(IE,J) >= ELOCAL) THEN
            XSUMMIT=.FALSE.
            GOTO 100
        ENDIF
        100   XSINK=.TRUE.
        IF (ELEVATION(I,JN) <= ELOCAL) THEN
            XSINK=.FALSE.
            GOTO 200
        ENDIF
        IF (ELEVATION(I,JS) <= ELOCAL) THEN
            XSINK=.FALSE.
            GOTO 200
        ENDIF
        IF (ELEVATION(IW,JN) <= ELOCAL) THEN
            XSINK=.FALSE.
            GOTO 200
        ENDIF
        IF (ELEVATION(IW,J) <= ELOCAL) THEN
            XSINK=.FALSE.
            GOTO 200
        ENDIF
        IF (ELEVATION(IW,JS) <= ELOCAL) THEN
            XSINK=.FALSE.
            GOTO 200
        ENDIF
        IF (ELEVATION(IE,JN) <= ELOCAL) THEN
            XSINK=.FALSE.
            GOTO 200
        ENDIF
        IF (ELEVATION(IE,J) <= ELOCAL) THEN
            XSINK=.FALSE.
            GOTO 200
        ENDIF
        IF (ELEVATION(IE,JS) <= ELOCAL) THEN
            XSINK=.FALSE.
            GOTO 200
        ENDIF
        200   XSADDLE=.TRUE.
        ICOUNT=0
        IIII=-1
        IF (ELEVATION(I,JN) > ELOCAL) THEN
            IIII=1
        ENDIF
        IF (IIII > 0) THEN
            IF (ELEVATION(IE,JN) < ELOCAL) THEN
                IIII=-1
                ICOUNT=ICOUNT+1
            ENDIF
        ELSE
            IF(ELEVATION(IE,JN) > ELOCAL) THEN
                IIII=1
                ICOUNT=ICOUNT+1
            ENDIF
        ENDIF
        IF (IIII > 0) THEN
            IF (ELEVATION(IE,J) < ELOCAL) THEN
                IIII=-1
                ICOUNT=ICOUNT+1
            ENDIF
        ELSE
            IF(ELEVATION(IE,J) > ELOCAL) THEN
                IIII=1
                ICOUNT=ICOUNT+1
            ENDIF
        ENDIF
        IF (IIII > 0) THEN
            IF (ELEVATION(IE,JS) < ELOCAL) THEN
                IIII=-1
                ICOUNT=ICOUNT+1
            ENDIF
        ELSE
            IF(ELEVATION(IE,JS) > ELOCAL) THEN
                IIII=1
                ICOUNT=ICOUNT+1
            ENDIF
        ENDIF
        IF (IIII > 0) THEN
            IF (ELEVATION(I,JS) < ELOCAL) THEN
                IIII=-1
                ICOUNT=ICOUNT+1
            ENDIF
        ELSE
            IF(ELEVATION(I,JS) > ELOCAL) THEN
                IIII=1
                ICOUNT=ICOUNT+1
            ENDIF
        ENDIF
        IF (IIII > 0) THEN
            IF (ELEVATION(IW,JS) < ELOCAL) THEN
                IIII=-1
                ICOUNT=ICOUNT+1
            ENDIF
        ELSE
            IF(ELEVATION(IW,JS) > ELOCAL) THEN
                IIII=1
                ICOUNT=ICOUNT+1
            ENDIF
        ENDIF
        IF (IIII > 0) THEN
            IF (ELEVATION(IW,J) < ELOCAL) THEN
                IIII=-1
                ICOUNT=ICOUNT+1
            ENDIF
        ELSE
            IF(ELEVATION(IW,J) > ELOCAL) THEN
                IIII=1
                ICOUNT=ICOUNT+1
            ENDIF
        ENDIF
        IF (IIII > 0) THEN
            IF (ELEVATION(IW,JN) < ELOCAL) THEN
                IIII=-1
                ICOUNT=ICOUNT+1
            ENDIF
        ELSE
            IF(ELEVATION(IW,JN) > ELOCAL) THEN
                IIII=1
                ICOUNT=ICOUNT+1
            ENDIF
        ENDIF
        IF (IIII > 0) THEN
            IF (ELEVATION(I,JN) < ELOCAL) THEN
                IIII=-1
                ICOUNT=ICOUNT+1
            ENDIF
        ELSE
            IF(ELEVATION(I,JN) > ELOCAL) THEN
                IIII=1
                ICOUNT=ICOUNT+1
            ENDIF
        ENDIF
        IF (ICOUNT < 3) XSADDLE=.FALSE.
        RETURN
    END !SUBROUTINE FIND_TOPOGRAPHIC_EXTREMA 
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DIVERGENCE5X5(I,J,DIVERGE)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,J
        INTEGER :: IW,IE,JN,JS,IWW,IEE,JNN,JSS
        REAL*4, INTENT(OUT) :: DIVERGE
        !    Calculates topographic diverence with a 5x5 kernel.  For non-periodic
        !      boundaries it returns zero for the first two rows/columns
                IW=I-1
                IE=I+1
                IWW=I-2
                IEE=I+2
                IF (IS_X_PERIODIC) THEN
                    IF (IW < 1) IW=IW+MX
                    IF (IE > MX) IE=IW-MX
                    IF (IWW < 1) IWW=IWW+MX
                    IF (IEE > MX) IEE=IEE-MX
                ELSE
                   IF ((IWW<1).OR.(IEE>MX)) THEN
                       DIVERGE=0.0
                       RETURN
                   ENDIF
                ENDIF
                JN=J-1
                JS=J+1
                JNN=J-2
                JSS=J+2
                IF (IS_Y_PERIODIC) THEN
                    IF (JN < 1) JN=JN+MX
                    IF (JS > MY) JS=JS-MX
                    IF (JNN < 1) JNN=JNN+MX
                    IF (JSS > MY) JSS=JSS-MX
                ELSE
                    IF ((JNN < 1).OR.(JSS>MX)) THEN
                        DIVERGE=0.0
                        RETURN
                    ENDIF
                ENDIF
                DIVERGE=(ELEVATION(I,JNN)+ELEVATION(IW,JN)+ELEVATION(IE,JN) &
                    + ELEVATION(IWW,J)+ELEVATION(IEE,J)+ELEVATION(IW,JS) &
                    + ELEVATION(IE,JS) + ELEVATION(I,JSS) &
                    + 2.0*(ELEVATION(I,JN)+ELEVATION(IW,J)+ELEVATION(IE,J) + ELEVATION(I,JS)) &
                    - 16.0*ELEVATION(I,J))/CELL_AREA
                RETURN
                END !   SUBROUTINE DIVERGENCE5X5    
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE CALCULATE_TOPO_DIVERGENCE()
        !      *********************************************************************
        !       Calculates and prints the moments of several morphometric variables:
        !               profile curvature
        !               planform curvature
        !               gradient divergence (topographic laplacian)
        !               ln(drainage area/gradient)
        !       Summarizes the average and variance of the following morphometric
        !         variables by category of support area (13 categories of area):
        !               profile curvature
        !               planform curvature
        !               gradient divergence
        !               ln(drainage area/gradient)
        !       Prints the correlation between the following variables:
        !               horizontal versus planform curvature
        !               gradient and divergence
        !               ln(drainage area/gradient) and divergence
        !       Prints the average channel and slope divergence and the range of divergence
        !       Prints the sink,summit, and saddle density
        !       Prints the relative frequencies of slope elements with the following
        !         classes:
        !               positive profile curvature & positive planform curvature
        !               positive profile curvature & negative planform curvature
        !               negative profile curvature & positive planform curvature
        !               negative profile curvature & negative planform curvature
        !
        !       primarily intended for steady-state simulations
        !   MODIFES: CFN, CFW, CFNE, CFNW
        !   CALLS:  MOMENTS, CALCULATE_MOMENTS, RESET_MOMENTS, PRINT_MOMENTS, FIND_TOPOGRAPHIC_EXTREMA
        !      *********************************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER I,J,NUMB1,NUMB2,K
        REAL(4) :: TEMP1,  RESULTANT
        REAL(4) :: CMAX,CMIN,CAVG,CSLOPEAVG,CCHANNELAVG
        REAL(4) :: NTOT,CSD
        REAL(4) :: Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10,Z11,Z12
        REAL(4) :: DD,EE,FF,GG,HH,DNM,SLP,DIV
        REAL(4) :: PROFCURV,PLANCURV,PROFPLANSUM,SLOPEDIVSUM
        REAL(4) :: LTBDIVSUM
        REAL(4) :: NSUMMIT,NSINK,NSADDLE,NPLPL,NPLMN,NMNPL,NMNMN
        REAL(4) :: DIVSUMS(5),PLANSUMS(5),PROFSUMS(5),SLOPESUMS(5)
        REAL(4) :: LTBSUMS(5),LTANB
        REAL(4) :: AVG,VARI,SKEW,KURT,TOTAL
        REAL(4) :: SUMXY,SUMX,SUMY,SUMXX,SUMYY,CORREL
        REAL(4) :: ALLTOT,SLPCAT(13),SLPSQCAT(13),DIVCAT(13),DIVSQCAT(13)
        REAL(4) :: PLANCAT(13),PLANSQCAT(13),PROFCAT(13),PROFSQCAT(13)
        REAL(4) :: LTBCAT(13),LTBSQCAT(13),AVGSLPCAT(13)
        REAL(4) :: NUMBCAT(13),LBOUND(13)
        LOGICAL XSUMMIT,XSINK,XSADDLE,XSKIP
        !     **********************************************************************
        !      The matrices hold the statistical sums for calculating average, variance,
        !         skewness, and kurtosis of various local properties.  the variables for
        !         which these moments are calculated are:
        !             divsums - slope divergence
        !             plansums - planimetric curvature
        !             profsums - profile curvature
        !             slopesums - local slope gradient (gradient to lowest of surrounding 8
        !                         points
        !             lbtsums - ln(local gradient/drainage area)
        !  CALLS: RESET_MOMENTS
        !     **********************************************************************
        CALL RESET_MOMENTS(DIVSUMS)
        CALL RESET_MOMENTS(PLANSUMS)
        CALL RESET_MOMENTS(PROFSUMS)
        CALL RESET_MOMENTS(SLOPESUMS)
        CALL RESET_MOMENTS(LTBSUMS)
        !     **********************************************************************
        !      The following sums are calculated for each of the thirteen ranges of drainage
        !         area:
        !            slpcat - sum of local slope gradient (gradient from
        !                     given point to lowest of surrounding 8 points)
        !            slpsqcat - sum of squares of local gradient
        !            divcat - sum of slope divergences
        !            divsqcat - sum of squares of slope divergences
        !            plancat - sum of planimetric curvatures
        !            plansqcat - sum of squares of planimentric curvatures
        !            profcat - sum of profile curvatures
        !            profsqcat - sum of squares of profile curvatures
        !            ltbcat  - sum of ln(drainage area/slope)
        !            ltbsqcat - sum of squares of ln(drainage area/sloep)
        !            avgslpcat - sum of average slope gradient (slope of trend surface fitted
        !                     through 9 points surrounding the given point
        !            numbcat - number of points in each range of drainage area
        !            ltbdivsum - sum of cross products of ln(area/gradient) and divergence
        !     **********************************************************************
        DO  I=1,13
            SLPCAT(I)=0.0
            SLPSQCAT(I)=0.0
            DIVCAT(I)=0.0
            DIVSQCAT(I)=0.0
            PLANCAT(I)=0.0
            PLANSQCAT(I)=0.0
            PROFCAT(I)=0.0
            PROFSQCAT(I)=0.0
            LTBCAT(I)=0.0
            LTBSQCAT(I)=0.0
            AVGSLPCAT(I)=0.0
            NUMBCAT(I)=0.0
        ENDDO
        !     **********************************************************************
        !      The values of lbound are the lower values of drainage area for the 13 classes
        !         of drainage area range for which various average properties are
        !         calculated.
        !      For example, the range of areas for the first four categories are:
        !         category     range of areas
        !            1               1
        !            2              2-3
        !            3              4-9
        !            4             10-19
        !     **********************************************************************
        LBOUND( 1)=1.0*CELL_SIZE*CELL_SIZE
        LBOUND( 2)=2.0*CELL_SIZE*CELL_SIZE
        LBOUND( 3)=4.0*CELL_SIZE*CELL_SIZE
        LBOUND( 4)=10.0*CELL_SIZE*CELL_SIZE
        LBOUND( 5)=20.0*CELL_SIZE*CELL_SIZE
        LBOUND( 6)=40.0*CELL_SIZE*CELL_SIZE
        LBOUND( 7)=100.0*CELL_SIZE*CELL_SIZE
        LBOUND( 8)=200.0*CELL_SIZE*CELL_SIZE
        LBOUND( 9)=400.0*CELL_SIZE*CELL_SIZE
        LBOUND(10)=1000.0*CELL_SIZE*CELL_SIZE
        LBOUND(11)=2000.0*CELL_SIZE*CELL_SIZE
        LBOUND(12)=4000.0*CELL_SIZE*CELL_SIZE
        LBOUND(13)=10000.0*CELL_SIZE*CELL_SIZE
        !     **********************************************************************
        !      Other variables:
        !         profplansum - sum of cross products of profile and planform curvature
        !         slopedivsum - sum of cross products of slope divergence and gradient
        !         nsummit - number of summits (points for which all surrounding points are
        !                   lower)
        !         nsink  - number of sinkholes (points for which all surrounding points are
        !                  higher)
        !         nsaddle - number of cols - points where profile is both convex and concave
        !                   in different directions
        !         alltot - total number of points
        !         nplpl - number of points with positive planform curvature and positive
        !                 profile curvature
        !         nplmn - number of points with positive planform curvature and negative
        !                 profile curvature
        !         nmnpl - number of points with negative planform curvature and positive
        !                 profile curvature
        !         nmnmn - number of points with both planform and profile curvature being
        !                 negative
        !         cmax -  maximum gradient divergence
        !         cmin -  minimum gradient divergence
        !         cavg -  sum of gradient divergences
        !         csd  -  sum of squares of gradient divergences
        !         cslopeavg - sum of divergences for locations with no fluvial erosion
        !         cchannelavg - sum of divergences for locations with fluvial erosion
        !         numb1 - number of locations that have fluvial erosion
        !         numb2 - number of locations lacking fluvial erosion
        !     **********************************************************************
        PROFPLANSUM=0.0
        SLOPEDIVSUM=0.0
        LTBDIVSUM=0.0
        ALLTOT=0.0
        NSUMMIT=0.0
        NSADDLE=0.0
        NSINK=0.0
        NPLPL=0.0
        NPLMN=0.0
        NMNPL=0.0
        NMNMN=0.0
        CMAX = -1.0E+12
        CMIN = 1.0E+12
        CAVG = 0.0
        CSD = 0.0
        CSLOPEAVG = 0.0
        CCHANNELAVG = 0.0
        NUMB1 = 0
        NUMB2 = 0
        !     **********************************************************************
        !      Cycle through all points in the interior of the matrix
        !     **********************************************************************
        L100: DO  J=2,MY-1
            M100: DO  I=2,MX-1
                !     **********************************************************************
                !      Various properties of the given locality are calculated using formulas based
                !         upon fitting a surface to the 9 points surrounding and including the given
                !         point - see zevenbergen and thorne (1987, earth surface processes and
                !         landforms, 12, 47-56)
                !      Local variables include
                !         slp - average slope gradient (stored in matrix cfnw)
                !         div - local gradient divergence (negative of this value stored in matrix cfw)
                !         profcurv - local profile curvature (stored in matrix cfn)
                !         plancurv - local planform curvature (stored in matrix cfne)
                !     **********************************************************************
                Z1=ELEVATION(I-1,J-1)
                Z2=ELEVATION(I,J-1)
                Z3=ELEVATION(I+1,J-1)
                Z4=ELEVATION(I-1,J)
                Z5=ELEVATION(I,J)
                Z6=ELEVATION(I+1,J)
                Z7=ELEVATION(I-1,J+1)
                Z8=ELEVATION(I,J+1)
                Z9=ELEVATION(I+1,J+1)
                DD=(Z1+Z3+Z4+Z6+Z7+Z9)/(6.0*CELL_SIZE)-(Z2+Z5+Z8)/(3.0*CELL_SIZE)
                EE=(Z1+Z2+Z3+Z7+Z8+Z9)/(6.0*CELL_SIZE)-(Z4+Z5+Z6)/(3.0*CELL_SIZE)
                FF=(Z1-Z3-Z7+Z9)/(4.0*CELL_SIZE)
                GG=(-Z1+Z3-Z4+Z6-Z7+Z9)/(6.0*CELL_SIZE)
                HH=(-Z1-Z2-Z3+Z7+Z8+Z9)/(6.0*CELL_SIZE)
                DNM=GG*GG+HH*HH
                SLP=SQRT(DNM)
                !          div=-(z2+z4+z6+z8-4*z5)*2.0/3.0-(z1+z3+z7+z9-4.0*z5)/6.0
                DIV=-(2.0*Z1-Z2+2.0*Z3-Z4-4.0*Z5-Z6+2.0*Z7-Z8+2.0*Z9)/(3.0*CELL_SIZE)
                IF (SLP > 0.0) THEN
                    PROFCURV=-2.0*(DD*GG*GG+EE*HH*HH+FF*GG*HH)/DNM
                ELSE
                    PROFCURV=0.0
                ENDIF
                IF (SLP > 0.0) THEN
                    PLANCURV=-2.0*(DD*HH*HH+EE*GG*GG-FF*GG*HH)/DNM
                ELSE
                    PLANCURV=0.0
                ENDIF
                ALLTOT=ALLTOT+1.0
                PROFPLANSUM=PROFPLANSUM+PROFCURV*PLANCURV
                SLOPEDIVSUM=SLOPEDIVSUM+D8_GRADIENT(I,J)*DIV
                !     **********************************************************************
                !      Calculate moments of variables
                !     **********************************************************************
                CALL MOMENTS(SLOPESUMS,SLP)
                CALL MOMENTS(DIVSUMS,DIV)
                CALL MOMENTS(PLANSUMS,PLANCURV)
                CALL MOMENTS(PROFSUMS,PROFCURV)
                !     **********************************************************************
                !      Calculate catagories of slope based upon planform and profile curvature
                !     **********************************************************************
                IF ((PROFCURV >= 0.0).AND.(PLANCURV >= 0.0)) &
                NPLPL=NPLPL+1.0
                IF ((PROFCURV >= 0.0).AND.(PLANCURV < 0.0)) &
                NPLMN=NPLMN+1.0
                IF ((PROFCURV < 0.0).AND.(PLANCURV >= 0.0)) &
                NMNPL=NMNPL+1.0
                IF ((PROFCURV < 0.0).AND.(PLANCURV < 0.0)) &
                NMNMN=NMNMN+1.0
                !     **********************************************************************
                !      Find sinks, summits, and saddles
                !     **********************************************************************
                CALL FIND_TOPOGRAPHIC_EXTREMA(I,J,XSUMMIT,XSINK,XSADDLE)
                IF (XSUMMIT) NSUMMIT=NSUMMIT+1.0
                IF (XSINK) NSINK=NSINK+1.0
                IF (XSADDLE) NSADDLE=NSADDLE+1.0
                !     **********************************************************************
                !      Calculate divergence statistics
                !     **********************************************************************
                RESULTANT=DIV
                IF (RESULTANT < CMIN) CMIN = RESULTANT
                IF (RESULTANT > CMAX) CMAX = RESULTANT
                CAVG = CAVG + RESULTANT
                CSD = CSD + RESULTANT*RESULTANT
                IF (ERODE_CHANNEL(I,J) < 0.0) THEN
                    CCHANNELAVG = CCHANNELAVG+RESULTANT
                    NUMB1 = NUMB1 + 1
                ELSE
                    CSLOPEAVG = CSLOPEAVG + RESULTANT
                    NUMB2 = NUMB2 + 1
                ENDIF
                !     **********************************************************************
                !      Store local values of variables for later use
                !     **********************************************************************
                CFW(I,J) = -DIV
                CFN(I,J)=PROFCURV
                CFNE(I,J)=PLANCURV
                CFNW(I,J)=SLP
                !     **********************************************************************
                !      Calculate moments of ln(drainage area/local gradient)
                !     **********************************************************************
                IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                    LTANB=LOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
                ELSE
                    LTANB=0.0
                ENDIF
                CALL MOMENTS(LTBSUMS,LTANB)
                LTBDIVSUM = LTBDIVSUM+LTANB*DIV
                !     **********************************************************************
                !      Do the analysis of variation of variables with size of contributing area
                !     **********************************************************************
                XSKIP=.FALSE.
                DO K=1,12
                    IF (DRAINAGE_AREA(I,J) < LBOUND(K+1)) THEN
                        NUMBCAT(K)=NUMBCAT(K)+1
                        DIVCAT(K)=DIVCAT(K)+DIV
                        DIVSQCAT(K)=DIVSQCAT(K)+DIV*DIV
                        SLPCAT(K)=SLPCAT(K)+D8_GRADIENT(I,J)
                        SLPSQCAT(K)=SLPSQCAT(K)+D8_GRADIENT(I,J)**2
                        PROFCAT(K)=PROFCAT(K)+PROFCURV
                        PROFSQCAT(K)=PROFSQCAT(K)+PROFCURV*PROFCURV
                        PLANCAT(K)=PLANCAT(K)+PLANCURV
                        PLANSQCAT(K)=PLANSQCAT(K)+PLANCURV*PLANCURV
                        LTBCAT(K)=LTBCAT(K)+LTANB
                        LTBSQCAT(K)=LTBSQCAT(K)+LTANB*LTANB
                        AVGSLPCAT(K)=AVGSLPCAT(K)+SLP
                        XSKIP=.TRUE.
                    ENDIF
                ENDDO
                IF (.NOT.XSKIP) THEN
                    NUMBCAT(13)=NUMBCAT(13)+1
                    DIVCAT(13)=DIVCAT(13)+DIV
                    DIVSQCAT(13)=DIVSQCAT(13)+DIV*DIV
                    SLPCAT(13)=SLPCAT(13)+D8_GRADIENT(I,J)
                    SLPSQCAT(13)=SLPSQCAT(13)+D8_GRADIENT(I,J)**2
                    PROFCAT(13)=PROFCAT(13)+PROFCURV
                    PROFSQCAT(13)=PROFSQCAT(13)+PROFCURV*PROFCURV
                    PLANCAT(13)=PLANCAT(13)+PLANCURV
                    PLANSQCAT(13)=PLANSQCAT(13)+PLANCURV*PLANCURV
                    LTBCAT(13)=LTBCAT(13)+LTANB
                    LTBSQCAT(13)=LTBSQCAT(13)+LTANB*LTANB
                    AVGSLPCAT(13)=AVGSLPCAT(13)+SLP
                ENDIF
            ENDDO M100
        ENDDO L100
        !     **********************************************************************
        !      The above analysis is repeated for boundary locations (edges and top) with
        !         allowances for behavior at boundaries
        !     **********************************************************************
        IF (.NOT.IS_Y_PERIODIC) THEN
            L110: DO  I=2,MX-1
                J=1
                Z1=ELEVATION(I+1,J+1)
                Z2=ELEVATION(I,J+1)
                Z3=ELEVATION(I-1,J+1)
                Z4=ELEVATION(I-1,J)
                Z5=ELEVATION(I,J)
                Z6=ELEVATION(I+1,J)
                Z7=ELEVATION(I-1,J+1)
                Z8=ELEVATION(I,J+1)
                Z9=ELEVATION(I+1,J+1)
                DD=(Z1+Z3+Z4+Z6+Z7+Z9)/6.0-(Z2+Z5+Z8)/3.0
                EE=(Z1+Z2+Z3+Z7+Z8+Z9)/6.0-(Z4+Z5+Z6)/3.0
                FF=(Z1-Z3-Z7+Z9)/4.0
                GG=(-Z1+Z3-Z4+Z6-Z7+Z9)/6.0
                HH=(-Z1-Z2-Z3+Z7+Z8+Z9)/6.0
                DNM=GG*GG+HH*HH
                SLP=SQRT(DNM)
                DIV=-(Z2+Z4+Z6+Z8-4*Z5)*2.0/3.0-(Z1+Z3+Z7+Z9-4.0*Z5)/6.0
                IF (SLP > 0.0) THEN
                    PROFCURV=-2.0*(DD*GG*GG+EE*HH*HH+FF*GG*HH)/DNM
                ELSE
                    PROFCURV=0.0
                ENDIF
                IF (SLP > 0.0) THEN
                    PLANCURV=-2.0*(DD*HH*HH+EE*GG*GG-FF*GG*HH)/DNM
                ELSE
                    PLANCURV=0.0
                ENDIF
                ALLTOT=ALLTOT+1.0
                PROFPLANSUM=PROFPLANSUM+PROFCURV*PLANCURV
                SLOPEDIVSUM=SLOPEDIVSUM+D8_GRADIENT(I,J)*DIV
                CALL MOMENTS(SLOPESUMS,SLP)
                CALL MOMENTS(DIVSUMS,DIV)
                CALL MOMENTS(PLANSUMS,PLANCURV)
                CALL MOMENTS(PROFSUMS,PROFCURV)
                IF ((PROFCURV >= 0.0).AND.(PLANCURV >= 0.0)) &
                NPLPL=NPLPL+1.0
                IF ((PROFCURV >= 0.0).AND.(PLANCURV < 0.0)) &
                NPLMN=NPLMN+1.0
                IF ((PROFCURV < 0.0).AND.(PLANCURV >= 0.0)) &
                NMNPL=NMNPL+1.0
                IF ((PROFCURV < 0.0).AND.(PLANCURV < 0.0)) &
                NMNMN=NMNMN+1.0
                CALL FIND_TOPOGRAPHIC_EXTREMA(I,J,XSUMMIT,XSINK,XSADDLE)
                IF (XSUMMIT) NSUMMIT=NSUMMIT+1.0
                IF (XSINK) NSINK=NSINK+1.0
                IF (XSADDLE) NSADDLE=NSADDLE+1.0
                RESULTANT=DIV
                IF (RESULTANT < CMIN) CMIN = RESULTANT
                IF (RESULTANT > CMAX) CMAX = RESULTANT
                CAVG = CAVG + RESULTANT
                CSD = CSD + RESULTANT*RESULTANT
                IF (ERODE_CHANNEL(I,J) < 0.0) THEN
                    CCHANNELAVG = CCHANNELAVG+RESULTANT
                    NUMB1 = NUMB1 + 1
                ELSE
                    CSLOPEAVG = CSLOPEAVG + RESULTANT
                    NUMB2 = NUMB2 + 1
                ENDIF
                CFW(I,J) = -DIV
                CFN(I,J)=PROFCURV
                CFNE(I,J)=PLANCURV
                CFNW(I,J)=SLP
                IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                    LTANB=LOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
                ELSE
                    LTANB=0.0
                ENDIF
                CALL MOMENTS(LTBSUMS,LTANB)
                LTBDIVSUM = LTBDIVSUM+LTANB*DIV
                XSKIP=.FALSE.
                DO  K=1,12
                    IF (DRAINAGE_AREA(I,J) < LBOUND(K+1)) THEN
                        NUMBCAT(K)=NUMBCAT(K)+1
                        DIVCAT(K)=DIVCAT(K)+DIV
                        DIVSQCAT(K)=DIVSQCAT(K)+DIV*DIV
                        SLPCAT(K)=SLPCAT(K)+D8_GRADIENT(I,J)
                        SLPSQCAT(K)=SLPSQCAT(K)+D8_GRADIENT(I,J)**2
                        PROFCAT(K)=PROFCAT(K)+PROFCURV
                        PROFSQCAT(K)=PROFSQCAT(K)+PROFCURV*PROFCURV
                        PLANCAT(K)=PLANCAT(K)+PLANCURV
                        PLANSQCAT(K)=PLANSQCAT(K)+PLANCURV*PLANCURV
                        LTBCAT(K)=LTBCAT(K)+LTANB
                        LTBSQCAT(K)=LTBSQCAT(K)+LTANB*LTANB
                        AVGSLPCAT(K)=AVGSLPCAT(K)+SLP
                        XSKIP=.TRUE.
                    ENDIF
                ENDDO
                IF (.NOT.XSKIP) THEN
                    NUMBCAT(13)=NUMBCAT(13)+1
                    DIVCAT(13)=DIVCAT(13)+DIV
                    DIVSQCAT(13)=DIVSQCAT(13)+DIV*DIV
                    SLPCAT(13)=SLPCAT(13)+D8_GRADIENT(I,J)
                    SLPSQCAT(13)=SLPSQCAT(13)+D8_GRADIENT(I,J)**2
                    PROFCAT(13)=PROFCAT(13)+PROFCURV
                    PROFSQCAT(13)=PROFSQCAT(13)+PROFCURV*PROFCURV
                    PLANCAT(13)=PLANCAT(13)+PLANCURV
                    PLANSQCAT(13)=PLANSQCAT(13)+PLANCURV*PLANCURV
                    LTBCAT(13)=LTBCAT(13)+LTANB
                    LTBSQCAT(13)=LTBSQCAT(13)+LTANB*LTANB
                    AVGSLPCAT(13)=AVGSLPCAT(13)+SLP
                ENDIF
                J=MY
                Z1=ELEVATION(I-1,J-1)
                Z2=ELEVATION(I,J-1)
                Z3=ELEVATION(I+1,J-1)
                Z4=ELEVATION(I-1,J)
                Z5=ELEVATION(I,J)
                Z6=ELEVATION(I+1,J)
                Z7=ELEVATION(I+1,J-1)
                Z8=ELEVATION(I,J-1)
                Z9=ELEVATION(I-1,J-1)
                DD=(Z1+Z3+Z4+Z6+Z7+Z9)/6.0-(Z2+Z5+Z8)/3.0
                EE=(Z1+Z2+Z3+Z7+Z8+Z9)/6.0-(Z4+Z5+Z6)/3.0
                FF=(Z1-Z3-Z7+Z9)/4.0
                GG=(-Z1+Z3-Z4+Z6-Z7+Z9)/6.0
                HH=(-Z1-Z2-Z3+Z7+Z8+Z9)/6.0
                DNM=GG*GG+HH*HH
                SLP=SQRT(DNM)
                DIV=-(Z2+Z4+Z6+Z8-4*Z5)*2.0/3.0-(Z1+Z3+Z7+Z9-4.0*Z5)/6.0
                IF (SLP > 0.0) THEN
                    PROFCURV=-2.0*(DD*GG*GG+EE*HH*HH+FF*GG*HH)/DNM
                ELSE
                    PROFCURV=0.0
                ENDIF
                IF (SLP > 0.0) THEN
                    PLANCURV=-2.0*(DD*HH*HH+EE*GG*GG-FF*GG*HH)/DNM
                ELSE
                    PLANCURV=0.0
                ENDIF
                ALLTOT=ALLTOT+1.0
                PROFPLANSUM=PROFPLANSUM+PROFCURV*PLANCURV
                SLOPEDIVSUM=SLOPEDIVSUM+D8_GRADIENT(I,J)*DIV
                CALL MOMENTS(SLOPESUMS,SLP)
                CALL MOMENTS(DIVSUMS,DIV)
                CALL MOMENTS(PLANSUMS,PLANCURV)
                CALL MOMENTS(PROFSUMS,PROFCURV)
                IF ((PROFCURV >= 0.0).AND.(PLANCURV >= 0.0)) &
                NPLPL=NPLPL+1.0
                IF ((PROFCURV >= 0.0).AND.(PLANCURV < 0.0))&
                NPLMN=NPLMN+1.0
                IF ((PROFCURV < 0.0).AND.(PLANCURV >= 0.0)) &
                NMNPL=NMNPL+1.0
                IF ((PROFCURV < 0.0).AND.(PLANCURV < 0.0)) &
                NMNMN=NMNMN+1.0
                CALL FIND_TOPOGRAPHIC_EXTREMA(I,J,XSUMMIT,XSINK,XSADDLE)
                IF (XSUMMIT) NSUMMIT=NSUMMIT+1.0
                IF (XSINK) NSINK=NSINK+1.0
                IF (XSADDLE) NSADDLE=NSADDLE+1.0
                RESULTANT=DIV
                IF (RESULTANT < CMIN) CMIN = RESULTANT
                IF (RESULTANT > CMAX) CMAX = RESULTANT
                CAVG = CAVG + RESULTANT
                CSD = CSD + RESULTANT*RESULTANT
                IF (ERODE_CHANNEL(I,J) < 0.0) THEN
                    CCHANNELAVG = CCHANNELAVG+RESULTANT
                    NUMB1 = NUMB1 + 1
                ELSE
                    CSLOPEAVG = CSLOPEAVG + RESULTANT
                    NUMB2 = NUMB2 + 1
                ENDIF
                CFW(I,J) = -DIV
                CFN(I,J)=PROFCURV
                CFNE(I,J)=PLANCURV
                CFNW(I,J)=SLP
                IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                    LTANB=LOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
                ELSE
                    LTANB=0.0
                ENDIF
                CALL MOMENTS(LTBSUMS,LTANB)
                LTBDIVSUM = LTBDIVSUM+LTANB*DIV
                XSKIP=.FALSE.
                DO  K=1,12
                    IF (DRAINAGE_AREA(I,J) < LBOUND(K+1)) THEN
                        NUMBCAT(K)=NUMBCAT(K)+1
                        DIVCAT(K)=DIVCAT(K)+DIV
                        DIVSQCAT(K)=DIVSQCAT(K)+DIV*DIV
                        SLPCAT(K)=SLPCAT(K)+D8_GRADIENT(I,J)
                        SLPSQCAT(K)=SLPSQCAT(K)+D8_GRADIENT(I,J)**2
                        PROFCAT(K)=PROFCAT(K)+PROFCURV
                        PROFSQCAT(K)=PROFSQCAT(K)+PROFCURV*PROFCURV
                        PLANCAT(K)=PLANCAT(K)+PLANCURV
                        PLANSQCAT(K)=PLANSQCAT(K)+PLANCURV*PLANCURV
                        LTBCAT(K)=LTBCAT(K)+LTANB
                        LTBSQCAT(K)=LTBSQCAT(K)+LTANB*LTANB
                        AVGSLPCAT(K)=AVGSLPCAT(K)+SLP
                        XSKIP=.TRUE.
                    ENDIF
                ENDDO
                IF (.NOT.XSKIP) THEN
                    NUMBCAT(13)=NUMBCAT(13)+1
                    DIVCAT(13)=DIVCAT(13)+DIV
                    DIVSQCAT(13)=DIVSQCAT(13)+DIV*DIV
                    SLPCAT(13)=SLPCAT(13)+D8_GRADIENT(I,J)
                    SLPSQCAT(13)=SLPSQCAT(13)+D8_GRADIENT(I,J)**2
                    PROFCAT(13)=PROFCAT(13)+PROFCURV
                    PROFSQCAT(13)=PROFSQCAT(13)+PROFCURV*PROFCURV
                    PLANCAT(13)=PLANCAT(13)+PLANCURV
                    PLANSQCAT(13)=PLANSQCAT(13)+PLANCURV*PLANCURV
                    LTBCAT(13)=LTBCAT(13)+LTANB
                    LTBSQCAT(13)=LTBSQCAT(13)+LTANB*LTANB
                    AVGSLPCAT(13)=AVGSLPCAT(13)+SLP
                ENDIF
            ENDDO L110
            L120: DO  J=2,MY-1
                I=1
                Z1=ELEVATION(MX,J-1)
                Z2=ELEVATION(I,J-1)
                Z3=ELEVATION(I+1,J-1)
                Z4=ELEVATION(MX,J)
                Z5=ELEVATION(I,J)
                Z6=ELEVATION(I+1,J)
                Z7=ELEVATION(MX,J+1)
                Z8=ELEVATION(I,J+1)
                Z9=ELEVATION(I+1,J+1)
                DD=(Z1+Z3+Z4+Z6+Z7+Z9)/6.0-(Z2+Z5+Z8)/3.0
                EE=(Z1+Z2+Z3+Z7+Z8+Z9)/6.0-(Z4+Z5+Z6)/3.0
                FF=(Z1-Z3-Z7+Z9)/4.0
                GG=(-Z1+Z3-Z4+Z6-Z7+Z9)/6.0
                HH=(-Z1-Z2-Z3+Z7+Z8+Z9)/6.0
                DNM=GG*GG+HH*HH
                SLP=SQRT(DNM)
                DIV=-(Z2+Z4+Z6+Z8-4*Z5)*2.0/3.0-(Z1+Z3+Z7+Z9-4.0*Z5)/6.0
                IF (SLP > 0.0) THEN
                    PROFCURV=-2.0*(DD*GG*GG+EE*HH*HH+FF*GG*HH)/DNM
                ELSE
                    PROFCURV=0.0
                ENDIF
                IF (SLP > 0.0) THEN
                    PLANCURV=-2.0*(DD*HH*HH+EE*GG*GG-FF*GG*HH)/DNM
                ELSE
                    PLANCURV=0.0
                ENDIF
                ALLTOT=ALLTOT+1.0
                PROFPLANSUM=PROFPLANSUM+PROFCURV*PLANCURV
                SLOPEDIVSUM=SLOPEDIVSUM+D8_GRADIENT(I,J)*DIV
                CALL MOMENTS(SLOPESUMS,SLP)
                CALL MOMENTS(DIVSUMS,DIV)
                CALL MOMENTS(PLANSUMS,PLANCURV)
                CALL MOMENTS(PROFSUMS,PROFCURV)
                IF ((PROFCURV >= 0.0).AND.(PLANCURV >= 0.0)) &
                NPLPL=NPLPL+1.0
                IF ((PROFCURV >= 0.0).AND.(PLANCURV < 0.0)) &
                NPLMN=NPLMN+1.0
                IF ((PROFCURV < 0.0).AND.(PLANCURV >= 0.0))&
                NMNPL=NMNPL+1.0
                IF ((PROFCURV < 0.0).AND.(PLANCURV < 0.0)) &
                NMNMN=NMNMN+1.0
                CALL FIND_TOPOGRAPHIC_EXTREMA(I,J,XSUMMIT,XSINK,XSADDLE)
                IF (XSUMMIT) NSUMMIT=NSUMMIT+1.0
                IF (XSINK) NSINK=NSINK+1.0
                IF (XSADDLE) NSADDLE=NSADDLE+1.0
                RESULTANT=DIV
                IF (RESULTANT < CMIN) CMIN = RESULTANT
                IF (RESULTANT > CMAX) CMAX = RESULTANT
                CAVG = CAVG + RESULTANT
                CSD = CSD + RESULTANT*RESULTANT
                IF (ERODE_CHANNEL(I,J) < 0.0) THEN
                    CCHANNELAVG = CCHANNELAVG+RESULTANT
                    NUMB1 = NUMB1 + 1
                ELSE
                    CSLOPEAVG = CSLOPEAVG + RESULTANT
                    NUMB2 = NUMB2 + 1
                ENDIF
                CFW(I,J) = -DIV
                CFN(I,J)=PROFCURV
                CFNE(I,J)=PLANCURV
                CFNW(I,J)=SLP
                IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                    LTANB=LOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
                ELSE
                    LTANB=0.0
                ENDIF
                CALL MOMENTS(LTBSUMS,LTANB)
                LTBDIVSUM = LTBDIVSUM+LTANB*DIV
                XSKIP=.FALSE.
                DO K=1,12
                    IF (DRAINAGE_AREA(I,J) < LBOUND(K+1)) THEN
                        NUMBCAT(K)=NUMBCAT(K)+1
                        DIVCAT(K)=DIVCAT(K)+DIV
                        DIVSQCAT(K)=DIVSQCAT(K)+DIV*DIV
                        SLPCAT(K)=SLPCAT(K)+D8_GRADIENT(I,J)
                        SLPSQCAT(K)=SLPSQCAT(K)+D8_GRADIENT(I,J)**2
                        PROFCAT(K)=PROFCAT(K)+PROFCURV
                        PROFSQCAT(K)=PROFSQCAT(K)+PROFCURV*PROFCURV
                        PLANCAT(K)=PLANCAT(K)+PLANCURV
                        PLANSQCAT(K)=PLANSQCAT(K)+PLANCURV*PLANCURV
                        LTBCAT(K)=LTBCAT(K)+LTANB
                        LTBSQCAT(K)=LTBSQCAT(K)+LTANB*LTANB
                        AVGSLPCAT(K)=AVGSLPCAT(K)+SLP
                        XSKIP=.TRUE. 
                    ENDIF
                ENDDO
                IF (.NOT.XSKIP) THEN
                    NUMBCAT(13)=NUMBCAT(13)+1
                    DIVCAT(13)=DIVCAT(13)+DIV
                    DIVSQCAT(13)=DIVSQCAT(13)+DIV*DIV
                    SLPCAT(13)=SLPCAT(13)+D8_GRADIENT(I,J)
                    SLPSQCAT(13)=SLPSQCAT(13)+D8_GRADIENT(I,J)**2
                    PROFCAT(13)=PROFCAT(13)+PROFCURV
                    PROFSQCAT(13)=PROFSQCAT(13)+PROFCURV*PROFCURV
                    PLANCAT(13)=PLANCAT(13)+PLANCURV
                    PLANSQCAT(13)=PLANSQCAT(13)+PLANCURV*PLANCURV
                    LTBCAT(13)=LTBCAT(13)+LTANB
                    LTBSQCAT(13)=LTBSQCAT(13)+LTANB*LTANB
                    AVGSLPCAT(13)=AVGSLPCAT(13)+SLP
                ENDIF
                I=MX
                Z1=ELEVATION(I-1,J-1)
                Z2=ELEVATION(I,J-1)
                Z3=ELEVATION(1,J-1)
                Z4=ELEVATION(I-1,J)
                Z5=ELEVATION(I,J)
                Z6=ELEVATION(1,J)
                Z7=ELEVATION(I-1,J+1)
                Z8=ELEVATION(I,J+1)
                Z9=ELEVATION(1,J+1)
                DD=(Z1+Z3+Z4+Z6+Z7+Z9)/6.0-(Z2+Z5+Z8)/3.0
                EE=(Z1+Z2+Z3+Z7+Z8+Z9)/6.0-(Z4+Z5+Z6)/3.0
                FF=(Z1-Z3-Z7+Z9)/4.0
                GG=(-Z1+Z3-Z4+Z6-Z7+Z9)/6.0
                HH=(-Z1-Z2-Z3+Z7+Z8+Z9)/6.0
                DNM=GG*GG+HH*HH
                SLP=SQRT(DNM)
                DIV=-(Z2+Z4+Z6+Z8-4*Z5)*2.0/3.0-(Z1+Z3+Z7+Z9-4.0*Z5)/6.0
                IF (SLP > 0.0) THEN
                    PROFCURV=-2.0*(DD*GG*GG+EE*HH*HH+FF*GG*HH)/DNM
                ELSE
                    PROFCURV=0.0
                ENDIF
                IF (SLP > 0.0) THEN
                    PLANCURV=-2.0*(DD*HH*HH+EE*GG*GG-FF*GG*HH)/DNM
                ELSE
                    PLANCURV=0.0
                ENDIF
                ALLTOT=ALLTOT+1.0
                PROFPLANSUM=PROFPLANSUM+PROFCURV*PLANCURV
                SLOPEDIVSUM=SLOPEDIVSUM+D8_GRADIENT(I,J)*DIV
                CALL MOMENTS(SLOPESUMS,SLP)
                CALL MOMENTS(DIVSUMS,DIV)
                CALL MOMENTS(PLANSUMS,PLANCURV)
                CALL MOMENTS(PROFSUMS,PROFCURV)
                IF ((PROFCURV >= 0.0).AND.(PLANCURV >= 0.0)) &
                NPLPL=NPLPL+1.0
                IF ((PROFCURV >= 0.0).AND.(PLANCURV < 0.0)) &
                NPLMN=NPLMN+1.0
                IF ((PROFCURV < 0.0).AND.(PLANCURV >= 0.0)) &
                NMNPL=NMNPL+1.0
                IF ((PROFCURV < 0.0).AND.(PLANCURV < 0.0)) &
                NMNMN=NMNMN+1.0
                CALL FIND_TOPOGRAPHIC_EXTREMA(I,J,XSUMMIT,XSINK,XSADDLE)
                IF (XSUMMIT) NSUMMIT=NSUMMIT+1.0
                IF (XSINK) NSINK=NSINK+1.0
                IF (XSADDLE) NSADDLE=NSADDLE+1.0
                RESULTANT=DIV
                IF (RESULTANT < CMIN) CMIN = RESULTANT
                IF (RESULTANT > CMAX) CMAX = RESULTANT
                CAVG = CAVG + RESULTANT
                CSD = CSD + RESULTANT*RESULTANT
                IF (ERODE_CHANNEL(I,J) < 0.0) THEN
                    CCHANNELAVG = CCHANNELAVG+RESULTANT
                    NUMB1 = NUMB1 + 1
                ELSE
                    CSLOPEAVG = CSLOPEAVG + RESULTANT
                    NUMB2 = NUMB2 + 1
                ENDIF
                CFW(I,J) = -DIV
                CFN(I,J)=PROFCURV
                CFNE(I,J)=PLANCURV
                CFNW(I,J)=SLP
                IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                    LTANB=LOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
                ELSE
                    LTANB=0.0
                ENDIF
                CALL MOMENTS(LTBSUMS,LTANB)
                LTBDIVSUM = LTBDIVSUM+LTANB*DIV
                XSKIP=.FALSE.
                DO  K=1,12
                    IF (DRAINAGE_AREA(I,J) < LBOUND(K+1)) THEN
                        NUMBCAT(K)=NUMBCAT(K)+1
                        DIVCAT(K)=DIVCAT(K)+DIV
                        DIVSQCAT(K)=DIVSQCAT(K)+DIV*DIV
                        SLPCAT(K)=SLPCAT(K)+D8_GRADIENT(I,J)
                        SLPSQCAT(K)=SLPSQCAT(K)+D8_GRADIENT(I,J)**2
                        PROFCAT(K)=PROFCAT(K)+PROFCURV
                        PROFSQCAT(K)=PROFSQCAT(K)+PROFCURV*PROFCURV
                        PLANCAT(K)=PLANCAT(K)+PLANCURV
                        PLANSQCAT(K)=PLANSQCAT(K)+PLANCURV*PLANCURV
                        LTBCAT(K)=LTBCAT(K)+LTANB
                        LTBSQCAT(K)=LTBSQCAT(K)+LTANB*LTANB
                        AVGSLPCAT(K)=AVGSLPCAT(K)+SLP
                        XSKIP=.TRUE. 
                    ENDIF
                ENDDO
                IF (.NOT.XSKIP) THEN
                    NUMBCAT(13)=NUMBCAT(13)+1
                    DIVCAT(13)=DIVCAT(13)+DIV
                    DIVSQCAT(13)=DIVSQCAT(13)+DIV*DIV
                    SLPCAT(13)=SLPCAT(13)+D8_GRADIENT(I,J)
                    SLPSQCAT(13)=SLPSQCAT(13)+D8_GRADIENT(I,J)**2
                    PROFCAT(13)=PROFCAT(13)+PROFCURV
                    PROFSQCAT(13)=PROFSQCAT(13)+PROFCURV*PROFCURV
                    PLANCAT(13)=PLANCAT(13)+PLANCURV
                    PLANSQCAT(13)=PLANSQCAT(13)+PLANCURV*PLANCURV
                    LTBCAT(13)=LTBCAT(13)+LTANB
                    LTBSQCAT(13)=LTBSQCAT(13)+LTANB*LTANB
                    AVGSLPCAT(13)=AVGSLPCAT(13)+SLP
                ENDIF
            ENDDO L120
            I=1
            J=1
            Z1=ELEVATION(I+1,J+1)
            Z2=ELEVATION(I,J+1)
            Z3=ELEVATION(MX,J+1)
            Z4=ELEVATION(MX,J)
            Z5=ELEVATION(I,J)
            Z6=ELEVATION(I+1,J)
            Z7=ELEVATION(MX,J+1)
            Z8=ELEVATION(I,J+1)
            Z9=ELEVATION(I+1,J+1)
            DD=(Z1+Z3+Z4+Z6+Z7+Z9)/6.0-(Z2+Z5+Z8)/3.0
            EE=(Z1+Z2+Z3+Z7+Z8+Z9)/6.0-(Z4+Z5+Z6)/3.0
            FF=(Z1-Z3-Z7+Z9)/4.0
            GG=(-Z1+Z3-Z4+Z6-Z7+Z9)/6.0
            HH=(-Z1-Z2-Z3+Z7+Z8+Z9)/6.0
            DNM=GG*GG+HH*HH
            SLP=SQRT(DNM)
            DIV=-(Z2+Z4+Z6+Z8-4*Z5)*2.0/3.0-(Z1+Z3+Z7+Z9-4.0*Z5)/6.0
            IF (SLP > 0.0) THEN
                PROFCURV=-2.0*(DD*GG*GG+EE*HH*HH+FF*GG*HH)/DNM
            ELSE
                PROFCURV=0.0
            ENDIF
            IF (SLP > 0.0) THEN
                PLANCURV=-2.0*(DD*HH*HH+EE*GG*GG-FF*GG*HH)/DNM
            ELSE
                PLANCURV=0.0
            ENDIF
            ALLTOT=ALLTOT+1.0
            PROFPLANSUM=PROFPLANSUM+PROFCURV*PLANCURV
            SLOPEDIVSUM=SLOPEDIVSUM+SLP*DIV
            CALL MOMENTS(SLOPESUMS,SLP)
            CALL MOMENTS(DIVSUMS,DIV)
            CALL MOMENTS(PLANSUMS,PLANCURV)
            CALL MOMENTS(PROFSUMS,PROFCURV)
            IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                LTANB=LOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
            ELSE
                LTANB=0.0
            ENDIF
            CALL MOMENTS(LTBSUMS,LTANB)
            LTBDIVSUM = LTBDIVSUM+LTANB*DIV
            IF ((PROFCURV >= 0.0).AND.(PLANCURV >= 0.0)) &
            NPLPL=NPLPL+1.0
            IF ((PROFCURV >= 0.0).AND.(PLANCURV < 0.0)) &
            NPLMN=NPLMN+1.0
            IF ((PROFCURV < 0.0).AND.(PLANCURV >= 0.0))  &
            NMNPL=NMNPL+1.0
            IF ((PROFCURV < 0.0).AND.(PLANCURV < 0.0)) &
            NMNMN=NMNMN+1.0
            CALL FIND_TOPOGRAPHIC_EXTREMA(I,J,XSUMMIT,XSINK,XSADDLE)
            IF (XSUMMIT) NSUMMIT=NSUMMIT+1.0
            IF (XSINK) NSINK=NSINK+1.0
            IF (XSADDLE) NSADDLE=NSADDLE+1.0
            RESULTANT=DIV
            IF (RESULTANT < CMIN) CMIN = RESULTANT
            IF (RESULTANT > CMAX) CMAX = RESULTANT
            CAVG = CAVG + RESULTANT
            CSD = CSD + RESULTANT*RESULTANT
            IF (ERODE_CHANNEL(I,J) < 0.0) THEN
                CCHANNELAVG = CCHANNELAVG+RESULTANT
                NUMB1 = NUMB1 + 1
            ELSE
                CSLOPEAVG = CSLOPEAVG + RESULTANT
                NUMB2 = NUMB2 + 1
            ENDIF
            CFW(I,J) = -DIV
            CFN(I,J)=PROFCURV
            CFNE(I,J)=PLANCURV
            CFNW(I,J)=SLP
            I=MX
            J=1
            Z1=ELEVATION(1,J+1)
            Z2=ELEVATION(I,J+1)
            Z3=ELEVATION(I-1,J+1)
            Z4=ELEVATION(I-1,J)
            Z5=ELEVATION(I,J)
            Z6=ELEVATION(1,J)
            Z7=ELEVATION(I-1,J+1)
            Z8=ELEVATION(I,J+1)
            Z9=ELEVATION(1,J+1)
            DD=(Z1+Z3+Z4+Z6+Z7+Z9)/6.0-(Z2+Z5+Z8)/3.0
            EE=(Z1+Z2+Z3+Z7+Z8+Z9)/6.0-(Z4+Z5+Z6)/3.0
            FF=(Z1-Z3-Z7+Z9)/4.0
            GG=(-Z1+Z3-Z4+Z6-Z7+Z9)/6.0
            HH=(-Z1-Z2-Z3+Z7+Z8+Z9)/6.0
            DNM=GG*GG+HH*HH
            SLP=SQRT(DNM)
            DIV=-(Z2+Z4+Z6+Z8-4*Z5)*2.0/3.0-(Z1+Z3+Z7+Z9-4.0*Z5)/6.0
            IF (SLP > 0.0) THEN
                PROFCURV=-2.0*(DD*GG*GG+EE*HH*HH+FF*GG*HH)/DNM
            ELSE
                PROFCURV=0.0
            ENDIF
            IF (SLP > 0.0) THEN
                PLANCURV=-2.0*(DD*HH*HH+EE*GG*GG-FF*GG*HH)/DNM
            ELSE
                PLANCURV=0.0
            ENDIF
            ALLTOT=ALLTOT+1.0
            PROFPLANSUM=PROFPLANSUM+PROFCURV*PLANCURV
            SLOPEDIVSUM=SLOPEDIVSUM+SLP*DIV
            CALL MOMENTS(SLOPESUMS,SLP)
            CALL MOMENTS(DIVSUMS,DIV)
            CALL MOMENTS(PLANSUMS,PLANCURV)
            CALL MOMENTS(PROFSUMS,PROFCURV)
            IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                LTANB=LOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
            ELSE
                LTANB=0.0
            ENDIF
            CALL MOMENTS(LTBSUMS,LTANB)
            LTBDIVSUM = LTBDIVSUM+LTANB*DIV
            IF ((PROFCURV >= 0.0).AND.(PLANCURV >= 0.0)) &
            NPLPL=NPLPL+1.0
            IF ((PROFCURV >= 0.0).AND.(PLANCURV < 0.0)) &
            NPLMN=NPLMN+1.0
            IF ((PROFCURV < 0.0).AND.(PLANCURV >= 0.0)) &
            NMNPL=NMNPL+1.0
            IF ((PROFCURV < 0.0).AND.(PLANCURV < 0.0)) &
            NMNMN=NMNMN+1.0
            CALL FIND_TOPOGRAPHIC_EXTREMA(I,J,XSUMMIT,XSINK,XSADDLE)
            IF (XSUMMIT) NSUMMIT=NSUMMIT+1.0
            IF (XSINK) NSINK=NSINK+1.0
            IF (XSADDLE) NSADDLE=NSADDLE+1.0
            RESULTANT=DIV
            IF (RESULTANT < CMIN) CMIN = RESULTANT
            IF (RESULTANT > CMAX) CMAX = RESULTANT
            CAVG = CAVG + RESULTANT
            CSD = CSD + RESULTANT*RESULTANT
            IF (ERODE_CHANNEL(I,J) < 0.0) THEN
                CCHANNELAVG = CCHANNELAVG+RESULTANT
                NUMB1 = NUMB1 + 1
            ELSE
                CSLOPEAVG = CSLOPEAVG + RESULTANT
                NUMB2 = NUMB2 + 1
            ENDIF
            CFW(I,J) = -DIV
            CFN(I,J)=PROFCURV
            CFNE(I,J)=PLANCURV
            CFNW(I,J)=SLP
            I=1
            J=MY
            Z1=ELEVATION(MX,J-1)
            Z2=ELEVATION(I,J-1)
            Z3=ELEVATION(I+1,J-1)
            Z4=ELEVATION(MX,J)
            Z5=ELEVATION(I,J)
            Z6=ELEVATION(I+1,J)
            Z7=ELEVATION(I+1,J-1)
            Z8=ELEVATION(I,J-1)
            Z9=ELEVATION(MX,J-1)
            DD=(Z1+Z3+Z4+Z6+Z7+Z9)/6.0-(Z2+Z5+Z8)/3.0
            EE=(Z1+Z2+Z3+Z7+Z8+Z9)/6.0-(Z4+Z5+Z6)/3.0
            FF=(Z1-Z3-Z7+Z9)/4.0
            GG=(-Z1+Z3-Z4+Z6-Z7+Z9)/6.0
            HH=(-Z1-Z2-Z3+Z7+Z8+Z9)/6.0
            DNM=GG*GG+HH*HH
            SLP=SQRT(DNM)
            DIV=-(Z2+Z4+Z6+Z8-4*Z5)*2.0/3.0-(Z1+Z3+Z7+Z9-4.0*Z5)/6.0
            IF (SLP > 0.0) THEN
                PROFCURV=-2.0*(DD*GG*GG+EE*HH*HH+FF*GG*HH)/DNM
            ELSE
                PROFCURV=0.0
            ENDIF
            IF (SLP > 0.0) THEN
                PLANCURV=-2.0*(DD*HH*HH+EE*GG*GG-FF*GG*HH)/DNM
            ELSE
                PLANCURV=0.0
            ENDIF
            ALLTOT=ALLTOT+1.0
            PROFPLANSUM=PROFPLANSUM+PROFCURV*PLANCURV
            SLOPEDIVSUM=SLOPEDIVSUM+SLP*DIV
            CALL MOMENTS(SLOPESUMS,SLP)
            CALL MOMENTS(DIVSUMS,DIV)
            CALL MOMENTS(PLANSUMS,PLANCURV)
            CALL MOMENTS(PROFSUMS,PROFCURV)
            IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                LTANB=LOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
            ELSE
                LTANB=0.0
            ENDIF
            CALL MOMENTS(LTBSUMS,LTANB)
            LTBDIVSUM = LTBDIVSUM+LTANB*DIV
            IF ((PROFCURV >= 0.0).AND.(PLANCURV >= 0.0))  &
            NPLPL=NPLPL+1.0
            IF ((PROFCURV >= 0.0).AND.(PLANCURV < 0.0)) &
            NPLMN=NPLMN+1.0
            IF ((PROFCURV < 0.0).AND.(PLANCURV >= 0.0)) &
            NMNPL=NMNPL+1.0
            IF ((PROFCURV < 0.0).AND.(PLANCURV < 0.0)) &
            NMNMN=NMNMN+1.0
            CALL FIND_TOPOGRAPHIC_EXTREMA(I,J,XSUMMIT,XSINK,XSADDLE)
            IF (XSUMMIT) NSUMMIT=NSUMMIT+1.0
            IF (XSINK) NSINK=NSINK+1.0
            IF (XSADDLE) NSADDLE=NSADDLE+1.0
            RESULTANT=DIV
            IF (RESULTANT < CMIN) CMIN = RESULTANT
            IF (RESULTANT > CMAX) CMAX = RESULTANT
            CAVG = CAVG + RESULTANT
            CSD = CSD + RESULTANT*RESULTANT
            IF (ERODE_CHANNEL(I,J) < 0.0) THEN
                CCHANNELAVG = CCHANNELAVG+RESULTANT
                NUMB1 = NUMB1 + 1
            ELSE
                CSLOPEAVG = CSLOPEAVG + RESULTANT
                NUMB2 = NUMB2 + 1
            ENDIF
            CFW(I,J) = -DIV
            CFN(I,J)=PROFCURV
            CFNE(I,J)=PLANCURV
            CFNW(I,J)=SLP
            I=MX
            J=MY
            Z1=ELEVATION(I-1,J-1)
            Z2=ELEVATION(I,J-1)
            Z3=ELEVATION(1,J-1)
            Z4=ELEVATION(I-1,J)
            Z5=ELEVATION(I,J)
            Z6=ELEVATION(1,J)
            Z7=ELEVATION(1,J-1)
            Z8=ELEVATION(I,J-1)
            Z9=ELEVATION(I-1,J-1)
            DD=(Z1+Z3+Z4+Z6+Z7+Z9)/6.0-(Z2+Z5+Z8)/3.0
            EE=(Z1+Z2+Z3+Z7+Z8+Z9)/6.0-(Z4+Z5+Z6)/3.0
            FF=(Z1-Z3-Z7+Z9)/4.0
            GG=(-Z1+Z3-Z4+Z6-Z7+Z9)/6.0
            HH=(-Z1-Z2-Z3+Z7+Z8+Z9)/6.0
            DNM=GG*GG+HH*HH
            SLP=SQRT(DNM)
            DIV=-(Z2+Z4+Z6+Z8-4*Z5)*2.0/3.0-(Z1+Z3+Z7+Z9-4.0*Z5)/6.0
            IF (SLP > 0.0) THEN
                PROFCURV=-2.0*(DD*GG*GG+EE*HH*HH+FF*GG*HH)/DNM
            ELSE
                PROFCURV=0.0
            ENDIF
            IF (SLP > 0.0) THEN
                PLANCURV=-2.0*(DD*HH*HH+EE*GG*GG-FF*GG*HH)/DNM
            ELSE
                PLANCURV=0.0
            ENDIF
            ALLTOT=ALLTOT+1.0
            PROFPLANSUM=PROFPLANSUM+PROFCURV*PLANCURV
            SLOPEDIVSUM=SLOPEDIVSUM+SLP*DIV
            CALL MOMENTS(SLOPESUMS,SLP)
            CALL MOMENTS(DIVSUMS,DIV)
            CALL MOMENTS(PLANSUMS,PLANCURV)
            CALL MOMENTS(PROFSUMS,PROFCURV)
            IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                LTANB=LOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
            ELSE
                LTANB=0.0
            ENDIF
            CALL MOMENTS(LTBSUMS,LTANB)
            LTBDIVSUM = LTBDIVSUM+LTANB*DIV
            IF ((PROFCURV >= 0.0).AND.(PLANCURV >= 0.0)) &
            NPLPL=NPLPL+1.0
            IF ((PROFCURV >= 0.0).AND.(PLANCURV < 0.0)) &
            NPLMN=NPLMN+1.0
            IF ((PROFCURV < 0.0).AND.(PLANCURV >= 0.0)) &
            NMNPL=NMNPL+1.0
            IF ((PROFCURV < 0.0).AND.(PLANCURV < 0.0)) &
            NMNMN=NMNMN+1.0
            CALL FIND_TOPOGRAPHIC_EXTREMA(I,J,XSUMMIT,XSINK,XSADDLE)
            IF (XSUMMIT) NSUMMIT=NSUMMIT+1.0
            IF (XSINK) NSINK=NSINK+1.0
            IF (XSADDLE) NSADDLE=NSADDLE+1.0
            RESULTANT=DIV
            IF (RESULTANT < CMIN) CMIN = RESULTANT
            IF (RESULTANT > CMAX) CMAX = RESULTANT
            CAVG = CAVG + RESULTANT
            CSD = CSD + RESULTANT*RESULTANT
            IF (ERODE_CHANNEL(I,J) < 0.0) THEN
                CCHANNELAVG = CCHANNELAVG+RESULTANT
                NUMB1 = NUMB1 + 1
            ELSE
                CSLOPEAVG = CSLOPEAVG + RESULTANT
                NUMB2 = NUMB2 + 1
            ENDIF
            CFW(I,J) = -DIV
            CFN(I,J)=PROFCURV
            CFNE(I,J)=PLANCURV
            CFNW(I,J)=SLP
        ENDIF
        !     **********************************************************************
        !      All of the raw variables have now been calculated, and the statistics and
        !         correlations are computed and printed - See format statements for
        !         explanation
        !     **********************************************************************
        SUMXY=PROFPLANSUM
        SUMX=PROFSUMS(1)
        SUMY=PLANSUMS(1)
        SUMXX=PROFSUMS(2)
        SUMYY=PLANSUMS(2)
        TOTAL=PLANSUMS(5)
        IF (TOTAL > 1.0) THEN
            CORREL=(SUMXY-SUMX*SUMY/TOTAL)/(SQRT((SUMXX-SUMX*SUMX/TOTAL) &
            )*SQRT((SUMYY-SUMY*SUMY/TOTAL)))
        ELSE
            CORREL=0.0
        ENDIF
        WRITE(OUTHIST,750)CORREL,TOTAL
        WRITE(OUTSUMMARY,906)CORREL,TOTAL
        750 FORMAT(' CORRELATION BETWEEN HORIZ. AND DOWNSL. CURV.=',G12.5  &
        ,' N=',G12.5)
        SUMXY=SLOPEDIVSUM
        SUMX=SLOPESUMS(1)
        SUMY=DIVSUMS(1)
        SUMXX=SLOPESUMS(2)
        SUMYY=DIVSUMS(2)
        TOTAL=DIVSUMS(5)
        IF (TOTAL > 1.0) THEN
            CORREL=(SUMXY-SUMX*SUMY/TOTAL)/(SQRT((SUMXX-SUMX*SUMX/TOTAL) &
            )*SQRT((SUMYY-SUMY*SUMY/TOTAL)))
        ELSE
            CORREL=0.0
        ENDIF
        WRITE(OUTHIST,751)CORREL,TOTAL
        WRITE(OUTSUMMARY,906)CORREL,TOTAL
        751 FORMAT(' CORRELATION BETWEEN GRADIENT AND DIVERGENCE=',G12.5, &
        ' N=',G13.6)
        SUMXY=LTBDIVSUM
        SUMX=LTBSUMS(1)
        SUMY=DIVSUMS(1)
        SUMXX=LTBSUMS(2)
        SUMYY=DIVSUMS(2)
        TOTAL=DIVSUMS(5)
        IF (TOTAL > 1.0) THEN
            CORREL=(SUMXY-SUMX*SUMY/TOTAL)/(SQRT((SUMXX-SUMX*SUMX/TOTAL) &
            )*SQRT((SUMYY-SUMY*SUMY/TOTAL)))
        ELSE
            CORREL=0.0
        ENDIF
        WRITE(OUTHIST,745)CORREL,TOTAL
        WRITE(OUTSUMMARY,906)CORREL,TOTAL
        745 FORMAT(' CORRELATION BETWEEN LN(A/S) AND DIVERGENCE=',G12.5,&
        ' N=',G13.6)
        CALL CALCULATE_MOMENTS(SLOPESUMS,AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTHIST,740)
        740 FORMAT(' SLOPE GRADIENT MOMENTS')
        CALL PRINT_MOMENTS(AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTSUMMARY,912)AVG,VARI,SKEW,KURT,TOTAL
        CALL CALCULATE_MOMENTS(DIVSUMS,AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTHIST,741)
        741 FORMAT(' SLOPE DIVERGENCE MOMENTS')
        CALL PRINT_MOMENTS(AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTSUMMARY,912)AVG,VARI,SKEW,KURT,TOTAL
        CALL CALCULATE_MOMENTS(PLANSUMS,AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTHIST,742)
        742 FORMAT(' PLAN (HORIZONTAL) CURVATURE MOMENTS')
        CALL PRINT_MOMENTS(AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTSUMMARY,912)AVG,VARI,SKEW,KURT,TOTAL
        CALL CALCULATE_MOMENTS(PROFSUMS,AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTHIST,743)
        743 FORMAT(' VERTICAL (DOWNSLOPE) CURVATURE MOMENTS')
        CALL PRINT_MOMENTS(AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTSUMMARY,912)AVG,VARI,SKEW,KURT,TOTAL
        CALL CALCULATE_MOMENTS(LTBSUMS,AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTHIST,744)
        744 FORMAT(' LN(A/S) MOMENTS')
        CALL PRINT_MOMENTS(AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTSUMMARY,912)AVG,VARI,SKEW,KURT,TOTAL
        NTOT = (MX-2)*(MY-2)
        TEMP1=CAVG/NTOT
        WRITE(OUTHIST,500) TEMP1
        WRITE(OUTSUMMARY,907) TEMP1
        500 FORMAT( ' AVERAGE DIVERGENCE=',G12.5)
        TEMP1=SQRT((CSD-CAVG*CAVG/NTOT)/(NTOT-1.0))
        WRITE(OUTHIST,501) TEMP1
        WRITE(OUTSUMMARY,907) TEMP1
        501 FORMAT( ' S.D. OF DIVERGENCE=',G12.5)
        TEMP1=   CCHANNELAVG/NUMB1
        WRITE(OUTHIST,502) TEMP1
        WRITE(OUTSUMMARY,907) TEMP1
        502 FORMAT( ' AVERAGE CHANNEL DIVERGENCE=',G12.5)
        TEMP1=CSLOPEAVG/NUMB2
        WRITE(OUTHIST,503) TEMP1
        WRITE(OUTSUMMARY,907) TEMP1
        503 FORMAT( ' AVERAGE SLOPE DIVERGENCE=',G12.5)
        WRITE(OUTHIST,504) CMAX,CMIN
        WRITE(OUTSUMMARY,906) CMAX,CMIN
        504 FORMAT( ' MAXIMUM.AND.MINIMUM DIVERGENCE=',2(G12.5,' '),/)
        TEMP1=NSUMMIT/NTOT
        WRITE(OUTHIST,755) TEMP1
        WRITE(OUTSUMMARY,907) TEMP1
        755 FORMAT(' SUMMIT DENSITY=',E15.8)
        TEMP1=NSINK/NTOT
        WRITE(OUTHIST,756) TEMP1
        WRITE(OUTSUMMARY,907) TEMP1
        756 FORMAT(' SINK DENSITY=',E15.8)
        TEMP1=NSADDLE/NTOT
        WRITE(OUTHIST,757) TEMP1
        WRITE(OUTSUMMARY,907) TEMP1
        757 FORMAT(' SADDLE DENSITY=',E15.8)
        TEMP1=NPLPL/NTOT
        WRITE(OUTSUMMARY,907) TEMP1
        WRITE(OUTHIST,758) TEMP1
        758 FORMAT(' DENSITY OF PROF+ PLAN+ POINTS=',E15.8)
        TEMP1=NPLMN/NTOT
        WRITE(OUTHIST,759) TEMP1
        WRITE(OUTSUMMARY,907) TEMP1
        759 FORMAT(' DENSITY OF PROF+ PLAN- POINTS=',E15.8)
        TEMP1=NMNPL/NTOT
        WRITE(OUTHIST,760) TEMP1
        WRITE(OUTSUMMARY,907) TEMP1
        760 FORMAT(' DENSITY OF PROF- PLAN+ POINTS=',E15.8)
        TEMP1=NMNMN/NTOT
        WRITE(OUTHIST,761) TEMP1
        WRITE(OUTSUMMARY,907) TEMP1
        761 FORMAT(' DENSITY OF PROF- PLAN- POINTS=',E15.8)
        WRITE (OUTHIST,810)
        810 FORMAT(/,' SUMMARY OF SLOPE CHARACTERISTICS BY SUPPORT AREA' &
        ,/,' AREA LOWER BOUND, NUMBER OF OBS',/,  &
        ' AVG AND VAR OF SLOPE, DIVERG., HORIZ. CURV, PROFILE CURV',/)
        L811: DO K=1,13
            Z1=NUMBCAT(K)
            IF (Z1 > 1.0) THEN
                Z2=SLPCAT(K)/Z1
                Z3=(SLPSQCAT(K)-Z1*Z2*Z2)/(Z1-1.0)
                Z4=DIVCAT(K)/Z1
                Z5=(DIVSQCAT(K)-Z1*Z4*Z4)/(Z1-1.0)
                Z6=PLANCAT(K)/Z1
                Z7=(PLANSQCAT(K)-Z1*Z6*Z6)/(Z1-1.0)
                Z8=PROFCAT(K)/Z1
                Z9=(PROFSQCAT(K)-Z1*Z8*Z8)/(Z1-1.0)
                Z10=LTBCAT(K)/Z1
                Z11=(LTBSQCAT(K)-Z1*Z10*Z10)/(Z1-1.0)
                Z12=AVGSLPCAT(K)/Z1
            ELSE
                Z2=0.0
                Z3=0.0
                Z4=0.0
                Z5=0.0
                Z6=0.0
                Z7=0.0
                Z8=0.0
                Z9=0.0
                Z10=0.0
                Z11=0.0
                Z12=0.0
            ENDIF
            WRITE(OUTHIST,812)LBOUND(K),Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10,Z11 &
            ,Z12
            WRITE(OUTSUMMARY,911)LBOUND(K),Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9, &
            Z10,Z11,Z12
            812 FORMAT(E13.5,3E13.5,/,6E13.5,/,3E13.5)
        ENDDO L811
        CALL PRINT_VARIABLE_SUMMARIES()
        906 FORMAT(2(' ',E15.8))
        907 FORMAT(' ',E15.8)
        908 FORMAT(2(' ',I5),3(' ',E15.8))
        909 FORMAT(2(' ',I5),2(' ',E15.8))
        910 FORMAT(' ',E15.8,' ',I5)
        911 FORMAT(13(' ',E13.5))
        912 FORMAT(5(' ',E14.7))
        913 FORMAT(6(' ',E12.5))
        914 FORMAT(3(' ',E15.8))
        RETURN
    END ! SUBROUTINE CALCULATE_TOPO_DIVERGENCE
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PRINT_MORPHOMETRY()
        !      *********************************************************************
        !       This is the main routine for printing out various summary statistics
        !         for morphology in simulations.
        !       The following variables are printed out in this routine:
        !              average and maximum drainage area
        !              maximum and minimum elevation
        !              minimum, average and maximum sediment thickness
        !              average relative elevation
        !              average gradient
        !              maximum gradient
        !      *********************************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I, J, NCC,NS,NB,NR
        REAL(4) :: AM
        !          integer am
        REAL(4) :: EP, EM, EV, AT, AG,MG
        REAL(4) :: RTHICK,RT,RTMIN,RTMAX
        REAL(4) :: STHICK,ST,STMIN,STMAX
        REAL(4) ::  TEMP1,TEMP2,TEMP3
        EM =-1.0E+20
        MG =-1.0E+20
        EP =1.0E+20
        ST = 0.0
        STMAX =-1.0E+20
        STMIN =1.0E+20
        NB=0
        NR=0
        RT=0.0
        RTMIN=STMIN
        RTMAX=STMAX
        NCC=0
        NS=0
        EV = 0.0
        AT = 0.0
        AM = 0
        AG = 0.0
        DO  J=1,MY
            DO  I=1,MX
                IF (ELEVATION(I,J)  <  EP) EP = ELEVATION(I,J)
            ENDDO
        ENDDO
        DO  J=1,MYY
            DO  I=1,MX
                IF (IS_ROCK_SURFACE(I,J)) THEN
                    NB=NB+1
                ELSE
                    NR=NR+1
                    RTHICK=REGOLITH(I,J)
                    RT=RT+RTHICK
                    IF (RTHICK > RTMAX) RTMAX=RTHICK
                    IF (RTHICK < RTMIN) RTMIN=RTHICK
                ENDIF
                NS = NS+1
                STHICK=ELEVATION(I,J)-SEDIMENT_BASE(I,J)
                ST=ST+STHICK
                IF (STHICK > STMAX) STMAX=STHICK
                IF (STHICK < STMIN) STMIN=STHICK
                AT = AT + DRAINAGE_AREA(I,J)
                IF (DRAINAGE_AREA(I,J)  >  AM) AM = DRAINAGE_AREA(I,J)
                EV = EV + ELEVATION(I,J) - EP
                IF (ELEVATION(I,J)  >  EM) EM = ELEVATION(I,J)
                AG = AG + D8_GRADIENT(I,J)
                IF (D8_GRADIENT(I,J)   >  MG) MG = D8_GRADIENT(I,J)
            ENDDO
        ENDDO
        TEMP1 = AT/(MX*MY)
        WRITE(OUTHIST,504) TEMP1,AM
        WRITE(OUTSUMMARY,910) TEMP1,AM
        504 FORMAT( 'AVERAGE.AND.MAXIMUM AREA =',G12.5, ' , ',G12.5)
        WRITE(OUTHIST,505) EP,EM
        WRITE(OUTSUMMARY,906) EP,EM
        505 FORMAT( 'MINIMUM.AND.MAXIMUM ELEVATIONS = ',G12.5,' , ',G12.5)
        TEMP1 = EV/(MX*MY)
        WRITE(OUTHIST,506) TEMP1
        WRITE(OUTSUMMARY,907) TEMP1
        506 FORMAT( 'AVERAGE RELATIVE ELEVATION = ',G12.5)
        TEMP1=ST/(MX*(MYY))
        WRITE(OUTHIST,701) STMIN,TEMP1,STMAX
        WRITE(OUTSUMMARY,801) STMIN,TEMP1,STMAX
        701 FORMAT(' MINIMUM, AVERAGE, AND MAXIMUM SEDIMENT THICKNESS=',/, &
        3(G12.5,' '))
        TEMP1=NB
        TEMP2=NR
        TEMP3=TEMP1/(TEMP1+TEMP2)
        WRITE(OUTHIST,345) TEMP3
        345   FORMAT(' FRACTION OF BEDROCK SLOPES=',G12.5)
        WRITE(OUTSUMMARY,907) TEMP3
        IF (TEMP2 > 0.0) THEN
            TEMP3=RT/TEMP2
        ELSE
            TEMP3=0.0
        ENDIF
        WRITE(OUTHIST,346) RTMIN,TEMP3,RTMAX
        346   FORMAT(' MINIMUM, AVERAGE, AND MAXIMUM REGOLITH THICKNESS=',/, &
        3(G12.5,' '))
        TEMP1 = AG/((MX)*(MYY))
        WRITE(OUTHIST,507) TEMP1
        WRITE(OUTSUMMARY,907) TEMP1
        507 FORMAT( 'AVERAGE GRADIENT = ',G12.5)
        WRITE(OUTHIST,508) MG
        WRITE(OUTSUMMARY,907) MG
        508 FORMAT( 'MAXIMUM GRADIENT = ',G12.5)
        GRADAVERAGE = AG/(MX*(MYY))
        !          call debugit()
        906 FORMAT(2(' ',E15.8))
        907 FORMAT(' ',E15.8)
        908 FORMAT(2(' ',I5),3(' ',E15.8))
        909 FORMAT(2(' ',I5),2(' ',E15.8))
        910 FORMAT(' ',E15.8,' ',E15.8)
        911 FORMAT(' ',I6,' ',E15.8)
        801 FORMAT(3(' ',E15.8))
        RETURN
    END ! SUBROUTINE PRINT_MORPHOMETRY
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PRINT_SIMULATION_INFORMATION()
        !      *********************************************************************
        !       This is the main routine for printing out various summary statistics
        !       for rate processes in simulations.
        !       The following variables are printed out in this routine:
        !              current time and iteration number
        !              number of active cells (excluding non-eroding border cells)
        !              current time increment
        !              min, average, and maximum elevation changes
        !              average slope elevation change & no. of slope elements
        !              maximum and minimum slope elevation changes
        !              average raw channel elevation change & no. of slope elements
        !              maximum and mininum raw channel elevation changes
        !              minimum, average, and maximum adjusted channel elev changes
        !      *********************************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I, J, NCC,NS
        REAL(4) :: CA,CMAX,CMIN,SCMAX,SCMIN,SEDIMENT_VOLUME,SEDDIFFERENCE,ELDIFF
        REAL(4) :: DIFFEL,TOTERODE,TOTDEPOSIT,NTOTDEPOSIT,NSEDVOLUME
        REAL(4) ::  SA,SMAX,SMIN,SC,TEMP1,TEMP2,RCA,RCMAX,RCMIN
        REAL(4) :: AVG,VARI,SKEW,KURT,TOTAL
        REAL(4) :: XERODESTATS(5),CRATERSTATS(5),EOSTATS(5),LAVASTATS(5)
        !  MODIFIES: CUMULATIVE_ELEVATION_CHANGE, CUMULATIVE_EOLIAN_CHANGE 
        !            CUMULATIVE_LAVA_CHANGE. CUMULATIVE_CRATERING_CHANGE
        !            DEPOSITWORK, ERODEWORK, CRATERWORK, LAVAWORK
        !            SLOPEWORK, SLOPEGRAV, DEPOSITGRAV, ERODEGRAV, CRATERGRAV
        !  CALLS: RESET_MOMENTS, MOMENTS, CALCULATE_MOMENTS, PRINT_MOMENTS
        SEDIMENT_VOLUME=0.0
        SEDDIFFERENCE=0.0
        TOTERODE=0.0
        TOTDEPOSIT=0.0
        NTOTDEPOSIT=0.0
        NSEDVOLUME=0.0
        CALL RESET_MOMENTS(XERODESTATS)
        CALL RESET_MOMENTS(CRATERSTATS)
        CALL RESET_MOMENTS(EOSTATS)
        CALL RESET_MOMENTS(LAVASTATS)
        DO  I=1,MX
            DO  J=1,MY
                DIFFEL=ELEVATION(I,J)-INITIAL_ELEVATION(I,J)
                IF (DIFFEL > 0.0) THEN
                    NTOTDEPOSIT=NTOTDEPOSIT+1
                    TOTDEPOSIT=TOTDEPOSIT+DIFFEL
                ELSE
                    TOTERODE=TOTERODE+DIFFEL
                ENDIF
                ELDIFF=ELEVATION(I,J)-SEDIMENT_BASE(I,J)
                IF (ELDIFF > 0.0) THEN
                    NSEDVOLUME=NSEDVOLUME+1.0
                    SEDIMENT_VOLUME=SEDIMENT_VOLUME+ELDIFF
                ELSE
                    SEDDIFFERENCE=SEDDIFFERENCE+ELDIFF
                ENDIF
            ENDDO
        ENDDO
        WRITE(OUTHIST,1010) SEDIMENT_VOLUME,NSEDVOLUME,SEDDIFFERENCE
        1010  FORMAT(' TOTAL SEDIMENT VOLUME=',G12.5,' N DEPOSITS=',G12.5,&
        ' LAG IN SEDBASE RECALCULATION=',G12.5)
        WRITE(OUTHIST,1020) TOTDEPOSIT,NTOTDEPOSIT,TOTERODE
        1020  FORMAT(' TOTAL AGGRADATIONAL VOLUME=',G12.5,' N AGGRADE=', &
        G12.5,' TOTAL EROSIONAL VOLUME=',G12.5)
        CMAX =-1.0E+20
        CMIN =1.0E+20
        CA = 0.0
        RCMAX =-1.0E+20
        RCMIN =1.0E+20
        RCA = 0.0
        SMAX =-1.0E+20
        SMIN =1.0E+20
        SA = 0.0
        SC = 0.0
        SCMAX =-1.0E+20
        SCMIN =1.0E+20
        NCC=0
        NS=0
        DO  J=1,MYY
            DO  I=1,MX
                NS = NS+1
                SA = SA + ERODE_SLOPE(I,J)
                IF (ERODE_SLOPE(I,J) > SMAX) SMAX = ERODE_SLOPE(I,J)
                IF (ERODE_SLOPE(I,J) < SMIN) SMIN = ERODE_SLOPE(I,J)
                SC = SC + CFNW(I,J)
                IF (CFNW(I,J) > SCMAX) SCMAX=CFNW(I,J)
                IF (CFNW(I,J) < SCMIN) SCMIN=CFNW(I,J)
                RCA = RCA + CFNE(I,J)
                IF (CFNE(I,J) > RCMAX) RCMAX=CFNE(I,J)
                IF (CFNE(I,J) < RCMIN) RCMIN=CFNE(I,J)
                NCC = NCC+1
                CA = CA+ERODE_CHANNEL(I,J)
                IF (ERODE_CHANNEL(I,J) > CMAX) CMAX = ERODE_CHANNEL(I,J)
                IF (ERODE_CHANNEL(I,J) < CMIN) CMIN = ERODE_CHANNEL(I,J)
            ENDDO
        ENDDO
        WRITE(OUTHIST,502) TOTAL_ITERATIONS, PRESENT_TIME,BOUNDARY_LOWERING_RATE
        WRITE(OUTSUMMARY,911) TOTAL_ITERATIONS, PRESENT_TIME,BOUNDARY_LOWERING_RATE
        502 FORMAT( ' ITERATION=',I9,' TIME=',G12.5,' EROSION RATE=',G12.5)
        WRITE(OUTHIST,503)TIME_INCREMENT
        WRITE(OUTSUMMARY,907)TIME_INCREMENT
        503 FORMAT( 'TIME INCREMENT=',G12.5)
        TEMP1 =TIME_INCREMENT*SC/(MX*(MYY))
        SCMAX=SCMAX*TIME_INCREMENT
        SCMIN=SCMIN*TIME_INCREMENT
        WRITE(OUTHIST,509) SCMIN,TEMP1,SCMAX
        WRITE(OUTSUMMARY,801 )SCMIN,TEMP1,SCMAX
        509 FORMAT( 'MIN,AVG&MAX ELEVATION CHANGE =',/,3(' ',G12.5))
        IF (NS > 0) THEN
            TEMP1 =TIME_INCREMENT*SA/NS
        ELSE
            TEMP1 = 0.0
        ENDIF
        WRITE(OUTHIST,510) TEMP1,NS
        WRITE(OUTSUMMARY,910) TEMP1,NS
        510 FORMAT( 'AVERAGE SLOPE CHANGE=',G12.5,' NO. SLOPE ELEMENTS=',I5)
        TEMP1=TIME_INCREMENT*SMIN
        TEMP2=TIME_INCREMENT*SMAX
        WRITE(OUTHIST,513) TEMP1,TEMP2
        WRITE(OUTSUMMARY,906) TEMP1,TEMP2
        513 FORMAT( 'MINIMUM.AND.MAXIMUM SLOPE ELEVATION CHANGES=',  &
        G12.5,' ',G12.5)
        IF (NCC > 0) THEN
            TEMP1=TIME_INCREMENT*CA/NCC
        ELSE
            TEMP1=0.0
        ENDIF
        WRITE(OUTHIST,511) TEMP1,NCC
        WRITE(OUTSUMMARY,910)TEMP1,NCC
        511 FORMAT( 'AVERAGE RAW CHANNEL CHANGE='  &
        ,G12.5,' NO. CHANNEL ELEMENTS=',I5)
        TEMP1=TIME_INCREMENT*CMIN
        TEMP2=TIME_INCREMENT*CMAX
        WRITE(OUTHIST,512) TEMP1,TEMP2
        WRITE(OUTSUMMARY,906) TEMP1,TEMP2
        512 FORMAT( 'MIN & MAX RAW CHANNEL ELEVATION CHANGES=',  &
        G12.5,' ',G12.5)
        IF (NCC > 0) THEN
            TEMP1=TIME_INCREMENT*RCA/NCC
        ELSE
            TEMP1=0.0
        ENDIF
        WRITE(OUTHIST,711) TEMP1,NCC
        WRITE(OUTSUMMARY,910)TEMP1,NCC
        711 FORMAT( 'AVERAGE ADJUSTED CHANNEL CHANGE='  &
        ,G12.5,' NO. CHANNEL ELEMENTS=',I5)
        TEMP1=TIME_INCREMENT*RCMIN
        TEMP2=TIME_INCREMENT*RCMAX
        WRITE(OUTHIST,712) TEMP1,TEMP2
        WRITE(OUTSUMMARY,906) TEMP1,TEMP2
        CALL PRINT_RATE_STATISTICS()
        CALL PRINT_BEDROCK_STATISTICS()
        IF ((ELAPSEDTIME > 0.0).AND.FLUVIAL_AND_SLOPE_MODELING) THEN
            DO  I=1,MX
                DO  J=1,MY
                    CUMULATIVE_ELEVATION_CHANGE(I,J)=CUMULATIVE_ELEVATION_CHANGE(I,J)!/ELAPSEDTIME
                    IF (MODEL_EOLIAN_CHANGES) CUMULATIVE_EOLIAN_CHANGE(I,J)=CUMULATIVE_EOLIAN_CHANGE(I,J)!/ELAPSEDTIME
                    IF(MODEL_LAVA_FLOWS) CUMULATIVE_LAVA_CHANGE(I,J)=CUMULATIVE_LAVA_CHANGE(I,J)!/ELAPSEDTIME
                    CUMULATIVE_CRATERING_CHANGE(I,J)=CUMULATIVE_CRATERING_CHANGE(I,J)!/ELAPSEDTIME
                    IF (MODEL_EOLIAN_CHANGES) CALL MOMENTS(EOSTATS,CUMULATIVE_EOLIAN_CHANGE(I,J))
                    IF (MODEL_LAVA_FLOWS) CALL MOMENTS(LAVASTATS,CUMULATIVE_LAVA_CHANGE(I,J))
                    CALL MOMENTS(CRATERSTATS,CUMULATIVE_CRATERING_CHANGE(I,J))
                    CALL MOMENTS(XERODESTATS,CUMULATIVE_ELEVATION_CHANGE(I,J))
                ENDDO
            ENDDO
            CUM_DEPOSIT_WORK=CUM_DEPOSIT_WORK+DEPOSITWORK*CELL_SIZE
            CUM_ERODE_WORK=CUM_ERODE_WORK+ERODEWORK/CELL_SIZE
            CUM_CRATER_WORK=CUM_CRATER_WORK+CRATERWORK
            CUM_LAVA_WORK=CUM_LAVA_WORK+LAVAWORK
            CUM_SLOPE_WORK=CUM_SLOPE_WORK+SLOPEWORK
            CUM_FLOW_WORK=CUM_FLOW_WORK+FLOWWORK
            CUM_DEPOSIT_GRAV=CUM_DEPOSIT_GRAV+DEPOSITGRAV
            CUM_ERODE_GRAV=CUM_ERODE_GRAV+ERODEGRAV/CELL_AREA
            CUM_CRATER_GRAV=CUM_CRATER_GRAV+CRATERGRAV
            CUM_LAVA_GRAV=CUM_LAVA_GRAV+LAVAGRAV
            CUM_SLOPE_GRAV=CUM_SLOPE_GRAV+SLOPEGRAV
            CUM_FLOW_GRAV=CUM_FLOW_GRAV+FLOWGRAV
            CUM_ACCRETION=CUM_ACCRETION+ACCRETE_WORK              
            DEPOSITWORK=CELL_SIZE*DEPOSITWORK/(MX*MY*ELAPSEDTIME)
            ERODEWORK=ERODEWORK/(CELL_SIZE*MX*MY*ELAPSEDTIME)
            CRATERWORK=CRATERWORK/(MX*MY*ELAPSEDTIME)
            LAVAWORK=LAVAWORK/(MX*MY*ELAPSEDTIME)
            SLOPEWORK=SLOPEWORK/(MX*MY*ELAPSEDTIME)
            SLOPEGRAV=SLOPEGRAV/(MX*MY*ELAPSEDTIME)
            DEPOSITGRAV=DEPOSITGRAV/(MX*MY*ELAPSEDTIME)
            ERODEGRAV=ERODEGRAV/(CELL_AREA*MX*MY*ELAPSEDTIME)
            CRATERGRAV=CRATERGRAV/(MX*MY*ELAPSEDTIME)
            FLOWWORK=FLOWWORK/(MX*MY*ELAPSEDTIME)
            FLOWGRAV=FLOWGRAV/(MX*MY+ELAPSEDTIME)
            ACCRETE_WORK=ACCRETE_WORK/(MX*MY*ELAPSEDTIME)
        ENDIF
        WRITE(OUTHIST,1711) SLOPEWORK,ERODEWORK,DEPOSITWORK,CRATERWORK, &
        LAVAWORK,FLOWWORK,ELAPSEDTIME
        1711  FORMAT(' SLOPETRANS=',G12.5,' ERODETRANS=',G12.5,' DEPOSITTRANS=',&
        G12.5,/,' CRATERTRANS=',G12.5,' LAVATRANS=',G12.5, &
        ' FLOWTRANS=',G12.5,' ELAPSEDTIME=',&
        G12.5)
        WRITE(OUTHIST,11711) SLOPEGRAV,ERODEGRAV,DEPOSITGRAV,CRATERGRAV,FLOWGRAV
        11711 FORMAT(' SLOPEWORK=',G12.5,' ERODEWORK=',G12.5,' DEPOSITWORK=',&
        G12.5,/,' CRATERWORK=',G12.5,' FLOWWORK=',G12.5)
        CALL CALCULATE_MOMENTS(XERODESTATS,AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTHIST,740)
        740 FORMAT(' RUNOFF EROSION RATE MOMENTS')
        CALL PRINT_MOMENTS(AVG,VARI,SKEW,KURT,TOTAL)
        CALL CALCULATE_MOMENTS(CRATERSTATS,AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTHIST,1740)
        1740 FORMAT(' CRATER MODIFICATION RATE MOMENTS')
        CALL PRINT_MOMENTS(AVG,VARI,SKEW,KURT,TOTAL)
        CALL CALCULATE_MOMENTS(EOSTATS,AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTHIST,2740)
        2740 FORMAT(' EOLIAN MODIFICATION RATE MOMENTS')
        CALL PRINT_MOMENTS(AVG,VARI,SKEW,KURT,TOTAL)
        CALL CALCULATE_MOMENTS(LAVASTATS,AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTHIST,3740)
        3740 FORMAT(' LAVA MODIFICATION RATE MOMENTS')
        CALL PRINT_MOMENTS(AVG,VARI,SKEW,KURT,TOTAL)
        WRITE(OUTWORK,4740) ITERATION, PRESENT_TIME,ELAPSEDTIME,SLOPEWORK, &
            ERODEWORK,DEPOSITWORK,CRATERWORK,LAVAWORK,FLOWWORK
4740        FORMAT('iTER=',I9,' TIME=',G12.5,' ELAPSEDTIME=',G12.5,/, &
            ' SLOPETRANS=',G12.5,' ERODETRANS=',G12.5,' DEPOSITTRANS=', &
            G12.5,/,' CRATERTRANS=',G12.5,' LAVATRANS=',G12.5, &
            ' FLOWTRANS=',G12.5)
        WRITE(OUTWORK,4750) SLOPEGRAV,ERODEGRAV,DEPOSITGRAV,CRATERGRAV,FLOWGRAV &
            ,LAVAGRAV,ACCRETE_WORK
4750        FORMAT(' SLOPEWORK=',G12.5,' ERODEWORK=',G12.5,' DEPOSITWORK=',&
             G12.5,/,' CRATERWORK=',G12.5,' FLOWWORK=',G12.5,' LAVAWORK=',G12.5, &
            'ACCRETEWORK=',G12.5,/)
        WRITE(OUTWORK,4760) CUM_SLOPE_WORK,CUM_ERODE_WORK,CUM_DEPOSIT_WORK, &
            CUM_CRATER_WORK,CUM_LAVA_WORK,CUM_FLOW_WORK
4760        FORMAT('CUMULATIVE VALUES: SLOPETRANS=',G12.5,' ERODETRANS=',G12.5 &
            ' DEPOSITTRANS=',G12.5,/,' CRATERTRANS=',G12.5,' LAVATRANS=',G12.5 &
            ' FLOWTRANS=',G12.5)
        WRITE(OUTWORK,4770) CUM_SLOPE_GRAV,CUM_ERODE_GRAV,CUM_DEPOSIT_GRAV, &
            CUM_CRATER_GRAV,CUM_LAVA_GRAV,CUM_FLOW_GRAV,CUM_ACCRETION
4770        FORMAT(' SLOPEWORK=',G12.5,' ERODEWORK=',G12.5,' DEPOSITWORK=',&
             G12.5,/,' CRATERWORK=',G12.5,' LAVAWORK=',G12.5,' FLOWWORK=',G12.5, &
            'ACCRETEWORK=',G12.5,/,/)          
        DEPOSITWORK=0.0
        ERODEWORK=0.0
        CRATERWORK=0.0
        LAVAWORK=0.0
        SLOPEGRAV=0.0
        ERODEGRAV=0.0
        DEPOSITGRAV=0.0
        CRATERGRAV=0.0
        ELAPSEDTIME=0.0
        SLOPEWORK=0.0
        FLOWWORK=0.0
        FLOWGRAV=0.0
        ACCRETE_WORK=0.0
        712 FORMAT( 'MIN & MAX ADJUSTED CHANNEL EL. CHANGES=', &
        G12.5,' ',G12.5)
        906 FORMAT(2(' ',E15.8))
        907 FORMAT(' ',E15.8)
        908 FORMAT(2(' ',I5),3(' ',E15.8))
        909 FORMAT(2(' ',I5),2(' ',E15.8))
        910 FORMAT(' ',E15.8,' ',I5)
        911 FORMAT(' ',I9,' ',E15.8,' ',E15.8)
        801 FORMAT(3(' ',E15.8))
        RETURN
    END ! SUBROUTINE PRINT_SIMULATION_INFORMATION
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PRINT_VARIABLE_SUMMARIES()
        !      *********************************************************************
        !       Calculates and prints the correlation between y (j) location and the
        !         following variables:
        !              gradient
        !              slope erosion rate
        !              channel erosion rate
        !       Calculates and prints an anova analysis of x (i) and y(j) dependency
        !         of the following variables:
        !              gradient
        !              slope erosion rate
        !              channel erosion rate
        !       Calculates and prints a table of relief ratios versus x (i) and y (j)
        !       Calculates the percentiles of the following variables:
        !              drainage area
        !              elevation
        !              gradient
        !              gradient convergence
        !              profile curvature
        !              planform curvature
        !              ln(drainage area/gradient)
        !
        !       Primarily intended for steady-state basin simulations
        !  MODIFIES: SORTING_VECTOR
        !  CALLS: PERCENTILES
        !      *********************************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,K,II,JJ,INDX,JNDX,NNNNN
        REAL(4) :: XX,YY,SUMX,SUMY,SUMXY,TOTAL,SUMXX,SUMYY,CORREL,LTANB
        REAL(4) :: RELIEFRATIOS(6,6,3),REGRSUMS(6,6)
        REAL(4) :: SLOPEANOVA(6,6,3),ESLOPEANOVA(6,6,3),ECHANANOVA(6,6,3)
        DO  I=1,6
            DO  J=1,6
                REGRSUMS(I,J)=0.0
                DO  K=1,3
                    RELIEFRATIOS(I,J,K)=0.0
                    SLOPEANOVA(I,J,K)=0.0
                    ESLOPEANOVA(I,J,K)=0.0
                    ECHANANOVA(I,J,K)=0.0
                ENDDO
            ENDDO
        ENDDO
        DO  JJ=1,MY
            YY=DFLOAT(JJ)
            JNDX=(JJ-1)/20 + 1
            DO  II=1,MX
                XX=CFNW(II,JJ)
                INDX=(II-1)/20 + 1
                REGRSUMS(1,1)=REGRSUMS(1,1)+XX
                REGRSUMS(1,2)=REGRSUMS(1,2)+XX*XX
                REGRSUMS(1,3)=REGRSUMS(1,3)+YY
                REGRSUMS(1,4)=REGRSUMS(1,4)+YY*YY
                REGRSUMS(1,5)=REGRSUMS(1,5)+XX*YY
                REGRSUMS(1,6)=REGRSUMS(1,6)+1.0
                IF ((INDX < 7).AND.(JNDX < 7)) THEN
                    SLOPEANOVA(INDX,JNDX,1)=SLOPEANOVA(INDX,JNDX,1)+XX
                    SLOPEANOVA(INDX,JNDX,2)=SLOPEANOVA(INDX,JNDX,2)+XX*XX
                    SLOPEANOVA(INDX,JNDX,3)=SLOPEANOVA(INDX,JNDX,3)+1.0
                ENDIF
                XX=ERODE_SLOPE(II,JJ)
                REGRSUMS(2,1)=REGRSUMS(2,1)+XX
                REGRSUMS(2,2)=REGRSUMS(2,2)+XX*XX
                REGRSUMS(2,3)=REGRSUMS(2,3)+YY
                REGRSUMS(2,4)=REGRSUMS(2,4)+YY*YY
                REGRSUMS(2,5)=REGRSUMS(2,5)+XX*YY
                REGRSUMS(2,6)=REGRSUMS(2,6)+1.0
                IF ((INDX < 7).AND.(JNDX < 7)) THEN
                    ESLOPEANOVA(INDX,JNDX,1)=ESLOPEANOVA(INDX,JNDX,1)+XX
                    ESLOPEANOVA(INDX,JNDX,2)=ESLOPEANOVA(INDX,JNDX,2)+XX*XX
                    ESLOPEANOVA(INDX,JNDX,3)=ESLOPEANOVA(INDX,JNDX,3)+1.0
                ENDIF
                XX=ERODE_CHANNEL(II,JJ)
                REGRSUMS(3,1)=REGRSUMS(3,1)+XX
                REGRSUMS(3,2)=REGRSUMS(3,2)+XX*XX
                REGRSUMS(3,3)=REGRSUMS(3,3)+YY
                REGRSUMS(3,4)=REGRSUMS(3,4)+YY*YY
                REGRSUMS(3,5)=REGRSUMS(3,5)+XX*YY
                REGRSUMS(3,6)=REGRSUMS(3,6)+1.0
                IF ((INDX < 7).AND.(JNDX < 7)) THEN
                    ECHANANOVA(INDX,JNDX,1)=ECHANANOVA(INDX,JNDX,1)+XX
                    ECHANANOVA(INDX,JNDX,2)=ECHANANOVA(INDX,JNDX,2)+XX*XX
                    ECHANANOVA(INDX,JNDX,3)=ECHANANOVA(INDX,JNDX,3)+1.0
                ENDIF
            ENDDO
        ENDDO
        SUMXY=REGRSUMS(1,5)
        SUMX= REGRSUMS(1,1)
        SUMY= REGRSUMS(1,3)
        SUMXX=REGRSUMS(1,2)
        SUMYY=REGRSUMS(1,4)
        TOTAL=REGRSUMS(1,6)
        IF (TOTAL > 1.0) THEN
            CORREL=(SUMXY-SUMX*SUMY/TOTAL)/(SQRT((SUMXX-SUMX*SUMX/TOTAL) &
            )*SQRT((SUMYY-SUMY*SUMY/TOTAL)))
        ELSE
            CORREL=0.0
        ENDIF
        WRITE(OUTHIST,751)CORREL,TOTAL
        WRITE(OUTSUMMARY,906) CORREL,TOTAL
        751 FORMAT(' CORRELATION BETWEEN YPOS. AND SLOPE=',G12.5, &
        ' N=',G13.6)
        SUMXY=REGRSUMS(2,5)
        SUMX= REGRSUMS(2,1)
        SUMY= REGRSUMS(2,3)
        SUMXX=REGRSUMS(2,2)
        SUMYY=REGRSUMS(2,4)
        TOTAL=REGRSUMS(2,6)
        IF (TOTAL > 1.0) THEN
            CORREL=(SUMXY-SUMX*SUMY/TOTAL)/(SQRT((SUMXX-SUMX*SUMX/TOTAL)  &
            )*SQRT((SUMYY-SUMY*SUMY/TOTAL)))
        ELSE
            CORREL=0.0
        ENDIF
        WRITE(OUTSUMMARY,906) CORREL,TOTAL
        WRITE(OUTHIST,752)CORREL,TOTAL
        752 FORMAT(' CORRELATION BETWEEN YPOS. AND SLOPE EROSION=',G12.5,&
        ' N=',G13.6)
        SUMXY=REGRSUMS(3,5)
        SUMX= REGRSUMS(3,1)
        SUMY= REGRSUMS(3,3)
        SUMXX=REGRSUMS(3,2)
        SUMYY=REGRSUMS(3,4)
        TOTAL=REGRSUMS(3,6)
        IF (TOTAL > 1.0) THEN
            CORREL=(SUMXY-SUMX*SUMY/TOTAL)/(SQRT((SUMXX-SUMX*SUMX/TOTAL)  &
            )*SQRT((SUMYY-SUMY*SUMY/TOTAL)))
        ELSE
            CORREL=0.0
        ENDIF
        WRITE(OUTSUMMARY,906) CORREL,TOTAL
        WRITE(OUTHIST,753) CORREL,TOTAL
        753 FORMAT(' CORRELATION BETWEEN YPOS. AND CHAN. EROSION=',G12.5, &
        ' N=',G13.6)
        WRITE(OUTHIST,756)
        756 FORMAT(/' ANOVA TABLE FOR SLOPE [I,J,SUMX,SUMXX,N]'/)
        DO  I=1,6
            DO  J=1,6
                WRITE (OUTSUMMARY,908) I,J,(SLOPEANOVA(I,J,K),K=1,3)
                IF (SLOPEANOVA(I,J,3) == 0.0) CYCLE !    goto 200
                WRITE (OUTHIST,757) I,J,(SLOPEANOVA(I,J,K),K=1,3)
                757 FORMAT(2I5,3E13.6)
            ENDDO
        ENDDO
        WRITE(OUTHIST,758)
        758 FORMAT(/' ANOVA TABLE FOR SLOPE ERODE RATE [I,J,SUMX,SUMXX,N]'/)
        DO  I=1,6
            DO  J=1,6
                WRITE (OUTSUMMARY,908) I,J,(ESLOPEANOVA(I,J,K),K=1,3)
                IF (ESLOPEANOVA(I,J,3) == 0.0) CYCLE !    goto 202
                WRITE (OUTHIST,757) I,J,(ESLOPEANOVA(I,J,K),K=1,3)
            ENDDO
        ENDDO
        WRITE(OUTHIST,759)
        759 FORMAT(/' ANOVA TABLE FOR CHAN. ERODE RATE [I,J,SUMX,SUMXX,N]'/)
        DO  I=1,6
            DO  J=1,6
                WRITE (OUTSUMMARY,908) I,J,(ECHANANOVA(I,J,K),K=1,3)
                IF (ECHANANOVA(I,J,3) == 0.0) CYCLE !    goto 204
                WRITE (OUTHIST,757) I,J,(ECHANANOVA(I,J,K),K=1,3)
            ENDDO
        ENDDO
        DO  J=1,6
            DO  I=1,6
                RELIEFRATIOS(I,J,1)=-1.0E+30
                RELIEFRATIOS(I,J,2)=1.0E+30
                RELIEFRATIOS(I,J,3)=0.0
            ENDDO
        ENDDO
        DO  I=1,MX
            INDX=(I-1)/20 + 1
            DO  J=1,MY
                JNDX=(J-1)/20 + 1
                IF ((INDX < 7).AND.(JNDX < 7)) THEN
                    IF (ELEVATION(I,J) > RELIEFRATIOS(INDX,JNDX,1))  &
                    RELIEFRATIOS(INDX,JNDX,1)=ELEVATION(I,J)
                    IF (ELEVATION(I,J) < RELIEFRATIOS(INDX,JNDX,2)) &
                    RELIEFRATIOS(INDX,JNDX,2)=ELEVATION(I,J)
                    RELIEFRATIOS(INDX,JNDX,3)=RELIEFRATIOS(INDX,JNDX,3)+1.0
                ENDIF
            ENDDO
        ENDDO
        WRITE(OUTHIST,762)
        762 FORMAT(/,' TABLE OF RELIEF RATIOS [DELTA EL, CHAR. DIST.]',/)
        DO  I=1,6
            DO  J=1,6
                IF (RELIEFRATIOS(I,J,3) > 350.0) THEN
                    XX=RELIEFRATIOS(I,J,1)-RELIEFRATIOS(I,J,2)
                    YY=SQRT(RELIEFRATIOS(I,J,3))
                    WRITE(OUTHIST,760) I,J,XX,YY
                    760 FORMAT (2I5,2G13.6)
                ELSE
                    XX=0.0
                    YY=0.0
                ENDIF
                WRITE(OUTSUMMARY,909)I,J,XX,YY
            ENDDO
        ENDDO
        NNNNN=0
        DO  J=1,MY
            DO  I=1,MX
                NNNNN=NNNNN+1
                SORTING_VECTOR(NNNNN)=DRAINAGE_AREA(I,J)
            ENDDO
        ENDDO
        WRITE(OUTHIST,445)
        445 FORMAT(' 16,25,50,75,84 PERCENTILES OF AREA')
        IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
        NNNNN=0
        DO  J=1,MY
            DO  I=1,MX
                NNNNN=NNNNN+1
                SORTING_VECTOR(NNNNN)=ELEVATION(I,J)
            ENDDO
        ENDDO
        WRITE(OUTHIST,447)
        447 FORMAT(' 16,25,50,75,84 PERCENTILES OF ELEVATION')
        IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
        NNNNN=0
        DO  J=1,MY
            DO  I=1,MX
                NNNNN=NNNNN+1
                SORTING_VECTOR(NNNNN)=CFNW(I,J)
            ENDDO
        ENDDO
        WRITE(OUTHIST,449)
        449 FORMAT(' 16,25,50,75,84 PERCENTILES OF GRADIENT')
        IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
        NNNNN=0
        DO  J=1,MY
            DO  I=1,MX
                NNNNN=NNNNN+1
                SORTING_VECTOR(NNNNN)=CFW(I,J)
            ENDDO
        ENDDO
        WRITE(OUTHIST,451)
        451 FORMAT(' 16,25,50,75,84 PERCENTILES OF CONVERGENCE')
        IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
        NNNNN=0
        DO  J=1,MY
            DO  I=1,MX
                NNNNN=NNNNN+1
                SORTING_VECTOR(NNNNN)=CFN(I,J)
            ENDDO
        ENDDO
        WRITE(OUTHIST,453)
        453 FORMAT(' 16,25,50,75,84 PERCENTILES OF PROF CURVATURE')
        IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
        NNNNN=0
        DO  J=1,MY
            DO  I=1,MX
                NNNNN=NNNNN+1
                SORTING_VECTOR(NNNNN)=CFNE(I,J)
            ENDDO
        ENDDO
        WRITE(OUTHIST,455)
        455 FORMAT(' 16,25,50,75,84 PERCENTILES OF PLAN CURVATURE')
        IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
        NNNNN=0
        DO  J=1,MY
            DO  I=1,MX
                NNNNN=NNNNN+1
                IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                    LTANB=LOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
                ELSE
                    LTANB=0.0
                ENDIF
                SORTING_VECTOR(NNNNN)=LTANB
            ENDDO
        ENDDO
        WRITE(OUTHIST,457)
        457 FORMAT(' 16,25,50,75,84 PERCENTILES OF LN(A/S)')
        IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
        906 FORMAT(2(' ',E15.8))
        907 FORMAT(' ',E15.8)
        908 FORMAT(2(' ',I5),3(' ',E15.8))
        909 FORMAT(2(' ',I5),2(' ',E15.8))
        RETURN
    END ! SUBROUTINE PRINT_VARIABLE_SUMMARIES
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE CORRELATE(A,B)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL(4), INTENT(IN) :: A(IMMX,JMMX),B(IMMX,JMMX)
        INTEGER :: I,J
        REAL(4) :: NSUM,XSUM,YSUM,XSUMSQ,YSUMSQ,XYSUM,CORR,T1,T2
        NSUM=0.0
        XSUM=0.0
        YSUM=0.0
        XSUMSQ=0.0
        YSUMSQ=0.0
        XYSUM=0.0
        DO  I=1,MX
            DO  J=1,MYY
                NSUM=NSUM+1.0
                XSUM=XSUM+A(I,J)
                YSUM=YSUM+B(I,J)
                XYSUM=XYSUM+A(I,J)*B(I,J)
                XSUMSQ=XSUMSQ+A(I,J)**2
                YSUMSQ=YSUMSQ+B(I,J)**2
            ENDDO
        ENDDO
        T1=XSUMSQ-XSUM*XSUM/NSUM
        T2=YSUMSQ-YSUM*YSUM/NSUM
        IF ((T1 > 0.0).AND.(T2 > 0.0)) THEN
            CORR=(XYSUM-XSUM*YSUM/NSUM)/(SQRT((T1)  &
            )*SQRT((T2)))
        ELSE
            CORR=0.0
        ENDIF
        WRITE(OUTHIST,751)CORR,NSUM
        WRITE(OUTSUMMARY,906)CORR,NSUM
        751 FORMAT(' CORRELATION=',G12.5, &
        ' N=',G13.6)
        906   FORMAT(G12.5,' ',G12.5)
        RETURN
    END ! SUBROUTINE CORRELATE
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PRINT_BEDROCK_STATISTICS()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER I,J
        REAL(4) :: NSUM,XSUM,YSUM,XSUMSQ,YSUMSQ,XYSUM,CORR
        REAL(4) :: BESUM, BGSUM, BDSUM, BRSUM, BNSUM, RESUM, RGSUM, RDSUM, &
        RRSUM, RNSUM
        !  CALLS: CORRELATE
        WRITE(OUTHIST,100)
        100   FORMAT(' CORRELATION BETWEEN RESISTANCE AND ELEVATION')
        CALL CORRELATE(ELEVATION,RELATIVE_RESISTANCE)
        WRITE(OUTHIST,110)
        110   FORMAT(' CORRELATION BETWEEN RESISTANCE AND GRADIENT')
        CALL CORRELATE(D8_GRADIENT,RELATIVE_RESISTANCE)
        WRITE(OUTHIST,120)
        120   FORMAT(' CORRELATION BETWEEN RESISTANCE AND EROSION RATE')
        CALL CORRELATE(CFNW,RELATIVE_RESISTANCE)
        WRITE(OUTHIST,130)
        130   FORMAT(' CORRELATION BETWEEN RESISTANCE AND DIVERGENCE')
        CALL CORRELATE(DIVERGENCE,RELATIVE_RESISTANCE)
        WRITE(OUTHIST,140)
        140   FORMAT(' CORRELATION BETWEEN EROSION RATE AND ELEVATION')
        CALL CORRELATE(ELEVATION,CFNW)
        WRITE(OUTHIST,150)
        150   FORMAT(' CORRELATION BETWEEN EROSION RATE AND GRADIENT')
        CALL CORRELATE(D8_GRADIENT,CFNW)
        WRITE(OUTHIST,160)
        160   FORMAT(' CORRELATION BETWEEN EROSION RATE AND DIVERGENCE')
        CALL CORRELATE(DIVERGENCE,CFNW)
        NSUM=0.0
        XSUM=0.0
        YSUM=0.0
        XSUMSQ=0.0
        YSUMSQ=0.0
        XYSUM=0.0
        DO  I=1,MX
            DO  J=1,MYY
                IF (IS_ROCK_SURFACE(I,J).OR.IS_SEDIMENT_COVERED(I,J)) CYCLE !    go to 170
                NSUM=NSUM+1.0
                XSUM=XSUM+REGOLITH(I,J)
                YSUM=YSUM+ELEVATION(I,J)
                XYSUM=XYSUM+REGOLITH(I,J)*ELEVATION(I,J)
                XSUMSQ=XSUMSQ+REGOLITH(I,J)**2
                YSUMSQ=YSUMSQ+ELEVATION(I,J)**2
            ENDDO
        ENDDO
        2170  FORMAT(6G12.5)
        IF (NSUM > 3.0) THEN
            CORR=(XYSUM-XSUM*YSUM/NSUM)/(SQRT((XSUMSQ-XSUM*XSUM/NSUM) &
            )*SQRT((YSUMSQ-YSUM*YSUM/NSUM)))
        ELSE
            CORR=0.0
        ENDIF
        WRITE(OUTHIST,751)CORR,NSUM
        WRITE(OUTSUMMARY,906)CORR,NSUM
        751 FORMAT(' CORRELATION BETWEEN REGOLITH THICK. AND ELEV.=',G12.5, &
        ' N=',G13.6)
        906   FORMAT(G12.5,' ',G12.5)
        NSUM=0.0
        XSUM=0.0
        YSUM=0.0
        XSUMSQ=0.0
        YSUMSQ=0.0
        XYSUM=0.0
        DO  I=1,MX
            DO  J=1,MYY
                IF (IS_ROCK_SURFACE(I,J).OR.IS_SEDIMENT_COVERED(I,J)) CYCLE !    go to 180
                NSUM=NSUM+1.0
                XSUM=XSUM+REGOLITH(I,J)
                YSUM=YSUM+D8_GRADIENT(I,J)
                XYSUM=XYSUM+REGOLITH(I,J)*D8_GRADIENT(I,J)
                XSUMSQ=XSUMSQ+REGOLITH(I,J)**2
                YSUMSQ=YSUMSQ+D8_GRADIENT(I,J)**2
            ENDDO
        ENDDO
        IF (NSUM > 3.0) THEN
            CORR=(XYSUM-XSUM*YSUM/NSUM)/(SQRT((XSUMSQ-XSUM*XSUM/NSUM)  &
            )*SQRT((YSUMSQ-YSUM*YSUM/NSUM)))
        ELSE
            CORR=0.0
        ENDIF
        WRITE(OUTHIST,752)CORR,NSUM
        WRITE(OUTSUMMARY,906)CORR,NSUM
        752 FORMAT(' CORRELATION BETWEEN REGOLITH THICK. AND GRADIENT=',G12.5,  &
        ' N=',G13.6)
        NSUM=0.0
        XSUM=0.0
        YSUM=0.0
        XSUMSQ=0.0
        YSUMSQ=0.0
        XYSUM=0.0
        DO  I=1,MX
            DO  J=1,MYY
                IF (IS_ROCK_SURFACE(I,J).OR.IS_SEDIMENT_COVERED(I,J)) CYCLE !    go to 190
                NSUM=NSUM+1.0
                XSUM=XSUM+REGOLITH(I,J)
                YSUM=YSUM+DIVERGENCE(I,J)
                XYSUM=XYSUM+REGOLITH(I,J)*DIVERGENCE(I,J)
                XSUMSQ=XSUMSQ+REGOLITH(I,J)**2
                YSUMSQ=YSUMSQ+DIVERGENCE(I,J)**2
            ENDDO
        ENDDO
        IF (NSUM > 3.0) THEN
            CORR=(XYSUM-XSUM*YSUM/NSUM)/(SQRT((XSUMSQ-XSUM*XSUM/NSUM) &
            )*SQRT((YSUMSQ-YSUM*YSUM/NSUM)))
        ELSE
            CORR=0.0
        ENDIF
        WRITE(OUTHIST,753)CORR,NSUM
        WRITE(OUTSUMMARY,906)CORR,NSUM
        753 FORMAT(' CORRELATION BETWEEN REGOLITH THICK. AND DIV.=',G12.5, &
        ' N=',G13.6)
        NSUM=0.0
        XSUM=0.0
        YSUM=0.0
        XSUMSQ=0.0
        YSUMSQ=0.0
        XYSUM=0.0
        DO  I=1,MX
            DO  J=1,MYY
                IF (IS_ROCK_SURFACE(I,J).OR.IS_SEDIMENT_COVERED(I,J)) CYCLE !    go to 200
                NSUM=NSUM+1.0
                XSUM=XSUM+REGOLITH(I,J)
                YSUM=YSUM+CFNW(I,J)
                XYSUM=XYSUM+REGOLITH(I,J)*CFNW(I,J)
                XSUMSQ=XSUMSQ+REGOLITH(I,J)**2
                YSUMSQ=YSUMSQ+CFNW(I,J)**2
            ENDDO
        ENDDO
        IF (NSUM > 3.0) THEN
            CORR=(XYSUM-XSUM*YSUM/NSUM)/(SQRT((XSUMSQ-XSUM*XSUM/NSUM) &
            )*SQRT((YSUMSQ-YSUM*YSUM/NSUM)))
        ELSE
            CORR=0.0
        ENDIF
        WRITE(OUTHIST,754)CORR,NSUM
        WRITE(OUTSUMMARY,906)CORR,NSUM
        754 FORMAT(' CORRELATION BETWEEN REGOLITH THICK. AND EROSION' &
        ,' RATE=',G12.5, &
        ' N=',G13.6)
        BESUM=0.0
        BGSUM=0.0
        BDSUM=0.0
        BRSUM=0.0
        BNSUM=0.0
        RESUM=0.0
        RGSUM=0.0
        RDSUM=0.0
        RRSUM=0.0
        RNSUM=0.0
        DO  I=1,MX
            DO  J=1,MYY
                IF (IS_SEDIMENT_COVERED(I,J)) CYCLE !    goto 210
                IF (IS_ROCK_SURFACE(I,J)) THEN
                    BESUM=BESUM+ELEVATION(I,J)
                    BGSUM=BGSUM+D8_GRADIENT(I,J)
                    BDSUM=BDSUM+DIVERGENCE(I,J)
                    BRSUM=BRSUM+CFNW(I,J)
                    BNSUM=BNSUM+1.0
                ELSE
                    RESUM=RESUM+ELEVATION(I,J)
                    RGSUM=RGSUM+D8_GRADIENT(I,J)
                    RDSUM=RDSUM+DIVERGENCE(I,J)
                    RRSUM=RRSUM+CFNW(I,J)
                    RNSUM=RNSUM+1.0
                ENDIF
            ENDDO
        ENDDO
        IF (BNSUM > 0.0) THEN
            BESUM=BESUM/BNSUM
            BGSUM=BGSUM/BNSUM
            BDSUM=BDSUM/BNSUM
            BRSUM=BRSUM/BNSUM
        ENDIF
        IF (RNSUM > 0.0) THEN
            RESUM=RESUM/RNSUM
            RGSUM=RGSUM/RNSUM
            RDSUM=RDSUM/RNSUM
            RRSUM=RRSUM/RNSUM
        ENDIF
        WRITE(OUTHIST,755) BNSUM,BESUM,BGSUM,BDSUM,BRSUM
        755   FORMAT(' AVERAGES FOR ',G12.5,' BEDROCK SLOPE SEGMENTS',/, &
        '   ELEVATION=',G12.5,/, &
        '   GRADIENT=',G12.5,/,  &
        '   DIVERGENCE=',G12.5,/, &
        '   EROSION RATE=',G12.5)
        WRITE(OUTSUMMARY,756) BNSUM,BESUM,BGSUM,BDSUM,BRSUM
        756   FORMAT(G12.5,/,G12.5,/,G12.5,/,G12.5,/,G12.5)
        WRITE(OUTHIST,759) RNSUM,RESUM,RGSUM,RDSUM,RRSUM
        759   FORMAT(' AVERAGES FOR ',G12.5,' REGOLITH SLOPE SEGMENTS',/, &
        '   ELEVATION=',G12.5,/, &
        '   GRADIENT=',G12.5,/, &
        '   DIVERGENCE=',G12.5,/, &
        '   EROSION RATE=',G12.5)
        WRITE(OUTSUMMARY,756) RNSUM,RESUM,RGSUM,RDSUM,RRSUM
        RETURN
    END ! SUBROUTINE PRINT_BEDROCK_STATISTICS
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PRINT_RATE_STATISTICS()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,NNNNN
        REAL(4) :: EMIN,EAVG,EMAX,TOTN
        !  MODIFIES: SORTING_VECTOR
        !  CALLS: PERCENTILES
        WRITE(OUTHIST,200)
        200   FORMAT(' STATISTICS OF TOTAL EROSION RATE')
        EMIN=1.0E+25
        EMAX=-EMIN
        EAVG=0.0
        DO  I=1,MX
            DO  J=1,MYY
                EAVG=EAVG+CFNW(I,J)
                IF (CFNW(I,J) > EMAX) EMAX=CFNW(I,J)
                IF (CFNW(I,J) < EMIN) EMIN=CFNW(I,J)
            ENDDO
        ENDDO
        EAVG=EAVG/(MX*(MYY))
        WRITE(OUTHIST,110) ITERATION,EMIN,EAVG,EMAX
        110   FORMAT(' I=',I9,' EMIN=',G12.5,' EAVG=',G12.5,' EMAX=',G12.5)
        NNNNN=0
        DO  J=1,MY
            DO  I=1,MX
                NNNNN=NNNNN+1
                SORTING_VECTOR(NNNNN)=CFNW(I,J)
            ENDDO
        ENDDO
        WRITE(OUTHIST,455)
        455 FORMAT(' 16,25,50,75,84 PERCENTILES OF EROSION RATE')
        IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
        WRITE(OUTHIST,1200)
        1200   FORMAT(' STATISTICS OF CHANNEL EROSION RATE')
        EMIN=1.0E+25
        EMAX=-EMIN
        EAVG=0.0
        DO  I=1,MX
            DO  J=1,MYY
                EAVG=EAVG+CFNE(I,J)
                IF (CFNE(I,J) > EMAX) EMAX=CFNE(I,J)
                IF (CFNE(I,J) < EMIN) EMIN=CFNE(I,J)
            ENDDO
        ENDDO
        EAVG=EAVG/(MX*(MYY))
        WRITE(OUTHIST,1110) ITERATION,EMIN,EAVG,EMAX
        1110   FORMAT(' I=',I9,' EMIN=',G12.5,' EAVG=',G12.5,' EMAX=',G12.5)
        NNNNN=0
        DO  J=1,MY
            DO  I=1,MX
                NNNNN=NNNNN+1
                SORTING_VECTOR(NNNNN)=CFNE(I,J)
            ENDDO
        ENDDO
        WRITE(OUTHIST,1455)
        1455 FORMAT(' 16,25,50,75,84 PERCENTILES OF EROSION RATE')
        IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
        WRITE(OUTHIST,2200)
        2200   FORMAT(' STATISTICS OF NON-ALLUVIAL TOTAL EROSION RATE')
        EMIN=1.0E+25
        EMAX=-EMIN
        EAVG=0.0
        TOTN=0.0
        DO I=1,MX
            DO  J=1,MYY
                IF (.NOT.IS_SEDIMENT_COVERED(I,J)) THEN
                    EAVG=EAVG+CFNW(I,J)
                    IF (CFNW(I,J) > EMAX) EMAX=CFNW(I,J)
                    IF (CFNW(I,J) < EMIN) EMIN=CFNW(I,J)
                    TOTN=TOTN+1
                ENDIF
            ENDDO
        ENDDO
        IF (TOTN > 0.0) EAVG=EAVG/TOTN
        WRITE(OUTHIST,2110) ITERATION,EMIN,EAVG,EMAX
        2110   FORMAT(' I=',I9,' EMIN=',G12.5,' EAVG=',G12.5,' EMAX=',G12.5)
        NNNNN=0
        DO  J=1,MY
            DO  I=1,MX
                IF (.NOT.IS_SEDIMENT_COVERED(I,J)) THEN
                    NNNNN=NNNNN+1
                    SORTING_VECTOR(NNNNN)=CFNW(I,J)
                ENDIF
            ENDDO
        ENDDO
        WRITE(OUTHIST,2455)
        2455 FORMAT(' 16,25,50,75,84 PERCENTILES OF EROSION RATE')
        IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
        WRITE(OUTHIST,3200)
        3200   FORMAT(' STATISTICS OF NON-ALLUVIAL CHANNEL EROSION RATE')
        EMIN=1.0E+25
        EMAX=-EMIN
        EAVG=0.0
        TOTN=0.0
        DO  I=1,MX
            DO  J=1,MYY
                IF (.NOT.IS_SEDIMENT_COVERED(I,J)) THEN
                    EAVG=EAVG+CFNE(I,J)
                    IF (CFNE(I,J) > EMAX) EMAX=CFNE(I,J)
                    IF (CFNE(I,J) < EMIN) EMIN=CFNE(I,J)
                    TOTN=TOTN+1
                ENDIF
            ENDDO
        ENDDO
        IF (TOTN > 0.0) EAVG=EAVG/TOTN
        WRITE(OUTHIST,3110) ITERATION,EMIN,EAVG,EMAX
        3110   FORMAT(' I=',I9,' EMIN=',G12.5,' EAVG=',G12.5,' EMAX=',G12.5)
        NNNNN=0
        DO  J=1,MY
            DO  I=1,MX
                IF (.NOT.IS_SEDIMENT_COVERED(I,J)) THEN
                    NNNNN=NNNNN+1
                    SORTING_VECTOR(NNNNN)=CFNE(I,J)
                ENDIF
            ENDDO
        ENDDO
        WRITE(OUTHIST,455)
        3455 FORMAT(' 16,25,50,75,84 PERCENTILES OF EROSION RATE')
        IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
        RETURN
    END ! SUBROUTINE PRINT_RATE_STATISTICS
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE SLOPE_RUNOFF_CHARACTERISTICS()
    USE ERODE_GLOBALS
    IMPLICIT NONE
    INTEGER :: I,J
    REAL(4) :: NTOTAL,N0,N001,N01,N1,N10,T0,T01,T1,T2,T5,T8,T9,TEMP1
    NTOTAL=0.0
    N0=0.0
    N001=0.0
    N01=0.0
    N1=0.0
    N10=0.0
    T0=0.0
    T01=0.0
    T1=0.0
    T2=0.0
    T5=0.0
    T8=0.8
    T9=0.0
    DO J=1,MY
        DO I=1,MX
            NTOTAL=NTOTAL+1.0
            IF (D8_GRADIENT(I,J)==0.0) N0=N0+1.0
            IF (D8_GRADIENT(I,J)>0.001) N001=N001+1.0
            IF (D8_GRADIENT(I,J)>0.01) N01=N01+1.0
            IF (D8_GRADIENT(I,J)>0.1) N1=N1+1.0
            IF (D8_GRADIENT(I,J)>1.0) N10=N10+1.0
            IF (D8_GRADIENT(I,J)>0.0) THEN
               TEMP1=1.0/(1.0+EXP(-SLOPE_RUNOFF_SCALE_FACTOR*(ALOG(D8_GRADIENT(I,J))-MEDIAN_SLOPE)))
            ELSE
               TEMP1=0.0
            ENDIF
            IF (TEMP1==0.0) T0=T0+1.0
            IF (TEMP1>0.01) T01=T01+1.0
            IF (TEMP1>0.1) T1=T1+1.0
            IF (TEMP1>0.2) T2=T2+1.0
            IF (TEMP1>0.5) T5=T5+1.0
            IF (TEMP1>0.8) T8=T8+1.0
            IF (TEMP1>0.9) T0=T9+1.0
        ENDDO
    ENDDO
    WRITE(OUTHIST,100)  EXP(MEDIAN_SLOPE), SLOPE_RUNOFF_SCALE_FACTOR
100 FORMAT('STATISTICS FOR SLOPE-DEPENDENT RUNOFF GENERATION',/, &
        '  GRADIENT GIVING MEDIAN RUNOFF=',G13.6,' SLOPE RUNOFF SCALE FACTOR=',G13.6)
    WRITE(OUTHIST,200) NTOTAL,N0/NTOTAL,N001/NTOTAL,N01/NTOTAL,N1/NTOTAL,N10/NTOTAL
200 FORMAT('TABLE OF GRADIENT CHARACTERISTICS:'/,'TOTAL SLOPES=',G13.6,/,'S=0 =',G13.6,/,'S>0.001 =',G13.6,/ &
        ,'S>0.01 =',G12.5,/,'S>0.1 =',G13.6,/,'S>1.0 =',G13.6)
    WRITE(OUTHIST,300) T0/NTOTAL,T01/NTOTAL,T1/NTOTAL,T2/NTOTAL,T5/NTOTAL,T8/NTOTAL,T9/NTOTAL
300 FORMAT('TABLE OF FRACTIONAL RUNOFF VALUES:',/,'NO RUNOFF =',G13.6,/,'F>0.01 =',G13.6,/,'F>0.1 =',G13.6,/ &
        ,'F>0.2 =',G13.6,/,'F>0.5 =',G13.6,/,'F>0.8 =',G13.6,/,'F>0.9 =',G13.6)
    RETURN
    END ! SUBROUTINE SLOPE_RUNOFF_CHARACTERISTICS
    
        