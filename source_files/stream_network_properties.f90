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
    !     ###################################################################################################
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE SUMMARIZE_CHANNELS()
        USE ERODE_GLOBALS, TOPO_CONVERGENCE => CFW, IS_CHANNEL => IDO
        IMPLICIT NONE
        !      *********************************************************************
        !       This routine defines the channel network for a sequence of values of the
        !         threshold value of critical convergence for the definition of channel
        !         sources.
        !         It alse calculates the drainage density and associated basin
        !         statistics.  It is intended primarily for steady-state
        !                network simulations.
        !         For each value of critical convergence (ranging from 0.0 through 1.6),
        !         the program prints the following statistics:
        !            percentiles of source drainage area
        !            percentiles of source gradient
        !            percentiles of source profile curvature
        !            percentiles of source planform curvature
        !            percentiles of source ln(drainage area/gradient)
        !            correlation between y (j) location and source occurrence
        !               (a test for steady state)
        !            anova table for source density versus x (i) and y (j) location
        !            number of sources
        !            average and variance of support area
        !            average and variance of source gradient
        !            average and variance of source planforn curvature
        !            average and variance of source profile curvature
        !            drainage density based on zero convergence
        !            drainage density based on occurrence of fluvial erosion
        !            drainage density based on critical convergence
        !   MODIFIES: IS_CHANNEL, SORTING_VECTOR
        !   CALLS: find_stream_network, percentiles
        !      *********************************************************************
        INTEGER :: I,J,K,NCOUNT,TCOUNT,DCOUNT,TCOUNTS(100)
        INTEGER :: IDD,I_OLD,J_OLD,I_NEW,J_NEW,LOCAL_DIRECTION,KKK,INDX,JNDX,NNNNN,NVALUES
        REAL(4) :: CRITVALUE,TEMP,TEMP1,CVALUES(100)
        REAL(4) :: TVALUES,PROFSUM,PROFSQSUM,PLANSUM,PLANSQSUM,ARSUM,ARSQSUM
        REAL(4) :: SLPSUM,SLPSQSUM,THRESHOLDS(6)
        REAL(4) :: REGRSUMS(6),AREASRCANOVA(6,6,3)
        REAL(4) :: SUMX,SUMXY,SUMY,SUMYY,SUMXX,TOTAL,CORREL,XX,YY,LTANB
        !      *********************************************************************
        !       The critical convergence values are defined below
        !      *********************************************************************
        THRESHOLDS(1)=0.0
        THRESHOLDS(2)=0.1
        THRESHOLDS(3)=0.2
        THRESHOLDS(4)=0.4
        THRESHOLDS(5)=0.8
        THRESHOLDS(6)=1.6
        !      *********************************************************************
        !       This defines threshold values normalized by average gradient
        !      *********************************************************************
        DO  I = 1, NCRITS+1
            CVALUES(I) = MAXCRIT*GRADAVERAGE*(I-1)/NCRITS
            TCOUNTS(I) = 0
        ENDDO
        !      *********************************************************************
        !       This counts the total length of channel for each value of critical
        !         convergence
        !      *********************************************************************
        DO  J = 1, MYY
            DO  I = 1, MX
                IF (ERODE_CHANNEL(I,J) < 0.0)  THEN
                    DO  K = 1, NCRITS+1
                        IF (TOPO_CONVERGENCE(I,J) >= CVALUES(K)) TCOUNTS(K) = TCOUNTS(K)+1
                    ENDDO
                ENDIF
            ENDDO
        ENDDO
        WRITE(OUTHIST,505)
        505 FORMAT(/,' TABLE OF CRENULATION-BASED DRAINAGE DENSITIES',/)
        !      *********************************************************************
        !       The routine cycles through each value of the critical threshold and prints
        !         out the drainage density that results for each of these values
        !      *********************************************************************
        DO  K =1,   NCRITS+1
            TEMP =FLOAT(TCOUNTS(K)+MX)/(MX*MY)
            TEMP1=CVALUES(K)/GRADAVERAGE
            WRITE(OUTHIST,506) TEMP1,TEMP
            WRITE(OUTSUMMARY,906)TEMP1,TEMP
            506 FORMAT(' CRIT. VALUE=',G12.5,' D.D.=',G12.5)
        ENDDO
        !      *********************************************************************
        !       We cycle through all the critical values, defining the channel network and
        !         calculating various statistics
        !      *********************************************************************
        L987: DO  KKK=1,6
            !      *********************************************************************
            !       This initializes various statistical parameters
            !      *********************************************************************
            TVALUES=0.0
            PROFSUM=0.0
            PROFSQSUM=0.0
            PLANSUM=0.0
            PLANSQSUM=0.0
            ARSUM=0.0
            ARSQSUM=0.0
            SLPSUM=0.0
            SLPSQSUM=0.0
            DO  I=1,6
                REGRSUMS(I)=0.0
                DO  J=1,6
                    DO  K=1,3
                        AREASRCANOVA(I,J,K)=0.0
                    ENDDO
                ENDDO
            ENDDO
            !      *********************************************************************
            !       Calculate the critical divergence, 'critvalue', for each threshold value
            !      *********************************************************************
            CRITVALUE = THRESHOLDS(KKK)*GRADAVERAGE
            WRITE(OUTHIST,986) THRESHOLDS(KKK),CRITVALUE
            WRITE(OUTSUMMARY,906) THRESHOLDS(KKK),CRITVALUE
            986 FORMAT(/,' CRITICAL CONVERGE. FOR DRAINAGE NETWORK=',2G15.7,/)
            DCOUNT = 0
            NCOUNT = 0
            TCOUNT = 0
            NVALUES = 0
            IS_CHANNEL=0
            !      *********************************************************************
            !       Here we cycle through all matrix points to determine locations with
            !         divergence values greater than critvalue.  local divergence values are in
            !         the TOPO_CONVERGENCE matrix.  The matrix IS_CHANNEL records locations that are channels as
            !         follows:
            !             IS_CHANNEL=0 for non-channels
            !             IS_CHANNEL=1 for channels
            !         The criterion is that fluvial erosion has occurred during the last
            !         iteration (erosion<0) and that the divergence exceeds the critical value
            !         ncount records the total number of channel locations
            !      *********************************************************************
            DO  J = 1, MYY
                DO  I = 1, MX
                    IS_CHANNEL(I,J)=0
                    IF ((ERODE_CHANNEL(I,J) < 0.0).AND.(TOPO_CONVERGENCE(I,J) >= CRITVALUE)) THEN
                        NCOUNT = NCOUNT+1
                        IS_CHANNEL(I,J)=1
                    ENDIF
                    IF (ERODE_CHANNEL(I,J) < 0.0) TCOUNT = TCOUNT + 1
                    IF (TOPO_CONVERGENCE(I,J) > 0.0)  DCOUNT = DCOUNT+1
                ENDDO
            ENDDO
            !      *********************************************************************
            !       IS_CHANNEL=0 for the flow exit boundary at j=my
            !      *********************************************************************
            IF (.NOT.IS_Y_PERIODIC) THEN
                DO  I=1,MX
                    IS_CHANNEL(I,MYY)=0
                ENDDO
            ENDIF
            IF (DO_FLOW_BOUNDARIES) THEN
                DO J=1,MYY
                    IS_CHANNEL(1,J)=0
                    IS_CHANNEL(MX,J)=0
                ENDDO
                DO I=1,MX
                    IS_CHANNEL(I,1)=0
                    IS_CHANNEL(I,MYY)=0
                ENDDO
            ENDIF
            !      *********************************************************************
            !       This is an 'enforcement' loop that makes sure that channels continue
            !         continuously downstream, until they either reach an exit at j=my or and
            !         enclosed depression - It starts at defined channel locations and proceeds
            !         downstream, altering any locations that were not defined as channels by
            !         the convergence criterion (i.e., they have IS_CHANNEL=0) and force them to be
            !         considered to be channels (IS_CHANNEL=1)
            !      *********************************************************************
            L210: DO  J=1,MYY
                M210: DO  I=1,MX
                    IF (IS_CHANNEL(I,J) /= 1) CYCLE M210
                    I_OLD=I
                    J_OLD=J
                    L220: DO
                        LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                        IF (LOCAL_DIRECTION <= 1) THEN
                            CYCLE M210
                        ENDIF
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
                        IF (IS_Y_PERIODIC) THEN
                        ELSE
                            IF (J_OLD == MY) CYCLE M210
                        ENDIF
                        IF (DO_FLOW_BOUNDARIES) THEN
                            IF (FLOW_DIRECTION(I_OLD,J_OLD)<=1) CYCLE M210
                        ENDIF
                        IF (IS_CHANNEL(I_OLD,J_OLD) == 0) THEN
                            IS_CHANNEL(I_OLD,J_OLD)=1
                            NCOUNT=NCOUNT+1
                        ELSE
                            CYCLE M210
                        ENDIF
                    ENDDO L220
                ENDDO M210
            ENDDO L210
            !      *********************************************************************
            !       We only do the following at times when printed output is requested and only
            !         for the value of critical convergence which is the default value
            !       This writes to an output file the i,j location and flow direction (id) for
            !         each channel segment
            !      *********************************************************************
            IF ((KKK == KWRITE).AND.(WRITE_OUT_CHANNELS)) THEN
                WRITE(OUTCHAN,500) TOTAL_ITERATIONS,NCOUNT
                500 FORMAT (2I5)
                DO  I = 1, MX
                    DO  J = 1, MYY
                        IDD=FLOW_DIRECTION(I,J)
                        IF (IDD <= 1) IDD=0
                        IF(IS_CHANNEL(I,J) == 1)  &
                        WRITE(OUTCHAN,501) I,J,IDD
                        501 FORMAT (I3,' ',I3,' ',I1)
                    ENDDO
                ENDDO
            ENDIF
            !      *********************************************************************
            !       This cycles through the matrix, and identifies those channel segments which
            !         originate from sources (no upstream contribution). at the close of this
            !         loop channel segments at sources have IS_CHANNEL=1 and those downstream have
            !         IS_CHANNEL=2.
            !      *********************************************************************
            DO  J=1,MYY
                DO  I=1,MX
                    IF (IS_CHANNEL(I,J) == 0) CYCLE
                    LOCAL_DIRECTION = FLOW_DIRECTION(I,J)
                    IF (LOCAL_DIRECTION <= 1) THEN
                        CYCLE
                    ENDIF
                    J_NEW=J+DOWNSTREAM(LOCAL_DIRECTION,2)
                    I_NEW=I+DOWNSTREAM(LOCAL_DIRECTION,1)
                    IF (IS_X_PERIODIC) THEN
                        IF (I_NEW < 1) I_NEW = MX
                        IF (I_NEW > MX) I_NEW=1
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (J_NEW < 1) J_NEW=MY
                        IF (J_NEW > MY) J_NEW=1
                    ENDIF
                    IF (IS_CHANNEL(I_NEW,J_NEW) == 1) IS_CHANNEL(I_NEW,J_NEW)=2
                ENDDO
            ENDDO
            !      *********************************************************************
            !       This routine looks for enclosed drainage basins (indicated by the variable
            !         'enclosed' and makes sure that the outlet of the enclosed basin at iout,
            !         jout has IS_CHANNEL=0
            !      *********************************************************************
            DO  K=1,TOTAL_BASINS
                IF (COMPLETE_RUNOFF) THEN
                    IF (ENCLOSED(K)) THEN
                        IS_CHANNEL(I_OUTFLOW(K),J_OUTFLOW(K))=0
                    ENDIF
                ENDIF
            ENDDO
            !      *********************************************************************
            !       Call the subroutine 'stream' to calculate a variety of morphometric
            !         parameters for interior and exterior links, first-order drainage basins,
            !         and the simulated basin as a whole.
            !      *********************************************************************
            IF (.NOT.IS_Y_PERIODIC) CALL FIND_STREAM_NETWORK()
            IF (.NOT.DO_FLOW_BOUNDARIES) CALL FIND_STREAM_NETWORK()
            !      *********************************************************************
            !       This calculates sums and sums of squares for channel segments for several
            !         variables:
            !             nvalues - (i) the total number of stream segments
            !             tvalues - (r) the total number of stream segments
            !             slpsum - (r) the sum of the average gradient
            !                       (defined as the slope of a planar
            !                       trend surface fitted to the 9 points surrounding the given
            !                       point - see summariz.f:  the gradient values have been
            !                       stored in the matrix cfnw
            !             slpsqsum - (r) sum of squares of average gradient
            !             arsum - (r) sum of the contributing drainage areas (from matrix ar)
            !             arsqsum - (r) sum of squares of contribution drainage area
            !             plansum - (r) sum of the planform curvature (from matrix cfne)
            !             plansqsum - (r) sum of squares of the planform curvature
            !             profsum - (r) sum of the profile curvature (from matrix cfn)
            !             profsqsum - (r) sum of squares of profile curvature
            !      *********************************************************************
            DO  J=1,MY
                DO  I=1,MX
                    IF (IS_CHANNEL(I,J) /= 1) CYCLE
                    NVALUES=NVALUES+1
                    TVALUES=TVALUES+1.0
                    SLPSUM=SLPSUM+CFNW(I,J)
                    SLPSQSUM=SLPSQSUM+CFNW(I,J)**2
                    ARSUM=ARSUM+DRAINAGE_AREA(I,J)
                    ARSQSUM=ARSQSUM+DRAINAGE_AREA(I,J)**2
                    PLANSUM=PLANSUM+CFNE(I,J)
                    PLANSQSUM=PLANSQSUM+CFNE(I,J)**2
                    PROFSUM=PROFSUM+CFN(I,J)
                    PROFSQSUM=PROFSQSUM+CFN(I,J)**2
                ENDDO
            ENDDO
            !      *********************************************************************
            !       For larger values of critical convergence write out some detailed
            !         information to the 'outsource' file of drainage area, local gradient,
            !         average gradient, planform curvature, profile curvature, and ltanb (the
            !         natural logarithm of drainage area divided by local gradient)
            !      *********************************************************************
            IF (THRESHOLDS(KKK) > 0.2) THEN
                WRITE(OUTSOURCE,408) NVALUES, THRESHOLDS(KKK)
                408 FORMAT(I6,',',G12.5)
                DO  J=1,MYY
                    DO  I=1,MX
                        IF (IS_CHANNEL(I,J) /= 1) CYCLE
                        IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                            LTANB=LOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
                        ELSE
                            LTANB=0.0
                        ENDIF
                        WRITE(OUTSOURCE,409) DRAINAGE_AREA(I,J),D8_GRADIENT(I,J),CFNW(I,J),CFN(I,J), &
                        CFNE(I,J),LTANB
                        409 FORMAT(' ',6G12.5)
                    ENDDO
                ENDDO
            ENDIF
            !      *********************************************************************
            !       Calculates coefficients necessary to regress contributing draiange area of
            !         channel sources with location relative to the flow exit at j=my (mostly used to
            !         tell if the basin is close to steady state) as well for as an analysis of
            !         variance with respect to location within the basin.
            !         The variable vect is a vector of all drainage areas of channel sources
            !      *********************************************************************
            NNNNN=0
            DO  J=1,MYY
                YY=FLOAT(J)
                JNDX=(J-1)/20 + 1
                DO  I=1,MX
                    IF (IS_CHANNEL(I,J) /= 1) CYCLE
                    XX=DRAINAGE_AREA(I,J)
                    NNNNN=NNNNN+1
                    SORTING_VECTOR(NNNNN)=DRAINAGE_AREA(I,J)
                    INDX=(I-1)/20 + 1
                    REGRSUMS(1)=REGRSUMS(1)+XX
                    REGRSUMS(2)=REGRSUMS(2)+XX*XX
                    REGRSUMS(3)=REGRSUMS(3)+YY
                    REGRSUMS(4)=REGRSUMS(4)+YY*YY
                    REGRSUMS(5)=REGRSUMS(5)+XX*YY
                    REGRSUMS(6)=REGRSUMS(6)+1.0
                    IF ((INDX < 7).AND.(JNDX < 7)) THEN
                        AREASRCANOVA(INDX,JNDX,1)=AREASRCANOVA(INDX,JNDX,1)+XX
                        AREASRCANOVA(INDX,JNDX,2)=AREASRCANOVA(INDX,JNDX,2)+XX*XX
                        AREASRCANOVA(INDX,JNDX,3)=AREASRCANOVA(INDX,JNDX,3)+1.0
                    ENDIF
                ENDDO
            ENDDO
            !      321 continue
            !      *********************************************************************
            !       Call the subroutine percentiles to calculate and print out the 16th, 25th,
            !         50th, 75th and 84th percentiles of source drainage areas
            !       The following loops print percentiles for average gradient at sources,
            !         profile curvature, and planform curvature at channel sources
            !      *********************************************************************
            WRITE(OUTHIST,775)
            775 FORMAT(' 16,25,50,75,85 PERCENTILES OF SOURCE AREA')
            IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
            NNNNN=0
            DO  J=1,MYY
                DO  I=1,MX
                    IF (IS_CHANNEL(I,J) /= 1) CYCLE
                    NNNNN=NNNNN+1
                    SORTING_VECTOR(NNNNN)=CFNW(I,J)
                ENDDO
            ENDDO
            WRITE(OUTHIST,328)
            328 FORMAT(' 16,25,50,75,84 PERCENTILES OF SOURCE GRADIENT')
            IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
            NNNNN=0
            DO  J=1,MYY
                DO  I=1,MX
                    IF (IS_CHANNEL(I,J) /= 1) CYCLE
                    NNNNN=NNNNN+1
                    SORTING_VECTOR(NNNNN)=CFN(I,J)
                ENDDO
            ENDDO
            WRITE(OUTHIST,330)
            330 FORMAT(' 16,25,50,75,84 PRCNTLS SOURCE PROF CURV')
            IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
            NNNNN=0
            DO  J=1,MYY
                DO  I=1,MX
                    IF (IS_CHANNEL(I,J) /= 1) CYCLE
                    NNNNN=NNNNN+1
                    SORTING_VECTOR(NNNNN)=CFNE(I,J)
                ENDDO
            ENDDO
            WRITE(OUTHIST,332)
            332 FORMAT(' 16,25,50,75,84 PRCNTLS SOURCE PLAN CURV')
            IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
            NNNNN=0
            DO  J=1,MYY
                DO  I=1,MX
                    IF (IS_CHANNEL(I,J) /= 1) CYCLE
                    NNNNN=NNNNN+1
                    IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                        LTANB=LOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
                    ELSE
                        LTANB=0.0
                    ENDIF
                    SORTING_VECTOR(NNNNN)=LTANB
                ENDDO
            ENDDO
            WRITE(OUTHIST,334)
            334 FORMAT(' 16,25,50,75,84 PRCNTLS SOURCE LN(A/S)')
            IF (NNNNN >= 0) CALL PERCENTILES(NNNNN)
            !      *********************************************************************
            !       The following calculate the correlation between y position and contributing
            !         area at sources as well as an anova table for both x and y position versus
            !         source areas
            !      *********************************************************************
            SUMXY=REGRSUMS(5)
            SUMX= REGRSUMS(1)
            SUMY= REGRSUMS(3)
            SUMXX=REGRSUMS(2)
            SUMYY=REGRSUMS(4)
            TOTAL=REGRSUMS(6)
            IF (TOTAL > 1.0) THEN
                CORREL=(SUMXY-SUMX*SUMY/TOTAL)/(SQRT((SUMXX-SUMX*SUMX/TOTAL) &
                )*SQRT((SUMYY-SUMY*SUMY/TOTAL)))
            ELSE
                CORREL=0.0
            ENDIF
            WRITE(OUTHIST,751)CORREL,TOTAL
            WRITE(OUTSUMMARY,906) CORREL,TOTAL
            751 FORMAT(' CORRELATION BETWEEN YPOS. AND SOURCE AREA=',G12.5, &
            ' N=',G13.6)
            WRITE(OUTHIST,756)
            756 FORMAT(/' ANOVA TABLE FOR SOURCE AREA[I,J,SUMX,SUMXX,N]'/)
            DO  I=1,6
                DO  J=1,6
                    WRITE (OUTSUMMARY,908) I,J,(AREASRCANOVA(I,J,K),K=1,3)
                    IF (AREASRCANOVA(I,J,3) == 0.0) CYCLE !    goto 322
                    WRITE (OUTHIST,757) I,J,(AREASRCANOVA(I,J,K),K=1,3)
                    757 FORMAT(2I5,3E13.6)
                ENDDO
            ENDDO
            !      *********************************************************************
            !       We now print out total number of source area and average and variance of the
            !         following parameters for channel sources
            !             support area
            !             source average gradient
            !             planform curvature
            !             profile curvature
            !      *********************************************************************
            WRITE(OUTHIST,630) TVALUES
            WRITE(OUTSUMMARY,907) TVALUES
            630 FORMAT (' NO. OF SOURCES=',G13.6)
            IF (TVALUES > 1.0) THEN
                TEMP=ARSUM/TVALUES
                TEMP1=(ARSQSUM-TVALUES*TEMP*TEMP)/(TVALUES-1.0)
            ELSE
                TEMP=0.0
                TEMP1=0.0
            ENDIF
            WRITE (OUTHIST,632) TEMP,TEMP1
            WRITE (OUTSUMMARY,906) TEMP,TEMP1
            632 FORMAT (' AVG. AND VAR. OF SUPPORT AREA=',2E15.6)
            IF (TVALUES > 1.0) THEN
                TEMP=SLPSUM/TVALUES
                TEMP1=(SLPSQSUM-TVALUES*TEMP*TEMP)/(TVALUES-1.0)
            ELSE
                TEMP=0.0
                TEMP1=0.0
            ENDIF
            WRITE (OUTHIST,634) TEMP,TEMP1
            WRITE (OUTSUMMARY,906) TEMP,TEMP1
            634 FORMAT (' AVG. AND VAR. OF SOURCE GRADIENT=',2E15.6)
            IF (TVALUES > 1.0) THEN
                TEMP=PLANSUM/TVALUES
                TEMP1=(PLANSQSUM-TVALUES*TEMP*TEMP)/(TVALUES-1.0)
            ELSE
                TEMP=0.0
                TEMP1=0.0
            ENDIF
            WRITE (OUTHIST,636) TEMP,TEMP1
            WRITE (OUTSUMMARY,906) TEMP,TEMP1
            636 FORMAT (' AVG. AND VAR. OF SOURCE PLAN CURVATURE=',2E15.6)
            IF (TVALUES > 1.0) THEN
                TEMP=PROFSUM/TVALUES
                TEMP1=(PROFSQSUM-TVALUES*TEMP*TEMP)/(TVALUES-1.0)
            ELSE
                TEMP=0.0
                TEMP1=0.0
            ENDIF
            WRITE (OUTHIST,638) TEMP,TEMP1
            WRITE (OUTSUMMARY,906) TEMP,TEMP1
            638 FORMAT (' AVG. AND VAR. OF SOURCE PROFILE CURVATURE',2E15.6)
            !      *********************************************************************
            !       Prints out drainage desity based upon convergence - i.e., the fraction of
            !         points for which the topography is convergence
            !      *********************************************************************
            TEMP = FLOAT(DCOUNT+MX)/((MX)*(MY))
            WRITE(OUTHIST,502) TEMP
            WRITE (OUTSUMMARY,907) TEMP
            502 FORMAT(' DRAINAGE DENSITY BASED ON CONVERGENCE=',G12.5)
            !      *********************************************************************
            !       Prints out drainage density based upon whether there is any fluvial erosion
            !         or not
            !      *********************************************************************
            TEMP = FLOAT(TCOUNT+MX)/((MX)*(MY))
            WRITE(OUTHIST,503)TEMP
            WRITE (OUTSUMMARY,907) TEMP
            503 FORMAT(' DRAINAGE DENSITY BASED ON FLUVIAL EROSION=',G12.5)
            !      *********************************************************************
            !       Prints out drainage density based on the critical value of slope divergence
            !      *********************************************************************
            TEMP = FLOAT(NCOUNT+MX)/(MX*MY)
            WRITE(OUTHIST,504) TEMP
            WRITE (OUTSUMMARY,907) TEMP
            504 FORMAT(' DRAINAGE DENSITY BASED ON CRENULATIONS=',G12.5)
        ENDDO L987
        906 FORMAT(2(' ',E15.8))
        907 FORMAT(' ',E15.8)
        908 FORMAT(2(' ',I5),3(' ',E15.8))
        RETURN
    END !  SUBROUTINE SUMMARIZE_CHANNELS
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PERCENTILES(N)
        USE ERODE_GLOBALS
        !      *********************************************************************
        !       Calculates the 16h, 25th, 50th, 75th and 84th percentiles of a vector
        !          of real numbers - sorting_vector(n).  Values must be entered into sorting_vector(n)
        !          prior to call.  Also prints these values. It first sorts the vector and
        !         then finds the percentiles at the specific locations relative to the total
        !         vector length
        !      *********************************************************************
        INTEGER I,J,L,IR,N,I16,I25,I50,I75,I84
        REAL(4) :: RRA,V16,V25,V50,V75,V84
        !      *********************************************************************
        !       Bubble sort
        !      *********************************************************************
        IF (N > 1) THEN
            L=N/2+1
            IR=N
            L10: DO
                IF (L > 1) THEN
                    L=L-1
                    RRA=SORTING_VECTOR(L)
                ELSE
                    RRA=SORTING_VECTOR(IR)
                    SORTING_VECTOR(IR)=SORTING_VECTOR(1)
                    IR=IR-1
                    IF (IR == 1) THEN
                        SORTING_VECTOR(1)=RRA
                        EXIT L10
                    ENDIF
                ENDIF
                I=L
                J=L+L
                L20: DO
                    IF(J <= IR) THEN
                        IF (J < IR) THEN
                            IF (SORTING_VECTOR(J) < SORTING_VECTOR(J+1)) J=J+1
                        ENDIF
                        IF (RRA < SORTING_VECTOR(J)) THEN
                            SORTING_VECTOR(I)=SORTING_VECTOR(J)
                            I=J
                            J=J+J
                        ELSE
                            J=IR+1
                        ENDIF
                        CYCLE L20
                    ENDIF
                    EXIT L20
                ENDDO L20
                SORTING_VECTOR(I)=RRA
            ENDDO L10
            !      *********************************************************************
            !       Define locations for the specified percentiles and then print
            !      *********************************************************************
            I16=INT(0.16*FLOAT(N+1))
            I25=(N+1)/4
            I50=(N+1)/2
            I75=3*(N+1)/4
            I84=INT(0.84*FLOAT(N+1))
            IF (I16 > 0)THEN
                V16=SORTING_VECTOR(I16)
            ELSE
                V16=0.0
            ENDIF
            IF(I25 > 0) THEN
                V25=SORTING_VECTOR(I25)
            ELSE
                V25=0.0
            ENDIF
            IF(I50 > 0) THEN
                V50=SORTING_VECTOR(I50)
            ELSE
                V50=0.0
            ENDIF
            IF (I75 > 0) THEN
                V75=SORTING_VECTOR(I75)
            ELSE
                V75 = 0.0
            ENDIF
            IF (I84 > 0) THEN
                V84 = SORTING_VECTOR(I84)
            ELSE
                V84=0.0
            ENDIF
        ELSE
            V16=0.0
            V25=0.0
            V50=0.0
            V75=0.0
            V84=0.0
            I16=0
            I25=0
            I50=0
            I75=0
            I84=0
        ENDIF
        WRITE (OUTHIST,100) I16,V16,I25,V25,I50,V50,I75,V75,I84,V84
        WRITE(OUTSUMMARY,100)I16,V16,I25,V25,I50,V50,I75,V75,I84,V84
        100 FORMAT(3(' ',I5,' ',E12.5),/,2(' ',I5,' ',E12.5))
        RETURN
    END ! SUBROUTINE PERCENTILES
