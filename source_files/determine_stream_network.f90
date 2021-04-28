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
    SUBROUTINE FIND_STREAM_NETWORK()
        !      *********************************************************************
        !       This subroutice calculates and prints information about stream
        !         network geometry -- it is intended primarily for steady-state
        !         drainage basin simulations.  The subroutine prints the following
        !         morphometric indices:
        !             For external (first-order) drainage basins
        !                 correlation between link area ratio and link gradient
        !              Forrelation between link density and link gradient
        !                 average and variance of link length
        !                 average and variance of link density
        !                 average and variance of link dimensionless area
        !                 average and variance of link area ratio
        !                 average and variance of link gradient
        !             For internal links
        !                 correlation of link length and magnitude
        !             For the whole drainage network
        !                 shreve kappa
        !                 form factor
        !                 fraction of t-s links
        !                 flowpath indirectness
        !                 exterior-interior link-length ratio
        !      *********************************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER(4),ALLOCATABLE, DIMENSION(:) :: ORDER,MAGNITUDE,IILOC,JJLOC
        INTEGER(4),ALLOCATABLE, DIMENSION(:) :: TOSEGMENT
        LOGICAL(4),ALLOCATABLE, DIMENSION(:) :: SOURCE,PROCESSED
        INTEGER I,J,K,L,KSEGS,NINDIR,PRESORDER,NEXTERN,NINTERN,THISMAG
        INTEGER NSOURCE,NSAMP,I_NEW,J_NEW,LOCAL_DIRECTION,IERROR
        REAL(4) :: ATOP,ABOT
        REAL(4) :: XDIST,INDSUM,INDSSQ,XSTART,YSTART,XSUM,XEND,YEND
        REAL(4) :: XXXX,XEXTSUM,XEXTSSQ,XMAGSUM,XMAGSSQ,XINTSUM,XINTSSQ
        REAL(4) :: NTS,LDENSUM,LDENSSQ,CORREL,TOTAL,SUMX,SUMY,SUMXY
        REAL(4) :: SUMXX,SUMYY,KAPPA,FORMFAC,TCHANLEN,ARRATSUM,ARRATSSQ
        REAL(4) :: RATARSUM,RATARSSQ,LINKSLP,LDENSLP,SLPSUM,SLPSSQ
        REAL(4) :: ETOP,GTOP,EBOT,GBOT,LRATIO,LTOT,TEMP,TEMP1,ARRAT,XINTMAG
        REAL(4) :: LINKDEN,ARRATSLP
        !      *********************************************************************
        !       Several variables are defined for each matrix location from which a
        !         channel eminates - these variables are vectors, with each entry
        !         corresponding to one matrix location (in this description, each portion of
        !         a stream of length 1 or sqrt(2) extending from a matrix location is called
        !         a "segment", as compared to a stream "link", which extends from a source
        !         to the first junction, or between junctions (that is, a link is composed
        !         of one or more segments).
        !             IILOC - (i) the i-location of the matrix location
        !             JJLOC - (i) the j-location of the matrix location
        !             tosegment - the location in the vector of the matrix location lying
        !                           immediately downstream (locations at the matrix edge
        !                           have tosegment = 0)
        !             source - (l) whether the location is a source (head of an exterior
        !                          link or first order stream)
        !             magnitude - (i) the magnitude of the stream at the location (total
        !                             number of sources upstream from the location,
        !                             including the link starting from that matrix location
        !             order - (i) the strahler order of the stream segment eminating from
        !                             the matrix location
        !             processed - (l) an index of whether statistical information has
        !                             already been collected for this stream segment
        !  MODIFIES: IDO, 
        !      *********************************************************************
        !      *********************************************************************
        !       First, zero out the segment variables
        !      *********************************************************************
        ALLOCATE(ORDER(LMMX),MAGNITUDE(LMMX),IILOC(LMMX),JJLOC(LMMX),STAT=IERROR)
		IF (IERROR > 0) CALL ALLOCERROR()
        ALLOCATE(TOSEGMENT(LMMX),SOURCE(LMMX),PROCESSED(LMMX),STAT=IERROR)
		IF (IERROR > 0) CALL ALLOCERROR()
        DO  L=1,LMMX
            ORDER(L)=0
            MAGNITUDE(L)=0
            SOURCE(L)=.FALSE.
            TOSEGMENT(L)=0
            PROCESSED(L)=.FALSE.
            IILOC(L)=0
            JJLOC(L)=0
        ENDDO
        !      *********************************************************************
        !        This sets up the vector of stream locations.  We look through all matrix
        !           points for those which are occupied by streams (ido(i,j)>0) and add that
        !           point to our vector, recording the i- and j- locations (IILOC,JJLOC)
        !           also, if the matrix location is a stream source (ido(i,j)=1) then we set
        !           source to true.
        !        ksegs is the total number of stream segments
        !      *********************************************************************
        K=0
        DO  J=1,MYY
            DO  I=1,MX
                IF (IDO(I,J) > 0) THEN
                    K=K+1
                    IILOC(K)=I
                    JJLOC(K)=J
                    IF (IDO(I,J) == 1) THEN
                        SOURCE(K)=.TRUE.
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        KSEGS=K
        !      *********************************************************************
        !        This determines the vector index of the matrix location of the segment
        !         lying immediately downstream (tosegment).  It does so by traveling
        !         downstream within the matrix to determine the i,j coordinates of the next
        !         point downstream.  Then the vector is searched for the index which has the
        !         same value of IILOC and JJLOC.  The first 7 lines of the do loop are
        !         standard code for moving downstream within the matrix, using the flow
        !         direction matrix (id) and the relative displacement matrix (im), also
        !         taking into account the periodic boundaries at i=1 and i=mx and that my is
        !         the flow exit location.
        !      *********************************************************************
        DO  K=1,KSEGS
            LOCAL_DIRECTION = FLOW_DIRECTION(IILOC(K),JJLOC(K))
            IF (LOCAL_DIRECTION < 0) CYCLE 
            J_NEW=JJLOC(K)+DOWNSTREAM(LOCAL_DIRECTION,2)
            I_NEW=IILOC(K)+DOWNSTREAM(LOCAL_DIRECTION,1)
            IF (IS_X_PERIODIC) THEN
                IF (I_NEW < 1) I_NEW = MX
                IF (I_NEW > MX) I_NEW=1
            ENDIF
            IF (IS_Y_PERIODIC) THEN
                IF (J_NEW < 1) J_NEW=MY
                IF (J_NEW > MY) J_NEW=1
            ELSE
                IF (J_NEW == MY) CYCLE
            ENDIF
            DO  L=1,KSEGS
                IF ((IILOC(L) == I_NEW).AND.(JJLOC(L) == J_NEW)) THEN
                    TOSEGMENT(K)=L
                    EXIT
                ENDIF
            ENDDO
        ENDDO
        !      *********************************************************************
        !        The following do-loop progresses downstream from each source, determining
        !           the magnitude of each segment and the flowpath length from source to
        !           exit.  the following variables are used
        !        xsum - flowpath length from source to exit
        !        xend, yend - location of the flow exit
        !        xstart,ystart - location of the source
        !        xdist - the sum of all flowpath lenghs for all sources
        !        xxxx - temporary variable recording the flowpath indirectness (a
        !               'sinuosity' measure defined as the ratio of the flowpath length to
        !               the straight-line distance from source to flow exit
        !        indsum - the sum of flowpath indirectness for each source
        !        indssq - the sum of squares if flowpath indirectness for each source
        !        nindir - the total number of flowpath indirectness measurements (equals the
        !                   the total number of sources
        !      *********************************************************************
        NINDIR=0
        XDIST=0.0
        INDSUM=0.0
        INDSSQ=0.0
        DO  K=1,KSEGS
            IF (SOURCE(K)) THEN
                XSTART=IILOC(K)
                YSTART=JJLOC(K)
                XSUM=0.0
                L=K
                L115: DO
                    MAGNITUDE(L)=MAGNITUDE(L)+1
                    IF (FLOW_DIRECTION(IILOC(L),JJLOC(L)) > 5) THEN
                        XSUM=XSUM+SQRTOFTWO
                    ELSE
                        XSUM=XSUM+1.0
                    ENDIF
                    XEND=IILOC(L)
                    YEND=JJLOC(L)
                    L=TOSEGMENT(L)
                    !      *********************************************************************
                    !        We continue downstream until we reach a flow exit (tosegment=0)
                    !      *********************************************************************
                    IF (L <= 0) EXIT L115 
                ENDDO L115
                NINDIR=NINDIR+1
                XDIST=XDIST+XSUM
                IF ((XEND /= XSTART).OR.(YEND /= YSTART)) THEN
                    XXXX=XSUM/SQRT((XEND-XSTART)**2+(YEND-YSTART)**2)
                    INDSUM=INDSUM+XXXX
                    INDSSQ=INDSSQ+XXXX*XXXX
                ELSE
                    INDSUM=INDSUM+1.0
                    INDSSQ=INDSSQ+1.0
                ENDIF
            ENDIF
        ENDDO
        !      *********************************************************************
        !        The next do loop determines the strahler stream order of the stream network
        !         by going downstream from each source, and incrementing the order of the
        !         next stream segment only if its existing order is less than the current
        !         stream segment.  It also sets the stream order of the next stream segment
        !         to equal the present stream segment if the next stream segment so far has
        !         a lesser order
        !      *********************************************************************
        DO  K=1,KSEGS
            IF (SOURCE(K)) THEN
                PRESORDER = 1
                ORDER(K)=1
                L=K
                L117: DO
                    IF (L <= 0) EXIT L117
                    L=TOSEGMENT(L)
                    IF (L > 0) THEN
                        IF (ORDER(L) == PRESORDER) THEN
                            PRESORDER = PRESORDER+1
                        ENDIF
                        IF (ORDER(L) < PRESORDER) THEN
                            ORDER(L)=PRESORDER
                        ELSE
                            EXIT L117
                        ENDIF
                    ENDIF
                ENDDO L117
            ENDIF
        ENDDO
        !      *********************************************************************
        !        This subroutine calculates length statistics for exterior and
        !         interior stream links.  the variables are
        !            nextern - the number of stream sources
        !            xsum - the length of an individual stream link
        !            xextsum - the total length of individual exterior stream links
        !            xextssq - the sum of squares of exterior stream link lengths
        !            xintsum - the total length of individual interior stream links
        !            xintssq - the sum of squares of interior stream link lengths
        !            xmagsum - the sum of magnitudes of interior links
        !            xmagssq - the sum of squares of interior stream magnitudes
        !            xintmag - the sum of the cross products of stream magnitude and
        !                      interior link length
        !            ninter  - the total number of interior links
        !         The vector 'processed' makes sure that we process interior links
        !         only once
        !      *********************************************************************
        NEXTERN=0
        NINTERN=0
        XEXTSUM=0.0
        XEXTSSQ=0.0
        XMAGSUM=0.0
        XMAGSSQ=0.0
        XINTSUM=0.0
        XINTSSQ=0.0
        XINTMAG=0.0
        L120: DO  K=1,KSEGS
            IF (SOURCE(K)) THEN
                NEXTERN=NEXTERN+1
                L=K
                XSUM=0.0
                L125: DO
                    IF (FLOW_DIRECTION(IILOC(L),JJLOC(L)) > 5) THEN
                        XSUM=XSUM+SQRTOFTWO
                    ELSE
                        XSUM=XSUM+1.0
                    ENDIF
                    L=TOSEGMENT(L)
                    IF (L == 0) GOTO 1120
                    IF (MAGNITUDE(L) > 1) THEN
                        XEXTSUM=XEXTSUM+XSUM
                        XEXTSSQ=XEXTSSQ+XSUM*XSUM
                        EXIT L125
                    ENDIF
                ENDDO L125
                L130: DO
                    IF (.NOT.PROCESSED(L)) THEN
                        PROCESSED(L)=.TRUE.
                        NINTERN=NINTERN+1
                        THISMAG=MAGNITUDE(L)
                        XSUM=0.0
                        L135: DO
                            IF (FLOW_DIRECTION(IILOC(L),JJLOC(L)) > 5) THEN
                                XSUM=XSUM+SQRTOFTWO
                            ELSE
                                XSUM=XSUM+1.0
                            ENDIF
                            L=TOSEGMENT(L)
                            IF (L == 0) GOTO 1120
                            IF (MAGNITUDE(L) > THISMAG) THEN
                                XINTSUM=XINTSUM+XSUM
                                XINTSSQ=XINTSSQ+XSUM*XSUM
                                XINTMAG=XINTMAG+XSUM*THISMAG
                                XMAGSUM=XMAGSUM+THISMAG
                                XMAGSSQ=XMAGSSQ+THISMAG*THISMAG
                                GOTO 1130
                            ELSE
                                PROCESSED(L)=.TRUE.
                            ENDIF
                        ENDDO L135
                    ENDIF
                    GOTO 1120
                    1130 CONTINUE
                ENDDO L130
                  
            ENDIF
            1120 CONTINUE
        ENDDO L120
        !      *********************************************************************
        !        This loop calculates a number of statistics on first order streams.
        !         for first order basins we calculate
        !         the following:
        !            atop - the contributing area to the stream source
        !            abot - the total contributing area to the first order basin
        !            etop - the elevation of the stream at its source
        !            ebot - the elevation of the mouth of the first order basin
        !            gtop - the stream gradient at its source
        !            gbot - the stream gradient at the mouth of the first order basin
        !            xsum - the length of the first order stream
        !            nsource - the total number of first order streams that join other first
        !                      order streams at their mouth
        !            nts     - the total number of first order streams that join higher
        !                      streams greater than first order at their mouth
        !         These variables are used to define the following derived variables for
        !         first order streams:
        !            linkden - the 'link density' the ratio of the square of the 1st order
        !                      channel length to the total basin area
        !            arrat   - the 'area ratio' is the ratio of the difference in
        !                      contributing area between the mouth and source of the basin
        !                      the square of the 1st order stream length
        !            linkslp - the average gradient of the stream link
        !         These variables are then summed over all sources into the following
        !         variables:
        !            nsamp - the total number of first order basins samples
        !            ldensum - sum of link densities
        !            ldesssq - sum of squares of link density
        !            arratsum - sum of area ratios
        !            arratssq - sum of squares of area ratios
        !            ratarsum - sum of the ratio of contributing area at the basin mouth to
        !                       that at the basin source
        !            ldenslp  - the cross product of link density and link gradient
        !            arratslp - the cross product of area ratio and link gradient
        !            slpsum  - the sum of average link gradient
        !            slpssq  - the sum of squares of average link gradient
        !      *********************************************************************
        NSOURCE=0.0
        NTS=0.0
        NSAMP=0
        LDENSUM=0.0
        LDENSSQ=0.0
        ARRATSUM=0.0
        ARRATSSQ=0.0
        RATARSUM=0.0
        RATARSSQ=0.0
        LDENSLP=0.0
        ARRATSLP=0.0
        SLPSUM=0.0
        SLPSSQ=0.0
        L131: DO  K=1,KSEGS
            IF (SOURCE(K)) THEN
                L=K
                ATOP = DRAINAGE_AREA(IILOC(K),JJLOC(K))
                ETOP = ELEVATION(IILOC(K),JJLOC(K))
                GTOP = D8_GRADIENT(IILOC(K),JJLOC(K))
                XSUM=0.0
                L132: DO
                    IF (FLOW_DIRECTION(IILOC(L),JJLOC(L)) > 5) THEN
                        XSUM=XSUM+SQRTOFTWO
                    ELSE
                        XSUM=XSUM+1.0
                    ENDIF
                    ABOT = DRAINAGE_AREA(IILOC(L),JJLOC(L))
                    EBOT = ELEVATION(IILOC(L),JJLOC(L))
                    GBOT = D8_GRADIENT(IILOC(L),JJLOC(L))
                    L=TOSEGMENT(L)
                    IF (L == 0) EXIT L132
                    IF (MAGNITUDE(L) == 1) CYCLE
                    EXIT L132
                ENDDO L132
                IF (L /= 0) THEN
                    IF (MAGNITUDE(L) == 2) THEN
                        NSOURCE=NSOURCE+1.0
                    ELSE
                        NTS=NTS+1.0
                    ENDIF
                ENDIF
                IF (ABOT == ATOP) ABOT=ATOP+1.0
                LINKDEN=XSUM*XSUM/ABOT
                ARRAT=(ABOT-ATOP)/XSUM**2
                IF (EBOT == ETOP) THEN
                    LINKSLP = GBOT
                ELSE
                    LINKSLP=(ETOP-EBOT)/XSUM
                ENDIF
                NSAMP=NSAMP+1
                LDENSUM=LDENSUM+LINKDEN
                LDENSSQ=LDENSSQ+LINKDEN*LINKDEN
                ARRATSUM=ARRATSUM+ARRAT
                ARRATSSQ=ARRATSSQ+ARRAT*ARRAT
                RATARSUM=RATARSUM+ABOT/ATOP
                RATARSSQ=RATARSSQ+(ABOT/ATOP)**2
                LDENSLP=LDENSLP+LINKDEN*LINKSLP
                ARRATSLP=ARRATSLP+ARRAT*LINKSLP
                SLPSUM=SLPSUM+LINKSLP
                SLPSSQ=SLPSSQ+LINKSLP*LINKSLP
            ENDIF
        ENDDO L131
        !      *********************************************************************
        !       The following loop just caculates the total length of channels in the
        !         network (tchanlen)
        !      *********************************************************************
        TCHANLEN=0.0
        DO  I=1,MX
            DO  J=1,MYY
                IF (IDO(I,J) > 0) THEN
                    IF (FLOW_DIRECTION(I,J) > 5) THEN
                        TCHANLEN=TCHANLEN+SQRTOFTWO
                    ELSE
                        TCHANLEN=TCHANLEN+1.0
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        !      *********************************************************************
        !        This calculates shreve's kappa: the ratio of the total number of stream
        !         links times total area divided by the square to total channel length
        !      *********************************************************************
        IF (TCHANLEN > 0.0) THEN
            KAPPA=(2*NEXTERN-1)*(MX*(MYY))/TCHANLEN**2
        ELSE
            KAPPA=0.0
        ENDIF
        !      *********************************************************************
        !        This calculates the 'form factor' the ratio of total channel length to the
        !         square of the average pathlength from source to exit
        !      *********************************************************************
        IF ((XDIST > 0.0).AND.(XDIST > 0.0)) THEN
            FORMFAC=TCHANLEN/((XDIST/NINDIR)**2)
        ELSE
            FORMFAC=0.0
        ENDIF
        !      *********************************************************************
        !        The rest of the subroutine calculate average and standard deviation of the
        !         variables defined above, as well as correlation coefficients between some
        !         variables, and prints out the variables.  See the output statements for
        !         individual descriptions
        !      *********************************************************************
        SUMXY=ARRATSLP
        SUMX=ARRATSUM
        SUMY=SLPSUM
        SUMXX=ARRATSSQ
        SUMYY=SLPSSQ
        TOTAL=NSAMP
        IF (TOTAL > 1.0) THEN
            CORREL=(SUMXY-SUMX*SUMY/TOTAL)/(SQRT((SUMXX-SUMX*SUMX/TOTAL) &
            )*SQRT((SUMYY-SUMY*SUMY/TOTAL)))
        ELSE
            CORREL=0.0
        ENDIF
        WRITE(OUTHIST,700)CORREL,TOTAL
        WRITE(OUTSUMMARY,906)CORREL,TOTAL
        700 FORMAT(' CORR. - EXT. AREA RATIO AND SLOPE=', &
        G12.5,' N=',G13.6)
        SUMXY=LDENSLP
        SUMX=LDENSUM
        SUMY=SLPSUM
        SUMXX=LDENSSQ
        SUMYY=SLPSSQ
        TOTAL=NSAMP
        IF (TOTAL > 1.0) THEN
            CORREL=(SUMXY-SUMX*SUMY/TOTAL)/(SQRT((SUMXX-SUMX*SUMX/TOTAL) &
            )*SQRT((SUMYY-SUMY*SUMY/TOTAL)))
        ELSE
            CORREL=0.0
        ENDIF
        WRITE(OUTHIST,705)CORREL,TOTAL
        WRITE(OUTSUMMARY,906)CORREL,TOTAL
        705 FORMAT(' CORR. - EXT. LINK DENSITY AND SLOPE=', &
        G12.5,' N=',G13.6)
        SUMXY=XINTMAG
        SUMX=XMAGSUM
        SUMY=XINTSUM
        SUMXX=XMAGSSQ
        SUMYY=XINTSSQ
        TOTAL=NINTERN
        IF (TOTAL > 1.0) THEN
            CORREL=(SUMXY-SUMX*SUMY/TOTAL)/(SQRT((SUMXX-SUMX*SUMX/TOTAL) &
            )*SQRT((SUMYY-SUMY*SUMY/TOTAL)))
        ELSE
            CORREL=0.0
        ENDIF
        WRITE(OUTHIST,710)CORREL,TOTAL
        WRITE(OUTSUMMARY,906)CORREL,TOTAL
        710 FORMAT(' CORR. - INT. LINK LENGTH AND MAG.=',  &
        G12.5,' N=',G13.6)
        WRITE(OUTHIST,715) KAPPA,FORMFAC
        715 FORMAT(' SHREVE KAPPA=',G12.5,' FORMFACTOR=',G12.5)
        WRITE(OUTSUMMARY,906) KAPPA,FORMFAC
        IF ((NTS+NSOURCE) > 0.0) THEN
            NTS = NTS/(NSOURCE+NTS)
        ELSE
            NTS = 0.0
        ENDIF
        WRITE(OUTHIST,720) NTS
        720 FORMAT(' FRACTION OF T-S LINKS=',G12.5)
        WRITE(OUTSUMMARY,907) NTS
        IF ((NEXTERN > 0.0).AND.(XINTSUM > 0.0)) THEN
            LRATIO=(XEXTSUM/NEXTERN)/(XINTSUM/(NEXTERN-1))
        ELSE
            LRATIO=0.0
        ENDIF
        LTOT=NINDIR
        IF (LTOT > 1) THEN
            TEMP=INDSUM/LTOT
            TEMP1=(INDSSQ-LTOT*TEMP*TEMP)/(LTOT-1.0)
        ELSE
            TEMP=0.0
            TEMP1=0.0
        ENDIF
        WRITE (OUTHIST,725) TEMP,TEMP1,LTOT
        WRITE (OUTSUMMARY,908) TEMP,TEMP1,LTOT
        725 FORMAT (' AV VAR & N   FLOWPATH INDIR.=',3E15.6)
        LTOT=NEXTERN
        IF (LTOT > 1) THEN
            TEMP=XEXTSUM/LTOT
            TEMP1=(XEXTSSQ-LTOT*TEMP*TEMP)/(LTOT-1.0)
        ELSE
            TEMP=0.0
            TEMP1=0.0
        ENDIF
        WRITE (OUTHIST,730) TEMP,TEMP1,LTOT
        WRITE (OUTSUMMARY,908) TEMP,TEMP1,LTOT
        730 FORMAT (' AV VAR & N   EXT. LINK LENGTH=',3E15.6)
        LTOT=NEXTERN-1
        IF (LTOT > 1) THEN
            TEMP=XINTSUM/LTOT
            TEMP1=(XINTSSQ-LTOT*TEMP*TEMP)/(LTOT-1.0)
        ELSE
            TEMP=0.0
            TEMP1=0.0
        ENDIF
        WRITE (OUTHIST,735) TEMP,TEMP1,LTOT
        WRITE (OUTSUMMARY,908) TEMP,TEMP1,LTOT
        735 FORMAT (' AV VAR & N   INT. LINK LENGTH=',3E15.6)
        WRITE (OUTHIST,740) LRATIO
        740 FORMAT(' EXT.-INT. LINK LENGTH RATIO=',G12.5)
        WRITE (OUTSUMMARY,907) LRATIO
        LTOT=NSAMP
        IF (LTOT > 1) THEN
            TEMP=LDENSUM/LTOT
            TEMP1=(LDENSSQ-LTOT*TEMP*TEMP)/(LTOT-1.0)
        ELSE
            TEMP=0.0
            TEMP1=0.0
        ENDIF
        WRITE (OUTHIST,745) TEMP,TEMP1,LTOT
        WRITE (OUTSUMMARY,908) TEMP,TEMP1,LTOT
        745 FORMAT (' AV VAR & N   1ST ORD. LINK DENSITY=',3E15.6)
        LTOT=NSAMP
        IF (LTOT > 1) THEN
            TEMP=ARRATSUM/LTOT
            TEMP1=(ARRATSSQ-LTOT*TEMP*TEMP)/(LTOT-1.0)
        ELSE
            TEMP=0.0
            TEMP1=0.0
        ENDIF
        WRITE (OUTHIST,750) TEMP,TEMP1,LTOT
        WRITE (OUTSUMMARY,908) TEMP,TEMP1,LTOT
        750 FORMAT (' AV VAR & N   1ST ORD. DIM. AREA =',3E15.6)
        LTOT=NSAMP
        IF (LTOT > 1) THEN
            TEMP=RATARSUM/LTOT
            TEMP1=(RATARSSQ-LTOT*TEMP*TEMP)/(LTOT-1.0)
        ELSE
            TEMP=0.0
            TEMP1=0.0
        ENDIF
        WRITE (OUTHIST,755) TEMP,TEMP1,LTOT
        WRITE (OUTSUMMARY,908) TEMP,TEMP1,LTOT
        755 FORMAT (' AV VAR & N   1ST ORD. AREA RATIO=',3E15.6)
        LTOT=NSAMP
        IF (LTOT > 1) THEN
            TEMP=SLPSUM/LTOT
            TEMP1=(SLPSSQ-LTOT*TEMP*TEMP)/(LTOT-1.0)
        ELSE
            TEMP=0.0
            TEMP1=0.0
        ENDIF
        WRITE (OUTHIST,760) TEMP,TEMP1,LTOT
        WRITE (OUTSUMMARY,908) TEMP,TEMP1,LTOT
        DEALLOCATE(ORDER,MAGNITUDE,IILOC,JJLOC,TOSEGMENT)
        DEALLOCATE(SOURCE,PROCESSED)
        760 FORMAT (' AV VAR & N   1ST ORD. LINK SLOPE=',3E15.6)
        906 FORMAT(2(' ',E15.8))
        907 FORMAT(' ',E15.8)
        908 FORMAT(3(' ',E15.8))
        RETURN
    END ! SUBROUTINE FIND_STREAM_NETWORK()

