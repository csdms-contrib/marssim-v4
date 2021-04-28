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
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    SUBROUTINE FIND_DOWNSTREAM_LOCATION(II,JI,IIO,JJO,LOCAL_DIRECTION,IS_OK)
        USE ERODE_GLOBALS
        USE SEDROUTE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: II,JI
        INTEGER, INTENT (OUT) :: IIO,JJO
        INTEGER, INTENT (OUT) :: LOCAL_DIRECTION
        INTEGER :: INDEPRESS
        INTEGER :: I1,I2,LLL,KKK,IOX,JOX,KKKUSE,NNN
        INTEGER :: IDEP
        REAL(4) ::  FF,EL1,EL2,MINGRAD
        REAL(4) ::  XGRAD(9)
        LOGICAL :: KOK(9),ISCIRCLE,IS_OK
        LOGICAL :: ROUNDABOUT,ATNEXTPOINT
        INTEGER :: IMUST,IROUND,IBACK,IDONE,ICIRCLE,IQIN,IQNEXT,IDEPR,ISUB,IENC
        !     **************************************************************
        !      This subroutine find the next location downstream for sediment routing
        !  MODIFIES: IS_DEPRESSION,DEFAULT_SEDIMENT_YIELD, DEFAULT_DRAINAGE_AREA
        !            ISTART, JSTART, IEND, JEND, SEDIMENT_CARRYOVER, DO_BACKTRACK
        !     **************************************************************
        ROUNDABOUT=.TRUE.
        IS_OK=.TRUE.
        IQIN=FLOW_DIRECTION(II,JI)
        IF (IS_DEPRESSION) THEN
            INDEPRESS=1
            GO TO 100
        ELSE
            INDEPRESS=0
            IF (FLOW_DIRECTION(II,JI) > 0) THEN
                LOCAL_DIRECTION=FLOW_DIRECTION(II,JI)
                IQNEXT=LOCAL_DIRECTION
            ELSE
                LOCAL_DIRECTION=-FLOW_DIRECTION(II,JI)
                IQNEXT=LOCAL_DIRECTION
                IS_DEPRESSION=.TRUE.
                IF (.NOT.COMPLETE_RUNOFF) MUST_SEARCH=.TRUE.
                DEFAULT_SEDIMENT_YIELD=SEDIMENT_CARRYOVER+SEDIMENT_YIELD(II,JI)
                DEFAULT_DRAINAGE_AREA=MAXIMUM_DISCHARGE(II,JI)
                IF (COMPLETE_RUNOFF) THEN
                    ISTART=II
                    JSTART=JI
                    INOW=II
                    JNOW=JI
                    NNN=BASIN_NUMBER(II,JI)
                    IEND=I_OUTFLOW(NNN)
                    JEND=J_OUTFLOW(NNN)
                    IDMAX=ABS(ISTART-IEND)+ABS(JSTART-JEND)+2
                    ICOUNT=1
                    IF (IS_X_PERIODIC) THEN
                        IF (ABS(IEND-ISTART) > (MX/2)) THEN
                            IF (ISTART < (MX/2)) THEN
                                IEND=IEND-MX
                            ELSE
                                IEND=IEND+MX
                            ENDIF
                         !   WRITE(*,8224)
                            8224 FORMAT('#')
                        ENDIF
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (ABS(JEND-JSTART) > (MY/2)) THEN
                            IF (JSTART < (MY/2)) THEN
                                JEND=JEND-MY
                            ELSE
                                JEND=JEND+MY
                            ENDIF
                         !   WRITE(*,8224)
                        ENDIF
                    ENDIF
                    IF ((ABS(ISTART-IEND) < 2).AND.(ABS(JSTART-JEND)     &
                     < 2))  THEN
                        IIO=IEND
                        JJO=JEND
                        IS_DEPRESSION=.FALSE.
                        SEDIMENT_CARRYOVER=DEFAULT_SEDIMENT_YIELD
                        LOCAL_DIRECTION=1
                        IF (IS_X_PERIODIC) THEN
                            IF (IIO < 1) IIO=MX
                            IF (IIO > MX) IIO=1
                        ENDIF
                        IF (IS_Y_PERIODIC) THEN
                            IF (JJO < 1) JJO=MY
                            IF (JJO > MY) JJO=1
                        ENDIF
                        GO TO 500
                    ELSE
                        IDX=IEND-ISTART
                        JDY=JEND-JSTART
                        IXX=ABS(IDX)
                        JYY=ABS(JDY)
                        INC1=MAX(IXX,JYY)
                        IXT=0
                        JYT=0
                        IIII=1
                        GO TO 100
                    ENDIF
                ELSE
                    JJO=JI+DOWNSTREAM(LOCAL_DIRECTION,2)
                    IIO=II+DOWNSTREAM(LOCAL_DIRECTION,1)
                    IF (IS_X_PERIODIC) THEN
                        IF (IIO < 1) IIO=MX
                        IF (IIO > MX) IIO=1
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JJO < 1) JJO=MY
                        IF (JJO > MY) JJO=1
                    ENDIF
                    IS_DEPRESSION=.FALSE.
                    SEDIMENT_CARRYOVER=DEFAULT_SEDIMENT_YIELD
                    ISCIRCLE=.FALSE.
                    IF (ROUNDABOUT) THEN
                        L1588: DO LLL=1,K
                            IF ((IIO == I_LOCATION(LLL)).AND.(JJO == J_LOCATION(LLL))) &
                            DO_BACKTRACK=.TRUE.
                        ENDDO L1588
                        IF (DONE_ONCE(IIO,JJO)) DO_BACKTRACK=.TRUE.
                        IF (DO_BACKTRACK) THEN
                            FF=1.0
                            EL1=ELEVATION(IIO,JJO)
                            L1589: DO KKK=2,9
                                KOK(KKK)=.TRUE.
                                XGRAD(KKK)=1.0E+25
                                IOX=II+DOWNSTREAM(KKK,1)
                                JOX=JI+DOWNSTREAM(KKK,2)
                                IF (KKK > 5) FF=0.7071
                                IF (IS_X_PERIODIC) THEN
                                    IF (IOX > MX) IOX=1
                                    IF (IOX < 1) IOX=MX
                                ELSE
                                    IF ((IOX < 1).OR.(IOX > MX)) THEN
                                        KOK(KKK)=.FALSE.
                                        CYCLE L1589
                                    ENDIF
                                ENDIF
                                IF (IS_Y_PERIODIC) THEN
                                    IF (JOX > MY) JOX=1
                                    IF (JOX < 1) JOX=MY
                                ELSE
                                    IF ((JOX < 1).OR.(JOX > MY)) THEN
                                        KOK(KKK)=.FALSE.
                                        CYCLE L1589
                                    ENDIF
                                ENDIF
                                EL2=ELEVATION(IOX,JOX)
                                DO LLL=1,K
                                    IF ((IOX == I_LOCATION(LLL)).AND.(JOX == J_LOCATION(LLL))) &
                                    KOK(KKK)=.FALSE.
                                    IF (DONE_ONCE(I_LOCATION(LLL),J_LOCATION(LLL))) KOK(KKK)=.FALSE.
                                ENDDO
                                XGRAD(KKK)=(EL1-EL2)*FF/CELL_SIZE
                            ENDDO L1589
                            MINGRAD=1.0E+25
                            KKKUSE=0
                            L1591: DO KKK=2,9
                                IF (KOK(KKK)) THEN
                                    IF (XGRAD(KKK) < MINGRAD) THEN
                                        MINGRAD=XGRAD(KKK)
                                        KKKUSE=KKK
                                    ENDIF
                                ENDIF
                            ENDDO L1591
                            IF (KKKUSE > 0) THEN
                                LOCAL_DIRECTION=KKKUSE
                                IIO=II+DOWNSTREAM(KKKUSE,1)
                                JJO=JI+DOWNSTREAM(KKKUSE,2)
                                IF (IS_X_PERIODIC) THEN
                                    IF (IIO > MX) IIO=1
                                    IF (IIO < 1) IIO=MX
                                ENDIF
                                IF (IS_Y_PERIODIC) THEN
                                    IF (JJO > MY) JJO=1
                                    IF (JJO < 1) JJO=MY
                                ENDIF
                                ISCIRCLE=.FALSE.
                            ELSE
                                ISCIRCLE=.TRUE.
                            ENDIF
                        ENDIF
                        IF (.NOT.ISCIRCLE) THEN
                            DO_BACKTRACK=.FALSE.
                        ELSE
                            DO_BACKTRACK=.TRUE.
                        ENDIF
                    ENDIF
                    GOTO 500
                ENDIF
            ENDIF
            JJO=JI+DOWNSTREAM(LOCAL_DIRECTION,2)
            IIO=II+DOWNSTREAM(LOCAL_DIRECTION,1)
            IF (IS_X_PERIODIC) THEN
                IF (IIO < 1) IIO=MX
                IF (IIO > MX) IIO=1
            ENDIF
            IF (IS_Y_PERIODIC) THEN
                IF (JJO < 1) JJO=MY
                IF (JJO > MY) JJO=1
            ENDIF
            ISCIRCLE=.FALSE.
            IF (ROUNDABOUT) THEN
                IF (MUST_SEARCH) THEN
                    L588: DO LLL=1,K
                        IF ((IIO == I_LOCATION(LLL)).AND.(JJO == J_LOCATION(LLL))) &
                        DO_BACKTRACK=.TRUE.
                    ENDDO L588
                    IF (DONE_ONCE(IIO,JJO)) DO_BACKTRACK=.TRUE.
                    IF (DO_BACKTRACK) THEN
                        FF=1.0
                        EL1=ELEVATION(IIO,JJO)
                        L589: DO KKK=2,9
                            KOK(KKK)=.TRUE.
                            XGRAD(KKK)=1.0E+25
                            IOX=II+DOWNSTREAM(KKK,1)
                            JOX=JI+DOWNSTREAM(KKK,2)
                            IF (KKK > 5) FF=0.7071
                            IF (IS_X_PERIODIC) THEN
                                IF (IOX > MX) IOX=1
                                IF (IOX < 1) IOX=MX
                            ELSE
                                IF ((IOX < 1).OR.(IOX > MX)) THEN
                                    KOK(KKK)=.FALSE.
                                    CYCLE L589
                                ENDIF
                            ENDIF
                            IF (IS_Y_PERIODIC) THEN
                                IF (JOX > MY) JOX=1
                                IF (JOX < 1) JOX=MY
                            ELSE
                                IF ((JOX < 1).OR.(JOX > MY)) THEN
                                    KOK(KKK)=.FALSE.
                                    CYCLE L589
                                ENDIF
                            ENDIF
                            EL2=ELEVATION(IOX,JOX)
                            L590: DO LLL=1,K
                                IF ((IOX == I_LOCATION(LLL)).AND.(JOX == J_LOCATION(LLL)))&
                                KOK(KKK)=.FALSE.
                                IF (DONE_ONCE(I_LOCATION(LLL),J_LOCATION(LLL))) KOK(KKK)=.FALSE.
                            ENDDO L590
                            XGRAD(KKK)=(EL1-EL2)*FF/CELL_SIZE
                        ENDDO L589
                        MINGRAD=1.0E+25
                        KKKUSE=0
                        L591: DO KKK=2,9
                            IF (KOK(KKK)) THEN
                                IF (XGRAD(KKK) < MINGRAD) THEN
                                    MINGRAD=XGRAD(KKK)
                                    KKKUSE=KKK
                                ENDIF
                            ENDIF
                        ENDDO L591
                        IF (KKKUSE > 0) THEN
                            LOCAL_DIRECTION=KKKUSE
                            IIO=II+DOWNSTREAM(KKKUSE,1)
                            JJO=JI+DOWNSTREAM(KKKUSE,2)
                            IF (IS_X_PERIODIC) THEN
                                IF (IIO > MX) IIO=1
                                IF (IIO < 1) IIO=MX
                            ENDIF
                            IF (IS_Y_PERIODIC) THEN
                                IF (JJO > MY) JJO=1
                                IF (JJO < 1) JJO=MY
                            ENDIF
                            ISCIRCLE=.FALSE.
                        ELSE
                            ISCIRCLE=.TRUE.
                        ENDIF
                    ENDIF
                    IF (.NOT.ISCIRCLE) THEN
                        DO_BACKTRACK=.FALSE.
                    ELSE
                        DO_BACKTRACK=.TRUE.
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
        IF ((II == IIO).AND.(JI == JJO)) THEN
            IF (IS_DEPRESSION) THEN
                IDEP=1
            ELSE
                IDEP=0
            ENDIF
            IF (ENCLOSED(NNN)) THEN
                IENC=1
            ELSE
                IENC=0
            ENDIF
            IF (SUBMERGED(II,JI)) THEN
                ISUB=1
            ELSE
                ISUB=0
            ENDIF
            IF (WRITE_SEDIMENT_DIAGNOSTICS) THEN
            OPEN(OUTSEDDEBUG,FILE='SEDTEMP.DAT',POSITION='APPEND')
            WRITE(OUTSEDDEBUG,707) II,JI,FLOW_DIRECTION(IIO,JJO),NNN,ISTART,JSTART, &
            IEND,JEND,INOW,JNOW!,DONE_ONCE(IIO,JJO)
            WRITE(*,1708) IENC,ISUB,ELEVATION(II,JI),I_OUTFLOW(NNN),J_OUTFLOW(NNN), &
            BASIN_DRAINAGE_AREA(NNN),LAKE_OUTLET_ELEVATION(NNN),LAKE_SURFACE_ELEVATION(NNN)
            707 FORMAT(' NO NEW NEXTPOINT I=',I5,' J=',I5,' ID=',I5,/,    &
            ' NNN=',I5,' IS=',I5,' JS=',I5,' IE=',I5,' JE=',I5, &
            ' INOW=',I5,' JNOW=',I5)!,' DONEONCE=',I5)
            1708 FORMAT(' ENCL=',I5,' SUBM=',I5,' EL=',G12.5,/, &
            ' IOUT=',I5,' JOUT=',I5,' DR_AREA=',G12.5,' LAKE_OUT=',G12.5,' LAKE_EL=',G12.5)
            IS_OK=.FALSE.
            !STOP
            CLOSE(OUTSEDDEBUG)
            ENDIF
        ENDIF
        RETURN
        100   CONTINUE
		IF (COMPLETE_RUNOFF) THEN
        ATNEXTPOINT=.FALSE.
        IXT=IXT+IXX
        JYT=JYT+JYY
        IF (IXT > INC1) THEN
            ATNEXTPOINT=.TRUE.
            IXT=IXT-INC1
            INOW=INOW+SIGN(1,IDX)
        ENDIF
        IF (JYT > INC1) THEN
            ATNEXTPOINT=.TRUE.
            JYT=JYT-INC1
            JNOW=JNOW+SIGN(1,JDY)
        ENDIF
        IF (ATNEXTPOINT) THEN
            IF ((INOW == ISTART).AND. &
            (JNOW == JSTART)) THEN
                IIII=IIII+1
                IF (IIII == INC1) THEN
                    JNOW=JEND
                    INOW=IEND
                    GOTO 102
                ENDIF
                GOTO 100
            ENDIF
            IF ((INOW == IEND).AND.(JNOW == JEND)) THEN
                GOTO 102
            ENDIF
            IIII=IIII+1
            IF (IIII > INC1) THEN
                INOW=IEND
                JNOW=JEND
            ENDIF
            GOTO 102
        ELSE
            IIII=IIII+1
            IF (IIII <= INC1) THEN
                GOTO 100
            ELSE
                INOW=IEND
                JNOW=JEND
            ENDIF
        ENDIF
        102   CONTINUE
        IIO=INOW
        JJO=JNOW
        IF ((IIO == IEND).AND.(JJO == JEND)) THEN
            IS_DEPRESSION=.FALSE.
            SEDIMENT_CARRYOVER=DEFAULT_SEDIMENT_YIELD
        ENDIF
        LOCAL_DIRECTION=1
        IF (IS_X_PERIODIC) THEN
            IF (IIO < 1) IIO=IIO+MX
            IF (IIO > MX) IIO=IIO-MX
        ENDIF
        IF (IS_Y_PERIODIC) THEN
            IF (JJO < 1) JJO=JJO+MY
            IF (JJO > MY) JJO=JJO-MY
        ENDIF
        ICOUNT=ICOUNT+1
        IF (ICOUNT > IDMAX) THEN
            IF (WRITE_SEDIMENT_DIAGNOSTICS) THEN
            OPEN(OUTSEDDEBUG,FILE='SEDTEMP.DAT',POSITION='APPEND')
            NNN=BASIN_NUMBER(ISTART,JSTART)
            WRITE(OUTSEDDEBUG,200) ISTART,JSTART,IEND,JEND,INOW,JNOW,II,JI,IIO,&
            JJO,FLOW_DIRECTION(ISTART,JSTART),NNN,I_OUTFLOW(NNN),J_OUTFLOW(NNN)
            200 FORMAT(' TROUBLE IN DEPRESSION ROUTING',/,' IST=',I5,    &
            ' JST=',I5,&
            ' IEND=',I5,' JEND=',I5,/,' II=',I5,' JJ=',I5, 'IO=',I5,&
            ' JJO=',I5,/,' IDSTART=',I5,'BASININ=',I5,' IOUT=',I5,&
            ' JOUT=',I5)
            NNN=BASIN_NUMBER(IEND,JEND)
            WRITE(OUTSEDDEBUG,210) FLOW_DIRECTION(IEND,JEND),FLOW_DIRECTION(INOW,JNOW),FLOW_DIRECTION(IIO,JJO),NNN
            210 FORMAT(' IDOUT=',I5,' IDNOW=',I5,' ID=',I5,' BASINOUT=',I5)
            CALL WRITE_DEBUG_DATA(IIO,JJO)
            CLOSE(OUTSEDDEBUG)
            ENDIF
            IS_OK=.FALSE.
            RETURN
           ! STOP
        		
        ENDIF
        ELSE
           IF (WRITE_SEDIMENT_DIAGNOSTICS) THEN
           OPEN(OUTSEDDEBUG,FILE='SEDTEMP.DAT',POSITION='APPEND')
           WRITE(OUTSEDDEBUG,1200) ISTART,JSTART,IEND,JEND,INOW,JNOW,II,JI,IIO,&
            JJO,FLOW_DIRECTION(ISTART,JSTART),NNN,I_OUTFLOW(NNN),J_OUTFLOW(NNN)
            1200 FORMAT(' NO NEW LOCATION FOUND',/,' IST=',I5,    &
            ' JST=',I5,&
            ' IEND=',I5,' JEND=',I5,/,' II=',I5,' JJ=',I5, 'IO=',I5,&
            ' JJO=',I5,/,' IDSTART=',I5,'BASININ=',I5,' IOUT=',I5,&
            ' JOUT=',I5)
            NNN=BASIN_NUMBER(IEND,JEND)
            WRITE(OUTSEDDEBUG,1210) FLOW_DIRECTION(IEND,JEND),FLOW_DIRECTION(INOW,JNOW),FLOW_DIRECTION(IIO,JJO),NNN
           1210 FORMAT(' IDOUT=',I5,' IDNOW=',I5,' ID=',I5,' BASINOUT=',I5)
            CLOSE (OUTSEDDEBUG)
            ENDIF
		   RETURN
		ENDIF
        500 CONTINUE
        IF ((IIO > MX).OR.(IIO < 1).OR.(JJO < 1).OR.(JJO > MY)) THEN
            IF (IS_DEPRESSION) THEN
                IDEPR=1
            ELSE
                IDEPR=0
            ENDIF
            IF (MUST_SEARCH) THEN
                IMUST=1
            ELSE
                IMUST=0
            ENDIF
            IF (ROUNDABOUT) THEN
                IROUND=1
            ELSE
                IROUND=0
            ENDIF
            IF (DO_BACKTRACK) THEN
                IBACK=1
            ELSE
                IBACK=0
            ENDIF
            IF (DONE_ONCE(II,JI))    THEN
                IDONE=1
            ELSE
                IDONE=0
            ENDIF
            IF (ISCIRCLE) THEN
                ICIRCLE=1
            ELSE
                ICIRCLE=0
            ENDIF
            IF (ENCLOSED(NNN)) THEN
                IENC=1
            ELSE
                IENC=0
            ENDIF
            IF (SUBMERGED(II,JI)) THEN
                ISUB=1
            ELSE
                ISUB=0
            ENDIF
            NNN=BASIN_NUMBER(IEND,JEND)
            I1=BASIN_NUMBER(II,JI)
            I2=BASIN_NUMBER(ISTART,JSTART)
            IF (WRITE_SEDIMENT_DIAGNOSTICS) THEN
            OPEN(OUTSEDDEBUG,FILE='SEDTEMP.DAT',POSITION='APPEND')
            WRITE(OUTSEDDEBUG,510) II,JI,IIO,JJO,LOCAL_DIRECTION,FLOW_DIRECTION(II,JI),ISTART,JSTART,&
            IEND,JEND,INOW,JNOW,IDEPR,DEFAULT_SEDIMENT_YIELD,SEDIMENT_CARRYOVER,&
            DEFAULT_DRAINAGE_AREA,IMUST,IROUND,IBACK,IDONE,ICIRCLE,INDEPRESS,&
            IQIN,IQNEXT
            510 FORMAT(' @@@@@ OUT OF BOUNDS @@@@ FROM II=',I5,' JJ=',I5,    &
            ' IIO=',I5,' JJO=',I5,' LOCAL_DIRECTION=',I5,' ID=',I5,/,' ISTART=',&
            I5,' JSTART=',I5,' IEND=',I5,' JEND=',I5,' INOW=',I5,&
            ' JNOW=',I5,&
            ' DEPR=',I5,/,' DEFSY=',G12.4,' SEDCARRY=',G12.5,&
            ' DEFAULT_DRAINAGE_AREA=',G12.5,/,' IMUST=',I5,' IROUND=',I5,&
            ' IBACK=',I5,' IDONE=',I5,' ICIRCLE=',I5,&
            ' INDEPRESS=',I5,' IQIN',I5,' IQNEXT',I5)
            WRITE(OUTSEDDEBUG,511) NNN,I1,I2,INCRX,INCRY,ICOUNT
            511 FORMAT(' BN(IEND,JEND)=',I8,' BN(II,JI)=',I8,&
            ' BN(ISTART,JSTART=',I8,/,' INCRX=',I5,' INCRY=',I5,&
            ' ICOUNT=',I5)
            CALL WRITE_DEBUG_DATA(IIO,JJO)
            CLOSE(OUTSEDDEBUG)
            ENDIF
            IS_OK=.FALSE.
            RETURN
            !STOP
        ENDIF
        IF ((II == IIO).AND.(JI == JJO)) THEN
            IF (IS_DEPRESSION) THEN
                IDEPR=1
            ELSE
                IDEPR=0
            ENDIF
            IF (MUST_SEARCH) THEN
                IMUST=1
            ELSE
                IMUST=0
            ENDIF
            IF (ROUNDABOUT) THEN
                IROUND=1
            ELSE
                IROUND=0
            ENDIF
            IF (DO_BACKTRACK) THEN
                IBACK=1
            ELSE
                IBACK=0
            ENDIF
            IF (DONE_ONCE(II,JI))    THEN
                IDONE=1
            ELSE
                IDONE=0
            ENDIF
            IF (ISCIRCLE) THEN
                ICIRCLE=1
            ELSE
                ICIRCLE=0
            ENDIF
            NNN=BASIN_NUMBER(IEND,JEND)
            I1=BASIN_NUMBER(II,JI)
            I2=BASIN_NUMBER(ISTART,JSTART)
            IF (ENCLOSED(NNN)) THEN
                IENC=1
            ELSE
                IENC=0
            ENDIF
            IF (SUBMERGED(II,JI)) THEN
                ISUB=1
            ELSE
                ISUB=0
            ENDIF
            IF (WRITE_SEDIMENT_DIAGNOSTICS) THEN
            OPEN(OUTSEDDEBUG,FILE='SEDTEMP.DAT',POSITION='APPEND')
            WRITE(*,707) II,IIO,FLOW_DIRECTION(IIO,JJO),NNN,ISTART,JSTART,&
            IEND,JEND,INOW,JNOW,DONE_ONCE(IIO,JJO)
            WRITE(*,708) II,JI,IIO,JJO,LOCAL_DIRECTION,FLOW_DIRECTION(II,JI),ISTART,JSTART,&
            IEND,JEND,INOW,JNOW,IDEPR,DEFAULT_SEDIMENT_YIELD,SEDIMENT_CARRYOVER, &
            DEFAULT_DRAINAGE_AREA,IMUST,IROUND,IBACK,IDONE,ICIRCLE,INDEPRESS,&
            IQIN,IQNEXT
            WRITE(OUTSEDDEBUG,511) NNN,I1,I2,INCRX,INCRY,ICOUNT
            WRITE(OUTHIST,707) II,IIO,FLOW_DIRECTION(IIO,JJO),NNN,ISTART,JSTART,&
            IEND,JEND,INOW,JNOW,DONE_ONCE(IIO,JJO)
            WRITE(OUTSEDDEBUG,708) II,JI,IIO,JJO,LOCAL_DIRECTION,FLOW_DIRECTION(II,JI),ISTART,JSTART,&
            IEND,JEND,INOW,JNOW,IDEPR,DEFAULT_SEDIMENT_YIELD,SEDIMENT_CARRYOVER, &
            DEFAULT_DRAINAGE_AREA,IMUST,IROUND,IBACK,IDONE,ICIRCLE,INDEPRESS,&
            IQIN,IQNEXT
            WRITE(OUTSEDDEBUG,511) NNN,I1,I2,INCRX,INCRY,ICOUNT
            708       FORMAT(' @@@@@ WENT NOWHERE @@@@ FROM II=',I5,' JJ=',I5,&
            ' IIO=',I5,' JJO=',I5,' LOCAL_DIRECTION=',I5,' ID=',I5,/,' ISTART=',&
            I5,' JSTART=',I5,' IEND=',I5,' JEND=',I5,' INOW=',I5,&
            ' JNOW=',I5,&
            ' DEPR=',I5,/,' DEFSY=',G12.4,' SEDCARRY=',G12.5,&
            ' DEFAULT_DRAINAGE_AREA=',G12.5,/,' IMUST=',I5,' IROUND=',I5,&
            ' IBACK=',I5,' IDONE=',I5,' ICIRCLE=',I5,&
            ' INDEPRESS=',I5,' IQIN',I5,' IQNEXT',I5)
            WRITE(*,1708) IENC,ISUB,ELEVATION(II,JI),I_OUTFLOW(NNN),J_OUTFLOW(NNN), &
            BASIN_DRAINAGE_AREA(NNN),LAKE_OUTLET_ELEVATION(NNN),LAKE_SURFACE_ELEVATION(NNN)
            CALL WRITE_DEBUG_DATA(II,JI)
            CLOSE(OUTSEDDEBUG)
            ENDIF
            IS_OK=.FALSE.
            !STOP
        ENDIF
        RETURN
	
    END ! SUBROUTINE FIND_DOWNSTREAM_LOCATION
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_DEBUG_DATA(II,JJ)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: II,JJ
        INTEGER :: I,J, ILOW, IHIGH, JLOW, JHIGH, ISUBM(5),IX
        ILOW=MAX(1,II-2)
        IHIGH=MIN(MX,II+2)
        JLOW=MAX(1,JJ-2)
        JHIGH=MIN(MY,JJ+2)
        WRITE(OUTHIST,100)
        100 FORMAT(' ID VALUES')
        DO J=JLOW,JHIGH
            WRITE(OUTSEDDEBUG,200) (FLOW_DIRECTION(I,J),I=ILOW,IHIGH)
            200 FORMAT(5(I2,' '))
        ENDDO
        WRITE(OUTSEDDEBUG,250)
        250 FORMAT(/,' SUBMERGE VALUES')
        DO J=JLOW,JHIGH
            ISUBM=0
            IX=1
            DO I=ILOW,IHIGH
               IF (SUBMERGED(I,J)) THEN
                   ISUBM(IX)=1
               ENDIF
               IX=IX+1
            ENDDO
            WRITE(OUTSEDDEBUG,200) (ISUBM(I),I=1,5)
        ENDDO
        WRITE(OUTSEDDEBUG,300)
        300 FORMAT(/,' ELEVATIONS')
        DO J=JLOW,JHIGH
            WRITE(OUTSEDDEBUG,350) (ELEVATION(I,J), I=ILOW,IHIGH)
            350 FORMAT(5(G12.5,' '))
        ENDDO
        WRITE(OUTSEDDEBUG,400)
        400 FORMAT(/,' GRADIENTS')
        DO J=JLOW,JHIGH
            WRITE(OUTSEDDEBUG,350) (D8_GRADIENT(I,J),I=ILOW,IHIGH)
        ENDDO    
        RETURN
    END ! SUBROUTINE WRITE_DEBUG_DATA
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE SEDIMENT_TRANSPORT_FLUX(BEDLOAD_FLUX,TRANSPORT_STAGE,GRADIENT,A1TERM,A2TERM)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        !        *********************************************************************
        !         This subroutine determines the potential sediment transport rate
        !           as a function of gradient and drainage area
        !        *********************************************************************
        REAL(4), INTENT(OUT):: BEDLOAD_FLUX,TRANSPORT_STAGE
        REAL(4), INTENT(IN) :: GRADIENT,A1TERM,A2TERM
        REAL(4):: XARGUMENT
        IF (GRADIENT > 0.0) THEN
            XARGUMENT=  (GRADIENT**SEDIMENT_GRADIENT_EXPONENT)*A2TERM*SEDIMENT_DISCHARGE_FACTOR  &
            -TRANSPORT_CRITICAL_DIM_SHEAR
            IF (XARGUMENT > 0.0) THEN
                BEDLOAD_FLUX=SEDIMENT_CONSTANT*A1TERM*  &
                XARGUMENT**SEDIMENT_TRANSPORT_EXPONENT
                TRANSPORT_STAGE=XARGUMENT/TRANSPORT_CRITICAL_DIM_SHEAR
            ELSE
                BEDLOAD_FLUX=0.0
                TRANSPORT_STAGE=0.0
            ENDIF
        ELSE
            BEDLOAD_FLUX=0.0
            TRANSPORT_STAGE=0.0
        ENDIF
        RETURN
    END ! SUBROUTINE SEDIMENT_TRANSPORT_FLUX
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE EQUILIBRIUM_SEDIMENT_GRADIENT(BEDLOAD_FLUX,GRADIENT,DISCHARGE_VALUE)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        !        *********************************************************************
        !         This subroutine determines the potential alluvial channel gradient
        !           as a function of sediment load and drainage area
        !        *********************************************************************
        REAL(4), INTENT(IN):: BEDLOAD_FLUX,DISCHARGE_VALUE
        REAL(4), INTENT(OUT) :: GRADIENT
        REAL(4) :: TERM1,DISCHARGE_TERM
        DISCHARGE_TERM=MAX(DISCHARGE_VALUE,0.1E00)*STEADY_DISCHARGE_FACTOR
        IF (BEDLOAD_FLUX > 0.0) THEN
            TERM1=(BEDLOAD_FLUX/(  &
            SEDIMENT_CONSTANT*(DISCHARGE_TERM)**SEDIMENT_1_EXPONENT))**  &
            (1.0/SEDIMENT_TRANSPORT_EXPONENT)
        ELSE
            TERM1=0.0
        ENDIF
        GRADIENT=((TERM1+TRANSPORT_CRITICAL_DIM_SHEAR)/(SEDIMENT_DISCHARGE_FACTOR*DISCHARGE_TERM  &
        **SEDIMENT_2_EXPONENT))** &
        (1.0/SEDIMENT_GRADIENT_EXPONENT)
        RETURN
    END ! SUBROUTINE EQUILIBRIUM_SEDIMENT_GRADIENT
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE SEDIMENT_FLUX_DIVERGENCE()
        USE ERODE_GLOBALS, SED_FLUX_DIVERGENCE=>CFW, CHANNEL_STATE_CHANGE=>IDO
        IMPLICIT NONE
        INTEGER :: I,J,K,IX,JY,NNN,LOCAL_DIRECTION,I_NEW,J_NEW
        REAL(4) :: XQQQ(9),XQSUM,XQMAX,A1TERM,A2TERM,X_DISCHARGE,GRADIENT,BEDLOAD_FLUX, &
        XQUOTE,TRANSPORT_STAGE
        !        *********************************************************************
        !         This calculate alluvial channel sediment diffusion based on local
        !           finite-differencing - calculated sediment divergence is stored in matrix
        !           sed_flux_divergence(i,j)
        !  MODIFIES: SEDIMENT_FLUX, SED_FLUX_DIVERGENCE, 
        !  CALLS: SEDIMENT_TRANSPORT_FLUX
        !        *********************************************************************
        L75: DO  J = 1, MYY
            M75: DO  I = 1, MX
                NNN=BASIN_NUMBER(I,J)
                X_DISCHARGE=MAXIMUM_DISCHARGE(I,J)
                A1TERM=(X_DISCHARGE)**SEDIMENT_1_EXPONENT
                A2TERM=(X_DISCHARGE)**SEDIMENT_2_EXPONENT
                !        *********************************************************************
                !          Skip most everything if we are not doing finite difference sed diffusion
                !        *********************************************************************
                IF (DO_SEDIMENT_DIFFUSION) THEN
                    !        *********************************************************************
                    !          Look at sediment-covered areas first
                    !        *********************************************************************
                    IF (IS_SEDIMENT_COVERED(I,J)) THEN
                        !        *********************************************************************
                        !          Skip lake depressions - No sediment export
                        !        *********************************************************************
                        IF (SUBMERGED(I,J)) CYCLE M75
                        !        *********************************************************************
                        !          Determine sediment transport rate to all surrounding locations with
                        !               gradients lower than local point - need to keep track of maximum and
                        !               total transport rate
                        !        *********************************************************************
                        XQSUM = 0.0
                        XQMAX = 0.0
                        L229: DO K=2,9
                            XQQQ(K)=0.0
                        ENDDO L229
                        L228: DO K=2,9
                            IX=I+DOWNSTREAM(K,1)
                            JY=J+DOWNSTREAM(K,2)
                            IF (IS_X_PERIODIC) THEN
                                IF (IX < 1) IX=MX
                                IF (IX > MX) IX=1
                            ELSE
                                IF (IX < 1) CYCLE L228
                                IF (IX > MX) CYCLE L228
                            ENDIF
                            IF (DO_FLOW_BOUNDARIES) THEN
                                IF (FLOW_DIRECTION(I,J) == 1) THEN
                                    CYCLE L228
                                ENDIF
                            ELSE
                                IF (IS_Y_PERIODIC) THEN
                                    IF (JY < 1) JY=MY
                                    IF (JY > MY) JY=1
                                ELSE
                                    IF (JY < 1) CYCLE L228
                                ENDIF
                            ENDIF
                            GRADIENT=ELEVATION(I,J)-ELEVATION(IX,JY)/CELL_SIZE
                            IF (GRADIENT > 0.0) THEN
                                IF (K > 5) GRADIENT=GRADIENT/1.4142
                                CALL SEDIMENT_TRANSPORT_FLUX(BEDLOAD_FLUX,TRANSPORT_STAGE, &
                                GRADIENT,A1TERM,A2TERM)
                                IF (BEDLOAD_FLUX > XQMAX) XQMAX=BEDLOAD_FLUX
                                XQSUM=XQSUM+BEDLOAD_FLUX
                                XQQQ(K)=BEDLOAD_FLUX
                            ENDIF
                        ENDDO L228
                        !        *********************************************************************
                        !         The sediment transport rate is the maximum of the rates in the various
                        !               directions - The following scales transport rates to add up to the
                        !               maximum rate
                        !        *********************************************************************
                        SEDIMENT_FLUX(I,J)=XQMAX
                        IF (XQMAX > 0.0) THEN
                            XQUOTE=XQMAX/XQSUM
                            L55: DO  K=2,9
                                IX=I+DOWNSTREAM(K,1)
                                JY=J+DOWNSTREAM(K,2)
                                IF (IS_X_PERIODIC) THEN
                                    IF (IX < 1) IX=MX
                                    IF (IX > MX) IX=1
                                ELSE
                                    IF (IX < 1) CYCLE L55
                                    IF (IX > MX) CYCLE L55
                                ENDIF
                                IF (DO_FLOW_BOUNDARIES) THEN
                                    IF (FLOW_DIRECTION(IX,JY) == 1) THEN
                                        CYCLE L55
                                    ENDIF
                                ELSE
                                    IF (IS_Y_PERIODIC) THEN
                                        IF (JY < 1) JY=MY
                                        IF (JY > MY) JY=1
                                    ELSE
                                        IF (JY < 1) CYCLE L55
                                        IF ((JY == (MYY)).AND.  &
                                        (NO_FLUX_LOWER_BOUNDARY)) CYCLE L55
                                    ENDIF
                                ENDIF
                                SED_FLUX_DIVERGENCE(I,J)=SED_FLUX_DIVERGENCE(I,J)-XQQQ(K)*XQUOTE/CELL_SIZE
                                SED_FLUX_DIVERGENCE(IX,JY)=SED_FLUX_DIVERGENCE(IX,JY)+XQQQ(K)*XQUOTE/CELL_SIZE
                            ENDDO L55
                        ENDIF
                    ELSE
                        !        *********************************************************************
                        !         Non-sedimentary locations - Need to check for contributions to sedimentary
                        !               locations
                        !        *********************************************************************
                        CALL SEDIMENT_TRANSPORT_FLUX(SEDIMENT_FLUX(I,J),TRANSPORT_STAGE, &
                        D8_GRADIENT(I,J),A1TERM,A2TERM)
                        LOCAL_DIRECTION=FLOW_DIRECTION(I,J)
                        !        ********************************************************************
                        !         Skip lakes and borders and depressions
                        !         Convert to is_sediment_covered if it is a depression
                        !        ********************************************************************
                        IF (LOCAL_DIRECTION < 0) CYCLE M75
                        I_NEW=I+DOWNSTREAM(LOCAL_DIRECTION,1)
                        J_NEW=J+DOWNSTREAM(LOCAL_DIRECTION,2)
                        IF (IS_X_PERIODIC) THEN
                            IF (I_NEW > MX) I_NEW=1
                            IF (I_NEW < 1) I_NEW=MX
                        ELSE
                            IF (I_NEW > MX) CYCLE M75
                            IF (I_NEW < 1) CYCLE M75
                        ENDIF
                        IF (DO_FLOW_BOUNDARIES) THEN
                            IF (FLOW_DIRECTION(I_NEW,J_NEW) == 1) THEN
                                CYCLE M75
                            ENDIF
                        ELSE
                            IF (IS_Y_PERIODIC) THEN
                                IF (J_NEW > MY) J_NEW=1
                                IF (J_NEW < 1) J_NEW=MY
                            ELSE
                                IF (J_NEW == MY) CYCLE M75
                            ENDIF
                        ENDIF
                        IF (SUBMERGED(I,J)) CYCLE M75
                        !        ********************************************************************
                        !         Otherwise look to see if it feeds a sedimentary area
                        !         If it does, add in sediment contribution
                        !        ********************************************************************
                        IF (IS_SEDIMENT_COVERED(I_NEW,J_NEW)) THEN
                            SED_FLUX_DIVERGENCE(I_NEW,J_NEW)=SED_FLUX_DIVERGENCE(I_NEW,J_NEW)+SEDIMENT_YIELD(I,J)/CELL_SIZE
                            !        ********************************************************************
                            !         Also deposit sediment if the neighboring point is submerged
                            !         and if not is_sediment_covered convert to is_sediment_sedcovered
                            !        ********************************************************************
                        ELSE
                            IF (SUBMERGED(I_NEW,J_NEW)  &
                            .OR.(FLOW_DIRECTION(I_NEW,J_NEW) < 0)) THEN
                                SED_FLUX_DIVERGENCE(I_NEW,J_NEW)=SED_FLUX_DIVERGENCE(I_NEW,J_NEW)+SEDIMENT_YIELD(I,J)/CELL_SIZE
                                CHANNEL_STATE_CHANGE(I_NEW,J_NEW)=1
                                ICASE(1)=ICASE(1)+1
                                SEDIMENT_FLUX(I_NEW,J_NEW)=0.0
                            ENDIF
                        ENDIF
                    ENDIF
                ELSE
                    !        *********************************************************************
                    !         If we are not doing finite-difference alluvial sediment divergence, all this
                    !               routine does is calculate potential sediment transport rate
                    !        *********************************************************************
                    CALL SEDIMENT_TRANSPORT_FLUX(SEDIMENT_FLUX(I,J),TRANSPORT_STAGE,D8_GRADIENT(I,J),A1TERM,A2TERM)
                ENDIF
            ENDDO M75
        ENDDO L75
        RETURN
    END !SUBROUTINE SEDIMENT_FLUX_DIVERGENCE 
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE ROUTE_SEDIMENT(I,J,DELT,BACKGRADE,SEDSOURCE,S1,S2)
        USE ERODE_GLOBALS, CHANNEL_STATE_CHANGE=>IDO, SED_SURFACE_ELEVATION=>CFN
        USE SEDROUTE_GLOBALS
        USE SEDDEBUG_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN):: I,J
        REAL(4), INTENT(IN) :: DELT
        LOGICAL, INTENT(OUT):: BACKGRADE
        LOGICAL, INTENT(INOUT) :: SEDSOURCE
        REAL(4), INTENT(OUT) :: S1,S2
        REAL(4):: CRITELCHANGE
        PARAMETER (CRITELCHANGE=150.0)
        REAL(4) ::  EGRADSUM,GDIFFSUM,  &
        VOLSTART,VOLEND,SEDREMAIN, SEDTOUSE,  &
        REFWATERLEVEL,XWATERLEVEL,DISCHARGEMAX,SEDQMAX, &
        EDIFF ,INITIAL_VOLUME &
        ,elev_max, elev_ref
        INTEGER :: L,LOCAL_DIRECTION, KI,KJ
        INTEGER :: II,JJ,ITYYP, kkkk,kmaxlength
        LOGICAL :: IDOIT,KDOIT,CLOSELOOP
        LOGICAL :: FIXED,MUSTDOIT,IS_OK,redo_elevations
        LOGICAL :: PRINTSTUFF,REVISIT
        INTEGER :: KLOWER
        INTEGER :: I999,I130,I119,I567,I120,I220,IXXX,JXX,IQX
        INTEGER :: I1220,I903
        REAL(4) :: TEMPAREA,ELX,GRADX,EGRADX,SEDX
        REAL(4) :: DELEL,DELEX,TEMPP,BEDLOAD_FLUX,ABSMAXCHG,LOCCHANGE
        !!!!!!!!write(outhist,4833)
4833    format('sediment called')        
        !        *********************************************************************
        !         This routine routes sediment downstream assuming a source from a
        !           non-alluvial channel and that sediment is deposited at equilibrium
        !           alluvial channel gradient for the given sediment yield and drainage area
        !           (discharge surrogate) - There are special cases if we get to a depression
        !           or fixed-elevation exit.
        !           closeloop tells us if we get to a depression with no exit
        !           depression tells us if there is a depression along the flowpath
        !           backgrade tells us if we are bacgrading from a fixed downstream point
        !           (true) or foregrading from an upstream sediment source (false)
        !           fixed tells us if, when we are backgrading, whether it is from an internal
        !           point (false) or from a matrix edge (true)
        !  MODIFIES:  (ALL VARIABLES IN SEDROUTE_GLOBALS AND SEDDEBUG_GLOBALS), CFN, IDO, MAXIMUM_SEDIMENT_YIELD
        !             MAXIMUM_DISCHARGE
        !  CALLS: EQUILIBRIUM_SEDIMENT_VOLUME, FIND_DOWNSTREAM_LOCATION, CHECK_IF_CHANGE_FLOW_DIRECTION
        !        *********************************************************************
        REWIND(OUTSEDDEBUG)
        PRINTSTUFF=.FALSE.
        DO_BACKTRACK=.FALSE.
        REVISIT=.FALSE.
        redo_elevations=.true.
        kmaxlength=5
        II=I
        JJ=J
        I999=0
        I130=0
        I119=0
        I567=0
        I120=0
        I220=0
        I1220=0
        I903=0
        DONE_ONCE=.FALSE.
        !        *********************************************************************
        !         The total sediment volume to be routed is equal to the sediment yield times
        !           the time increment
        !        *********************************************************************
        IF (SEDIMENT_YIELD(II,JJ) <= 0.0) THEN
            S1=0.0
            S2=0.0
            BACKGRADE=.FALSE.
            SEDSOURCE=.FALSE.
            RETURN
        ENDIF
        SEDIMENT_VOLUME=DELT*SEDIMENT_YIELD(II,JJ)*SEDBIAS
        IF (COMPLETE_RUNOFF.AND.MODEL_PELAGIC_DEPOSITION) THEN
            SEDIMENT_VOLUME=SEDIMENT_VOLUME*(1.0-WASHLOAD_FRACTION)
            PELAGIC_SEDIMENT_VOLUME(BASIN_NUMBER(II,JJ))=PELAGIC_SEDIMENT_VOLUME(BASIN_NUMBER(II,JJ)) &
            +  SEDIMENT_VOLUME*WASHLOAD_FRACTION
            8046  FORMAT('V=',G12.5)
        ENDIF
        INITIAL_VOLUME=SEDIMENT_VOLUME
        IF (SEDIMENT_VOLUME <= 0.0) THEN
            S1=0.0
            S2=0.0
            BACKGRADE=.FALSE.
            SEDSOURCE=.FALSE.
            RETURN
        ENDIF
        SEDIMENT_CARRYOVER=0.0
        DEFAULT_DRAINAGE_AREA=0
        DEFAULT_SEDIMENT_YIELD=0.0
        VOLSTART=0.0
        VOLEND=0.0
        SEDREMAIN=0.0
        S1=SEDIMENT_YIELD(II,JJ)
        S2=0
        KDOIT=.FALSE.
        IDOIT=.FALSE.
        DISCHARGEMAX=0.0
        SEDQMAX=0.0
        MUSTDOIT=.FALSE.
        L999: DO
            CLOSELOOP=.FALSE.
            MAXCHG=0.0
            ABSMAXCHG=0.0
            MUST_SEARCH=.FALSE.
            I999=I999+1
            IF (I999 > RMMX) THEN
                RETURN
            ENDIF
            IS_DEPRESSION=.FALSE.
            FIXED=.TRUE.
            !        *********************************************************************
            !         is_bedrock_channel is false if the local site is alluvial channel, otherwise it is true
            !           (fixed) if it is a bedrock channel location
            !        *********************************************************************
            IS_BEDROCK_CHANNEL=.FALSE.
            !        *********************************************************************
            !         We may already have gone through this location with routing from another
            !           source - We will know this has been done if sed_surface_elevation(i,j) is less than the
            !           dummy value of 1.0e+24.  sed_surface_elevation(i,j) holds the new elevation of the alluvial
            !           channel after routing -- upstream_elevation is the elevation of the head of channel
            !           section to be routed through - It is generally a non-alluvial location
            !        *********************************************************************
            IF (SED_SURFACE_ELEVATION(II,JJ) < 1.0E24) THEN
                UPSTREAM_ELEVATION=SED_SURFACE_ELEVATION(II,JJ)
            ELSE
                UPSTREAM_ELEVATION=ELEVATION(II,JJ)
            ENDIF
            IF (COMPLETE_RUNOFF) THEN
                REFWATERLEVEL=AMAX1(OCEAN_ELEVATION,LAKE_SURFACE_ELEVATION(BASIN_NUMBER(II,JJ)))
            ELSE
                REFWATERLEVEL=OCEAN_ELEVATION
            ENDIF
            K=1
            !        *********************************************************************
            !         start_step_distance and alluvial_start_gradient are the distance and equilibrium alluvial gradient from
            !           location 0 (non-alluvial head) and the first alluvial location (index 1)
            !           actual_start_gradient is the present gradient from location 0 to location 1
            !        *********************************************************************
            START_STEP_DISTANCE=CELL_SIZE
            !        *********************************************************************
            !          Standard code for moving to next downstream location
            !        *********************************************************************
            I_OLD=II
            J_OLD=JJ
            CALL FIND_DOWNSTREAM_LOCATION(I_OLD,J_OLD,I_NEW,J_NEW,LOCAL_DIRECTION,IS_OK)
            IF (.NOT.IS_OK) THEN
               WRITE(OUTHIST,8376) I,J
               RETURN
            ENDIF
            8376 FORMAT(' RETURNING, CANNOT ROUTE, I=',I5,' J=',I5)
            I_LOCATION(K)=I_NEW
            J_LOCATION(K)=J_NEW
            IF (COMPLETE_RUNOFF) THEN
                WATER_LEVEL(K)=AMAX1(OCEAN_ELEVATION,LAKE_SURFACE_ELEVATION(BASIN_NUMBER(I_NEW,J_NEW)))
            ELSE
                WATER_LEVEL(K)=OCEAN_ELEVATION
            ENDIF
            !        *********************************************************************
            !          sed_surface_elevation(ii,jj) is the new elevation - When we are to modify it, we need to
            !           check whether it has already been modified by a previous pass downstream
            !           (i.e., if sed_surface_elevation(ii,jj)<1.e24 - If it has been modified, then the present
            !           elevation is sed_surface_elevation(ii,jj), otherwise it is elevation(ii,jj))
            !           starting_elevation(k) are the initial elevations along the flowpath downstream
            !           k is the index of position along the flowpath
            !           i_location(k) and j_location(k) are the matrix locations of the downstream locations
            !        *********************************************************************
            IF (SED_SURFACE_ELEVATION(I_NEW,J_NEW) < 1.0E24) THEN
                ACTUAL_START_GRADIENT=(UPSTREAM_ELEVATION-SED_SURFACE_ELEVATION(I_NEW,J_NEW))/CELL_SIZE
                STARTING_ELEVATION(K)=SED_SURFACE_ELEVATION(I_NEW,J_NEW)
            ELSE
                ACTUAL_START_GRADIENT=(UPSTREAM_ELEVATION-ELEVATION(I_NEW,J_NEW))/CELL_SIZE
                STARTING_ELEVATION(K)=ELEVATION(I_NEW,J_NEW)
            ENDIF
            !        *********************************************************************
            !          We have to take into account whether the next location is along a matrix
            !           dirction or matrix diagonal when we compute the gradient
            !        *********************************************************************
            IF (LOCAL_DIRECTION > 5) THEN
                ACTUAL_START_GRADIENT=ACTUAL_START_GRADIENT/1.4142
                START_STEP_DISTANCE=1.4142*CELL_SIZE
            ENDIF
            DELEL = WATER_LEVEL(K) - STARTING_ELEVATION(K)
            IF (DELEL > 0.0) THEN
                DELEX = DELEL / DELTA_FORESET_GRADIENT
                IF (DELEX < START_STEP_DISTANCE) THEN
                    SEDQMAX=AMAX1(SEDQMAX,MAXIMUM_SEDIMENT_YIELD(II,JJ))
                    DISCHARGEMAX=AMAX1(DISCHARGEMAX,MAXIMUM_DISCHARGE(II,JJ))
                    IF (SEDQMAX > MAXIMUM_SEDIMENT_YIELD(II,JJ)) MAXIMUM_SEDIMENT_YIELD(II,JJ)=SEDQMAX
                    IF (DISCHARGEMAX > MAXIMUM_DISCHARGE(II,JJ)) MAXIMUM_DISCHARGE(II,JJ)=DISCHARGEMAX
                    BEDLOAD_FLUX=SEDQMAX
                    CALL EQUILIBRIUM_SEDIMENT_GRADIENT(BEDLOAD_FLUX,ALLUVIAL_START_GRADIENT,DISCHARGEMAX)
                    TEMPP=ALLUVIAL_START_GRADIENT
                    ALLUVIAL_GRADIENT(K)=ALLUVIAL_START_GRADIENT
                    ALLUVIAL_START_GRADIENT = ALLUVIAL_START_GRADIENT + (DELTA_FORESET_GRADIENT   &
                    - ALLUVIAL_START_GRADIENT)*DELEL/(DELTA_FORESET_GRADIENT*START_STEP_DISTANCE)
                    IF (ALLUVIAL_START_GRADIENT < 0.0) THEN
                        WRITE(OUTHIST,12158) UPSTREAM_ELEVATION,STARTING_ELEVATION(K)  &
                        ,REFWATERLEVEL,DELEX,START_STEP_DISTANCE &
                        ,DELEL,DELTA_FORESET_GRADIENT,TEMPP,ALLUVIAL_START_GRADIENT
                    ENDIF
                ELSE
                    ALLUVIAL_START_GRADIENT = DELTA_FORESET_GRADIENT
                ENDIF
            ELSE
                SEDQMAX=AMAX1(SEDQMAX,MAXIMUM_SEDIMENT_YIELD(II,JJ))
                DISCHARGEMAX=AMAX1(DISCHARGEMAX,MAXIMUM_DISCHARGE(II,JJ))
                IF (SEDQMAX > MAXIMUM_SEDIMENT_YIELD(II,JJ)) MAXIMUM_SEDIMENT_YIELD(II,JJ)=SEDQMAX
                IF (DISCHARGEMAX > MAXIMUM_DISCHARGE(II,JJ)) MAXIMUM_DISCHARGE(II,JJ)=DISCHARGEMAX
                BEDLOAD_FLUX=SEDQMAX
                CALL EQUILIBRIUM_SEDIMENT_GRADIENT(BEDLOAD_FLUX,ALLUVIAL_START_GRADIENT,DISCHARGEMAX)
            ENDIF
            IF (.NOT.IS_SEDIMENT_COVERED(I_NEW,J_NEW)) THEN
                IF (CHANNEL_STATE_CHANGE(I_NEW,J_NEW) == 0) IS_BEDROCK_CHANNEL(K)=.TRUE.
            ENDIF
            IF (DO_FLOW_BOUNDARIES) THEN
                IF (FLOW_DIRECTION(I_NEW,J_NEW) == 1) THEN
                    IS_BEDROCK_CHANNEL(K)=.TRUE.
                ENDIF
            ELSE
                IF (IS_Y_PERIODIC) THEN
                ELSE
                    IF ((J_NEW == MY).AND.(.NOT.NO_FLUX_LOWER_BOUNDARY)) IS_BEDROCK_CHANNEL(K)=.TRUE.
                ENDIF
            ENDIF
            !        *********************************************************************
            !         This is the main loop - Basically we continue downstream until we have
            !           deposited all of the input sediment (sediment_volume) or until we reach an
            !           undrained is_depression, or a bedrock channel that cannot be buried by the
            !           amount of alluvium to be deposited - if we reach a fixed downstream
            !           location before depositing all alluvium, we backgrade from that fixed
            !           location (label 120)
            !        *********************************************************************
            L130: DO
                I130=I130+1
                IF (I130 > RMMX) THEN
                    RETURN
                ENDIF
                IF (K > RMMX) THEN
                    RETURN
                ENDIF
                !        *********************************************************************
                !         The following few lines are the basic procedure for finding the next
                !           downstream location and keeping track of whether there are depressions or
                !           fixed bedrock locations.
                !           i_location(k) and j_location(k) are the i,j locations of the downstream points
                !        *********************************************************************
                IF (DO_FLOW_BOUNDARIES) THEN
                    IF (FLOW_DIRECTION(I_NEW,J_NEW) == 1) THEN
                        IF (I120 > 2*RMMX) THEN
                        ENDIF
                        GOTO 120
                    ENDIF
                ELSE
                    IF (IS_Y_PERIODIC) THEN
                    ELSE
                        IF (J_NEW == MY) THEN
                            IF (I120 > 2*RMMX) THEN
                            ENDIF
                            IF (.NOT.NO_FLUX_LOWER_BOUNDARY) GOTO 120
                        ENDIF
                    ENDIF
                ENDIF
                I_OLD=I_NEW
                J_OLD=J_NEW
                SEDQMAX=AMAX1(MAXIMUM_SEDIMENT_YIELD(I_OLD,J_OLD),SEDQMAX)
                BEDLOAD_FLUX=SEDQMAX
                DISCHARGEMAX=AMAX1(MAXIMUM_DISCHARGE(I_OLD,J_OLD),DISCHARGEMAX)
                IF (SEDQMAX > MAXIMUM_SEDIMENT_YIELD(I_OLD,J_OLD)) MAXIMUM_SEDIMENT_YIELD(I_OLD,J_OLD)=SEDQMAX
                IF (DISCHARGEMAX > MAXIMUM_DISCHARGE(I_OLD,J_OLD)) MAXIMUM_DISCHARGE(I_OLD,J_OLD)=DISCHARGEMAX
                CALL FIND_DOWNSTREAM_LOCATION(I_OLD,J_OLD,I_NEW,J_NEW,LOCAL_DIRECTION,IS_OK)
                IF (.NOT.IS_OK) THEN
                   WRITE(OUTHIST,8376) I,J
                   RETURN
                ENDIF
                
                I_LOCATION(K+1)=I_NEW
                J_LOCATION(K+1)=J_NEW
                !        *********************************************************************
                !         starting_elevation(k) are the initial elevations of the downstream locations
                !        *********************************************************************
                IF (COMPLETE_RUNOFF) THEN
                    WATER_LEVEL(K+1)=AMAX1(OCEAN_ELEVATION,LAKE_SURFACE_ELEVATION(BASIN_NUMBER(I_NEW,J_NEW)))
                ELSE
                    WATER_LEVEL(K+1)=OCEAN_ELEVATION
                ENDIF
                IF (SED_SURFACE_ELEVATION(I_NEW,J_NEW) < 1.0E24) THEN
                    EL2=SED_SURFACE_ELEVATION(I_NEW,J_NEW)
                    STARTING_ELEVATION(K+1)=SED_SURFACE_ELEVATION(I_NEW,J_NEW)
                ELSE
                    EL2=ELEVATION(I_NEW,J_NEW)
                    STARTING_ELEVATION(K+1)=ELEVATION(I_NEW,J_NEW)
                ENDIF
                IF (SED_SURFACE_ELEVATION(I_OLD,J_OLD) < 1.0E24) THEN
                    EL1=SED_SURFACE_ELEVATION(I_OLD,J_OLD)
                ELSE
                    EL1=ELEVATION(I_OLD,J_OLD)
                ENDIF
                !        *********************************************************************
                !          end_elevation is the elevation of the downstream non-alluviating location.
                !           actual_gradient(k) are the actual (initial) gradients of the downstream locations.
                !           step_distance(k) are the point-to-point distances between nodes (1 or sqrt(2)).
                !           alluvial_gradient(k) are the equilibrium alluvial gradients of the downstream points.
                !        *********************************************************************
                END_ELEVATION=EL2
                ACTUAL_GRADIENT(K)=(EL1-EL2)/CELL_SIZE
                IF (.NOT.IS_SEDIMENT_COVERED(I_NEW,J_NEW)) THEN
                    IF (CHANNEL_STATE_CHANGE(I_NEW,J_NEW) == 0) IS_BEDROCK_CHANNEL(K+1)=.TRUE.
                ENDIF
                IF (DO_FLOW_BOUNDARIES) THEN
                    IF (FLOW_DIRECTION(I_NEW,J_NEW) == 1) THEN
                        IS_BEDROCK_CHANNEL(K)=.TRUE.
                    ENDIF
                ELSE
                    IF (IS_Y_PERIODIC) THEN
                    ELSE
                        IF ((J_NEW == MY).AND.(.NOT.NO_FLUX_LOWER_BOUNDARY)) IS_BEDROCK_CHANNEL(K)=.TRUE.
                    ENDIF
                ENDIF
                STEP_DISTANCE(K)=CELL_SIZE
                IF (LOCAL_DIRECTION > 5) THEN
                    ACTUAL_GRADIENT(K)=ACTUAL_GRADIENT(K)/1.4142
                    STEP_DISTANCE(K)=1.4142*CELL_SIZE
                ENDIF
                IF (IS_DEPRESSION) THEN
                    DELEL = WATER_LEVEL(K+1) - EL2
                    IF (DELEL > 0.0) THEN
                        DELEX = DELEL / DELTA_FORESET_GRADIENT
                        IF (DELEX < STEP_DISTANCE(K)) THEN
                            CALL EQUILIBRIUM_SEDIMENT_GRADIENT(BEDLOAD_FLUX,ALLUVIAL_GRADIENT(K),DISCHARGEMAX)
                            TEMPP=ALLUVIAL_GRADIENT(K)
                            ALLUVIAL_GRADIENT(K) = ALLUVIAL_GRADIENT(K) + &
                            (DELTA_FORESET_GRADIENT - ALLUVIAL_GRADIENT(K))*DELEL/(DELTA_FORESET_GRADIENT*STEP_DISTANCE(K))
                            IF (ALLUVIAL_GRADIENT(K) < 0.0) THEN
                                WRITE(OUTHIST,12158) EL1,EL2,WATER_LEVEL(K),DELEX, &
                                STEP_DISTANCE(K),DELEL,DELTA_FORESET_GRADIENT,TEMPP,ALLUVIAL_GRADIENT(K)
                                12158  FORMAT(' NEG EGRAD, EL1=',G12.5,' EL2=',G12.5,' WATER_LEVEL=', &
                                G12.5,' DELEX=',G12.5,' STEP_DISTANCE=',G12.5,' DELEL=',G12.5,/, &
                                ' DELTA_FORESET_GRADIENT=',G12.5,' OLD EGRAD=',G12.5,' NEW EGRAD=',G12.5)
                            ENDIF
                        ELSE
                            ALLUVIAL_GRADIENT(K) = DELTA_FORESET_GRADIENT
                        ENDIF
                    ELSE
                        CALL EQUILIBRIUM_SEDIMENT_GRADIENT(BEDLOAD_FLUX,ALLUVIAL_GRADIENT(K),DISCHARGEMAX)
                    ENDIF
                ELSE
                    SEDTOUSE=SEDIMENT_YIELD(I_OLD,J_OLD)+SEDIMENT_CARRYOVER
                    TEMPAREA=MAXIMUM_DISCHARGE(I_OLD,J_OLD)
                    SEDQMAX=AMAX1(SEDQMAX,MAXIMUM_SEDIMENT_YIELD(I_OLD,J_OLD))
                    BEDLOAD_FLUX=SEDQMAX
                    DISCHARGEMAX=AMAX1(MAXIMUM_DISCHARGE(I_OLD,J_OLD),DISCHARGEMAX)
                    IF (SEDQMAX > MAXIMUM_SEDIMENT_YIELD(II,JJ)) MAXIMUM_SEDIMENT_YIELD(II,JJ)=SEDQMAX
                    IF (DISCHARGEMAX > MAXIMUM_DISCHARGE(II,JJ)) MAXIMUM_DISCHARGE(II,JJ)=DISCHARGEMAX
                    IF (TEMPAREA < DEFAULT_DRAINAGE_AREA) TEMPAREA=DEFAULT_DRAINAGE_AREA
                    DELEL = WATER_LEVEL(K+1) - EL2
                    IF (DELEL > 0.0) THEN
                        DELEX = DELEL / DELTA_FORESET_GRADIENT
                        IF (DELEX < STEP_DISTANCE(K)) THEN
                            CALL EQUILIBRIUM_SEDIMENT_GRADIENT(BEDLOAD_FLUX,ALLUVIAL_GRADIENT(K),DISCHARGEMAX)
                            TEMPP=ALLUVIAL_GRADIENT(K)
                            ALLUVIAL_GRADIENT(K) = ALLUVIAL_GRADIENT(K) +  &
                            (DELTA_FORESET_GRADIENT - ALLUVIAL_GRADIENT(K))*DELEL/(DELTA_FORESET_GRADIENT*STEP_DISTANCE(K))
                            IF (ALLUVIAL_GRADIENT(K) < 0.0) THEN
                                WRITE(OUTHIST,12158) EL1,EL2,WATER_LEVEL(K),DELEX, &
                                STEP_DISTANCE(K),DELEL,DELTA_FORESET_GRADIENT,TEMPP,ALLUVIAL_GRADIENT(K)
                            ENDIF
                        ELSE
                            ALLUVIAL_GRADIENT(K) = DELTA_FORESET_GRADIENT
                        ENDIF
                    ELSE
                        CALL EQUILIBRIUM_SEDIMENT_GRADIENT(BEDLOAD_FLUX,ALLUVIAL_GRADIENT(K),DISCHARGEMAX)
                    ENDIF
                ENDIF
                IF (SEDSOURCE) THEN
                    GDIFFSUM=0.0
                    IF (K > 1) THEN
                        DO  L=1,K-1
                            GDIFFSUM=GDIFFSUM+(L+1)*(ACTUAL_GRADIENT(L)-ALLUVIAL_GRADIENT(L))*STEP_DISTANCE(L)
                        ENDDO
                    ENDIF
                    PROVISIONAL_END_GRADIENT=ACTUAL_GRADIENT(K)+(SEDIMENT_VOLUME/CELL_AREA+GDIFFSUM+ &
                    (ACTUAL_START_GRADIENT-ALLUVIAL_START_GRADIENT)*START_STEP_DISTANCE  &
                    )/((K+1)*STEP_DISTANCE(K))
                    IF (.NOT.DO_BACKTRACK) THEN
                        IF (PROVISIONAL_END_GRADIENT > ALLUVIAL_GRADIENT(K)) THEN
                            IF (IS_Y_PERIODIC) THEN
                                K=K+1
                                CYCLE L130
                            ELSE
                                IF ((J_NEW == MY).AND.(NO_FLUX_LOWER_BOUNDARY)) THEN
                                ELSE
                                    K=K+1
                                    CYCLE L130
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDIF
                    NEW_ELEVATION(K)=STARTING_ELEVATION(K+1)+PROVISIONAL_END_GRADIENT*STEP_DISTANCE(K)
                    IF (.NOT.DO_BACKTRACK) THEN
                        IF (IS_BEDROCK_CHANNEL(K)) THEN
                            IF (NEW_ELEVATION(K) <= ELEVATION(I_LOCATION(K),J_LOCATION(K))) THEN
                                KKK= K
                                GOTO 119
                            ENDIF
                        ELSE
                            IF (NEW_ELEVATION(K) < SEDIMENT_BASE(I_LOCATION(K),J_LOCATION(K))) THEN
                                KKK = K
                                CHANNEL_STATE_CHANGE(I_LOCATION(K),J_LOCATION(K))=2
                                ICASE(2)=ICASE(2)-1
                                GOTO 119
                            ENDIF
                        ENDIF
                    ENDIF
                    IF (K > 1) THEN
                        L11160: DO  L=K-1,1,-1
                            KI=I_LOCATION(L)
                            KJ=J_LOCATION(L)
                            NEW_ELEVATION(L)=NEW_ELEVATION(L+1)+ALLUVIAL_GRADIENT(L)*STEP_DISTANCE(L)
                            IF (.NOT.DO_BACKTRACK) THEN
                                IF (IS_BEDROCK_CHANNEL(L)) THEN
                                    IF (NEW_ELEVATION(L) <= ELEVATION(KI,KJ)) THEN
                                        KKK=L
                                        IF (NEW_ELEVATION(L) < ELOWEST) MUSTDOIT=.TRUE.
                                        GOTO 119
                                    ENDIF
                                ELSE
                                    IF (NEW_ELEVATION(L) < SEDIMENT_BASE(KI,KJ)) THEN
                                        KKK=L
                                        CHANNEL_STATE_CHANGE(KI,KJ)=2
                                        ICASE(3)=ICASE(3)-1
                                        IF (NEW_ELEVATION(L) < ELOWEST) MUSTDOIT=.TRUE.
                                        GOTO 119
                                    ENDIF
                                ENDIF
                            ENDIF
                        ENDDO L11160
                    ENDIF
                    SED_SURFACE_ELEVATION(II,JJ)=NEW_ELEVATION(1)+ALLUVIAL_START_GRADIENT*START_STEP_DISTANCE
                    LOCCHANGE=SED_SURFACE_ELEVATION(II,JJ)-ELEVATION(II,JJ)
                    IF (ABS(LOCCHANGE) > ABSMAXCHG) THEN
                        ABSMAXCHG=ABS(LOCCHANGE)
                        MAXCHG=LOCCHANGE
                        IMAX=II
                        JMAX=JJ
                    ENDIF
                    PROVISIONAL_START_GRADIENT=ALLUVIAL_START_GRADIENT
                ELSE
                    !        *********************************************************************
                    !         egradsum and gdiffsum are summation terms used to calculate the new
                    !               gradients of the furthest upstream and furthest downstream segments
                    !               being considered
                    !        *********************************************************************
                    EGRADSUM=0.0
                    GDIFFSUM=0.0
                    IF (K > 1) THEN
                        DO  L=1,K-1
                            GDIFFSUM=GDIFFSUM+(K-L)*(ALLUVIAL_GRADIENT(L)-ACTUAL_GRADIENT(L))*STEP_DISTANCE(L)
                            EGRADSUM=EGRADSUM+ALLUVIAL_GRADIENT(L)*STEP_DISTANCE(L)
                        ENDDO
                    ENDIF
                    !        *********************************************************************
                    !          provisional_start_gradient is the (provisional) new gradient of the the upstream segment
                    !           (between locations 0 to 1)
                    !          provisional_end_gradient is the (provisional) new gradient of the last downstream segment
                    !           (between locations k and k+1)
                    !          new_elevation(l) are the (provisional) new elevations of the alluvial channel nodes
                    !        *********************************************************************
                    PROVISIONAL_START_GRADIENT=ACTUAL_START_GRADIENT-(SEDIMENT_VOLUME/CELL_AREA+GDIFFSUM)/(K*START_STEP_DISTANCE)
                    PROVISIONAL_END_GRADIENT=(UPSTREAM_ELEVATION-END_ELEVATION- &
                    PROVISIONAL_START_GRADIENT*START_STEP_DISTANCE-EGRADSUM)/STEP_DISTANCE(K)
                    NEW_ELEVATION(1)=UPSTREAM_ELEVATION-PROVISIONAL_START_GRADIENT*START_STEP_DISTANCE
                    IF (.NOT.DO_BACKTRACK) THEN
                        IF (IS_BEDROCK_CHANNEL(1)) THEN
                            IF (NEW_ELEVATION(1) <= ELEVATION(I_LOCATION(1),J_LOCATION(1))) THEN
                                KKK= 1
                                GOTO 119
                            ENDIF
                        ELSE
                            IF (NEW_ELEVATION(1) < SEDIMENT_BASE(I_LOCATION(1),J_LOCATION(1))) THEN
                                KKK = 1
                                CHANNEL_STATE_CHANGE(I_LOCATION(1),J_LOCATION(1))=2
                                ICASE(2)=ICASE(2)-1
                                GOTO 119
                            ENDIF
                        ENDIF
                    ENDIF
                    IF (K > 1) THEN
                        L160: DO  L=2,K
                            KI=I_LOCATION(L)
                            KJ=J_LOCATION(L)
                            NEW_ELEVATION(L)=NEW_ELEVATION(L-1)   &
                            -ALLUVIAL_GRADIENT(L-1)*STEP_DISTANCE(L-1)
                            !        *********************************************************************
                            !         The following code sends us to label 119 if there are fixed, non alluvial
                            !                   channels along the path that cannot be buried with sediment
                            !        *********************************************************************
                            IF (.NOT.DO_BACKTRACK) THEN
                                IF (IS_BEDROCK_CHANNEL(L)) THEN
                                    IF (NEW_ELEVATION(L) <= ELEVATION(KI,KJ)) THEN
                                        KKK=L
                                        IF (NEW_ELEVATION(L) < ELOWEST) MUSTDOIT=.TRUE.
                                        GOTO 119
                                    ENDIF
                                ELSE
                                    IF (NEW_ELEVATION(L) < SEDIMENT_BASE(KI,KJ)) THEN
                                        KKK=L
                                        CHANNEL_STATE_CHANGE(KI,KJ)=2
                                        ICASE(3)=ICASE(3)-1
                                        IF (NEW_ELEVATION(L) < ELOWEST) MUSTDOIT=.TRUE.
                                        GOTO 119
                                    ENDIF
                                ENDIF
                            ENDIF
                        ENDDO L160
                    ENDIF
                    !        *********************************************************************
                    !         If the calculated value of provisional_end_gradient is larger than the equilibrium alluvial
                    !           gradient for that segment, then we have not exhausted the deposition of
                    !           sediment supplied from upstream, so we need to progress downstream one
                    !           more node and redo the calculations
                    !         However, if it is a closed loop, we stop and accept the new elevations
                    !        *********************************************************************
                    IF (.NOT.DO_BACKTRACK) THEN
                        IF (PROVISIONAL_END_GRADIENT > ALLUVIAL_GRADIENT(K)) THEN
                            IF (IS_Y_PERIODIC) THEN
                                K=K+1
                                CYCLE L130
                            ELSE
                                IF ((J_NEW == MY).AND.(NO_FLUX_LOWER_BOUNDARY)) THEN
                                    EXIT L130
                                ELSE
                                    K=K+1
                                    CYCLE L130
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDIF
                ENDIF
                EXIT L130
            ENDDO L130
            !        *********************************************************************
            !         If we are here, we have exhausted the supply of sediment for upstream, so we
            !           set the new alluvial channel elevations, sed_surface_elevation(ii,jj) to the calculated
            !           equilibrium values (new_elevation(l).
            !         Also, if we bury a non-alluvial segment of channel, we need to convert it
            !           to alluvial and to set the base of the sedimentary deposit for future
            !           reference and to determine the new equilibrium gradient
            !        *********************************************************************
            !       ```````````````````````````````````````
            if (do_backtrack) then
                if(redo_elevations.and.(k<=kmaxlength)) then
                IF (SED_SURFACE_ELEVATION(II,JJ) < 1.0E24) THEN
                   UPSTREAM_ELEVATION=SED_SURFACE_ELEVATION(II,JJ)
                ELSE
                   UPSTREAM_ELEVATION=ELEVATION(II,JJ)
                ENDIF
                
                elev_max=upstream_elevation
                if (k>1) then
                    do kkkk=1,k
                        IF (SED_SURFACE_ELEVATION(i_location(kkkk),j_location(kkkk)) < 1.0E24) THEN
                            elev_ref=sed_surface_elevation(i_location(kkkk),j_location(kkkk))
                        ELSE
                            elev_ref=elevation(i_location(kkkk),j_location(kkkk))
                        endif
                        if (elev_ref>elev_max) elev_max=elev_ref
                    enddo
                endif
                if (.NOT.IS_SEDIMENT_COVERED(ii,jj)) then
                      channel_state_change(ii,jj)=1
                      sediment_base(ii,jj)=upstream_elevation
                      sed_surface_elevation(ii,jj)=elev_max+0.01
                endif
                if (k>1) then
                   do kkkk=1,k
                       if (.not.is_sediment_covered(i_location(kkkk),j_location(kkkk))) then
                          channel_state_change(i_location(kkkk),j_location(kkkk))=1
                          sediment_base(i_location(kkkk),j_location(kkkk))=elevation(i_location(kkkk),j_location(kkkk))
                       endif
                       sed_surface_elevation(i_location(kkkk),j_location(kkkk))=sed_surface_elevation(ii,jj)-kkkk*0.01
                   enddo
                endif
                endif
                return
            endif
            SED_SURFACE_ELEVATION(I_LOCATION(1),J_LOCATION(1))=NEW_ELEVATION(1)
            LOCCHANGE=SED_SURFACE_ELEVATION(I_LOCATION(1),J_LOCATION(1))-ELEVATION(I_LOCATION(1),J_LOCATION(1))
            IF (ABS(LOCCHANGE) > ABSMAXCHG) THEN
                ABSMAXCHG=ABS(LOCCHANGE)
                MAXCHG=LOCCHANGE
                IMAX=I_LOCATION(1)
                JMAX=J_LOCATION(1)
            ENDIF
            IF (NEW_ELEVATION(1) < ELOWEST) MUSTDOIT=.TRUE.
            IF (DO_ALLUVIAL_REROUTING) CALL CHECK_IF_CHANGE_FLOW_DIRECTION(I_LOCATION(1),J_LOCATION(1))
            IF (DONE_ONCE(I_LOCATION(1),J_LOCATION(1))) THEN
            ELSE
                DONE_ONCE(I_LOCATION(1),J_LOCATION(1))=.TRUE.
            ENDIF
            KI=I_LOCATION(1)
            KJ=J_LOCATION(1)
            IF (.NOT.IS_SEDIMENT_COVERED(KI,KJ)) THEN
                IF (ELEVATION(KI,KJ) < SED_SURFACE_ELEVATION(KI,KJ)) THEN
                    CHANNEL_STATE_CHANGE(KI,KJ)=1
                    ICASE(4)=ICASE(4)+1
                ENDIF
            ENDIF
            VOLSTART=STARTING_ELEVATION(1)*CELL_AREA
            VOLEND=NEW_ELEVATION(1)*CELL_AREA
            IF (K > 1) THEN
                L445: DO  L=2,K
                    KI=I_LOCATION(L)
                    KJ=J_LOCATION(L)
                    SED_SURFACE_ELEVATION(KI,KJ)=NEW_ELEVATION(L)
                    LOCCHANGE=SED_SURFACE_ELEVATION(KI,KJ)-ELEVATION(KI,KJ)
                    IF (ABS(LOCCHANGE) > ABSMAXCHG) THEN
                        ABSMAXCHG=ABS(LOCCHANGE)
                        MAXCHG=LOCCHANGE
                        IMAX=KI
                        JMAX=KJ
                    ENDIF
                    IF (NEW_ELEVATION(L) < ELOWEST) MUSTDOIT=.TRUE.
                    DEPOSITWORK=DEPOSITWORK+(NEW_ELEVATION(L)-STARTING_ELEVATION(L))* &
                    SQRT(REAL((I-KI)**2+(J-KJ)**2))
                    DEPOSITGRAV=DEPOSITGRAV+(NEW_ELEVATION(L)-STARTING_ELEVATION(L))* &
                    (ELEVATION(I,J)-0.5*(NEW_ELEVATION(L)+STARTING_ELEVATION(L)))
                    IF (DO_ALLUVIAL_REROUTING) CALL CHECK_IF_CHANGE_FLOW_DIRECTION(KI,KJ)
                    IF (DONE_ONCE(KI,KJ)) THEN
                    ELSE
                        DONE_ONCE(KI,KJ)=.TRUE.
                    ENDIF
                    VOLSTART=VOLSTART+STARTING_ELEVATION(L)*CELL_AREA
                    VOLEND=VOLEND+NEW_ELEVATION(L)*CELL_AREA
                    IF (.NOT.IS_SEDIMENT_COVERED(KI,KJ)) THEN
                        IF (ELEVATION(KI,KJ) < SED_SURFACE_ELEVATION(KI,KJ)) THEN
                            CHANNEL_STATE_CHANGE(KI,KJ)=1
                            ICASE(5)=ICASE(5)+1
                        ENDIF
                    ENDIF
                ENDDO L445
            ENDIF
            BACKGRADE=.FALSE.
            TEMPP=VOLEND-VOLSTART
            !        *********************************************************************
            !         In some cases, particularly if there is a downstream depression, the
            !           calculated gradient for the segment between nodes 0 and 1 may be
            !           calculated to be negative. If that is the case, we convert the upstream
            !           point 0 to alluvial and deposit some sediment there
            !        *********************************************************************
            IF (PROVISIONAL_START_GRADIENT < ALLUVIAL_START_GRADIENT) THEN
                I220=I220+1
                IF (I220 > 10*RMMX) THEN
                    RETURN
                ENDIF
                IF (PROVISIONAL_START_GRADIENT < 0.0) THEN
                ENDIF
                CHANNEL_STATE_CHANGE(II,JJ)=1
                IF (HAVE_INFLUENT_RIVERS.AND.IS_INFLUENT_RIVER_LOCATION(II,JJ)) THEN
                    IF (K > 1) SED_SURFACE_ELEVATION(II,JJ)=NEW_ELEVATION(1)+ALLUVIAL_START_GRADIENT*START_STEP_DISTANCE
                    LOCCHANGE=SED_SURFACE_ELEVATION(II,JJ)-ELEVATION(II,JJ)
                    IF (ABS(LOCCHANGE) > ABSMAXCHG) THEN
                        ABSMAXCHG=ABS(LOCCHANGE)
                        MAXCHG=LOCCHANGE
                        IMAX=II
                        JMAX=JJ
                    ENDIF

                ENDIF
                ICASE(6)=ICASE(6)+1
            ENDIF
            !      ******************************************************
            !      ******************************************************
            !        Some debugging
            DO L=1,K-1
                TEMPP=(NEW_ELEVATION(L)-NEW_ELEVATION(L+1))/CELL_SIZE
                IF (TEMPP < 0.0) THEN
                    WRITE(OUTHIST,3882) I_LOCATION(L),J_LOCATION(L),TEMPP,NEW_ELEVATION(L),NEW_ELEVATION(L+1),WATER_LEVEL(L)
                    3882  FORMAT(' NEG SED GRAD A, I=',I5,' J=',I5,' S=',G12.5,/, &
                    ' E=',G12.5,' EN=',G12.5, ' EW=',G12.5)
                ENDIF
            ENDDO
            !      ******************************************************
            !      ******************************************************
            IF (ABSMAXCHG >= CRITELCHANGE) THEN
                ITYYP=1
                CALL  PRINT_SEDIMENT_DIAGNOSTICS(IMAX,JMAX,II,JJ,ITYYP,SEDSOURCE,ABSMAXCHG,DELT)
            ENDIF
            !        *********************************************************************
            !         Print out some data if we are debugging this subroutine
            !        *********************************************************************
            1220  CONTINUE
            I1220=I1220+1
            IF (I1220 > RMMX) THEN
                RETURN
            ENDIF
            S2=S2+(VOLEND-VOLSTART)/DELT
            903   CONTINUE
            I903=I903+1
            IF (I903 > RMMX) THEN
                RETURN
            ENDIF
            IF (DO_FLOW_BOUNDARIES) THEN
                IF (FLOW_DIRECTION(I_NEW,J_NEW) == 1) THEN
                    RETURN
                ENDIF
            ELSE
                IF (IS_Y_PERIODIC) THEN
                ELSE
                    IF (J_NEW == MY) THEN
                        IF (PRINTSTUFF) THEN
                            RETURN
                        ELSE
                            RETURN
                        ENDIF
                    ENDIF
                ENDIF
            ENDIF
            SEDREMAIN = SEDIMENT_VOLUME-VOLEND+VOLSTART
            IF (IS_SEDIMENT_COVERED(I_NEW,J_NEW).OR.(CHANNEL_STATE_CHANGE(I_NEW,J_NEW) == 1)) THEN
                IF (ABS(SEDREMAIN) > (0.01*SEDIMENT_VOLUME)) THEN
                    IF (BACKGRADE) THEN
                        IF (CHANNEL_STATE_CHANGE(I_NEW,J_NEW) == 1) THEN
                            CHANNEL_STATE_CHANGE(I_NEW,J_NEW)=0
                            ICASE(7)=ICASE(7)-1
                        ELSE
                        ENDIF
                    ENDIF
                ENDIF
                IF (BACKGRADE) THEN
                    IF (CHANNEL_STATE_CHANGE(I_NEW,J_NEW) == 1) THEN
                        CHANNEL_STATE_CHANGE(I_NEW,J_NEW)=0
                        ICASE(8)=ICASE(8)-1
                    ELSE
                    ENDIF
                ENDIF
            ENDIF
            IF (SED_SURFACE_ELEVATION(I_NEW,J_NEW) < 1.0E24) THEN
                EL2=SED_SURFACE_ELEVATION(I_NEW,J_NEW)
            ELSE
                EL2=ELEVATION(I_NEW,J_NEW)
            ENDIF
            IF (BACKGRADE.AND.(SEDREMAIN > 0.0)) THEN
                IS_DEPRESSION=.FALSE.
                I567=0
                L567: DO
                    I_OLD=I_NEW
                    I567=I567+1
                    IF (I567 > 10*IMMX) THEN
                        RETURN
                    ENDIF
                    J_OLD=J_NEW
                    CALL FIND_DOWNSTREAM_LOCATION(I_OLD,J_OLD,I_NEW,J_NEW,LOCAL_DIRECTION,IS_OK)
                    IF (.NOT.IS_OK) THEN
                       WRITE(OUTHIST,8376) I,J
                       RETURN
                    ENDIF
                    
                    IF (SED_SURFACE_ELEVATION(I_NEW,J_NEW) < 1.0E24) THEN
                        EL2=SED_SURFACE_ELEVATION(I_NEW,J_NEW)
                    ELSE
                        EL2=ELEVATION(I_NEW,J_NEW)
                    ENDIF
                    !      ******************************************************
                    !      ******************************************************
                    !        Some debugging
                    DO L=1,K-1
                        TEMPP=(NEW_ELEVATION(L)-NEW_ELEVATION(L+1))/CELL_SIZE
                        IF (TEMPP < 0.0) THEN
                            WRITE(*,3883) I_LOCATION(L),J_LOCATION(L),TEMPP,NEW_ELEVATION(L),NEW_ELEVATION(L+1),WATER_LEVEL(L)
                            3883  FORMAT(' NEG SED GRAD B, I=',I5,' J=',I5,' S=',G12.5,/, &
                            ' E=',G12.5,' EN=',G12.5, ' EW=',G12.5)
                        ENDIF
                    ENDDO
                    !      ******************************************************
                    !      ******************************************************
                    IF (DO_FLOW_BOUNDARIES) THEN
                        IF (FLOW_DIRECTION(I_NEW,J_NEW) == 1) THEN
                            RETURN
                        ENDIF
                    ELSE
                        IF (IS_Y_PERIODIC) THEN
                        ELSE
                            IF (J_NEW == MY) THEN
                                IF (PRINTSTUFF) THEN
                                    RETURN
                                ELSE
                                    RETURN
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDIF
                    CALL FIND_DOWNSTREAM_LOCATION(I_NEW,J_NEW,IXXX,JXX,IQX,IS_OK)
                    IF (.NOT.IS_OK) THEN
                        WRITE(OUTHIST,8376) I,J
                        RETURN
                    ENDIF
                    
                    IF (SED_SURFACE_ELEVATION(IXXX,JXX) < 1.0E24) THEN
                        ELX=SED_SURFACE_ELEVATION(IXXX,JXX)
                    ELSE
                        ELX=ELEVATION(IXXX,JXX)
                    ENDIF
                    IF (COMPLETE_RUNOFF) THEN
                        XWATERLEVEL=AMAX1(OCEAN_ELEVATION, &
                        LAKE_SURFACE_ELEVATION(BASIN_NUMBER(IXXX,JXX)))
                    ELSE
                        XWATERLEVEL=OCEAN_ELEVATION
                    ENDIF
                    GRADX=(EL2-ELX)/CELL_SIZE
                    IF (IQX > 5) GRADX=GRADX/1.4142
                    SEDX=SEDREMAIN/(DELT)
                    BEDLOAD_FLUX=AMAX1(SEDQMAX,MAXIMUM_SEDIMENT_YIELD(I_NEW,J_NEW))
                    DISCHARGEMAX=AMAX1(DISCHARGEMAX,MAXIMUM_DISCHARGE(I_NEW,J_NEW))
                    IF (SEDQMAX > MAXIMUM_SEDIMENT_YIELD(I_NEW,J_NEW)) MAXIMUM_SEDIMENT_YIELD(I_NEW,J_NEW)=SEDQMAX
                    IF (DISCHARGEMAX > MAXIMUM_DISCHARGE(I_NEW,J_NEW)) MAXIMUM_DISCHARGE(I_NEW,J_NEW)=DISCHARGEMAX
                    IF (EL2 < XWATERLEVEL) THEN
                        EGRADX = DELTA_FORESET_GRADIENT
                    ELSE
                        CALL EQUILIBRIUM_SEDIMENT_GRADIENT(BEDLOAD_FLUX,EGRADX,DISCHARGEMAX)
                    ENDIF
                    IF (IS_SEDIMENT_COVERED(I_NEW,J_NEW).OR.(GRADX < EGRADX).OR.(FLOW_DIRECTION(I_NEW, &
                    J_NEW) < 0.0)) THEN
                        II=I_OLD
                        JJ=J_OLD
                        SEDIMENT_VOLUME=SEDREMAIN
                        DEFAULT_SEDIMENT_YIELD=0.0
                        IF (I999 > RMMX) THEN
                            IF (IS_SEDIMENT_COVERED(I_NEW,J_NEW)) THEN
                            ELSE
                            ENDIF
                        ENDIF
                        CYCLE L999
                    ENDIF
                    CYCLE L567
                ENDDO L567
            ENDIF
            !      ******************************************************
            !      ******************************************************
            !        Some debugging
            DO L=1,K-1
                TEMPP=(NEW_ELEVATION(L)-NEW_ELEVATION(L+1))/CELL_SIZE
                IF (TEMPP < 0.0) THEN
                    WRITE(*,3884) I_LOCATION(L),J_LOCATION(L),TEMPP,NEW_ELEVATION(L),NEW_ELEVATION(L+1),WATER_LEVEL(L)
                    3884  FORMAT(' NEG SED GRAD C, I=',I5,' J=',I5,' S=',G12.5,/, &
                    ' E=',G12.5,' EN=',G12.5, ' EW=',G12.5)
                ENDIF
            ENDDO
            !      ******************************************************
            !      ******************************************************
            IF (PRINTSTUFF) THEN
                RETURN
            ELSE
                RETURN
            ENDIF

            !        *********************************************************************
            !           We come here if we have a fixed downstream non-alluvial point that cannot
            !           be buried - We will backgrade from that point
            !        *********************************************************************
            119   CONTINUE
            K=KKK
            I119=I119+1
            IF (I119 > 10*RMMX) THEN
                RETURN 
            ENDIF
            FIXED=.FALSE.
            I_NEW=I_LOCATION(K)
            J_NEW=J_LOCATION(K)
            !        *********************************************************************
            !          We end up here if we are backgrading from a fixed non-alluvial location
            !           we work upstream, making each elevation correspond to the product of
            !           elevation of the next point downstream plus the product of distance and
            !           gradient of the intervening segment
            !        *********************************************************************
            120  CONTINUE
            END_ELEVATION=ELEVATION(I_NEW,J_NEW)
            IF (K == 1) THEN
                VOLEND=0.0
                VOLSTART=0.0
                BACKGRADE=.TRUE.
                GOTO 903
            ENDIF
            I120=I120+1
            IF (I120 > 10*RMMX) THEN
                RETURN
            ENDIF
            SEDIMENT_CARRYOVER=0.0
            DEFAULT_DRAINAGE_AREA=0
            !        *********************************************************************
            !          volstart is the initial volume beneath the downstream flowpath
            !          volend is the final volume beneath the downstream flowpath
            !          provisional_elevation(k) holds the provisional new elevations backgraded from location k
            !        *********************************************************************
            KI=I_LOCATION(K-1)
            KJ=J_LOCATION(K-1)
            VOLSTART=STARTING_ELEVATION(K-1)*CELL_AREA
            KLOWER=1
            PROVISIONAL_ELEVATION(K-1)=END_ELEVATION+ALLUVIAL_GRADIENT(K-1)*STEP_DISTANCE(K-1)
            !        *********************************************************************
            !         If we submerged a non-alluvial location, convert it to alluvial, and note the
            !           elevation of the base of the deposit for later reference - We also need to
            !           compute its equilibrium gradient
            !        *********************************************************************
            IF (.NOT.IS_SEDIMENT_COVERED(KI,KJ)) THEN
                IF (ELEVATION(KI,KJ) > PROVISIONAL_ELEVATION(K-1)) THEN
                    PROVISIONAL_ELEVATION(K-1)=ELEVATION(KI,KJ)
                    KLOWER = K-1
                ENDIF
            ENDIF
            IF (PROVISIONAL_ELEVATION(K-1) < SEDIMENT_BASE(KI,KJ)) THEN
                PROVISIONAL_ELEVATION(K-1)=SEDIMENT_BASE(KI,KJ)
                KLOWER = K-1
            ENDIF
            VOLEND=PROVISIONAL_ELEVATION(K-1)*CELL_AREA
            IF (K > 2) THEN
                L200: DO  L=K-2,1,-1
                    KI=I_LOCATION(L)
                    KJ=J_LOCATION(L)
                    PROVISIONAL_ELEVATION(L)=PROVISIONAL_ELEVATION(L+1)+ALLUVIAL_GRADIENT(L)*STEP_DISTANCE(L)
                    IF (.NOT.IS_SEDIMENT_COVERED(KI,KJ)) THEN
                        IF (ELEVATION(KI,KJ) > PROVISIONAL_ELEVATION(L)) THEN
                            PROVISIONAL_ELEVATION(L)=ELEVATION(KI,KJ)
                            IF(KLOWER == 1) KLOWER=L
                        ENDIF
                    ENDIF
                    IF (PROVISIONAL_ELEVATION(L) < SEDIMENT_BASE(KI,KJ)) THEN
                        PROVISIONAL_ELEVATION(L)=SEDIMENT_BASE(KI,KJ)
                        IF (KLOWER == 1) KLOWER = L
                    ENDIF
                    IF(KLOWER == 1 ) THEN
                        VOLEND=VOLEND+PROVISIONAL_ELEVATION(L)*CELL_AREA
                        VOLSTART=VOLSTART+STARTING_ELEVATION(L)*CELL_AREA
                    ENDIF
                ENDDO L200
            ENDIF
            !        *********************************************************************
            !         We need to check that the new deposited volume is not greater than the
            !           available volume for deposition
            !         This should not occur, because we come here when we are limited by
            !           a downstream fixed location and we are generally not using all available
            !           sediment volume
            !         Assuming that there is no problem, then we make sed_surface_elevation(ii,jj) equal to the
            !           provisionan new elevations (provisional_elevation(k), and if there are bedrock locations
            !           thathat are buried, we set channel_state_change(ii,jj) equal to 1 so that they will be
            !           converted to alluvial
            !        *********************************************************************
            IF ((VOLEND-VOLSTART) > SEDIMENT_VOLUME) THEN
                CUMULEXCESS=CUMULEXCESS+VOLEND-VOLSTART-SEDIMENT_VOLUME
            ENDIF
            KI=I_LOCATION(K-1)
            KJ=J_LOCATION(K-1)
            SED_SURFACE_ELEVATION(KI,KJ)=PROVISIONAL_ELEVATION(K-1)
            LOCCHANGE=SED_SURFACE_ELEVATION(KI,KJ)-ELEVATION(KI,KJ)
            IF (ABS(LOCCHANGE) > ABSMAXCHG) THEN
                ABSMAXCHG=ABS(LOCCHANGE)
                MAXCHG=LOCCHANGE
                IMAX=KI
                JMAX=KJ
            ENDIF
            IF (DONE_ONCE(KI,KJ)) THEN
            ELSE
                DONE_ONCE(KI,KJ)=.TRUE.
            ENDIF
            IF (.NOT.IS_SEDIMENT_COVERED(KI,KJ)) THEN
                IF (ELEVATION(KI,KJ) > PROVISIONAL_ELEVATION(K-1)) THEN
                    SED_SURFACE_ELEVATION(KI,KJ)=ELEVATION(KI,KJ)
                    CHANNEL_STATE_CHANGE(KI,KJ)=0
                    ICASE(9)=ICASE(9)-1
                ELSE
                    CHANNEL_STATE_CHANGE(KI,KJ)=1
                    ICASE(10)=ICASE(10)+1
                    BEDLOAD_FLUX=MAXIMUM_SEDIMENT_YIELD(KI,KJ)
                    CALL EQUILIBRIUM_SEDIMENT_GRADIENT(BEDLOAD_FLUX,EQUILIBRIUM_GRADIENT(KI,KJ),MAXIMUM_DISCHARGE(KI,KJ))
                ENDIF
            ENDIF
            IF (SED_SURFACE_ELEVATION(KI,KJ) < ELOWEST) MUSTDOIT=.TRUE.
            IF (DO_ALLUVIAL_REROUTING) CALL CHECK_IF_CHANGE_FLOW_DIRECTION(KI,KJ)
            IF (K > 2) THEN
                L202: DO  L=K-2,KLOWER,-1
                    KI=I_LOCATION(L)
                    KJ=J_LOCATION(L)
                    SED_SURFACE_ELEVATION(KI,KJ)=PROVISIONAL_ELEVATION(L)
                    LOCCHANGE=SED_SURFACE_ELEVATION(KI,KJ)-ELEVATION(KI,KJ)
                    IF (ABS(LOCCHANGE) > ABSMAXCHG) THEN
                        ABSMAXCHG=ABS(LOCCHANGE)
                        MAXCHG=LOCCHANGE
                        IMAX=KI
                        JMAX=KJ
                    ENDIF
                    IF (DONE_ONCE(KI,KJ)) THEN
                    ELSE
                        DONE_ONCE(KI,KJ)=.TRUE.
                    ENDIF
                    IF (.NOT.IS_SEDIMENT_COVERED(KI,KJ)) THEN
                        IF (ELEVATION(KI,KJ) > PROVISIONAL_ELEVATION(L)) THEN
                            SED_SURFACE_ELEVATION(KI,KJ)=ELEVATION(KI,KJ)
                            IF (ABS(LOCCHANGE) > ABSMAXCHG) THEN
                                ABSMAXCHG=ABS(LOCCHANGE)
                                MAXCHG=LOCCHANGE
                                IMAX=KI
                                JMAX=KJ
                            ENDIF
                            CHANNEL_STATE_CHANGE(KI,KJ)=0
                            ICASE(11)=ICASE(11)-1
                        ELSE
                            CHANNEL_STATE_CHANGE(KI,KJ)=1
                            ICASE(12)=ICASE(12)+1
                            BEDLOAD_FLUX=MAXIMUM_SEDIMENT_YIELD(KI,KJ)
                            CALL EQUILIBRIUM_SEDIMENT_GRADIENT(BEDLOAD_FLUX, &
                            EQUILIBRIUM_GRADIENT(KI,KJ),MAXIMUM_DISCHARGE(KI,KJ))
                        ENDIF
                    ENDIF
                    IF (SED_SURFACE_ELEVATION(KI,KJ) < ELOWEST) MUSTDOIT=.TRUE.
                    IF (DO_ALLUVIAL_REROUTING) CALL CHECK_IF_CHANGE_FLOW_DIRECTION(KI,KJ)
                ENDDO L202
            ENDIF
            IF (ABSMAXCHG >= CRITELCHANGE) THEN
                ITYYP=2
                CALL  PRINT_SEDIMENT_DIAGNOSTICS(IMAX,JMAX,II,JJ,ITYYP,SEDSOURCE,ABSMAXCHG,DELT)
            ENDIF
            IF ((VOLEND-VOLSTART) > (1.1*SEDIMENT_VOLUME)) THEN
                !        *********************************************************************
                !         We should only be here if there is a problem
                !        *********************************************************************
                EDIFF=VOLEND-VOLSTART
            ENDIF
            BACKGRADE=.TRUE.
            TEMPP=VOLEND-VOLSTART
            GOTO 1220
            EXIT L999
            ENDDO L999

        !      ******************************************************
        !      ******************************************************
        !        Some debugging
        DO L=1,K-1
            TEMPP=(NEW_ELEVATION(L)-NEW_ELEVATION(L+1))/CELL_SIZE
            IF (TEMPP < 0.0) THEN
                WRITE(*,3885) I_LOCATION(L),J_LOCATION(L),TEMPP,NEW_ELEVATION(L),NEW_ELEVATION(L+1),WATER_LEVEL(L)
                3885  FORMAT(' NEG SED GRAD D, I=',I5,' J=',I5,' S=',G12.5,/, &
                ' E=',G12.5,' EN=',G12.5, ' EW=',G12.5)
            ENDIF
        ENDDO
        !      ******************************************************
        !      ******************************************************
        RETURN
   
    END ! SUBROUTINE ROUTE_SEDIMENT! 
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     SUBROUTINE SMOOTHSED
        USE ERODE_GLOBALS, CHANNEL_STATE_CHANGE=>IDO, WORKING_MATRIX=>CFNE,  &
        SED_FLUX_DIVERGENCE=>CFW, SED_SURFACE_ELEVATION=>CFN
        IMPLICIT NONE
        INTEGER I,J,TOTALWEIGHT,IE,IW,JN,JS
        LOGICAL USENW,USENE,USESW,USESE,USEW,USEE,USEN,USES
        LOGICAL USEIDO
        REAL(4) :: ENW,ENE,ESE,ESW,EE,EW,EN,ES,EDIFF,ECENTER,WORKING_MATRIXSUM
        USEIDO=.TRUE.
        WORKING_MATRIX=1.0E+25
        WORKING_MATRIXSUM=0.0
        DO  J=1,MYY
            DO  I=1,MX
                IF (SEDIMENT_DEPOSITED(I,J)) THEN
                    TOTALWEIGHT=0
                    USENW=.FALSE.
                    USENE=.FALSE.
                    USESW=.FALSE.
                    USESE=.FALSE.
                    USEW=.FALSE.
                    USEE=.FALSE.
                    USEN=.FALSE.
                    USES=.FALSE.
                    IE=I+1
                    IW=I-1
                    IF (IS_X_PERIODIC) THEN
                        IF (IE > MX) IE=1
                        IF (IW < 1) IW=MX
                    ELSE
                        IF (IE > MX) IE=MX-1
                        IF (IW < 1) IW=2
                    ENDIF
                    JN=J-1
                    JS=J+1
                    IF (IS_Y_PERIODIC) THEN
                        IF (JN < 1)JN=MY
                        IF (JS > MY) JS=1
                    ELSE
                        IF (JN < 1) JN=2
                    ENDIF
                    IF (((IS_SEDIMENT_COVERED(IW,JN)).OR.   &
                    (USEIDO.AND.(CHANNEL_STATE_CHANGE(IW,JN) > 0)))  &
                    .AND.(SEDIMENT_DEPOSITED(IW,JN))) THEN
                        IF (SED_SURFACE_ELEVATION(IW,JN) < 1.0E+24) THEN
                            ENW=SED_SURFACE_ELEVATION(IW,JN)
                        ELSE
                            ENW=ELEVATION(IW,JN)
                        ENDIF
                        TOTALWEIGHT=TOTALWEIGHT+1
                        USENW=.TRUE.
                    ENDIF
                    IF ((IS_Y_PERIODIC).OR.(JS < MY)) THEN
                        IF (((IS_SEDIMENT_COVERED(IW,JS)).OR.  &
                        (USEIDO.AND.(CHANNEL_STATE_CHANGE(IW,JS) > 0))) &
                        .AND.(SEDIMENT_DEPOSITED(IW,JS))) THEN
                            IF (SED_SURFACE_ELEVATION(IW,JS) < 1.0E+24) THEN
                                ESW=SED_SURFACE_ELEVATION(IW,JS)
                            ELSE
                                ESW=ELEVATION(IW,JS)
                            ENDIF
                            TOTALWEIGHT=TOTALWEIGHT+1
                            USESW=.TRUE.
                        ENDIF
                    ENDIF
                    IF (((IS_SEDIMENT_COVERED(IE,JN)).OR. &
                    (USEIDO.AND.(CHANNEL_STATE_CHANGE(IE,JN) > 0))) &
                    .AND.(SEDIMENT_DEPOSITED(IE,JN))) THEN
                        IF (SED_SURFACE_ELEVATION(IE,JN) < 1.0E+24) THEN
                            ENE=SED_SURFACE_ELEVATION(IE,JN)
                        ELSE
                            ENE=ELEVATION(IE,JN)
                        ENDIF
                        TOTALWEIGHT=TOTALWEIGHT+1
                        USENE=.TRUE.
                    ENDIF
                    IF ((IS_Y_PERIODIC).OR.(JS < MY)) THEN
                        IF (((IS_SEDIMENT_COVERED(IE,JS)).OR. &
                        (USEIDO.AND.(CHANNEL_STATE_CHANGE(IE,JS) > 0))) &
                        .AND.(SEDIMENT_DEPOSITED(IE,JS))) THEN
                            IF (SED_SURFACE_ELEVATION(IE,JS) < 1.0E+24) THEN
                                ESE=SED_SURFACE_ELEVATION(IE,JS)
                            ELSE
                                ESE=ELEVATION(IE,JS)
                            ENDIF
                            TOTALWEIGHT=TOTALWEIGHT+1
                            USESE=.TRUE.
                        ENDIF
                    ENDIF
                    IF (((IS_SEDIMENT_COVERED(I,JN)).OR. &
                    (USEIDO.AND.(CHANNEL_STATE_CHANGE(I,JN) > 0))) &
                    .AND.(SEDIMENT_DEPOSITED(I,JN))) THEN
                        IF (SED_SURFACE_ELEVATION(I,JN) < 1.0E+24) THEN
                            EN=SED_SURFACE_ELEVATION(I,JN)
                        ELSE
                            EN=ELEVATION(I,JN)
                        ENDIF
                        TOTALWEIGHT=TOTALWEIGHT+1
                        USEN=.TRUE.
                    ENDIF
                    IF ((IS_Y_PERIODIC).OR.(JS < MY)) THEN
                        IF (((IS_SEDIMENT_COVERED(I,JS)).OR.  &
                        (USEIDO.AND.(CHANNEL_STATE_CHANGE(I,JS) > 0))) &
                        .AND.(SEDIMENT_DEPOSITED(I,JS))) THEN
                            IF (SED_SURFACE_ELEVATION(I,JS) < 1.0E+24) THEN
                                ES=SED_SURFACE_ELEVATION(I,JS)
                            ELSE
                                ES=ELEVATION(I,JS)
                            ENDIF
                            TOTALWEIGHT=TOTALWEIGHT+1
                            USES=.TRUE.
                        ENDIF
                    ENDIF
                    IF (((IS_SEDIMENT_COVERED(IW,J)).OR.   &
                    (USEIDO.AND.(CHANNEL_STATE_CHANGE(IW,J) > 0))) &
                    .AND.(SEDIMENT_DEPOSITED(IW,J))) THEN
                        IF (SED_SURFACE_ELEVATION(IW,J) < 1.0E+24) THEN
                            EW=SED_SURFACE_ELEVATION(IW,J)
                        ELSE
                            EW=ELEVATION(IW,J)
                        ENDIF
                        TOTALWEIGHT=TOTALWEIGHT+1
                        USEW=.TRUE.
                    ENDIF
                    IF (((IS_SEDIMENT_COVERED(IE,J)).OR. &
                    (USEIDO.AND.(CHANNEL_STATE_CHANGE(IE,J) > 0))) &
                    .AND.(SEDIMENT_DEPOSITED(IE,J))) THEN
                        IF (SED_SURFACE_ELEVATION(IE,J) < 1.0E+24) THEN
                            EE=SED_SURFACE_ELEVATION(IE,J)
                        ELSE
                            EE=ELEVATION(IE,J)
                        ENDIF
                        TOTALWEIGHT=TOTALWEIGHT+1
                        USEE=.TRUE.
                    ENDIF
                    IF (TOTALWEIGHT > 0) THEN
                        IF (WORKING_MATRIX(I,J) > 1.0E+24) WORKING_MATRIX(I,J)=0.0
                        IF (SED_SURFACE_ELEVATION(I,J) < 1.0E+24) THEN
                            ECENTER=SED_SURFACE_ELEVATION(I,J)
                        ELSE
                            ECENTER=ELEVATION(I,J)
                        ENDIF
                        IF (USENW) THEN
                            EDIFF=(ECENTER-ENW)
                            WORKING_MATRIX(I,J)=WORKING_MATRIX(I,J)-EDIFF
                            IF(WORKING_MATRIX(IW,JN) > 1.0E+24) THEN
                                WORKING_MATRIX(IW,JN)=EDIFF
                            ELSE
                                WORKING_MATRIX(IW,JN)=WORKING_MATRIX(IW,JN)+EDIFF
                            ENDIF
                        ENDIF
                        IF (USENE)THEN
                            EDIFF=(ECENTER-ENE)
                            WORKING_MATRIX(I,J)=WORKING_MATRIX(I,J)-EDIFF
                            IF(WORKING_MATRIX(IE,JN) > 1.0E+24) THEN
                                WORKING_MATRIX(IE,JN)=EDIFF
                            ELSE
                                WORKING_MATRIX(IE,JN)=WORKING_MATRIX(IE,JN)+EDIFF
                            ENDIF
                        ENDIF
                        IF (USESE)THEN
                            EDIFF=(ECENTER-ESE)
                            WORKING_MATRIX(I,J)=WORKING_MATRIX(I,J)-EDIFF
                            IF(WORKING_MATRIX(IE,JS) > 1.0E+24) THEN
                                WORKING_MATRIX(IE,JS)=EDIFF
                            ELSE
                                WORKING_MATRIX(IE,JS)=WORKING_MATRIX(IE,JS)+EDIFF
                            ENDIF
                        ENDIF
                        IF (USESW)THEN
                            EDIFF=(ECENTER-ESW)
                            WORKING_MATRIX(I,J)=WORKING_MATRIX(I,J)-EDIFF
                            IF(WORKING_MATRIX(IW,JS) > 1.0E+24) THEN
                                WORKING_MATRIX(IW,JS)=EDIFF
                            ELSE
                                WORKING_MATRIX(IW,JS)=WORKING_MATRIX(IW,JS)+EDIFF
                            ENDIF
                        ENDIF
                        IF (USES)THEN
                            EDIFF=4.0*(ECENTER-ES)
                            WORKING_MATRIX(I,J)=WORKING_MATRIX(I,J)-EDIFF
                            IF(WORKING_MATRIX(I,JS) > 1.0E+24) THEN
                                WORKING_MATRIX(I,JS)=EDIFF
                            ELSE
                                WORKING_MATRIX(I,JS)=WORKING_MATRIX(I,JS)+EDIFF
                            ENDIF
                        ENDIF
                        IF (USEN)THEN
                            EDIFF=4.0*(ECENTER-EN)
                            WORKING_MATRIX(I,J)=WORKING_MATRIX(I,J)-EDIFF
                            IF(WORKING_MATRIX(I,JN) > 1.0E+24) THEN
                                WORKING_MATRIX(I,JN)=EDIFF
                            ELSE
                                WORKING_MATRIX(I,JN)=WORKING_MATRIX(I,JN)+EDIFF
                            ENDIF
                        ENDIF
                        IF (USEW)THEN
                            EDIFF=4.0*(ECENTER-EW)
                            WORKING_MATRIX(I,J)=WORKING_MATRIX(I,J)-EDIFF
                            IF(WORKING_MATRIX(IW,J) > 1.0E+24) THEN
                                WORKING_MATRIX(IW,J)=EDIFF
                            ELSE
                                WORKING_MATRIX(IW,J)=WORKING_MATRIX(IW,J)+EDIFF
                            ENDIF
                        ENDIF
                        IF (USEE)THEN
                            EDIFF=4.0*(ECENTER-EE)
                            WORKING_MATRIX(I,J)=WORKING_MATRIX(I,J)-EDIFF
                            IF(WORKING_MATRIX(IE,J) > 1.0E+24) THEN
                                WORKING_MATRIX(IE,J)=EDIFF
                            ELSE
                                WORKING_MATRIX(IE,J)=WORKING_MATRIX(IE,J)+EDIFF
                            ENDIF
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        DO  J=1,MYY
            DO  I=1,MX
                IF (WORKING_MATRIX(I,J) < 1.0E+24) THEN
                    IF (SED_SURFACE_ELEVATION(I,J) < 1.0E+24) THEN
                        SED_SURFACE_ELEVATION(I,J)=SED_SURFACE_ELEVATION(I,J)+WORKING_MATRIX(I,J)*ALLUVIUM_SMOOTHING_FACTOR
                    ELSE
                        IF (SED_FLUX_DIVERGENCE(I,J) < 0.0) THEN
                            SED_SURFACE_ELEVATION(I,J)=ELEVATION(I,J)+SED_FLUX_DIVERGENCE(I,J)*TIME_INCREMENT &
                            +WORKING_MATRIX(I,J)*ALLUVIUM_SMOOTHING_FACTOR
                        ELSE
                            SED_SURFACE_ELEVATION(I,J)=ELEVATION(I,J)+WORKING_MATRIX(I,J)*ALLUVIUM_SMOOTHING_FACTOR
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        RETURN
    END ! SUBROUTINE SMOOTHSED
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE CHECK_IF_CHANGE_FLOW_DIRECTION(I,J)
        USE ERODE_GLOBALS, SED_SURFACE_ELEVATION=>CFN
        IMPLICIT NONE
        INTEGER, INTENT(IN):: I,J
        INTEGER :: II,JJ,KK,IADD,JADD,IX,JX,IPROV,JPROV
        REAL(4) :: EL1,EL2,FF,PROVGRAD
    !  MODIFIES:  D8_GRADIENT, FLOW_DIRECTION
        IF (FLOW_DIRECTION(I,J)==1) RETURN        
        L401: DO  IX=I-1,I+1
            II=IX
            IF (IS_X_PERIODIC) THEN
                IF (II < 1) II=II+MX
                IF (II > MX) II=II-MX
            ELSE
                IF (II < 1) CYCLE L401
                IF (II > MX) CYCLE L401
            ENDIF
            L400: DO  JX=J-1,J+1
                JJ=JX
                IF (IS_Y_PERIODIC) THEN
                    IF (JJ < 1) JJ=JJ+MY
                    IF (JJ > MY) JJ=JJ-MY
                ELSE
                    IF (JJ < 1) CYCLE L400
                    IF (JJ >= MY) CYCLE L400
                ENDIF
                D8_GRADIENT(II,JJ)=-1.0E+24
            ENDDO L400
        ENDDO L401
        L100: DO KK=2,9
            FF=1.0
            IF (KK > 5) FF = 0.7071068
            IADD = DOWNSTREAM(KK,1)
            JADD = DOWNSTREAM(KK,2)
            L200: DO JX=J-1,J+1
                JJ=JX
                IF (IS_Y_PERIODIC) THEN
                    IF (JJ < 1) JJ=JJ+MY
                    IF (JJ > MY) JJ=JJ-MY
                ELSE
                    IF (JJ < 1) CYCLE L200
                    IF (JJ >= MY) CYCLE L200
                ENDIF
                JPROV =JJ+JADD
                IF (IS_Y_PERIODIC) THEN
                    IF (JPROV > MY) JPROV=JPROV-MY
                    IF (JPROV < 1) JPROV=JPROV+MY
                ELSE
                    IF (JPROV > MY) CYCLE L200
                    IF (JPROV < 1) CYCLE L200
                ENDIF
                L201: DO IX=I-1,I+1
                    II=IX
                    IPROV =II+IADD
                    IF (IS_X_PERIODIC) THEN
                        IF (IPROV > MX) IPROV=IPROV-MX
                        IF (IPROV < 1) IPROV=IPROV+MX
                        IF (II > MX) II=II-MX
                        IF (II < 1) II=II+MX
                    ELSE
                        IF (IPROV > MX) CYCLE L201
                        IF (IPROV < 1) CYCLE L201
                        IF (II > MX) CYCLE L201
                        IF (II < 1) CYCLE L201
                    ENDIF
                    IF (SED_SURFACE_ELEVATION(II,JJ) < 1.0E+24) THEN
                        EL1=SED_SURFACE_ELEVATION(II,JJ)
                    ELSE
                        EL1=ELEVATION(II,JJ)
                    ENDIF
                    IF (SED_SURFACE_ELEVATION(IPROV,JPROV) < 1.0E+24) THEN
                        EL2=SED_SURFACE_ELEVATION(IPROV,JPROV)
                    ELSE
                        EL2=ELEVATION(IPROV,JPROV)
                    ENDIF
                    PROVGRAD=(EL1-EL2)*FF/CELL_SIZE
                    IF (PROVGRAD > D8_GRADIENT(II,JJ)) THEN
                        D8_GRADIENT(II,JJ)=PROVGRAD
!!!! TESTING BUGFIX                        
                        IF (FLOW_DIRECTION(II,JJ)/=1) FLOW_DIRECTION(II,JJ)=KK
                    ENDIF
                ENDDO L201
            ENDDO L200
        ENDDO L100
        L300: DO  IX=I-1,I+1
            II=IX
            IF (IS_X_PERIODIC) THEN
                IF (II < 1) II=II+MX
                IF (II > MX) II=II-MX
            ELSE
                IF (II < 1) CYCLE L300
                IF (II > MX) CYCLE L300
            ENDIF
            L301: DO JX=J-1,J+1
                JJ=JX
                IF (IS_Y_PERIODIC) THEN
                    IF (JJ < 1) JJ=JJ+MY
                    IF (JJ > MY) JJ=JJ-MY
                ELSE
                    IF (JJ < 1) CYCLE L301
                    IF (JJ >= MY) CYCLE L301
                ENDIF
                IF (D8_GRADIENT(II,JJ) <= 0.0) THEN
!!!! TESTING BUGFIX                
                     IF (FLOW_DIRECTION(II,JJ)/=1) FLOW_DIRECTION(II,JJ)=-FLOW_DIRECTION(II,JJ)
                    D8_GRADIENT(II,JJ)=0.0
                ENDIF
            ENDDO L301
        ENDDO L300
        RETURN
    END !  SUBROUTINE CHECK_IF_CHANGE_FLOW_DIRECTION
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PRINT_SEDIMENT_DIAGNOSTICS(I,J,IX,JX,ITYYP,SEDSOURCE,MAXCHANGE,DELT)
        USE ERODE_GLOBALS, CHANNEL_STATE_CHANGE=>IDO, SED_SURFACE_ELEVATION=>CFN
        USE SEDROUTE_GLOBALS
        USE SEDDEBUG_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,J, ITYYP
        INTEGER :: L,IX, JX
        INTEGER :: ITEMP1,ITEMP2,ITEMP3,ITEMP4,ITEMP5
        INTEGER :: IDEPR,   &
        I_LAVA,IBACK
        LOGICAL, INTENT(IN) :: SEDSOURCE
        REAL(4) :: XWATER,MAXCHANGE,DELT
        IF (WRITE_SEDIMENT_DIAGNOSTICS) THEN
        !RETURN
        OPEN(OUTSEDDEBUG,FILE='SEDTEMP.DAT',POSITION='APPEND')
        IF (IS_DEPRESSION) THEN
            IDEPR=1
        ELSE
            IDEPR=0
        ENDIF
        IF (SEDSOURCE) THEN
            I_LAVA=1
        ELSE
            I_LAVA=0
        ENDIF
        IF (DO_BACKTRACK) THEN
            IBACK=1
        ELSE
            IBACK=0
        ENDIF
        WRITE(OUTSEDDEBUG,8119) ITYYP,MAXCHANGE,I,J,K,ITERATION,I_NEW,J_NEW,DELT
        8119     FORMAT(           &
        'ITYPE=',I5,' MXCHG=',G12.5,' I=',I5,' J=',I5,/,' K=',I5,' ITER=',I9, &
        ' I_NEW=',I5,' J_NEW=',I5,' DELT=',G12.5)
        WRITE(OUTSEDDEBUG,8709) UPSTREAM_ELEVATION,END_ELEVATION,EL1,EL2,ALLUVIAL_START_GRADIENT, &
        PROVISIONAL_START_GRADIENT,PROVISIONAL_END_GRADIENT,SEDIMENT_VOLUME
        8709           FORMAT(' UPSTR_EL=',G12.5,' END_EL=',G12.5,  &
        ' EL1=',   &
        G12.5,' EL2=',G12.5,/,' ALL_STRT_GRD=',G12.5,' PROV_STRT_GRD=' &
        ,G12.5, &
        ' PROV_END_GRD=',G12.5,' SED_VOL=',G12.5)
        WRITE(OUTSEDDEBUG,18119) IDEPR,I_LAVA,IBACK
        18119 FORMAT(' DEPR=',I5,' SEDSRC=',I5,' BKTRK=',I5)
        WRITE(OUTSEDDEBUG,1447)
1447    FORMAT('  K  I   J  ID  CHG  STRT_EL  SED_SRF_EL     EL      '&
        ,'  ACT_GRD   ', &
        '   ALVL GRD','   SED_FLX',  &
        '  PROV_EL   ',    &
        'NEW_EL     WTR_LEV    DISCH     SDBSE    ', &
        '  SED_YLD   BC ONCE SC  SB  RX     BN  LAKE_EL')
         IF (K > 0) THEN
            ITEMP1=0
            ITEMP2=0
            IF (IS_SEDIMENT_COVERED(IX,JX)) THEN
                ITEMP3=1
            ELSE
                ITEMP3=0
            ENDIF
            IF (SUBMERGED(IX,JX)) THEN
                ITEMP4=1
            ELSE
                ITEMP4=0
            ENDIF
            IF (IS_ROCK_SURFACE(IX,JX)) THEN
                ITEMP5=1
            ELSE
                ITEMP5=0
            ENDIF
            L=0
            IF (COMPLETE_RUNOFF) THEN
                XWATER=AMAX1(OCEAN_ELEVATION,LAKE_SURFACE_ELEVATION(BASIN_NUMBER(IX,JX)))
            ELSE
                XWATER=OCEAN_ELEVATION
            ENDIF
            WRITE(OUTSEDDEBUG,1449) L,IX,JX,FLOW_DIRECTION(IX,  &
            JX),CHANNEL_STATE_CHANGE(IX,JX),UPSTREAM_ELEVATION,  &
            SED_SURFACE_ELEVATION(IX,JX),ELEVATION(IX,JX),    &
            ACTUAL_START_GRADIENT,ALLUVIAL_START_GRADIENT,SEDIMENT_FLUX(IX,JX),  &
            UPSTREAM_ELEVATION,UPSTREAM_ELEVATION,XWATER,MAXIMUM_DISCHARGE(IX, &
            JX),SEDIMENT_BASE(IX,JX),    &
            SEDIMENT_YIELD(IX,JX),ITEMP1,ITEMP2, &
            ITEMP3,ITEMP4,ITEMP5, &
            BASIN_NUMBER(IX,JX),LAKE_SURFACE_ELEVATION(BASIN_NUMBER(IX,JX))
            DO  L=1,K
                IF (IS_BEDROCK_CHANNEL(L)) THEN
                    ITEMP1=1
                ELSE
                    ITEMP1=0
                ENDIF
                IF (DONE_ONCE(I_LOCATION(L),J_LOCATION(L))) THEN
                    ITEMP2=1
                ELSE
                    ITEMP2=0
                ENDIF
                IF (IS_SEDIMENT_COVERED(I_LOCATION(L),J_LOCATION(L))) THEN
                    ITEMP3=1
                ELSE
                    ITEMP3=0
                ENDIF
                IF (SUBMERGED(I_LOCATION(L),J_LOCATION(L))) THEN
                    ITEMP4=1
                ELSE
                    ITEMP4=0
                ENDIF
                IF (IS_ROCK_SURFACE(I_LOCATION(L),J_LOCATION(L))) THEN
                    ITEMP5=1
                ELSE
                    ITEMP5=0
                ENDIF
                WRITE(OUTSEDDEBUG,1449) L,I_LOCATION(L),J_LOCATION(L),FLOW_DIRECTION(I_LOCATION(L), &
                J_LOCATION(L)),CHANNEL_STATE_CHANGE(I_LOCATION(L),J_LOCATION(L)),STARTING_ELEVATION(L),  &
                SED_SURFACE_ELEVATION(I_LOCATION(L),J_LOCATION(L)),ELEVATION(I_LOCATION(L),J_LOCATION(L)), &
                ACTUAL_GRADIENT(L),ALLUVIAL_GRADIENT(L),SEDIMENT_FLUX(I_LOCATION(L),J_LOCATION(L)), &
                PROVISIONAL_ELEVATION(L),NEW_ELEVATION(L),WATER_LEVEL(L),MAXIMUM_DISCHARGE(I_LOCATION(L), &
                J_LOCATION(L)),SEDIMENT_BASE(I_LOCATION(L),J_LOCATION(L)),  &
                SEDIMENT_YIELD(I_LOCATION(L),J_LOCATION(L)),ITEMP1,ITEMP2, &
                ITEMP3,ITEMP4,ITEMP5, &
                BASIN_NUMBER(I_LOCATION(L),J_LOCATION(L))  &
                ,LAKE_SURFACE_ELEVATION(BASIN_NUMBER(I_LOCATION(L),J_LOCATION(L)))
                1449   FORMAT(5(I3,' '),12(G11.4),' ',5(I1,'   '),' ',I5,' ',G11.4)
            ENDDO
            ITEMP1=0
            ITEMP2=0
            IF (IS_SEDIMENT_COVERED(I_NEW,J_NEW)) THEN
                ITEMP3=1
            ELSE
                ITEMP3=0
            ENDIF
            IF (SUBMERGED(I_NEW,J_NEW)) THEN
                ITEMP4=1
            ELSE
                ITEMP4=0
            ENDIF
            IF (IS_ROCK_SURFACE(I_NEW,J_NEW)) THEN
                ITEMP5=1
            ELSE
                ITEMP5=0
            ENDIF
            L=K+1
            IF (COMPLETE_RUNOFF) THEN
                XWATER=AMAX1(OCEAN_ELEVATION,LAKE_SURFACE_ELEVATION(BASIN_NUMBER(I_NEW,J_NEW)))
            ELSE
                XWATER=OCEAN_ELEVATION
            ENDIF
            WRITE(OUTSEDDEBUG,1449) L,I_NEW,J_NEW,FLOW_DIRECTION(I_NEW,  &
            J_NEW),CHANNEL_STATE_CHANGE(I_NEW,J_NEW),END_ELEVATION, &
            SED_SURFACE_ELEVATION(I_NEW,J_NEW),ELEVATION(I_NEW,J_NEW),  &
            PROVISIONAL_END_GRADIENT,PROVISIONAL_END_GRADIENT,SEDIMENT_FLUX(I_NEW,J_NEW),  &
            END_ELEVATION,END_ELEVATION,XWATER,MAXIMUM_DISCHARGE(I_NEW, &
            J_NEW),SEDIMENT_BASE(I_NEW,J_NEW), &
            SEDIMENT_YIELD(I_NEW,J_NEW),ITEMP1,ITEMP2,&
            ITEMP3,ITEMP4,ITEMP5,&
            BASIN_NUMBER(I_NEW,J_NEW)  &
            ,LAKE_SURFACE_ELEVATION(BASIN_NUMBER(I_NEW,J_NEW))
        ENDIF
        CLOSE(OUTSEDDEBUG)
        ENDIF
        RETURN
    END ! UBROUTINE PRINT_SEDIMENT_DIAGNOSTICS
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
