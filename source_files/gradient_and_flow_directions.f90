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
    SUBROUTINE GRADIENT_AND_FLOW_DIRECTION()
        USE ERODE_GLOBALS, PAST_FLOW_DIRECTION=>IDO, PAST_GRADIENT=>CFNW
        IMPLICIT NONE
        REAL(4) :: PROVGRAD,FF,GRADIENT_RATIO,GRADOLD,GRADAVG,NGRAD, RRAND
        INTEGER :: I,J,KK,IADD,JADD,IPROV,JPROV,IIIDO,IIGRAD,JJGRAD
        LOGICAL :: CHANGE_FLOW_DIRECTION
        EXTERNAL CHANGE_FLOW_DIRECTION, RRAND
        !      *******************************************************************
        !       Calculates maximum downslope gradient and drainage directions
        !         Uses d8 flow direction modeling or alternatively a random gradient determination
        !   MODIFIES: IDO, CFNW, D8_GRADIENT,FLOW_DIRECTION, IS_ROCK_SURFACE, REGOLITH
        !   CALLS: RRAND, CHANGE_FLOW_DIRECTION
        !      *******************************************************************
        !      *******************************************************************
        !       maxgradient is maximum gradient in matrix
        !      *******************************************************************
        MAXGRADIENT = -1.0E+10
        GRADAVG=0.0
        NGRAD=0.0
        !      *******************************************************************
        !       Set up initial values
        !      *******************************************************************
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            PAST_GRADIENT=D8_GRADIENT
        ENDIF
        PAST_FLOW_DIRECTION=FLOW_DIRECTION
        D8_GRADIENT=-1.0E+5
        !      *******************************************************************
        !       Cycle through all directions and all matrix points, calculating
        !       the maximum gradient and the direction of the maximum gradient
        !      *******************************************************************
        DO  KK=2,5
            IADD = DOWNSTREAM(KK,1)
            JADD = DOWNSTREAM(KK,2)
            DO  J=2,MY-1
                JPROV =J+JADD
                DO  I=2,MX-1
                    IPROV =I+IADD
                    PROVGRAD=(ELEVATION(I,J)-ELEVATION(IPROV,JPROV))
                    IF (PROVGRAD > D8_GRADIENT(I,J)) THEN
                        D8_GRADIENT(I,J)=PROVGRAD
                        FLOW_DIRECTION(I,J)=KK
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
        FF = 0.7071068
        DO  KK=6,9
            IF (USE_RANDOM_FLOWROUTING) THEN
                FF=1.0/(2.0-RRAND())
            ENDIF
            IADD = DOWNSTREAM(KK,1)
            JADD = DOWNSTREAM(KK,2)
            DO  J=2,MY-1
                JPROV =J+JADD
                DO  I=2,MX-1
                    IPROV =I+IADD
                    PROVGRAD=(ELEVATION(I,J)-ELEVATION(IPROV,JPROV))*FF
                    IF (PROVGRAD > D8_GRADIENT(I,J)) THEN
                        D8_GRADIENT(I,J)=PROVGRAD
                        FLOW_DIRECTION(I,J)=KK
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
        !      *******************************************************************
        !       Do the same for the drainage divide at the top of the matrix
        !      *******************************************************************
       ! IF (.NOT.DO_FLOW_BOUNDARIES) THEN
            J = 1
            DO  KK=2,5
                IADD = DOWNSTREAM(KK,1)
                JADD = DOWNSTREAM(KK,2)
                JPROV =J+JADD
                IF (IS_Y_PERIODIC) THEN
                    IF (JPROV < 1) JPROV=MY
                ENDIF
                DO  I=2,MX-1
                    IPROV =I+IADD
                    IF (JPROV > 0) THEN
                        PROVGRAD=(ELEVATION(I,J)-ELEVATION(IPROV,JPROV))
                        IF (PROVGRAD > D8_GRADIENT(I,J)) THEN
                            D8_GRADIENT(I,J)=PROVGRAD
                            FLOW_DIRECTION(I,J)=KK
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
            FF = 0.7071068
            DO  KK=6,9
                IF (USE_RANDOM_FLOWROUTING) THEN
                    FF=1.0/(2.0-RRAND())
                ENDIF
                IADD = DOWNSTREAM(KK,1)
                JADD = DOWNSTREAM(KK,2)
                JPROV =J+JADD
                IF (IS_Y_PERIODIC) THEN
                    IF (JPROV < 1) JPROV=MY
                ENDIF
                DO  I=2,MX-1
                    IPROV =I+IADD
                    IF (JPROV > 0) THEN
                        PROVGRAD=(ELEVATION(I,J)-ELEVATION(IPROV,JPROV))*FF
                        IF (PROVGRAD > D8_GRADIENT(I,J)) THEN
                            D8_GRADIENT(I,J)=PROVGRAD
                            FLOW_DIRECTION(I,J)=KK
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
            !      *******************************************************************
            !       Do the same at the bottom of the matrix
            !      *******************************************************************
            IF (IS_Y_PERIODIC) THEN
                J = MY
                DO  KK=2,5
                    IADD = DOWNSTREAM(KK,1)
                    JADD = DOWNSTREAM(KK,2)
                    JPROV =J+JADD
                    IF (JPROV > MY) JPROV=1
                    L504: DO  I=1,MX
                        IPROV =I+IADD
                        IF (IS_X_PERIODIC) THEN
                            IF (IPROV < 1) IPROV=MX
                            IF (IPROV > MX) IPROV=1
                        ELSE
                            IF (IPROV < 1) CYCLE L504
                            IF (IPROV > MX) CYCLE L504
                        ENDIF
                        IF (JPROV > 0) THEN
                            PROVGRAD=(ELEVATION(I,J)-ELEVATION(IPROV,JPROV))
                            IF (PROVGRAD > D8_GRADIENT(I,J)) THEN
                                D8_GRADIENT(I,J)=PROVGRAD
                                FLOW_DIRECTION(I,J)=KK
                            ENDIF
                        ENDIF
                    ENDDO L504
                ENDDO
                FF = 0.7071068
                L533: DO  KK=6,9
                    IF (USE_RANDOM_FLOWROUTING) THEN
                        FF=1.0/(2.0-RRAND())
                    ENDIF
                    IADD = DOWNSTREAM(KK,1)
                    JADD = DOWNSTREAM(KK,2)
                    JPROV =J+JADD
                    IF (JPROV > MY) JPROV=1
                    L534: DO  I=1,MX
                        IPROV =I+IADD
                        IF (IS_X_PERIODIC) THEN
                            IF (IPROV < 1) IPROV=MX
                            IF (IPROV > MX) IPROV=1
                        ELSE
                            IF (IPROV < 1) CYCLE L534
                            IF (IPROV > MX) CYCLE L534
                        ENDIF
                        IF (JPROV > 0) THEN
                            PROVGRAD=(ELEVATION(I,J)-ELEVATION(IPROV,JPROV))*FF
                            IF (PROVGRAD > D8_GRADIENT(I,J)) THEN
                                D8_GRADIENT(I,J)=PROVGRAD
                                FLOW_DIRECTION(I,J)=KK
                            ENDIF
                        ENDIF
                    ENDDO L534
                ENDDO    L533
            ENDIF

            !      *******************************************************************
            !       Do the same for the periodic left and right boundaries
            !      *******************************************************************
            I = 1
            DO  KK=2,5
                IADD = DOWNSTREAM(KK,1)
                JADD = DOWNSTREAM(KK,2)
                IPROV =I+IADD
                IF (IS_X_PERIODIC) THEN
                    IF (IPROV < 1) IPROV = MX
                ENDIF
                DO  J=1,MY-1
                    JPROV =J+JADD
                    IF (IS_Y_PERIODIC) THEN
                        IF (JPROV < 1) JPROV=MY
                    ENDIF
                    IF ((JPROV > 0).AND.(IPROV > 0)) THEN
                        PROVGRAD=(ELEVATION(I,J)-ELEVATION(IPROV,JPROV))
                        IF (PROVGRAD > D8_GRADIENT(I,J)) THEN
                            D8_GRADIENT(I,J)=PROVGRAD
                            FLOW_DIRECTION(I,J)=KK
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
            FF = 0.7071068
            DO  KK=6,9
                IF (USE_RANDOM_FLOWROUTING) THEN
                    FF=1.0/(2.0-RRAND())
                ENDIF
                IADD = DOWNSTREAM(KK,1)
                JADD = DOWNSTREAM(KK,2)
                IPROV =I+IADD
                IF (IS_X_PERIODIC) THEN
                    IF (IPROV < 1) IPROV = MX
                ENDIF
                DO  J=1,MY-1
                    JPROV =J+JADD
                    IF (IS_Y_PERIODIC) THEN
                        IF (JPROV < 1) JPROV=MY
                    ENDIF
                    IF ((JPROV > 0).AND.(IPROV > 0)) THEN
                        PROVGRAD=(ELEVATION(I,J)-ELEVATION(IPROV,JPROV))*FF
                        IF (PROVGRAD > D8_GRADIENT(I,J)) THEN
                            D8_GRADIENT(I,J)=PROVGRAD
                            FLOW_DIRECTION(I,J)=KK
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
            I = MX
            DO  KK=2,5
                IADD = DOWNSTREAM(KK,1)
                JADD = DOWNSTREAM(KK,2)
                IPROV =I+IADD
                IF (IS_X_PERIODIC) THEN
                    IF (IPROV > MX) IPROV = 1
                ENDIF
                DO  J=1,MY-1
                    JPROV =J+JADD
                    IF (IS_Y_PERIODIC) THEN
                        IF (JPROV < 1) JPROV=MY
                    ENDIF
                    IF ((JPROV > 0).AND.(IPROV <= MX)) THEN
                        PROVGRAD=(ELEVATION(I,J)-ELEVATION(IPROV,JPROV))
                        IF (PROVGRAD > D8_GRADIENT(I,J)) THEN
                            D8_GRADIENT(I,J)=PROVGRAD
                            FLOW_DIRECTION(I,J)=KK
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
            FF = 0.7071068
            DO  KK=6,9
                IF (USE_RANDOM_FLOWROUTING) THEN
                    FF=1.0/(2.0-RRAND())
                ENDIF
                IADD = DOWNSTREAM(KK,1)
                JADD = DOWNSTREAM(KK,2)
                IPROV =I+IADD
                IF (IS_X_PERIODIC) THEN
                    IF (IPROV > MX) IPROV = 1
                ENDIF
                DO  J=1,MY-1
                    JPROV =J+JADD
                    IF (IS_Y_PERIODIC) THEN
                        IF (JPROV < 1) JPROV=MY
                    ENDIF
                    IF ((JPROV > 0).AND.(IPROV <= MX)) THEN
                        PROVGRAD=(ELEVATION(I,J)-ELEVATION(IPROV,JPROV))*FF
                        IF (PROVGRAD > D8_GRADIENT(I,J)) THEN
                            D8_GRADIENT(I,J)=PROVGRAD
                            FLOW_DIRECTION(I,J)=KK
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
      !  ENDIF
        IF (USE_RANDOM_FLOWROUTING) THEN
            DO  I=1,MX
                DO  J=1,MYY
                    IF (FLOW_DIRECTION(I,J) > 5) THEN
                        IADD=DOWNSTREAM(FLOW_DIRECTION(I,J),1)
                        JADD=DOWNSTREAM(FLOW_DIRECTION(I,J),2)
                        IPROV=I+IADD
                        IF (IS_X_PERIODIC) THEN
                            IF (IPROV > MX) IPROV=1
                            IF (IPROV < 1) IPROV=MX
                        ENDIF
                        JPROV=J+JADD
                        IF (IS_Y_PERIODIC) THEN
                            IF (JPROV > MY) JPROV=1
                            IF (JPROV < 1) JPROV=MY
                        ENDIF
                        D8_GRADIENT(I,J)=(ELEVATION(I,J)- &
						     ELEVATION(IPROV,JPROV))*0.7071068
                    ENDIF
                ENDDO
            ENDDO
        ENDIF
        !      *******************************************************************
        !       Set gradient to zero if it is a depression (grad<=0.0) and also
        !       negate sign of flow_direction(i,j) to indicate a depression
        !      *******************************************************************
        DO  J=1,MY
            DO  I=1,MX
                IF (D8_GRADIENT(I,J) <= 0.0) THEN
                    FLOW_DIRECTION(I,J)=-FLOW_DIRECTION(I,J)
                    D8_GRADIENT(I,J)=0.0
                ENDIF
            ENDDO
        ENDDO
        !      *********************************************************************
        !        The id value of the exit border is 1
        !      *********************************************************************
            IF (.NOT.IS_Y_PERIODIC) THEN
                DO  I=1,MX
                    D8_GRADIENT(I,MY)=0.0
                    FLOW_DIRECTION(I,MY)=1
                ENDDO
            ENDIF
        IF (DO_SEDIMENT_TRANSPORT.AND.STICKY_SEDIMENT_ROUTING) THEN
            DO  J=1,MYY
                DO  I=1,MX
!!! TEST DEBUG                
                    IF (FLOW_DIRECTION(I,J)==1) CYCLE
                    IF (IS_SEDIMENT_COVERED(I,J).AND.(FLOW_DIRECTION(I,J) > 1)) THEN
                        IIIDO=PAST_FLOW_DIRECTION(I,J)
                        IF (IIIDO > 1) THEN
                            IF (FLOW_DIRECTION(I,J) /= IIIDO) THEN
                                GRADOLD=ELEVATION(I,J)-ELEVATION(I+DOWNSTREAM(IIIDO,1), &
								     J+DOWNSTREAM(IIIDO,2))
                                IF (GRADOLD > 0.0) THEN
                                    GRADIENT_RATIO=D8_GRADIENT(I,J)/GRADOLD
                                ELSE
                                    GRADIENT_RATIO=1.0E+25
                                ENDIF
                                IF (.NOT.CHANGE_FLOW_DIRECTION(GRADIENT_RATIO)) THEN
                                    FLOW_DIRECTION(I,J)=PAST_FLOW_DIRECTION(I,J)
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
        ENDIF		
        IF (DO_FLOW_BOUNDARIES) THEN  !! TEST DEC 9 2017
		    J=1
			DO I=1,MX
			    FLOW_DIRECTION(I,J)=1
			ENDDO
			J=MY
			DO I=1,MX
                FLOW_DIRECTION(I,J)=1
			    
			ENDDO
			I=1
			DO J=1,MY
                FLOW_DIRECTION(I,J)=1
                
			ENDDO
			I=MX
			DO J=1,MY
                FLOW_DIRECTION(I,J)=1
			    
			ENDDO
		ENDIF
        
        !      *******************************************************************
        !       Determine maximum gradient in matrix
        !      *******************************************************************
        DO  J=1,MY
            DO  I=1,MX
                IF (D8_GRADIENT(I,J) > 0.0) THEN
                    D8_GRADIENT(I,J)=D8_GRADIENT(I,J)/CELL_SIZE
                    GRADAVG=GRADAVG+D8_GRADIENT(I,J)
                    NGRAD=NGRAD+1.0
                ELSE
                    !     ********************************************************************
                    !      A depression is assumed to be a regolith-covered location
                    !     ********************************************************************
                    IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                       IS_ROCK_SURFACE(I,J)=.FALSE.
                       IF (REGOLITH(I,J) < 0.0) REGOLITH(I,J)=0.0
                    ENDIF
                ENDIF
                !      *******************************************************************
                !         Calculate the number of changes of drainage direction since
                !         the last iteration (numbidchange)
                !         and the total change in absolute and signed gradient
                !      *******************************************************************
                IF(FLOW_DIRECTION(I,J) /= PAST_FLOW_DIRECTION(I,J)) NUMBIDCHANGE=NUMBIDCHANGE+1.0
                SUMGRADCHANGE=SUMGRADCHANGE+D8_GRADIENT(I,J)-PAST_GRADIENT(I,J)
                ABSGRADCHANGE=ABSGRADCHANGE+ABS(D8_GRADIENT(I,J)-PAST_GRADIENT(I,J))
                IF (D8_GRADIENT(I,J) > MAXGRADIENT) THEN
                    MAXGRADIENT=D8_GRADIENT(I,J)
                    IIGRAD=I
                    JJGRAD=J
                ENDIF
                IF (D8_GRADIENT(I,J) > 10000.0) THEN
                    WRITE(OUTHIST,411) I,J
                    411           FORMAT(' GRADIENT>10000 AT I=',I5,' J=',I5)
                ENDIF
            ENDDO
        ENDDO
        GRADAVG=GRADAVG/NGRAD
        GRADCUT=CUT_RATIO*GRADAVG
        RETURN
    END !  SUBROUTINE GRADIENT_AND_FLOW_DIRECTION
