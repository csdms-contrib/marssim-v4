!     This is a program to generate a pseudo-fractal 2-D (X and Y) matrix for
!      input to the simulation (e.g., for initial elevations or rock resistance,
!      or hydraulic conductivity)  See the GETPARAMETERS subroutine for
!      further description
!     
      MODULE MATRIX_GLOBALS
      IMPLICIT NONE
      SAVE
      INTEGER :: MAXDIMENSION, MAXLEVEL
      real*4, ALLOCATABLE, DIMENSION(:,:) :: ELEVATION,ENEW
      real*4, ALLOCATABLE, DIMENSION(:) :: YAVG,LNEWELEVS,RNEWELEVS
      real*4, ALLOCATABLE, DIMENSION(:) :: TNEWELEVS,BNEWELEVS
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: ACTIVE, NEWACTIVE, DONE
      real*4 ::  RANDTEMP
      INTEGER, DIMENSION (2,4) ::  APLACE, BPLACE
      real*4, DIMENSION (2,2,4) :: WEIGHTS
      real*4, dimension(-2:2,-2:2) :: SMOOTH
      INTEGER :: ISHOW, KKK,INN, JNN,POINTS
      INTEGER :: RDUMMY
      real*4 :: FACTOR
      real*4 :: CONDUCTIVITY,ERROR, CHARDIMENSION
      INTEGER :: XREFLECT,CURLEVEL, INTERVAL,DOBLENDX,DODETREND
      INTEGER :: DOLEVELLOWER,YREFLECT,DOBLENDY,NO_BLENDS,BLEND_WIDTH
      real*4 :: MULFACTOR, ELAVERAGE
      real*4 :: SEED, ISEED1
      INTEGER :: TNUM,IACT, JACT, IST,IFN,JST,JFN,DOEXP
      real*4 :: AVGVAL,MAXVAL,MINVAL,REGIONALSLOPE
      LOGICAL :: USEEXP
      END MODULE MATRIX_GLOBALS
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real*4 FUNCTION NORMAL1(DUMMY)
	   USE MATRIX_GLOBALS
      IMPLICIT NONE
             INTEGER DUMMY,I,J
             real*4 X, TEMP,RAN3
             real*4 RAN,T1
              J = 1
              TEMP = 0.0
              DO 10 I = 1,48
                  CALL RANDOM_NUMBER(T1)
                   TEMP = TEMP + T1 -0.5
              NORMAL1 = TEMP / 2.0
10    CONTINUE
	 RETURN
	END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real*4 FUNCTION NORMAL2(DUMMY)
	USE MATRIX_GLOBALS
      IMPLICIT NONE
      INTEGER DUMMY, J
      real*4 V1, V2, SS,RAN3
      real*4 RAN,T1,T2
           J = 1
100       CONTINUE
                  CALL RANDOM_NUMBER(T1)
                  CALL RANDOM_NUMBER(T2)
      	V1 = 2.0*T1-1.0
              V2 = 2.0*T2-1.0
              SS = V1*V1 + V2*V2
           IF (SS.GE.1.0) GOTO 100
      SS = SQRT(-2.0*LOG(SS)/SS)
           NORMAL2 = V1*SS
      RETURN
      END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE SETUP
      USE MATRIX_GLOBALS      
      IMPLICIT NONE
      INTEGER I,J,K
      real*4 SUM,smoothsum
           APLACE(1,1) = -1
           APLACE(2,1) = -1
           APLACE(1,2) = -1
           APLACE(2,2) = 1
           APLACE(1,3) = 1
           APLACE(2,3) = -1
           APLACE(1,4) = 1
           APLACE(2,4) = 1
           BPLACE(1,1) = 0
           BPLACE(2,1) = 1
           BPLACE(1,2) = 0
           BPLACE(2,2) = -1
           BPLACE(1,3) = 1
           BPLACE(2,3) = 0
           BPLACE(1,4) = -1
           BPLACE(2,4) = 0
          WEIGHTS(1,1,1)=4.0/4.0
          WEIGHTS(1,1,2)=0.0/4.0
          WEIGHTS(1,1,3)=0.0/4.0
          WEIGHTS(1,1,4)=0.0/4.0
          WEIGHTS(2,1,1)=2.0/4.0
          WEIGHTS(2,1,2)=2.0/4.0
          WEIGHTS(2,1,3)=0.0/4.0
          WEIGHTS(2,1,4)=0.0/4.0
          WEIGHTS(1,2,1)=2.0/4.0
          WEIGHTS(1,2,2)=0.0/4.0
          WEIGHTS(1,2,3)=2.0/4.0
          WEIGHTS(1,2,4)=0.0/4.0
          WEIGHTS(2,2,1)=1.0/4.0
          WEIGHTS(2,2,2)=1.0/4.0
          WEIGHTS(2,2,3)=1.0/4.0
          WEIGHTS(2,2,4)=1.0/4.0
          DO 100 I=1, 2
          DO 100 J= 1, 2
              SUM = 0.0
              DO 110 K= 1,4
110              SUM = SUM+WEIGHTS(I,J,K)
              IF (SUM.NE.1.0) THEN 
                 WRITE(*,120) I,J
120              FORMAT(' BAD WEIGHTS AT I=',I5,' J=',I5)
                 STOP
              ENDIF
100       CONTINUE
         smooth(-2,-2)=-0.07428571
         smooth(-1,-2)=0.01142857
         smooth(0,-2)=0.04
         smooth(1,-2)=smooth(-1,-2)
         smooth(2,-2)=smooth(-2,-2)
         smooth(-2,-1)=smooth(-1,-2)
         smooth(-1,-1)=0.09714286
         smooth(0,-1)=0.12571429
         smooth(1,-1)=smooth(-1,-1)
         smooth(2,-1)=smooth(-2,-1)
         smooth(-2,0)=smooth(0,-2)
         smooth(-1,0)=smooth(0,-1)
         smooth(0,0)=0.15428571
         smooth(1,0)=smooth(-1,0)
         smooth(2,0)=smooth(-2,0)
         do j=1,2
             do i=-2,2
                 smooth(i,j)=smooth(i,-j)
             enddo
         enddo
         do j=-2,2
             write(*,170)(smooth(i,j),i=-2,2)
170          format (5(g13.6,' '))
             
         enddo
         smoothsum=0.0
         do i=-2,2
             do j=-2,2
                 smoothsum=smoothsum+smooth(i,j)
             enddo
         enddo
         write(*,150) smoothsum
150      format('sum of smoothing filter=',g13.6)
         !stop
         
      RETURN
      END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE GETPARAMETERS
      USE MATRIX_GLOBALS      
      IMPLICIT NONE
      INTEGER RNSEED, NRANDOM, IERROR
      INTEGER, DIMENSION(:), ALLOCATABLE :: RANDSEED
      READ(10,*) MAXDIMENSION, MAXLEVEL
      write(*,833) maxdimension, maxlevel
833   format(' maxdimension=',i5, ' maxlevel=',i5)
      ALLOCATE(ELEVATION(MAXDIMENSION,MAXDIMENSION), ENEW(MAXDIMENSION,MAXDIMENSION),STAT=IERROR)
      ALLOCATE(ACTIVE(MAXDIMENSION,MAXDIMENSION), STAT=IERROR)
      ALLOCATE(NEWACTIVE(MAXDIMENSION,MAXDIMENSION), STAT=IERROR)
      ALLOCATE(DONE(MAXDIMENSION,MAXDIMENSION), STAT=IERROR)
      ALLOCATE(YAVG(MAXDIMENSION), STAT=IERROR)
!		CONDUCTIVITY determines the magnitude of variation of the simulated values
      READ(10,*) CONDUCTIVITY
!        ERROR is also a multiplicative factor that determines
!        the magnitude of variation of the simulated values
!        which differs from CONDUCTIVITY in that if lognormal
!        values are generated it multiplies before the exponental
!        transformation is taken rather than after as with
!        CONDUCTIVITY
      READ(10,*) ERROR
!		CHARDIMENSION determines how rapidly the variation declines with scale
!		it is closely related to fractal dimension and should range between 0.0
!		(uncorrelated random noise) to 1.0 (rough at large scales, smooth at small
!		scales 
      READ(10,*) CHARDIMENSION
!        A random number generator seed
      READ(10,*) RNSEED
!		Set FACTOR to 1.0 to allow random variation of the corners and edges of the matrix
!			otherwise 0.0 (usually use 1.0)
! 
!		Several parameters are used when the X or Y boundaries (or both) are periodic, otherwise set
!			all of these to zero.
!		If IACT,JACT are equal to the maximum compiled dimension (256x256) then it is best
!		to set XREFLECT to 1 if the matrix is periodic in the X direction, and the same with
!		YREFLECT for Y periodicity.  If IACT and/or JACT is less than the maximum dimensions, then
!		set XREFLECT and/or YREFLECT to zero and use DOBLENDX and/or DOBLENDY to unity.
!		This accomplishes more or less the same thing by doing a seletive averaging of generated
!		variates on the two sides of the periodic boundary.
!        set DOLEVELLOWER to unity if the values on the lower boundary
!        are all at the same value (e.g., for simulations with elevations set to have a level lower boundary), 
!           otherwise (and usually) equals zero
      READ(10,*) FACTOR, XREFLECT, YREFLECT,DOLEVELLOWER
!		The actual X,and Y dimensions of the output - the maximum values are 256x256
      READ(10,*) IACT,JACT
!        if this is <>0.0 then a regional trend of the generated
!           values is imposed whose slope is governed by REGIONALSLOPE      
      READ(10,*) REGIONALSLOPE
!        See above for DOBLENDX and DOBLENDY
!        if DODETREND=1 then the simulated values are linearly detrended
!           otherwise (and usually) set to zero
!
!		Set doexp to unity if you want the output to be a lognormal distribution (e.g., to
!			simulate rock erodibility or hydraulic conductivity, which are presumably >= 0)
!           otherwise 0 for a normal
!			distribution (allowing negative values)

      READ(10,*) DOBLENDX,DOBLENDY, DODETREND,DOEXP
      IF (DOEXP.GT.0) THEN 
             USEEXP = .TRUE. 
      ELSE 
          USEEXP = .FALSE.
      ENDIF
      READ(10,*) NO_BLENDS, BLEND_WIDTH
      ALLOCATE(LNEWELEVS(BLEND_WIDTH),RNEWELEVS(BLEND_WIDTH),TNEWELEVS(BLEND_WIDTH),BNEWELEVS(BLEND_WIDTH), STAT=IERROR)
      IST = (MAXDIMENSION - IACT) / 2+1
      IFN = IST+IACT-1
      JST = (MAXDIMENSION - JACT) / 2+1
      JFN = JST+JACT-1
      WRITE(*,100) IST,IFN,JST,JFN
100   FORMAT(' ISTART=',I5,' IEND=',I5,' JSTART=',I5,' JEND=',I5)
      CALL RANDOM_SEED(SIZE=NRANDOM)
      ALLOCATE(RANDSEED(NRANDOM))
      RANDSEED=RNSEED
      CALL RANDOM_SEED(PUT=RANDSEED)
      DEALLOCATE (RANDSEED)
      RETURN
      END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE DETREND
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      INTEGER I,J,III,NNN 
      real*4	INTR,SLP,YYYY,YPRED,SY,SSY,SX,SSX,SXY
              YYYY = 0.0
      	DO 100 I = IST,IFN
      		III = I-IST+1
      		YAVG(III) = 0.0
      		DO 100 J = JST,JFN
      			YAVG(III) = YAVG(III)+ELEVATION(I,J)
      		YAVG(III) = YAVG(III)/(JFN-JST+1)
                      YYYY = YYYY + YAVG(III)
100       CONTINUE
              YYYY = YYYY/(IFN-IST+1)
      	SY = 0.0
      	SSY = 0.0
      	SX = 0.0
      	SSX = 0.0
      	SXY = 0.0
      	DO 150 I = IST,IFN
      		III = I-IST+1
      		SY = SY + YAVG(III)
      		SSY = SSY + YAVG(III)*YAVG(III)
      		SX = SX +III
      		SSX = SSX + III*III
      		SXY = SXY + YAVG(III)*III
150       CONTINUE
      	NNN = IFN-IST+1
      	SLP = (SXY-(SX*SY)/NNN)/(SSX-SX*SX/NNN)
      	INTR= (SY-SLP*SX)/NNN
       	DO 160 I = IST,IFN
      		III = I-IST+1
      		YPRED = III*SLP+INTR
      		DO 160 J = JST,JFN
      			ELEVATION(I,J) = ELEVATION(I,J)-YPRED+YYYY
160       CONTINUE
      RETURN
      END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE BLENDX
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      INTEGER :: I,J,IL,IR,JL,JR 
      INTEGER :: ILOC,JLOC,II,JJ,IREL,JU
      real*4 ::   XST,XL,XR,XFN,WGT,AA,BB
      real*4 :: FILTSUM,WGGT
      LOGICAL :: USE_OLD
      USE_OLD=.FALSE.
      !real*4 LNEWELEVS(),RNEWELEVS(29)
      IF(USE_OLD) THEN
            XST = IST
            XFN = IFN
          AA = 0.5/29.0
          BB = 0.5 - AA
          DO 100 J = JST,JFN
          DO 110 I = 1,29
            IL=IST+I-1
            IR=IFN-I+1
            WGT = I*AA + BB
            LNEWELEVS(I) = WGT*ELEVATION(IL,J)+ &
                 (1.0-WGT)*ELEVATION(IR,J)
            RNEWELEVS(I) = WGT*ELEVATION(IR,J)+ &
                 (1.0-WGT)*ELEVATION(IL,J)
!      { 
!            XL = IL
!            XR = IR
!            LNEWELEVS(I) = ((XL-XST)/58.0+0.5)*ELEVATION(IL,J)+
!                            ((XR-XFN)/58.0+0.5)*ELEVATION(IR,J)
!            RNEWELEVS(I) = ((XST-XL)/58.0+0.5)*ELEVATION(IL,J)+
!                            ((XFN-XR)/58.0+0.5)*ELEVATION(IR,J)}
110       CONTINUE
          DO 120 I = 1,29
            IL = IST+I-1
            IR = IFN-I+1
            ELEVATION(IL,J) = LNEWELEVS(I)
            ELEVATION(IR,J) = RNEWELEVS(I)
120       CONTINUE
100       CONTINUE
      ELSE
          
          ENEW=0.0
          DO J=JST,JFN
          DO I=IST,IST+BLEND_WIDTH
              FILTSUM=0.0
              !ENEW=0.0
              DO JJ=-2,2
                  IF (((JJ+J).GE.JST).AND.((JJ+J).LE.JFN)) THEN
                    JLOC=J+JJ
                     DO II=-2,2
                        ILOC=I+II
                        IREL=ILOC-IST
                        IF (ILOC.LT.IST) ILOC=IFN+IREL+1
                        FILTSUM=FILTSUM+SMOOTH(II,JJ)
                        ENEW(I,J)=ENEW(I,J)+ELEVATION(ILOC,JLOC)*SMOOTH(II,JJ)
                    ENDDO
                  ENDIF
              ENDDO
              ENEW(I,J)=ENEW(I,J)/FILTSUM
          ENDDO
          DO I=IFN-BLEND_WIDTH,IFN
              FILTSUM=0.0
              !ENEW=0.0
              DO JJ=-2,2
                  IF (((JJ+J).GE.JST).AND.((JJ+J).LE.JFN)) THEN
                    JLOC=J+JJ
                     DO II=-2,2
                        ILOC=I+II
                        IREL=ILOC-IFN
                        IF (ILOC.GT.IFN) ILOC=IST+IREL-1
                        FILTSUM=FILTSUM+SMOOTH(II,JJ)
                        ENEW(I,J)=ENEW(I,J)+ELEVATION(ILOC,JLOC)*SMOOTH(II,JJ)
                    ENDDO
                  ENDIF
              ENDDO
              ENEW(I,J)=ENEW(I,J)/FILTSUM              
          ENDDO
          ENDDO
           AA=-1.0/(BLEND_WIDTH)
          BB=1.0-AA*IST
          DO J=JST,JFN
              DO I=IST,IST+BLEND_WIDTH
                  WGGT=I*AA+BB
                  ELEVATION(I,J)=ELEVATION(I,J)*(1.0-WGGT)+WGGT*ENEW(I,J)
              ENDDO
          ENDDO
          AA=1.0/BLEND_WIDTH
          BB=1.0-AA*JFN
          DO J=JST,JFN
              DO I=IFN-BLEND_WIDTH,IFN
                  WGGT=I*AA+BB
                  ELEVATION(I,J)=ELEVATION(I,J)*(1.0-WGGT)+WGGT*ENEW(I,J)
              ENDDO
          ENDDO
      ENDIF
      RETURN
      END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE BLENDY
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      INTEGER I,J,JT,JB,IT,IB
      INTEGER :: II,ILOC, JJ, JLOC, JREL
      real*4    YST,YT,YB,YFN,WGT,AA,BB
      real*4 :: FILTSUM,WGGT
      LOGICAL :: USE_OLD
      USE_OLD=.FALSE.
      IF (USE_OLD) THEN
      !real*4 TNEWELEVS(29),BNEWELEVS(29)
            YST = JST
            YFN = JFN
          AA = 0.5/29.0
          BB = 0.5 - AA
          DO 100 I = IST,IFN
          DO 110 J = 1,29
            JT=JST+J-1
            JB=JFN-J+1
            WGT = J*AA + BB
            TNEWELEVS(J) = WGT*ELEVATION(I,JT)+  &
                 (1.0-WGT)*ELEVATION(I,JB)
            BNEWELEVS(J) = WGT*ELEVATION(I,JB)+  &
                 (1.0-WGT)*ELEVATION(I,JT)
!      { 
!            YT = JT
!            YB = JB
!            TNEWELEVS(J) = ((YT-YST)/58.0+0.5)*ELEVATION(JT,J)+
!                            ((YB-YFN)/58.0+0.5)*ELEVATION(JB,J)
!            BNEWELEVS(J) = ((YST-YT)/58.0+0.5)*ELEVATION(JT,J)+
!                            ((YFN-YB)/58.0+0.5)*ELEVATION(JB,J)}
110       CONTINUE
          DO 120 J = 1,29
            JT = JST+J-1
            JB = JFN-J+1
            ELEVATION(I,JT) = TNEWELEVS(J)
            ELEVATION(I,JB) = BNEWELEVS(J)
120       CONTINUE
100       CONTINUE
      ELSE
          ENEW=0.0
          DO I=IST,IFN
          DO J=JST,JST+BLEND_WIDTH
              FILTSUM=0.0
              !ENEW=0.0
              DO II=-2,2
                  IF (((II+I).GE.IST).AND.((II+I).LE.IFN)) THEN
                    ILOC=I+II
                     DO JJ=-2,2
                        JLOC=J+JJ
                        JREL=JLOC-JST
                        IF (JLOC.LT.JST) JLOC=JFN+JREL+1
                        FILTSUM=FILTSUM+SMOOTH(II,JJ)
                        ENEW(I,J)=ENEW(I,J)+ELEVATION(ILOC,JLOC)*SMOOTH(II,JJ)
                    ENDDO
                  ENDIF
              ENDDO
              ENEW(I,J)=ENEW(I,J)/FILTSUM
          ENDDO
          DO J=JFN-BLEND_WIDTH,JFN
              FILTSUM=0.0
              !ENEW=0.0
              DO II=-2,2
                  IF (((II+I).GE.IST).AND.((II+I).LE.IFN)) THEN
                    ILOC=I+II
                     DO JJ=-2,2
                        JLOC=J+JJ
                        JREL=JLOC-JFN
                        IF (JLOC.GT.JFN) JLOC=JST+JREL-1
                        FILTSUM=FILTSUM+SMOOTH(II,JJ)
                        ENEW(I,J)=ENEW(I,J)+ELEVATION(ILOC,JLOC)*SMOOTH(II,JJ)
                    ENDDO
                  ENDIF
              ENDDO
              ENEW(I,J)=ENEW(I,J)/FILTSUM              
          ENDDO
          ENDDO
           AA=-1.0/(BLEND_WIDTH)
          BB=1.0-AA*IST
          DO I=IST,IFN
              DO J=JST,JST+BLEND_WIDTH
                  WGGT=J*AA+BB
                  ELEVATION(I,J)=ELEVATION(I,J)*(1.0-WGGT)+WGGT*ENEW(I,J)
              ENDDO
          ENDDO
          AA=1.0/BLEND_WIDTH
          BB=1.0-AA*IFN
          DO I=IST,IFN
              DO J=JFN-BLEND_WIDTH,JFN
                  WGGT=J*AA+BB
                  ELEVATION(I,J)=ELEVATION(I,J)*(1.0-WGGT)+WGGT*ENEW(I,J)
              ENDDO
          ENDDO
      ENDIF  
      RETURN
      END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE TILT
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      INTEGER I,J
      real*4	 LOCALELEV,AA,BB
      	AA=-REGIONALSLOPE
      	BB = 0.0-AA*JFN
      	DO 100 J = JST,JFN
      	DO 100 I = IST,IFN
      		ELEVATION(I,J) = ELEVATION(I,J)+AA*J+BB
              IF (ELEVATION(I,J).LT.0.0) ELEVATION(I,J) = 0.0
100       CONTINUE
          RETURN
      END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PROGRAM MATRIX3
      USE MATRIX_GLOBALS
      IMPLICIT NONE
	INTEGER I,J,K,IE,IW,JN,JS
	real*4 NORMAL2,MAXGRAD,RANGE
      real*4 RAN,TEMPX
      INTEGER ILL,IHH,JLL,JHH
      real*4 ELMIN,ELMAX
              J = 1
           OPEN(10,FILE='matrix_2d.prm')
           OPEN(11,FILE='matrix_2d.out')
           OPEN(12,FILE='matrix_2d.res')
           call SETUP()
           CALL GETPARAMETERS()
          ! CALL SETUP()
           DO 100 I = 1,MAXDIMENSION
           DO 100 J = 1,MAXDIMENSION
           
               DONE(I,J) = .FALSE.
               ACTIVE(I,J) = .FALSE.
               NEWACTIVE(I,J) = .FALSE.
100       CONTINUE
!	    write(12,805)
805       format('ready to go')
!           INTERVAL = MAXDIMENSION
           INTERVAL = (MAXDIMENSION-1)/2
           ACTIVE(1,1) = .TRUE.
           IF (USEEXP) THEN
              ELEVATION(1,1) = CONDUCTIVITY* &
                     EXP(FACTOR*NORMAL2(RDUMMY)*ERROR)
           ELSE
              ELEVATION(1,1) = CONDUCTIVITY +  &
                     FACTOR*NORMAL2(RDUMMY)*ERROR
           ENDIF
           ACTIVE(1,MAXDIMENSION) = .TRUE.
           IF (YREFLECT.GT.0) THEN
              ELEVATION(1,MAXDIMENSION)=ELEVATION(1,1)
           ELSE
             IF (DOLEVELLOWER.GT.0) THEN
                ELEVATION(1,MAXDIMENSION) = CONDUCTIVITY
             ELSE
               IF (USEEXP) THEN
                  ELEVATION(1,MAXDIMENSION) = CONDUCTIVITY*EXP(FACTOR* &
                  NORMAL2(RDUMMY)*ERROR)
               ELSE
                  ELEVATION(1,MAXDIMENSION) = CONDUCTIVITY + &
                     FACTOR*NORMAL2(RDUMMY)*ERROR
               ENDIF
             ENDIF
           ENDIF
           ACTIVE(MAXDIMENSION,1) = .TRUE.
           IF (XREFLECT.EQ.0) THEN
            IF (USEEXP) THEN
              ELEVATION(MAXDIMENSION,1) = CONDUCTIVITY*EXP( &
                  FACTOR*NORMAL2(RDUMMY)*ERROR)
            ELSE
              ELEVATION(MAXDIMENSION,1) = CONDUCTIVITY + &
                     FACTOR*NORMAL2(RDUMMY)*ERROR     
            ENDIF
           ELSE
               ELEVATION(MAXDIMENSION,1) = ELEVATION(1,1)
           ENDIF
           ACTIVE(MAXDIMENSION,MAXDIMENSION) = .TRUE.
           IF (YREFLECT.GT.0) THEN
              ELEVATION(MAXDIMENSION,   &
                     MAXDIMENSION)=ELEVATION(MAXDIMENSION,1)
           ELSE
             IF (DOLEVELLOWER.GT.0) THEN
                ELEVATION(MAXDIMENSION,MAXDIMENSION) = CONDUCTIVITY
             ELSE
               IF (XREFLECT.EQ.0) THEN
                 IF (USEEXP) THEN
                   ELEVATION(MAXDIMENSION,MAXDIMENSION) =  &
                    CONDUCTIVITY*EXP( FACTOR*NORMAL2(RDUMMY)*ERROR)
                 ELSE
                   ELEVATION(MAXDIMENSION,MAXDIMENSION) = CONDUCTIVITY  &
                     +FACTOR*NORMAL2(RDUMMY)*ERROR
                 ENDIF
               ELSE
                 ELEVATION(MAXDIMENSION,MAXDIMENSION) =  &
                     ELEVATION(1,MAXDIMENSION)
               ENDIF
             ENDIF
           ENDIF
           CURLEVEL = 1
200        CONTINUE
              MULFACTOR = ERROR*EXP(-CURLEVEL*CHARDIMENSION)
              WRITE(*,600) CURLEVEL, MULFACTOR, INTERVAL
600           FORMAT(' CURLEVEL= ',I5,' MULFACTOR =',G15.5, &
                         ' INTERVAL=',I5)
              J = 1
250           CONTINUE
      	      IF ((J.EQ.MAXDIMENSION).AND.(YREFLECT.GT.0)) THEN
      	        DO 310 I = 2,MAXDIMENSION
      	          ACTIVE(I,J) = ACTIVE(I,1)
      	          ELEVATION(I,J) = ELEVATION(I,1)
310               CONTINUE
      	      ELSE
                DO 260 I = 1,MAXDIMENSION
                  IF (ACTIVE(I,J)) THEN
                      INN = I + INTERVAL
                      IF ((INN .GT. 1).AND. &
                             (INN .LT. MAXDIMENSION)) THEN   
                          IF (.NOT.ACTIVE(INN,J)) THEN
                             ELAVERAGE = (ELEVATION(INN-INTERVAL,J) + &
                                  ELEVATION(INN+INTERVAL,J))/2.0
      		              IF ((DOLEVELLOWER.GT.0).AND.  &
                            (J.EQ.MAXDIMENSION)) THEN
      		                 ELEVATION(INN,J) = CONDUCTIVITY
      		              ELSE
      		                IF (USEEXP) THEN
                                 ELEVATION(INN,J) = ELAVERAGE*   &
                                 EXP(NORMAL2(RDUMMY)*MULFACTOR)
      		                ELSE
      		                  ELEVATION(INN,J) = ELAVERAGE+   &
             			               MULFACTOR*NORMAL2(RDUMMY)
                              ENDIF
                            ENDIF
                            NEWACTIVE(INN,J) = .TRUE.
                            DONE(INN,J) = .TRUE.
                          ENDIF
                      ENDIF
                  ENDIF
260             CONTINUE
                ENDIF
                DO 270 I = 1,MAXDIMENSION
                      IF (NEWACTIVE(I,J)) ACTIVE(I,J) = .TRUE.
                      NEWACTIVE(I,J) = .FALSE.
270             CONTINUE
                J = J + MAXDIMENSION-1
              IF (J.LE.MAXDIMENSION) GOTO 250 
              I = 1
300           CONTINUE
      	      IF ((I.EQ.MAXDIMENSION).AND.(XREFLECT.GT.0)) THEN
      	        DO 311 J = 2,MAXDIMENSION
      	          ACTIVE(I,J) = ACTIVE(1,J)
      	          ELEVATION(I,J) = ELEVATION(1,J)
311               CONTINUE
      	      ELSE
                  DO 320 J = 1,MAXDIMENSION
                    IF (ACTIVE(I,J)) THEN
                      JNN = J + INTERVAL
                      IF ((JNN .GT. 1).AND.  &
                             (JNN .LT. MAXDIMENSION)) THEN
                          IF (.NOT.ACTIVE(I,JNN)) THEN
                             ELAVERAGE = (ELEVATION(I,JNN-INTERVAL) +  &
                                  ELEVATION(I,JNN+INTERVAL))/2.0
      		               IF (USEEXP) THEN
                               ELEVATION(I,JNN) = ELAVERAGE  &
                                  *EXP(NORMAL2(RDUMMY)*MULFACTOR)
      		               ELSE
      		                 ELEVATION(I,JNN) = ELAVERAGE +  &
      			              NORMAL2(RDUMMY)*MULFACTOR
                             ENDIF
                             NEWACTIVE(I,JNN) = .TRUE.
                             DONE(I,JNN) = .TRUE.
                          ENDIF
                      ENDIF
                    ENDIF
320               CONTINUE
                  DO 330 J = 1,MAXDIMENSION
                      IF (NEWACTIVE(I,J)) ACTIVE(I,J) = .TRUE.
                      NEWACTIVE(I,J) = .FALSE.
330               CONTINUE
      	      ENDIF
              I = I + MAXDIMENSION-1
              IF (I.LE.MAXDIMENSION) GOTO 300
              CURLEVEL = CURLEVEL + 1
              INTERVAL = INTERVAL / 2
           IF (INTERVAL.GT.0) GOTO 200
           CURLEVEL = 1
           INTERVAL = (MAXDIMENSION-1) / 2
           DO 340 I = 1,MAXDIMENSION
           DO 340 J = 1,MAXDIMENSION
                NEWACTIVE(I,J) = .FALSE.
                ACTIVE(I,J) = .FALSE.
340        CONTINUE
           ACTIVE(1,1) = .TRUE.
           ACTIVE(1,MAXDIMENSION) = .TRUE.
           ACTIVE(MAXDIMENSION,1) = .TRUE.
           ACTIVE(MAXDIMENSION,MAXDIMENSION) = .TRUE.
           ACTIVE((INTERVAL+1),MAXDIMENSION) = .TRUE.
           ACTIVE(1,(INTERVAL+1)) = .TRUE.
           ACTIVE(MAXDIMENSION,(INTERVAL+1)) = .TRUE.
           ACTIVE((INTERVAL+1),(INTERVAL+1)) = .TRUE.
           ACTIVE((INTERVAL+1),1) = .TRUE.
           ELAVERAGE = (ELEVATION(1,1)+ELEVATION(1,MAXDIMENSION)+  &
               ELEVATION(MAXDIMENSION,1)+   &
                     ELEVATION(MAXDIMENSION,MAXDIMENSION))/4.0
           MULFACTOR = ERROR*EXP(-CURLEVEL*CHARDIMENSION)
!      {     WRITE(*,)CURLEVEL= ',CURLEVEL,' MULFACTOR.EQ.',MULFACTOR:15:5)}
           IF (USEEXP) THEN
              ELEVATION((INTERVAL+1),(INTERVAL+1)) = ELAVERAGE &
                  *EXP(NORMAL2(RDUMMY)*MULFACTOR)
           ELSE
              ELEVATION((INTERVAL+1),(INTERVAL+1)) = ELAVERAGE +  &
      	    NORMAL2(RDUMMY)*MULFACTOR
           ENDIF
           DONE((INTERVAL+1),(INTERVAL+1)) = .TRUE.
           DO 370 J = 1,MAXDIMENSION
               DONE(1,J) = .TRUE.
               DONE(MAXDIMENSION,J) = .TRUE.
370        CONTINUE
           DO 360 I = 1,MAXDIMENSION
               DONE(I,1) = .TRUE.
               DONE(I,MAXDIMENSION) = .TRUE.
360        CONTINUE
           CURLEVEL = 2
           INTERVAL = INTERVAL / 2
350        CONTINUE
              MULFACTOR = ERROR*EXP(-CURLEVEL*CHARDIMENSION)
              DO 380 I = 1,MAXDIMENSION
              DO 380 J = 1,MAXDIMENSION
                 IF (ACTIVE(I,J)) THEN
                    DO 390 K = 1,4
                      INN = I + APLACE(1,K)*INTERVAL
                      JNN = J + APLACE(2,K)*INTERVAL
                      IF ((INN .GT.1).AND.(INN .LT. MAXDIMENSION)) THEN
                      IF ((JNN.GT.1).AND.(JNN .LT. MAXDIMENSION)) THEN
                          IF (.NOT.DONE(INN,JNN)) THEN
                             ELAVERAGE = (ELEVATION(INN+INTERVAL, &
                                 JNN+INTERVAL)    &
                               +ELEVATION(INN+INTERVAL,JNN-INTERVAL)  &
                               +ELEVATION(INN-INTERVAL,JNN+INTERVAL)  &
                               +ELEVATION(INN-INTERVAL,JNN-INTERVAL)  &
                               )/4.0
      		               IF (USEEXP) THEN
                               ELEVATION(INN,JNN) = &
                              ELAVERAGE*EXP(NORMAL2(RDUMMY)*MULFACTOR)
      			           ELSE
      			             ELEVATION(INN,JNN) =  &
             			             ELAVERAGE+NORMAL2(RDUMMY)*  &
                                         MULFACTOR
                             ENDIF
                                NEWACTIVE(INN,JNN) = .TRUE.
                                DONE(INN,JNN) = .TRUE.
                          ENDIF
                      ENDIF
                      ENDIF
390                 CONTINUE
                 ENDIF
380           CONTINUE
              DO 410 I = 1,MAXDIMENSION
              DO 410 J = 1,MAXDIMENSION
                 IF (ACTIVE(I,J)) THEN
                   DO 420 K = 1,4
                     INN = I + BPLACE(1,K)*INTERVAL
                     JNN = J + BPLACE(2,K)*INTERVAL
                     IF ((INN .GT.1).AND.(INN .LT. MAXDIMENSION)) THEN
                     IF ((JNN.GT.1).AND.(JNN .LT. MAXDIMENSION)) THEN
                        IF (.NOT.DONE(INN,JNN)) THEN
                            ELAVERAGE = (ELEVATION(INN+INTERVAL,JNN)  &
                                 +ELEVATION(INN-INTERVAL,JNN) &
                                 +ELEVATION(INN,JNN+INTERVAL) &
                                 +ELEVATION(INN,JNN-INTERVAL) &
                                 )/4.0
      			          IF (USEEXP) THEN   
                               ELEVATION(INN,JNN) =  &
                              ELAVERAGE*EXP(NORMAL2(RDUMMY)*MULFACTOR)
      			          ELSE
                              ELEVATION(INN,JNN) = ELAVERAGE +  &
                                     MULFACTOR*NORMAL2(RDUMMY)
                            ENDIF
                            NEWACTIVE(INN,JNN) = .TRUE.
                            DONE(INN,JNN) = .TRUE.
                        ENDIF
                     ENDIF
                     ENDIF
420                CONTINUE
                 ENDIF
410           CONTINUE
              DO 430 I = 1,MAXDIMENSION
              DO 430 J = 1,MAXDIMENSION
                  IF (NEWACTIVE(I,J)) ACTIVE(I,J) = .TRUE.
                  NEWACTIVE(I,J) = .FALSE.
430           CONTINUE
           CURLEVEL = CURLEVEL + 1               
           INTERVAL = INTERVAL / 2
           IF (INTERVAL.GT.0) GOTO 350
!           WRITE (12,538) ELEVATION(1,1),ELEVATION(1,MAXDIMENSION),
!     +            ELEVATION(MAXDIMENSION,1),ELEVATION(MAXDIMENSION,
!     +            MAXDIMENSION)
538   FORMAT(' corners ',4(G12.5,' '))
           DO 537 I=1,MAXDIMENSION
              TEMPX=ELEVATION(I,1)-ELEVATION(I,MAXDIMENSION)
!              WRITE(12,536) TEMPX
536   FORMAT(' tb ',G12.5)
537        CONTINUE
           DO 535 J=1,MAXDIMENSION
              TEMPX=ELEVATION(1,J)-ELEVATION(MAXDIMENSION,J)
!              WRITE(12,534) TEMPX
534   format(' lr ',G12.5) 
535        CONTINUE
           IF (DOBLENDX .GT.0) THEN
               DO I=1,NO_BLENDS
                   CALL BLENDX()
               ENDDO
           ENDIF
           IF (DOBLENDY .GT.0) THEN
               DO I=1,NO_BLENDS
                   CALL BLENDY()
               ENDDO
           ENDIF
           IF (DODETREND .GT.0) CALL DETREND()
           IF (REGIONALSLOPE.GT.0.0) CALL TILT()
      ELMIN=1.0E+25
      ELMAX=-ELMIN
      DO 422 I=1,MAXDIMENSION
      DO 422 J=1,MAXDIMENSION
	 IF (ELEVATION(I,J).LT.ELMIN) THEN
           ELMIN=ELEVATION(I,J)
           ILL=I
           JLL=J
         ENDIF
         IF (ELEVATION(I,J).GT.ELMAX) THEN
           ELMAX=ELEVATION(I,J)
           IHH=I
           JHH=J
         ENDIF
422   CONTINUE
      WRITE(*,423) ILL,JLL,IHH,JHH
423   FORMAT(' LOWEST POINT AT I=',I5, ' AND J=',I5,/,  &
         ' HIGHEST POINT AT I=',I5,' AND J=',I5)
!     HERE IS WHERE THE OUTPUT FILE IS WRITTEN -- IF YOU
!     WANT IT IN A DIFFERENT FORMAT, MODIFY THIS SECTION
           WRITE(11,610) IACT, JACT
610        FORMAT(2I6)
           AVGVAL = 0.0
           TNUM = 0
           MAXVAL = -1.0E+10
           MINVAL = 1.0E+10
           MAXGRAD=0.0
           DO 500 I = IST,IFN
              if (dolevellower.gt.0) then
                 elevation(i,jfn)=conductivity
              endif
               DO 500 J = JST,JFN
                    IE=I+1
                    IF (IE.GT.IFN) THEN
                      IF (XREFLECT.GT.0) THEN
                          IE=IST
                      ELSE
                          IE=I
                      ENDIF
                    ENDIF
                    IW=I-1
                    IF (IW.LT.IST) THEN
                      IF (XREFLECT.GT.0) THEN
                          IW=IFN
                      ELSE
                          IW=I
                      ENDIF
                    ENDIF
                    JS=J+1
                    IF (JS.GT.JFN) THEN
                      IF (YREFLECT.GT.0) THEN
                          JS=JST
                      ELSE
                          JS=J
                      ENDIF
                    ENDIF
                    JN=J-1
                    IF (JN.LT.JST) THEN
                      IF (YREFLECT.GT.0) THEN
                          JN=JFN
                      ELSE
                          JN=J
                      ENDIF
                    ENDIF
                    WRITE(11,620) ELEVATION(I,J)
620                 FORMAT(G12.5)
                    TNUM = TNUM+1
                    TEMPX=abs(ELEVATION(I,J)-ELEVATION(I,JS))
                    IF (TEMPX.GT.MAXGRAD) MAXGRAD=TEMPX
                    TEMPX=abs(ELEVATION(I,J)-ELEVATION(I,JN))
                    IF (TEMPX.GT.MAXGRAD) MAXGRAD=TEMPX
                    TEMPX=abs(ELEVATION(I,J)-ELEVATION(IE,J))
                    IF (TEMPX.GT.MAXGRAD) MAXGRAD=TEMPX
                    TEMPX=abs(ELEVATION(I,J)-ELEVATION(IW,J))
                    IF (TEMPX.GT.MAXGRAD) MAXGRAD=TEMPX
                    TEMPX=abs(ELEVATION(I,J)-ELEVATION(IE,JS))
                    IF (TEMPX.GT.MAXGRAD) MAXGRAD=TEMPX
                    TEMPX=abs(ELEVATION(I,J)-ELEVATION(IE,JN))
                    IF (TEMPX.GT.MAXGRAD) MAXGRAD=TEMPX
                    TEMPX=abs(ELEVATION(I,J)-ELEVATION(IW,JN))
                    IF (TEMPX.GT.MAXGRAD) MAXGRAD=TEMPX
                    TEMPX=abs(ELEVATION(I,J)-ELEVATION(IW,JS))
                    IF (TEMPX.GT.MAXGRAD) MAXGRAD=TEMPX
                    AVGVAL = AVGVAL+ELEVATION(I,J)
                    IF (ELEVATION(I,J).GT.MAXVAL)  &
                         MAXVAL =ELEVATION(I,J)
                    IF (ELEVATION(I,J).LT.MINVAL)  &
                         MINVAL =ELEVATION(I,J)
500       CONTINUE
          AVGVAL=AVGVAL/TNUM
          RANGE=MAXVAL-MINVAL
          WRITE(12,630) AVGVAL,MINVAL,MAXVAL,RANGE,MAXGRAD
630       FORMAT(' AVERAGE=',G12.5,' MIN=',G12.5,' MAX=',G12.5,/   &
                 ' RANGE=',G12.5,' MAXGRAD=',G12.5)
          close(11)
	    call topo()
          CALL SHADE()
!          CALL TOPO()
 !             CLOSE(11)
              CLOSE(12)
      DEALLOCATE(ELEVATION,ACTIVE,NEWACTIVE,DONE)
      STOP
      END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE TOPO()
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      real*4 RANGE
      INTEGER I,J,IMX,JMX
      BYTE BVAL
      INTEGER ITEMP
      CHARACTER ACHAR
      IMX=iact
      JMX=jact
      MAXVAL=-1.0E+25
      MINVAL=-MAXVAL
      DO 100 I=1,IMX
      DO 100 J=1,JMX
          IF (ELEVATION(I,J).GT.MAXVAL) MAXVAL=ELEVATION(I,J)
          IF (ELEVATION(I,J).LT.MINVAL) MINVAL=ELEVATION(I,J)
100   CONTINUE
      RANGE=MAXVAL-MINVAL
      OPEN(45,NAME='m3topo.raw')
      DO 120 J=1,JMX
      DO 120 I=1,IMX
      ITEMP=INT((ELEVATION(I,J)-MINVAL)*253.0/RANGE)+1
      IF (ITEMP.LT.0) ITEMP=0
      IF (ITEMP.GT.255) ITEMP=255
                    BVAL=ITEMP
                    ACHAR=CHAR(BVAL)
                    WRITE(44,130) ACHAR
      WRITE(45,130) ACHAR
130   FORMAT(A1,$)
120   CONTINUE
      CLOSE(45)
      RETURN
      END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE SHADE()
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      real*4 HALFPI,PI,RADIAN,NOMANGLE,DXX
	real*4 ETEMP(MAXDIMENSION,MAXDIMENSION)
      real*4 YDIV,SLX,SLY,SLOPE,SLMAX,SLSUN,ANGLE,SLOPEANGLE,NEWVALUE
      INTEGER I,J,IW,IE,JS,JN, SELECT,SLOPENUMB,IMX,JMX
      real*4	 NOMSLMAX,RRRR,RANGE,SLOPEAVG,XDIV
      INTEGER II,JJ,IISTART,IIEND,JJSTART,JJEND
      real*4 XVAL1,XVAL2,XVAL3,XVAL4,XAVG
      BYTE BVAL
	INTEGER ITEMP
	CHARACTER ACHAR
	IMX=iact-1
	JMX=jact-1
!
! ******  MAKES A SHADED RELIEF RAW IMAGE OF THE SIMULATED TERRAIN AND
!         WRITES IT OUT TO A FILE
!
	HALFPI =1.5708
      PI=3.141593
      RADIAN=0.01745329
      NOMANGLE=25.0
          SLMAX=-1.0E+25
          SLOPEAVG=0.0
          SLOPENUMB= 0
          DO 101 I = 2,IMX-1 
          DO 101 J = 2,JMX-1 
             IW=I-1
             IE=I+1
             JS=J+1
             JN=J-1
             SLX = ( ELEVATION(IE,J)- ELEVATION(IW,J))/2.0
             SLY = ( ELEVATION(I,JN)- ELEVATION(I,JS))/2.0
             SLOPE =SQRT(SLX**2+SLY**2)
             IF (SLOPE.GT.SLMAX) SLMAX = SLOPE
             SLOPENUMB = SLOPENUMB+1
             SLOPEAVG = SLOPEAVG+SLOPE
101        CONTINUE
      	 NOMSLMAX = SIN(HALFPI-NOMANGLE*RADIAN)/  &
             COS(HALFPI-NOMANGLE*RADIAN)
      	 DXX = SLMAX/NOMSLMAX
      	 SLMAX = NOMSLMAX
          SLOPEAVG = SLOPEAVG/(SLOPENUMB*DXX)
          SLSUN=HALFPI-atan(SLMAX)
          DO 102 I = 1,IMX 
          DO 102 J = 1,JMX 
             IW=I-1
             IE=I+1
             JS=J+1
             JN=J-1
!             IF (IW .LT.1) IW=JMX
!             IF (IE.GT.JMX) IE=1
             YDIV=2.0*DXX
             XDIV=2.0*DXX										             
		   IF (JS.GT.JMX) THEN
               YDIV=1.0*DXX
               JS=JMX
             ENDIF
             IF (JN.LT.1) THEN
               YDIV=1.0*DXX
               JN=1
             ENDIF
             IF (IW.LT.1) THEN
                XDIV=1.0*DXX
                IW=1
             ENDIF
             IF (IE.GT.IMX) THEN
                XDIV=1.0*DXX
                IE=IMX
             ENDIF
             SLX = ( ELEVATION(IE,J)- ELEVATION(IW,J))/XDIV
             SLY = ( ELEVATION(I,JS)- ELEVATION(I,JN))/YDIV
             SLOPE =sqrt(SLX**2+SLY**2)
             SLOPEANGLE= atan(SLOPE)
             IF ((SLY.EQ.0.0).OR. (SLX.EQ.0.0)) THEN
               IF (SLY.EQ.0.0) THEN
                 IF (SLX.GT.0.0) THEN
                     ANGLE =  HALFPI
                 ELSE
                     ANGLE = -HALFPI
                 ENDIF
               ELSE
                   IF (SLY.GT.0.0) THEN
                     ANGLE = 0.0
                   ELSE
                     ANGLE = PI
                   ENDIF
               ENDIF
             ELSE
                 ANGLE = atan ((SLX) / (SLY))
                 IF (((SLY) .LT. 0.0).AND.((SLX) .LT. 0.0)) THEN
                   ANGLE = ANGLE + PI
                 ELSE
                   IF ((SLY) .LT. 0.0) THEN
                      ANGLE = ANGLE + PI
                   ENDIF
                 ENDIF
              ENDIF
               ANGLE = ANGLE-HALFPI
               NEWVALUE =       &
                  (COS(SLSUN)*COS(SLOPEANGLE)+SIN(SLSUN)  &
                         *SIN(SLOPEANGLE)*   &
                      COS(ANGLE))
               IF (NEWVALUE.LT.0.0) NEWVALUE = 0.0
               ETEMP(I,J) = NEWVALUE
102       CONTINUE
         MAXVAL = -1.0E+25
         MINVAL = -MAXVAL
         DO 100 I = 1,IMX 
         DO 100 J = 1,JMX 
           RRRR = ETEMP(I,J)
           IF (RRRR.LT.MINVAL) MINVAL = RRRR
           IF (RRRR.GT. MAXVAL) MAXVAL = RRRR
100      CONTINUE
      RANGE=MAXVAL-MINVAL
      OPEN(44,NAME='m2dshade.raw')
      IF (YREFLECT.GT.0) THEN
          JJSTART=JMX-JMX/4
          JJEND=JMX
      DO 1120 J=JJSTART,JJEND
          DO 1121 JJ=1,2
             IF (XREFLECT.GT.0) THEN
                IISTART=IMX-IMX/4
                IIEND=IMX
             DO 1123 I=IISTART,IIEND
                 XVAL1=ETEMP(I,J)
                 IF (I.EQ.IIEND) THEN
                  XVAL2=ETEMP(1,J)
                  IF (J.EQ.JJEND) THEN
                      XVAL4=ETEMP(1,1)
                  ELSE
                      XVAL4=ETEMP(1,J+1)
                  ENDIF
                 ELSE
                  XVAL2=ETEMP(I+1,J)
                  IF (J.EQ.JJEND) THEN
                      XVAL4=ETEMP(I+1,1)
                  ELSE
                      XVAL4=ETEMP(I+1,J+1)
                  ENDIF
                 ENDIF
                 IF (J.EQ.JJEND) THEN
                  XVAL3=ETEMP(I,1)
                 ELSE
                  XVAL3=ETEMP(I,J+1)
                 ENDIF
                 DO 1124 II=1,2
                    XAVG=XVAL1*WEIGHTS(II,JJ,1)+  &
                         XVAL2*WEIGHTS(II,JJ,2)+ &
                         XVAL3*WEIGHTS(II,JJ,3)+  &
                         XVAL4*WEIGHTS(II,JJ,4)
                    ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                    IF (ITEMP.LT.0) ITEMP=0
                    IF (ITEMP.GT.255) ITEMP=255
                    BVAL=ITEMP
                    ACHAR=char(BVAL)
                    WRITE(44,130) ACHAR
1124              CONTINUE
1123          CONTINUE
             ENDIF
             DO 1125 I=1,IMX-1
                 XVAL1=ETEMP(I,J)
                  XVAL2=ETEMP(I+1,J)
                  IF (J.EQ.JJEND) THEN
                      XVAL4=ETEMP(I+1,1)
                      XVAL3=ETEMP(I,1)
                  ELSE
                      XVAL4=ETEMP(I+1,J+1)
                      XVAL3=ETEMP(I,J+1)
                  ENDIF
                 DO 1126 II=1,2
                    XAVG=XVAL1*WEIGHTS(II,JJ,1)+   &
                         XVAL2*WEIGHTS(II,JJ,2)+ &
                         XVAL3*WEIGHTS(II,JJ,3)+  &
                         XVAL4*WEIGHTS(II,JJ,4)
                    ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                    IF (ITEMP.LT.0) ITEMP=0
                    IF (ITEMP.GT.255) ITEMP=255
                    BVAL=ITEMP
                    ACHAR=char(BVAL)
                    WRITE(44,130) ACHAR
1126              CONTINUE
1125          CONTINUE
          IF (XREFLECT.GT.0) THEN
              IISTART=0
              IIEND=IMX/4
             DO 1127 I=IISTART,IIEND
                 XVAL2=ETEMP(I+1,J)
                 IF (I.EQ.IISTART) THEN
                  XVAL1=ETEMP(IMX,J)
                  IF (J.EQ.JJEND) THEN
                      XVAL3=ETEMP(IMX,1)
                  ELSE
                      XVAL3=ETEMP(IMX,J+1)
                  ENDIF
                 ELSE
                  XVAL1=ETEMP(I+1,J)
                  IF (J.EQ.JJEND) THEN
                      XVAL3=ETEMP(I,1)
                  ELSE
                      XVAL3=ETEMP(I,J+1)
                  ENDIF
                 ENDIF
                 IF (J.EQ.JJEND) THEN
                  XVAL4=ETEMP(I+1,1)
                 ELSE
                  XVAL4=ETEMP(I+1,J+1)
                 ENDIF
                 DO 1128 II=1,2
                    XAVG=XVAL1*WEIGHTS(II,JJ,1)+ &
                         XVAL2*WEIGHTS(II,JJ,2)+  &
                         XVAL3*WEIGHTS(II,JJ,3)+  &
                         XVAL4*WEIGHTS(II,JJ,4)
                    ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                    IF (ITEMP.LT.0) ITEMP=0
                    IF (ITEMP.GT.255) ITEMP=255
                    BVAL=ITEMP
                    ACHAR=char(BVAL)
                    WRITE(44,130) ACHAR
1128              CONTINUE
1127          CONTINUE
          ENDIF
1121       CONTINUE
1120   CONTINUE
      ENDIF
      DO 120 J=1,JMX-1
          DO 121 JJ=1,2
             IF (XREFLECT.GT.0) THEN
                IISTART=IMX-IMX/4
                IIEND=IMX
             DO 123 I=IISTART,IIEND
                 XVAL1=ETEMP(I,J)
                 IF (I.EQ.IIEND) THEN
                  XVAL2=ETEMP(1,J)
                  XVAL4=ETEMP(1,J+1)
                 ELSE
                  XVAL2=ETEMP(I+1,J)
                  XVAL4=ETEMP(I+1,J+1)
                 ENDIF
                 XVAL3=ETEMP(I,J+1)
                 DO 124 II=1,2
                    XAVG=XVAL1*WEIGHTS(II,JJ,1)+   &
                         XVAL2*WEIGHTS(II,JJ,2)+   &
                         XVAL3*WEIGHTS(II,JJ,3)+   &
                         XVAL4*WEIGHTS(II,JJ,4)
                    ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                    IF (ITEMP.LT.0) ITEMP=0
                    IF (ITEMP.GT.255) ITEMP=255
                    BVAL=ITEMP
                    ACHAR=char(BVAL)
                    WRITE(44,130) ACHAR
124              CONTINUE
123          CONTINUE
             ENDIF
             DO 125 I=1,IMX-1
                 XVAL1=ETEMP(I,J)
                  XVAL2=ETEMP(I+1,J)
                  XVAL4=ETEMP(I+1,J+1)
                 XVAL3=ETEMP(I,J+1)
                 DO 126 II=1,2
                    XAVG=XVAL1*WEIGHTS(II,JJ,1)+  &
                         XVAL2*WEIGHTS(II,JJ,2)+ &
                         XVAL3*WEIGHTS(II,JJ,3)+  &
                         XVAL4*WEIGHTS(II,JJ,4)
                    ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                    IF (ITEMP.LT.0) ITEMP=0
                    IF (ITEMP.GT.255) ITEMP=255
                    BVAL=ITEMP
                    ACHAR=char(BVAL)
                    WRITE(44,130) ACHAR
126              CONTINUE
125          CONTINUE
          IF (XREFLECT.GT.0) THEN
              IISTART=0
              IIEND=IMX/4
             DO 127 I=IISTART,IIEND
                 XVAL2=ETEMP(I+1,J)
                 IF (I.EQ.IISTART) THEN
                  XVAL1=ETEMP(IMX,J)
                  XVAL3=ETEMP(IMX,J+1)
                 ELSE
                  XVAL1=ETEMP(I+1,J)
                  XVAL3=ETEMP(I+1,J+1)
                 ENDIF
                 XVAL4=ETEMP(I+1,J+1)
                 DO 128 II=1,2
                    XAVG=XVAL1*WEIGHTS(II,JJ,1)+   &
                         XVAL2*WEIGHTS(II,JJ,2)+  &
                         XVAL3*WEIGHTS(II,JJ,3)+ &
                         XVAL4*WEIGHTS(II,JJ,4)
                    ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                    IF (ITEMP.LT.0) ITEMP=0
                    IF (ITEMP.GT.255) ITEMP=255
                    BVAL=ITEMP
                    ACHAR=char(BVAL)
                    WRITE(44,130) ACHAR
128              CONTINUE
127          CONTINUE
          ENDIF
121       CONTINUE
120   CONTINUE
      IF (YREFLECT.GT.0) THEN
          JJSTART=0
          JJEND=JMX/4
      DO 2120 J=JJSTART,JJEND
          DO 2121 JJ=1,2
             IF (XREFLECT.GT.0) THEN
                IISTART=IMX-IMX/4
                IIEND=IMX
             DO 2123 I=IISTART,IIEND
                 IF (J.EQ.JJSTART) THEN
                  XVAL1=ETEMP(I,JMX)
                 ELSE
                  XVAL1=ETEMP(I,J)
                 ENDIF
                 IF (I.EQ.IIEND) THEN
                  IF (J.EQ.JJSTART) THEN
                      XVAL2=ETEMP(1,JMX)
                  ELSE
                      XVAL2=ETEMP(1,J)
                  ENDIF
                  XVAL4=ETEMP(1,J+1)
                 ELSE
                  IF (J.EQ.JJSTART) THEN
                      XVAL2=ETEMP(I+1,JMX)
                  ELSE
                      XVAL2=ETEMP(I+1,J)
                  ENDIF
                  XVAL4=ETEMP(I+1,J+1)
                 ENDIF
                  XVAL3=ETEMP(I,J+1)
                 DO 2124 II=1,2
                    XAVG=XVAL1*WEIGHTS(II,JJ,1)+  &
                         XVAL2*WEIGHTS(II,JJ,2)+  &
                         XVAL3*WEIGHTS(II,JJ,3)+  &
                         XVAL4*WEIGHTS(II,JJ,4)
                    ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                    IF (ITEMP.LT.0) ITEMP=0
                    IF (ITEMP.GT.255) ITEMP=255
                    BVAL=ITEMP
                    ACHAR=char(BVAL)
                    WRITE(44,130) ACHAR
2124              CONTINUE
2123          CONTINUE
             ENDIF
             DO 2125 I=1,IMX-1
                  IF (J.EQ.JJSTART) THEN
                     XVAL1=ETEMP(I,JMX)
                     XVAL2=ETEMP(I+1,JMX)
                  ELSE
                     XVAL1=ETEMP(I,J)
                     XVAL2=ETEMP(I+1,J)
                  ENDIF
                  XVAL4=ETEMP(I+1,J+1)
                  XVAL3=ETEMP(I,J+1)
                 DO 2126 II=1,2
                    XAVG=XVAL1*WEIGHTS(II,JJ,1)+   &
                         XVAL2*WEIGHTS(II,JJ,2)+  &
                         XVAL3*WEIGHTS(II,JJ,3)+   &
                         XVAL4*WEIGHTS(II,JJ,4)
                    ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                    IF (ITEMP.LT.0) ITEMP=0
                    IF (ITEMP.GT.255) ITEMP=255
                    BVAL=ITEMP
                    ACHAR=char(BVAL)
                    WRITE(44,130) ACHAR
2126              CONTINUE
2125          CONTINUE
          IF (XREFLECT.GT.0) THEN
              IISTART=0
              IIEND=IMX/4
             DO 2127 I=IISTART,IIEND
                 IF (J.EQ.JJSTART) THEN
                  XVAL2=ETEMP(I+1,JMX)
                 ELSE
                  XVAL2=ETEMP(I+1,J)
                 ENDIF
                 IF (I.EQ.IISTART) THEN
                  IF (J.EQ.JJSTART) THEN
                      XVAL1=ETEMP(IMX,JMX)
                  ELSE
                      XVAL1=ETEMP(IMX,J)
                  ENDIF
                      XVAL3=ETEMP(IMX,J+1)
                 ELSE
                  IF (J.EQ.JJSTART) THEN
                      XVAL1=ETEMP(I,JMX)
                  ELSE
                      XVAL1=ETEMP(I,J)
                  ENDIF
                  XVAL3=ETEMP(I,J+1)
                 ENDIF
                 XVAL4=ETEMP(I+1,J+1)
                 DO 2128 II=1,2
                    XAVG=XVAL1*WEIGHTS(II,JJ,1)+  &
                         XVAL2*WEIGHTS(II,JJ,2)+ &
                         XVAL3*WEIGHTS(II,JJ,3)+  &
                         XVAL4*WEIGHTS(II,JJ,4)
                    ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                    IF (ITEMP.LT.0) ITEMP=0
                    IF (ITEMP.GT.255) ITEMP=255
                    BVAL=ITEMP
                    ACHAR=char(BVAL)
                    WRITE(44,130) ACHAR
2128              CONTINUE
2127          CONTINUE
          ENDIF
2121       CONTINUE
2120   CONTINUE
      ENDIF
130   FORMAT(A1,$)
      CLOSE(44)
      RETURN
      END
