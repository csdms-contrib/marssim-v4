      MODULE MATRIX_GLOBALS
      IMPLICIT NONE
      SAVE
      INTEGER :: MAXDIMENSION,MAXLEVEL
      REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: ELEVATION
      LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: ACTIVE
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: NEWACTIVE, DONE
      LOGICAL :: USEBLENDX,USEBLENDY,USEPERIODX,USEPERIODY,USEEXP
      REAL*4 :: EFACTOR, EDECAY, FACTOR, CONDUCTIVITY, ERROR, CHARDIMENSION
      REAL*4 :: MULFACTOR, ELAVERAGE, ELNUM
      REAL*4 :: XTOT, XSUM, XSSS, XOVER, XSTAVG, YTOT, YOVER
      REAL*4 :: YSUM, YSSS
      INTEGER :: CURLEVEL, INTERVAL, PARAMS, OUTFILE, DEBUG
      INTEGER :: IACT, JACT, KACT, IST, IFN, JST, JFN, KST, KFN
      INTEGER :: INN, JNN, KNN, KK, JJ, II,LL, IMM, JMM, KMM, RDUMMY
      INTEGER, DIMENSION(2,4) :: APLACE, BPLACE
      INTEGER, DIMENSION(3,26) :: CPLACE, DPLACE
      END MODULE MATRIX_GLOBALS
! *******************************************************************
! *******************************************************************     
      PROGRAM MATRIX_3D
!        This is a program to produce a 3-D pseudo-fractal
!        "cube".  See the subroutine GETPARAMETERS for info
!        on the simulation parameters
!        The results are output as a direct access file so
!        that the large resulting file does not have to
!        reside in memory.  For an example of the use of
!        such a cube to simulate rock resistance variability
!        see the subroutines FINDERODIBILITY and READERODIBILITY
!        in the file BOUNDARY.F
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      integer iisz,jjsz,kksz,irec,i,j,k
!	integer debug
       real*4 xxxx,normal2,XAVG
      params=20
      outfile=21
	   debug=22 
      open(PARAMS,file='matrix_3d.prm') 
      OPEN(outfile,FILE='matrix_3d.out', &
      ACCESS='DIRECT',FORM='unformatted',RECL=4)
	   open(debug,file='matrix_3d.dbg')
      call GETPARAMETERS() 
      call setup()
	  YOVER = 0.0 
	  YTOT = 0.0 
	  YSUM = 0.0 
	  YSSS = 0.0 
      DO I = 0 , MAXDIMENSION 
      DO J = 0 , MAXDIMENSION
	  DO K = 0 , MAXDIMENSION
         ACTIVE(I,J,K) = .FALSE. 
		ELEVATION(I,J,K) = -1.0E+25 
      ENDdo
      enddo
      enddo 
      INTERVAL = MAXDIMENSION /  2
       
      ACTIVE(0,0,0) = .TRUE. 
      ELEVATION(0,0,0) = NORMAL2(RDUMMY)*FACTOR*ERROR
       
      ACTIVE(0,0,MAXDIMENSION) = .TRUE. 
      ELEVATION(0,0,MAXDIMENSION) = NORMAL2(RDUMMY)*FACTOR*ERROR
       
      ACTIVE(0,MAXDIMENSION,0) = .TRUE.
      if (useperiody) then
         elevation(0,maxdimension,0)=elevation(0,0,0)
      else
         ELEVATION(0,MAXDIMENSION,0) = NORMAL2(RDUMMY)*FACTOR*ERROR
      endif
       
      ACTIVE(0,MAXDIMENSION,MAXDIMENSION) = .TRUE.
      if (useperiody) then
      ELEVATION(0,MAXDIMENSION,MAXDIMENSION) =  &
            elevation(0,0,maxdimension)
      else
      ELEVATION(0,MAXDIMENSION,MAXDIMENSION) =  &
            NORMAL2(RDUMMY)*FACTOR*ERROR
      endif
       
       ACTIVE(MAXDIMENSION,0,0) = .TRUE. 
      IF (USEPERIODX) THEN
         ELEVATION(MAXDIMENSION,0,0) =ELEVATION(0,0,0)
      ELSE
         ELEVATION(MAXDIMENSION,0,0) =  &
            NORMAL2(RDUMMY)*ERROR*FACTOR
      endif
       
      ACTIVE(MAXDIMENSION,0,MAXDIMENSION) = .TRUE. 
      IF (USEPERIODX) THEN
       ELEVATION(MAXDIMENSION,0,MAXDIMENSION) =  &
         ELEVATION(0,0,MAXDIMENSION)
      ELSE
       ELEVATION(MAXDIMENSION,0,MAXDIMENSION) = &
            NORMAL2(RDUMMY)*ERROR*FACTOR
      endif
      
     
      ACTIVE(MAXDIMENSION,MAXDIMENSION,0) = .TRUE. 
      IF (USEPERIODX) THEN
        ELEVATION(MAXDIMENSION,MAXDIMENSION,0) =  &
          ELEVATION(0,MAXDIMENSION,0)
      ELSE
         if(useperiody) then
            ELEVATION(MAXDIMENSION,MAXDIMENSION,0) = &
           elevation(maxdimension,0,0)
         else
           ELEVATION(MAXDIMENSION,MAXDIMENSION,0) = &
           NORMAL2(RDUMMY)*FACTOR*ERROR
         endif
      endif
       
      
       ACTIVE(MAXDIMENSION,MAXDIMENSION,MAXDIMENSION) = .TRUE.
      IF (USEPERIODX) THEN
       ELEVATION(MAXDIMENSION,MAXDIMENSION,MAXDIMENSION) =  &
          ELEVATION(0,MAXDIMENSION,MAXDIMENSION)
       ELSE
          if(useperiody) then
               ELEVATION(MAXDIMENSION,MAXDIMENSION,MAXDIMENSION) = &
              elevation(maxdimension,0,maxdimension)
          else
               ELEVATION(MAXDIMENSION,MAXDIMENSION,MAXDIMENSION) = &
                NORMAL2(RDUMMY)*ERROR*FACTOR
          endif
       endif
        
	  XSTAVG = (ELEVATION(0,0,0) &
      				+ELEVATION(MAXDIMENSION,MAXDIMENSION,MAXDIMENSION) &
      				+ELEVATION(0,MAXDIMENSION,MAXDIMENSION) &
      				+ELEVATION(MAXDIMENSION,0,MAXDIMENSION) &
      				+ELEVATION(MAXDIMENSION,MAXDIMENSION,0) &
      				+ELEVATION(MAXDIMENSION,0,0)  &
      				+ELEVATION(0,0,MAXDIMENSION) &
      				+ELEVATION(0,MAXDIMENSION,0))/8.0 
       CURLEVEL = 1 
100      CONTINUE
        MULFACTOR = ERROR*EXP(-CURLEVEL*CHARDIMENSION) 
        WRITE(*,110)curlevel,mulfactor,interval
110     format('CURLEVEL= ',i5,' MULFACTOR = ',g12.5,  &
              ' INTERVAL=',i5) 
		 call DOIEDGE(0,MAXDIMENSION) 
		 call DOIEDGE(0,0)
		 call DOJEDGE(0,MAXDIMENSION) 
		 call DOJEDGE(0,0) 
		  call DOKEDGE(0,0) 
          if (.not.useperiody) then
		   call DOIEDGE(MAXDIMENSION,MAXDIMENSION) 
		   call DOIEDGE(MAXDIMENSION,0)
          endif 
		  IF (.not. USEPERIODX) THEN 
		     call DOJEDGE(MAXDIMENSION,MAXDIMENSION)
		  endif 
		  IF (.NOT. USEPERIODX) THEN
		     call DOJEDGE(MAXDIMENSION,0)
		  endif 
		  if (.not.useperiody) then
              call DOKEDGE(0,MAXDIMENSION)
           endif 
		  IF (.NOT. USEPERIODX) THEN 
		     call DOKEDGE(MAXDIMENSION,0)
		  endif 
		  IF (.not.USEPERIODX) THEN
              if (.not.useperiody) then
		       call DOKEDGE(MAXDIMENSION,maxdimension)
              endif
		  endif 
        CURLEVEL = CURLEVEL + 1 
        INTERVAL = INTERVAL /  2 
      if (interval.gt.0) goto 100
		  IF (USEPERIODX) THEN
		     DO J = 0 , MAXDIMENSION
		        ELEVATION(MAXDIMENSION,J,MAXDIMENSION) &
      		 = ELEVATION(0,J,MAXDIMENSION)
	            active(maxdimension,j,maxdimension)=.true.
              enddo
           endif 
		  IF (USEPERIODX) THEN
		     DO J = 0 , MAXDIMENSION 
		        ELEVATION(MAXDIMENSION,J,0) &
      		 = ELEVATION(0,J,0)
	            active(maxdimension,j,0)=.true.
              enddo
           endif
		  IF (USEPERIODX) THEN
		     DO K = 0 , MAXDIMENSION
		        ELEVATION(MAXDIMENSION,0,K) &
      		 = ELEVATION(0,0,K)
	            active(maxdimension,0,k)=.true.
              enddo
           endif
           if (useperiody) then
              do i=0,maxdimension
                 elevation(i,maxdimension,0)= &
                elevation(i,0,0)
                 active(i,maxdimension,0)=.true.
              enddo
           endif
           if (useperiody) then
              do i=0,maxdimension
                 elevation(i,maxdimension,maxdimension)= &
                elevation(i,0,maxdimension) 
	           active(i,maxdimension,maxdimension)=.true.
              enddo
           endif
           if (useperiody) then
              do k=0,maxdimension
                 elevation(0,maxdimension,k)= &
                 elevation(0,0,k)
	           active(0,maxdimension,k)=.true.
              enddo
           endif
           if (useperiody) then
              do k=0,maxdimension
                 elevation(maxdimension,maxdimension,k)= &
                 elevation(maxdimension,0,k)
	           active(maxdimension,maxdimension,k)=.true.
              enddo
           endif
		  IF (USEPERIODX) THEN
		     DO K = 0 , MAXDIMENSION
		        ELEVATION(MAXDIMENSION,MAXDIMENSION,K) &
      		 = ELEVATION(0,MAXDIMENSION,K)
	            active(maxdimension,maxdimension,k)=.true.
              enddo
	     endif
	call checkedges()
	  call DOIJFACE(0) 
	  call DOIJFACE(MAXDIMENSION) 
	  call DOJKFACE(0) 
	  call DOIKFACE(0)
	  IF (USEPERIODX) THEN
	    DO J = 0 , MAXDIMENSION 
	    DO K = 0 , MAXDIMENSION
	        ELEVATION(MAXDIMENSION,J,K) = &
      	  ELEVATION(0,J,K)
	        active(maxdimension,j,k)=.true.
          enddo
          enddo
        ELSE
	    call DOJKFACE(MAXDIMENSION)
        endif 
        if (useperiody) then
           do i=0,maxdimension
              do k=0,maxdimension
                 elevation(i,maxdimension,k)=  &
                   elevation(i,0,k)
	           active(i,maxdimension,k)=.true.
              enddo
           enddo
        else
	     call DOIKFACE(MAXDIMENSION)
        endif
	  call checkfaces() 
      CURLEVEL = 1 
      INTERVAL = MAXDIMENSION /  2 
      ELAVERAGE = (ELEVATION(0,0,0)+ELEVATION(0,MAXDIMENSION,0)+  &
         ELEVATION(MAXDIMENSION,0,0)  &
             +ELEVATION(MAXDIMENSION,MAXDIMENSION,0)+ &
         ELEVATION(0,0,MAXDIMENSION) &
         +ELEVATION(0,MAXDIMENSION,MAXDIMENSION)+ &
         ELEVATION(MAXDIMENSION,0,MAXDIMENSION) &
             +ELEVATION(MAXDIMENSION,MAXDIMENSION,MAXDIMENSION))/8.0 
      MULFACTOR = ERROR*EXP(-CURLEVEL*CHARDIMENSION) 
	  ACTIVE(INTERVAL,INTERVAL,INTERVAL) = .TRUE. 
      ELEVATION(INTERVAL,INTERVAL,INTERVAL) = EFACTOR*ELAVERAGE  &
            +MULFACTOR*NORMAL2(RDUMMY)  
      CURLEVEL = 2 
      INTERVAL = INTERVAL /  2 
200     CONTINUE
         MULFACTOR = ERROR*EXP(-CURLEVEL*CHARDIMENSION) 
        WRITE(*,110)curlevel,mulfactor,interval
			KK = INTERVAL 
210			CONTINUE
			JJ = INTERVAL 
220			CONTINUE
			II = INTERVAL 
230			CONTINUE
			  ELAVERAGE = 0.0 
			  ELNUM = 0.0 
			  IF (.NOT. ACTIVE(II,JJ,KK)) THEN
               DO LL = 0 , 12
                 INN = II + CPLACE(1,(2*LL+1))*INTERVAL 
                 JNN = JJ + CPLACE(2,(2*LL+1))*INTERVAL 
					KNN = KK + CPLACE(3,(2*LL+1))*INTERVAL 
                 IMM = II + CPLACE(1,(2*LL+2))*INTERVAL 
                 JMM = JJ + CPLACE(2,(2*LL+2))*INTERVAL 
					KMM = KK + CPLACE(3,(2*LL+2))*INTERVAL 
					IF (ACTIVE(INN,JNN,KNN).AND.ACTIVE(IMM,JMM,KMM)) &
     					THEN
						ELAVERAGE = ELAVERAGE + ELEVATION(INN,JNN,KNN)  &
      					    +ELEVATION(IMM,JMM,KMM) 
						ELNUM = ELNUM + 2.0 
					ENDif 
			  ENDdo 
			  IF (ELNUM.gt.0.0) THEN
              ELEVATION(II,JJ,KK) = EFACTOR*  &
                  ELAVERAGE/ELNUM+MULFACTOR*NORMAL2(RDUMMY) 
				  ACTIVE(II,JJ,KK) = .TRUE. 
			  ENDif 
		     ENDif 
			II = II+2*INTERVAL 
			if (II.le.(MAXDIMENSION-INTERVAL)) goto 230 
			JJ = JJ+2*INTERVAL 
			if (JJ.le.(MAXDIMENSION-INTERVAL)) goto 220 
			KK = KK+2*INTERVAL 
			if (KK.le.(MAXDIMENSION-INTERVAL)) goto 210 

			KK = INTERVAL 
240			CONTINUE
			JJ = INTERVAL 
250			CONTINUE
			II = INTERVAL 
260			CONTINUE
			  ELAVERAGE = 0.0 
			  ELNUM = 0.0 
			  IF (.NOT. ACTIVE(II,JJ,KK)) THEN
              DO LL = 0 , 12
                INN = II + CPLACE(1,(2*LL+1))*INTERVAL 
                JNN = JJ + CPLACE(2,(2*LL+1))*INTERVAL 
					KNN = KK + CPLACE(3,(2*LL+1))*INTERVAL 
                IMM = II + CPLACE(1,(2*LL+2))*INTERVAL 
                JMM = JJ + CPLACE(2,(2*LL+2))*INTERVAL 
				 KMM = KK + CPLACE(3,(2*LL+2))*INTERVAL 
				 IF (ACTIVE(INN,JNN,KNN).AND. ACTIVE(IMM,JMM,KMM)) &
     				  THEN
						ELAVERAGE = ELAVERAGE + ELEVATION(INN,JNN,KNN) &
      					    +ELEVATION(IMM,JMM,KMM) 
						ELNUM = ELNUM + 2.0 
				 ENDif 
			  ENDdo 
			  IF (ELNUM.gt.0.0) THEN
              ELEVATION(II,JJ,KK) = EFACTOR* &
                  ELAVERAGE/ELNUM+MULFACTOR*NORMAL2(RDUMMY) 
				  ACTIVE(II,JJ,KK) = .TRUE. 
			  ENDif 
		     ENDif 
			II = II+INTERVAL 
			if (II.le.(MAXDIMENSION-INTERVAL)) goto 260 
			JJ = JJ+INTERVAL 
			if (JJ.le.(MAXDIMENSION-INTERVAL)) goto 250 
			KK = KK+INTERVAL 
			if ( KK.le.(MAXDIMENSION-INTERVAL)) goto 240 

         CURLEVEL = CURLEVEL + 1 
         INTERVAL = INTERVAL /  2 
       if (interval.gt.0) goto 200
	  XSUM = 0.0 
	  XTOT = 0.0 
	  XSSS = 0.0 
	  XOVER = 0.0 
	  DO K = 0 , MAXDIMENSION
	  DO I = 0 , MAXDIMENSION
	  DO J = 0 , MAXDIMENSION
	    if (elevation(i,j,k).lt.-1.0e+24) then
	       write(debug,356) i,j,k
356    format(3i5)
          endif
	    IF (ELEVATION(I,J,K).gt.XSTAVG) XOVER = XOVER+1.0 
	    XSUM = XSUM + ELEVATION(I,J,K) 
		 XTOT = XTOT + 1.0 
		 XSSS = XSSS + ELEVATION(I,J,K)**2 
	  ENDdo
        enddo
        enddo 
	  XSSS=SQRT((XSSS-XSUM**2/XTOT)/(XTOT-1)) 
	  XSUM=XSUM/XTOT
        xxxx=xover/xtot 
        WRITE(*,300)xsum,xsss,xxxx,xstavg
300     format('AVG=',G12.5,' STD=',G12.5, ' FRACT> =',  &
        G12.5,' xstavg=',g12.5)
      IF (USEBLENDx) call BLENDX()
      if (useblendy) call blendy()
	  XSUM = 0.0 
	  XTOT = 0.0
	  XSSS = 0.0 
      iisz=ifn-ist+1
      jjsz=jfn-jst+1
      kksz=kfn-kst+1
	WRITE(*,822) IISZ,IST,IFN,JJSZ,JST,JFN,KKSZ,KST,KFN
822   FORMAT('ISIZE=',I5,' IST=',I5,' IFN=',I5,/, &
      'JSIZE=',I5,' JST=',I5,' JFN=',I5,/, &
      'KSIZE=',I5,' KST=',I5,' KFN=',I5)
      do k=kst,kfn
         do i=ist,ifn
            do j=jst,jfn
               irec=j+(jjsz)*(i)+(iisz)*(jjsz)*(k)+1
               xxxx=elevation(i,j,k)
               if (useexp) xxxx=exp(xxxx)
	ELEVATION(I,J,K)=XXXX
	    XSUM = XSUM + xxxx
		xsss = xsss + xxxx*xxxx 
		 XTOT = XTOT + 1.0 
        !       write(outfile,rec=irec,err=333) xxxx
            enddo
         enddo
      enddo
	  XSSS=SQRT((XSSS-XSUM**2/XTOT)/(XTOT-1)) 
	  XSUM=XSUM/XTOT
        WRITE(*,336)xsum,xsss
336     format('EXP AVG=',G12.5,' STD=',G12.5)
        XAVG=XSUM
        XSUM=0.0
	  XTOT=0.0
	XSSS=0.0
      do k=kst,kfn
         do i=ist,ifn
            do j=jst,jfn
               irec=j+(jjsz)*(i)+(iisz)*(jjsz)*(k)+1
               xxxx=elevation(i,j,k)/XAVG
	    XSUM = XSUM + xxxx
		xsss = xsss + xxxx*xxxx 
		 XTOT = XTOT + 1.0 
!               write(outfile,rec=irec,err=333) xxxx
            enddo
         enddo
      enddo
	  XSSS=SQRT((XSSS-XSUM**2/XTOT)/(XTOT-1)) 
	  XSUM=XSUM/XTOT
        WRITE(*,337)xsum,xsss
337     format('OUTPUT AVG=',G12.5,' STD=',G12.5)

	goto 335
333   write(*,334) i,j,k,irec,ist,ifn,jst,jfn,kst,kfn,iisz,jjsz,kksz
334   format('record error: ',13i9)
!      DO I = IST , IFN
!         DO J = JST , JFN 
!				do K = KST , KFN
!      IREC=J+MY*(I-1)+MX*MY*(K-1)
!      READ(INRESIST,REC=IREC) RRRR
!               xxxx=elevation(i,j,k)-xsum/xsss
!               write(outfile,310) xxxx
!310            format(g8.5)
!               enddo
!         enddo
!      enddo
335   continue
        CLOSE(OUTFILE)
	  close(debug)
      DEALLOCATE(ELEVATION,ACTIVE,NEWACTIVE,DONE)
      stop 
      END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
      real*4 FUNCTION NORMAL1(DUMMY)
      use ifport
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      integer dummy
      integer I,J
      real*4  X, TEMP,xxxx
        J = 1 
        TEMP = 0.0 
        do I = 1,48
	       CALL RANDOM_NUMBER(XXXX) 
             TEMP = TEMP +  xxxx -0.5 
        enddo
        NORMAL1 = TEMP / 2.0 
      return
      end
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      real*4 FUNCTION NORMAL2(DUMMY)
      use ifport
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      integer dummy,j
      real*4 V1, V2, SS
	   real*4 temp1,temp2
      J = 1 
100   continue
      CALL RANDOM_NUMBER(TEMP1) 
	   CALL RANDOM_NUMBER(TEMP2)
      V1 = 2.0*temp1-1.0 
      V2 = 2.0*temp2-1.0 
      SS = V1*V1 + V2*V2 
      if (ss.ge.1.0) goto 100
      SS = SQRT(-2.0*Log(SS)/SS) 
      NORMAL2 = V1*SS 
      return
      END 
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      real*4 FUNCTION LOGNOR(MULLFACTOR)
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      real*4 mulLfactor
      real*4 XXXX,normal2
      XXXX = EXP(NORMAL2(RDUMMY)*MULFACTOR) 
      IF (XXXX.gt.1.0) YOVER = YOVER+1.0 
      YTOT = YTOT+1.0 
      YSUM = YSUM+XXXX 
      YSSS = YSSS+xxxx*xxxx 
      LOGNOR = XXXX 
      return
      end
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SETUP() 
      USE MATRIX_GLOBALS
      IMPLICIT NONE
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
      CPLACE(1,1) = -1 
      CPLACE(2,1) = -1 
 	  CPLACE(3,1) = -1 
      CPLACE(1,2) = 1 
      CPLACE(2,2) = 1 
 	  CPLACE(3,2) = 1 
      CPLACE(1,3) = -1 
      CPLACE(2,3) = 1 
 	  CPLACE(3,3) = -1 
      CPLACE(1,4) = 1 
      CPLACE(2,4) = -1 
 	  CPLACE(3,4) = 1  
      CPLACE(1,5) = 1 
      CPLACE(2,5) = -1 
 	  CPLACE(3,5) = -1 
      CPLACE(1,6) = -1 
      CPLACE(2,6) = 1 
 	  CPLACE(3,6) = 1 
      CPLACE(1,7) = 1 
      CPLACE(2,7) = 1 
 	  CPLACE(3,7) = -1 
      CPLACE(1,8) = -1 
      CPLACE(2,8) = -1 
 	  CPLACE(3,8) = 1  
      CPLACE(1,9) = 0 
      CPLACE(2,9) = 0 
 	  CPLACE(3,9) = -1 
      CPLACE(1,10) = 0 
      CPLACE(2,10) = 0 
	  CPLACE(3,10) = 1 
      CPLACE(1,11) = 0 
      CPLACE(2,11) = 1 
	  CPLACE(3,11) = 0 
      CPLACE(1,12) = 0 
      CPLACE(2,12) = -1 
	  CPLACE(3,12) = 0 
      CPLACE(1,13) = 1 
      CPLACE(2,13) = 0 
	  CPLACE(3,13) = 0 
      CPLACE(1,14) = -1 
      CPLACE(2,14) = 0 
	  CPLACE(3,14) = 0 
      CPLACE(1,15) = 1 
      CPLACE(2,15) = 1 
	  CPLACE(3,15) = 0 
      CPLACE(1,16) = -1 
      CPLACE(2,16) = -1 
	  CPLACE(3,16) = 0 
      CPLACE(1,17) = -1 
      CPLACE(2,17) = 1 
	  CPLACE(3,17) = 0 
      CPLACE(1,18) = 1 
      CPLACE(2,18) = -1 
	  CPLACE(3,18) = 0 
      CPLACE(1,19) = 0 
      CPLACE(2,19) = 1 
	  CPLACE(3,19) = 1 
      CPLACE(1,20) = 0 
      CPLACE(2,20) = -1 
	  CPLACE(3,20) = -1  
      CPLACE(1,21) = 0 
      CPLACE(2,21) = -1 
	  CPLACE(3,21) = 1 
      CPLACE(1,22) = 0 
      CPLACE(2,22) = 1 
	  CPLACE(3,22) = -1 
      CPLACE(1,23) = 1 
      CPLACE(2,23) = 0 
	  CPLACE(3,23) = -1 
      CPLACE(1,24) = -1 
      CPLACE(2,24) = 0 
	  CPLACE(3,24) = 1 
      CPLACE(1,25) = 1 
      CPLACE(2,25) = 0 
	  CPLACE(3,25) = 1  
      CPLACE(1,26) = -1 
      CPLACE(2,26) = 0 
	  CPLACE(3,26) = -1  
      DPLACE(1,1) = -1 
      DPLACE(2,1) = -1 
      DPLACE(1,2) = -1 
      DPLACE(2,2) = 1 
      DPLACE(1,3) = 1 
      DPLACE(2,3) = -1 
      DPLACE(1,4) = 1 
      DPLACE(2,4) = 1 
      DPLACE(1,5) = 0 
      DPLACE(2,5) = 1 
      DPLACE(1,6) = 0 
      DPLACE(2,6) = -1 
      DPLACE(1,7) = 1 
      DPLACE(2,7) = 0 
      DPLACE(1,8) = -1 
      DPLACE(2,8) = 0 
      return
      END 
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      subroutine GETPARAMETERS() 	
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      INTEGER, DIMENSION(:), ALLOCATABLE :: RANDSEED
      integer BLENDUSEx, blendusey,PERIODUSEx,periodusey
      integer expuse,IERROR
      INTEGER RNSEED, NRANDOM
      READ(PARAMS,*) MAXDIMENSION, MAXLEVEL
      write(*,833) maxdimension, maxlevel
833   format(' maxdimension=',i5, ' maxlevel=',i5)
      ALLOCATE(ELEVATION(0:MAXDIMENSION,0:MAXDIMENSION,0:MAXDIMENSION),STAT=IERROR)
      ALLOCATE(ACTIVE(0:MAXDIMENSION,0:MAXDIMENSION,0:MAXDIMENSION),STAT=IERROR)
      ALLOCATE(NEWACTIVE(0:MAXDIMENSION,0:MAXDIMENSION),STAT=IERROR)
      ALLOCATE(DONE(0:MAXDIMENSION,0:MAXDIMENSION),STAT=IERROR)
      READ(PARAMS,*)ERROR
!		ERROR determines the magnitude of variation of the simulated values
      READ(PARAMS,*)CHARDIMENSION
!      CHARDIMENSION = -CHARDIMENSION * Log(1.0/EXP(Log(2.0)/3.0))
!		CHARDIMENSION determines how rapidly the variation declines with scale
!		it is closely related to fractal dimension and should range between 0.0
!		(uncorrelated random noise) to 1.0 (rough at large scales, smooth at small
!		scales 
      WRITE(*,100) chardimension
100   format('NEW CHARACTERISIC DIMENSION=',g12.5)
!        A random number generator seed
      READ(PARAMS,*)RNSEED 
!        ISEED1 = rnSEED/4.6566128752458E-10
!      call srand(rdummy)
!		Set FACTOR to 1.0 to allow random variation of the corners and edges of the matrix
!			otherwise 0.0 (usually use 1.0) 
      READ(PARAMS,*) FACTOR
!		The actual X,Y,and Z dimensions of the output - the maximum values are 256x256x256
      READ(PARAMS,*)IACT,JACT,KACT
!		These parameters are used when the X or Y boundaries (or both) are periodic, otherwise set
!			all of these to zero.
!		If IACT,JACT are equal to the maximum compiled dimension (256x256) then it is best
!		to set PERIODUSEX to 1 if the matrix is periodic in the X direction, and the same with
!		PERIODUSEY for Y periodicity.  If IACT and/or JACT is less than the maximum dimensions, then
!		set PERIODUSEX and/or PERIODUSEY to zero and use BLENDUSEX and/or BLENDUSEY to unity.
!		This accomplishes more or less the same thing by doing a seletive averaging of generated
!		variates on the two sides of the periodic boundary.  NOTE:  the program assumes that the
!		matrix is never periodic in the Z direction 
      READ(PARAMS,*)BLENDUSEx,blendusey,PERIODUSEx,periodusey
!		I've forgotten what this is used for, just let it be a negative number.
      READ(PARAMS,*)EDECAY
!		Set expuse to unity if you want the output to be a lognormal distribution (e.g., to
!			simulate rock erodibility, which is presumably >= 0) otherwise 0 for a normal
!			distribution (allowing negative values)
      read(params,*)expuse
      IF (EDECAY.gt.0.0) THEN
       EFACTOR = 1.0 - EXP(-EDECAY*CHARDIMENSION)
      ELSE
       EFACTOR = 1.0 
      endif
      IF (BLENDUSEx.gt.0)then
          USEBLENDx = .TRUE. 
      ELSE 
         USEBLENDx = .FALSE. 
      endif
      IF (PERIODUSEx.gt.0) THEN 
         USEPERIODx = .TRUE. 
      ELSE 
         USEPERIODx = .FALSE. 
      endif
      IF (BLENDUSEy.gt.0)then
          USEBLENDy = .TRUE. 
      ELSE 
         USEBLENDy = .FALSE. 
      endif
      IF (PERIODUSEy.gt.0) THEN 
         USEPERIODy = .TRUE. 
      ELSE 
         USEPERIODy = .FALSE. 
      endif
      if (expuse.gt.0) then
         useexp=.true.
      else
         useexp=.false.
      endif
      IST = (MAXDIMENSION - IACT) / 2  
      IFN = IST+IACT-1 
      JST = (MAXDIMENSION - JACT) / 2  
      JFN = JST+JACT-1 
      KST =(MAXDIMENSION-KACT) / 2  
      KFN = KST+KACT-1 
      CALL RANDOM_SEED(SIZE=NRANDOM)
      ALLOCATE(RANDSEED(NRANDOM))
      RANDSEED=RNSEED
      CALL RANDOM_SEED(PUT=RANDSEED)
      DEALLOCATE (RANDSEED)
      return
      END 

!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE BLENDX
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      INTEGER I,J,k,IL,IR,JL,JR 
      REAL*8    XST,XL,XR,XFN,WGT,AA,BB
      REAL*8 LNEWELEVS(29),RNEWELEVS(29)
            XST = IST
            XFN = IFN
          AA = 0.5/29.0
          BB = 0.5 - AA
          do 100 k = kst,kfn
          DO 100 J = JST,JFN
          DO 110 I = 1,29
            IL=IST+I-1
            IR=IFN-I+1
            WGT = I*AA + BB
            LNEWELEVS(I) = WGT*ELEVATION(IL,J,k)+ &
                 (1.0-WGT)*ELEVATION(IR,J,k)
            RNEWELEVS(I) = WGT*ELEVATION(IR,J,k)+ &
                 (1.0-WGT)*ELEVATION(IL,J,k)
110       CONTINUE
          DO 120 I = 1,29
            IL = IST+I-1
            IR = IFN-I+1
            ELEVATION(IL,J,k) = LNEWELEVS(I)
            ELEVATION(IR,J,k) = RNEWELEVS(I)
120       CONTINUE
100       CONTINUE
      RETURN
      END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE BLENDY
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      INTEGER I,J,k,JT,JB,IT,IB 
      REAL*8    YST,YT,YB,YFN,WGT,AA,BB
      REAL*8 TNEWELEVS(29),BNEWELEVS(29)
            YST = JST
            YFN = JFN
          AA = 0.5/29.0
          BB = 0.5 - AA
          do 100 k = kst,kfn
          DO 100 I = IST,IFN
          DO 110 J = 1,29
            JT=JST+J-1
            JB=JFN-J+1
            WGT = J*AA + BB
            TNEWELEVS(J) = WGT*ELEVATION(I,JT,k)+  &
                 (1.0-WGT)*ELEVATION(I,JB,k)
            BNEWELEVS(J) = WGT*ELEVATION(I,JB,k)+  &
                 (1.0-WGT)*ELEVATION(I,JT,k)
110       CONTINUE
          DO 120 J = 1,29
            JT = JST+J-1
            JB = JFN-J+1
            ELEVATION(I,JT,k) = TNEWELEVS(J)
            ELEVATION(I,JB,k) = BNEWELEVS(J)
120       CONTINUE
100       CONTINUE
      RETURN
      END
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DOIEDGE(J,k) 
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      integer i,j,k
	real*4 normal2
!      write(*,50)
50    format('doiedge') 
      DO I = 0 , MAXDIMENSION
         DONE(I,J) = .FALSE. 
         NEWACTIVE(I,J) = .FALSE. 
      ENDDO 
        DO I = 0 , MAXDIMENSION
            IF (ACTIVE(I,J,K)) THEN
                INN = I + INTERVAL 
                IF ((INN .gt. 0).AND. (INN .lt. MAXDIMENSION)) THEN
                    IF (.NOT. ACTIVE(INN,J,K)) THEN
                       ELAVERAGE = (ELEVATION(INN-INTERVAL,J,K) + &
                          ELEVATION(INN+INTERVAL,J,K))/2.0 
                       ELEVATION(INN,J,K) = EFACTOR*ELAVERAGE+ &
                           MULFACTOR*NORMAL2(RDUMMY) 
                       NEWACTIVE(INN,J) = .TRUE. 
                       DONE(INN,J) = .TRUE. 
                    ENDIF
                    ENDIF 
            ENDIF
            ENDDO 
        DO I = 0 , MAXDIMENSION
                IF (NEWACTIVE(I,J)) ACTIVE(I,J,K) = .TRUE. 
                NEWACTIVE(I,J) = .FALSE. 
        ENDDO 
        RETURN
        END 
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DOJEDGE(I,K) 
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      integer i,j,k
	real*4 normal2 
!      write(*,50)
50    format('dojedge') 
      DO J = 0 , MAXDIMENSION
         DONE(J,K) = .FALSE. 
         NEWACTIVE(J,K) = .FALSE. 
      ENDDO 
        DO J = 0 , MAXDIMENSION
            IF (ACTIVE(I,J,K)) THEN
                JNN = J + INTERVAL 
                IF ((JNN .gt. 0).AND. (JNN .lt. MAXDIMENSION)) THEN
                    IF (.NOT. ACTIVE(I,JNN,K)) THEN
                       ELAVERAGE = (ELEVATION(I,JNN-INTERVAL,K) + &
                          ELEVATION(I,JNN+INTERVAL,K))/2.0 
                       ELEVATION(I,JNN,K) = EFACTOR*ELAVERAGE+ &
                           MULFACTOR*NORMAL2(RDUMMY) 
                       NEWACTIVE(JNN,K) = .TRUE. 
                       DONE(JNN,K) = .TRUE. 
                    ENDIF
                    ENDIF 
            ENDIF
            ENDDO 
        DO J = 0 , MAXDIMENSION
                IF (NEWACTIVE(J,K)) ACTIVE(I,J,K) = .TRUE. 
                NEWACTIVE(J,K) = .FALSE. 
        ENDDO 
        RETURN
        END 

!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine DOKEDGE(I,J) 
      USE MATRIX_GLOBALS
      IMPLICIT NONE
      integer i,j,k
	real*4 normal2 
!      write(*,50)
50    format('dokedge') 
      DO K = 0 , MAXDIMENSION
         DONE(I,K) = .FALSE. 
         NEWACTIVE(I,K) = .FALSE. 
      ENDDO 
        DO K = 0 , MAXDIMENSION
            IF (ACTIVE(I,J,K)) THEN
                KNN = K + INTERVAL 
                IF ((KNN .gt. 0).AND. (KNN .lt. MAXDIMENSION)) THEN
                    IF (.NOT. ACTIVE(I,J,KNN)) THEN
                       ELAVERAGE = (ELEVATION(I,J,KNN-INTERVAL) + &
                          ELEVATION(I,J,KNN+INTERVAL))/2.0 
                       ELEVATION(I,J,KNN) = EFACTOR*ELAVERAGE+ &
                           MULFACTOR*NORMAL2(RDUMMY) 
                       NEWACTIVE(I,KNN) = .TRUE. 
                       DONE(I,KNN) = .TRUE. 
                    ENDIF
                    ENDIF 
            ENDIF
            ENDDO 
        DO K = 0 , MAXDIMENSION
                IF (NEWACTIVE(I,K)) ACTIVE(I,J,K) = .TRUE. 
                NEWACTIVE(I,K) = .FALSE. 
        ENDDO 
        RETURN
        END 
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DOIJFACE(K) 
      USE MATRIX_GLOBALS
      IMPLICIT NONE
	INTEGER I,J,K
	real*4 normal2
!      write(*,50)
50    format('doijface') 
      DO I = 0 , MAXDIMENSION
      DO J = 0 , MAXDIMENSION
         DONE(I,J) = .FALSE. 
         NEWACTIVE(I,J) = .FALSE. 
      ENDDO 
	ENDDO
      CURLEVEL = 1 
      INTERVAL = MAXDIMENSION /  2 
      ELAVERAGE = (ELEVATION(0,0,K)+ELEVATION(0,MAXDIMENSION,K)+ &
         ELEVATION(MAXDIMENSION,0,K)  &
             +ELEVATION(MAXDIMENSION,MAXDIMENSION,K))/4.0 
      MULFACTOR = ERROR*EXP(-CURLEVEL*CHARDIMENSION) 
      ELEVATION(INTERVAL,INTERVAL,K) = EFACTOR*ELAVERAGE &
            +MULFACTOR*NORMAL2(RDUMMY)  
      DONE(INTERVAL,INTERVAL) = .TRUE. 
      ACTIVE(INTERVAL,INTERVAL,K)= .TRUE. 
      DO J = 0 , MAXDIMENSION
         DONE(0,J) = .TRUE. 
         DONE(MAXDIMENSION,J) = .TRUE. 
      ENDDO 
      DO I = 0 , MAXDIMENSION
         DONE(I,0) = .TRUE. 
         DONE(I,MAXDIMENSION) = .TRUE. 
      ENDDO 
      CURLEVEL = 2 
      INTERVAL = INTERVAL /  2 
100   CONTINUE
         MULFACTOR = ERROR*EXP(-CURLEVEL*CHARDIMENSION) 
			JJ = INTERVAL 
120			CONTINUE
			II = INTERVAL 
130			CONTINUE
			  ELAVERAGE = 0.0 
			  ELNUM = 0.0 
			  IF (.NOT. ACTIVE(II,JJ,K)) THEN
                DO LL = 1 , 4
                  INN = II + DPLACE(1,LL)*INTERVAL 
                  JNN = JJ + DPLACE(2,LL)*INTERVAL
!				write(debug,641)ii,jj,ll,inn,jnn,k
641    format(6i6) 
					IF (ACTIVE(INN,JNN,K)) THEN
						ELAVERAGE = ELAVERAGE + ELEVATION(INN,JNN,K) 
						ELNUM = ELNUM + 1.0 
					ENDIF 
			    ENDDO 
			    IF (ELNUM.gt.0.0) THEN
!				write(debug,641)ii,jj,k
                  ELEVATION(II,JJ,K) = EFACTOR* &
                  ELAVERAGE/ELNUM+MULFACTOR*NORMAL2(RDUMMY) 
				   ACTIVE(II,JJ,K) = .TRUE. 
			    ENDIF 
		     ENDIF 
			II = II+2*INTERVAL 
			IF (II.LE.(MAXDIMENSION-INTERVAL)) GOTO 130
			JJ = JJ+2*INTERVAL 
			IF (JJ.LE.(MAXDIMENSION-INTERVAL)) GOTO 120 

			JJ = INTERVAL 
140   		CONTINUE
			II = INTERVAL 
150			CONTINUE
			  ELAVERAGE = 0.0 
			  ELNUM = 0.0 
			  IF (.NOT. ACTIVE(II,JJ,K)) THEN
                DO LL = 5 , 8
                  INN = II + DPLACE(1,LL)*INTERVAL 
                  JNN = JJ + DPLACE(2,LL)*INTERVAL 
					IF (ACTIVE(INN,JNN,K)) THEN
						ELAVERAGE = ELAVERAGE + ELEVATION(INN,JNN,K) 
						ELNUM = ELNUM + 1.0 
					ENDIF 
			    ENDDO 
			    IF (ELNUM.gt.0.0) THEN
!				write(debug,641)ii,jj,k
                  ELEVATION(II,JJ,K) = EFACTOR*  &
                   ELAVERAGE/ELNUM+MULFACTOR*NORMAL2(RDUMMY) 
				   ACTIVE(II,JJ,K) = .TRUE. 
			    ENDIF 
		     ENDIF 
			II = II+INTERVAL 
			IF (II.LE.(MAXDIMENSION-INTERVAL)) GOTO 150 
			JJ = JJ+INTERVAL 
			IF (JJ.LE.(MAXDIMENSION-INTERVAL)) GOTO 140 

           CURLEVEL = CURLEVEL + 1 
           INTERVAL = INTERVAL /  2 
         IF (INTERVAL.GT.0) GOTO 100
      RETURN 
      END 
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DOIKFACE(J) 
      USE MATRIX_GLOBALS
      IMPLICIT NONE
	INTEGER I,J,K
	real*4 normal2
!      write(*,50)
50    format('doikface') 
      DO I = 0 , MAXDIMENSION
      DO K = 0 , MAXDIMENSION
         DONE(I,K) = .FALSE. 
         NEWACTIVE(I,K) = .FALSE. 
      ENDDO 
	ENDDO
      CURLEVEL = 1 
      INTERVAL = MAXDIMENSION /  2 
      ELAVERAGE = (ELEVATION(0,J,0)+ELEVATION(0,J,MAXDIMENSION)+ &
         ELEVATION(MAXDIMENSION,J,0) &
             +ELEVATION(MAXDIMENSION,J,MAXDIMENSION))/4.0 
      MULFACTOR = ERROR*EXP(-CURLEVEL*CHARDIMENSION) 
      ELEVATION(INTERVAL,J,INTERVAL) = EFACTOR*ELAVERAGE &
            +MULFACTOR*NORMAL2(RDUMMY)  
      DONE(INTERVAL,INTERVAL) = .TRUE. 
      ACTIVE(INTERVAL,J,INTERVAL)= .TRUE. 
      DO K = 0 , MAXDIMENSION
         DONE(0,K) = .TRUE. 
         DONE(MAXDIMENSION,K) = .TRUE. 
      ENDDO 
      DO I = 0 , MAXDIMENSION
         DONE(I,0) = .TRUE. 
         DONE(I,MAXDIMENSION) = .TRUE. 
      ENDDO 
      CURLEVEL = 2 
      INTERVAL = INTERVAL /  2 
100   CONTINUE
         MULFACTOR = ERROR*EXP(-CURLEVEL*CHARDIMENSION) 
			KK = INTERVAL 
120			CONTINUE
			II = INTERVAL 
130			CONTINUE
			  ELAVERAGE = 0.0 
			  ELNUM = 0.0 
			  IF (.NOT. ACTIVE(II,J,KK)) THEN
                DO LL = 1 , 4
                  INN = II + DPLACE(1,LL)*INTERVAL 
                  KNN = KK + DPLACE(2,LL)*INTERVAL 
					IF (ACTIVE(INN,J,KNN)) THEN
						ELAVERAGE = ELAVERAGE + ELEVATION(INN,J,KNN) 
						ELNUM = ELNUM + 1.0 
					ENDIF 
			    ENDDO 
			    IF (ELNUM.gt.0.0) THEN
                  ELEVATION(II,J,KK) = EFACTOR*  &
                  ELAVERAGE/ELNUM+MULFACTOR*NORMAL2(RDUMMY) 
				   ACTIVE(II,J,KK) = .TRUE. 
			    ENDIF 
		     ENDIF 
			II = II+2*INTERVAL 
			IF (II.LE.(MAXDIMENSION-INTERVAL)) GOTO 130
			KK = KK+2*INTERVAL 
			IF (KK.LE.(MAXDIMENSION-INTERVAL)) GOTO 120 

			KK = INTERVAL 
140   		CONTINUE
			II = INTERVAL 
150			CONTINUE
			  ELAVERAGE = 0.0 
			  ELNUM = 0.0 
			  IF (.NOT. ACTIVE(II,JJ,KK)) THEN
                DO LL = 5 , 8
                  INN = II + DPLACE(1,LL)*INTERVAL 
                  KNN = KK + DPLACE(2,LL)*INTERVAL 
					IF (ACTIVE(INN,J,KNN)) THEN
						ELAVERAGE = ELAVERAGE + ELEVATION(INN,J,KNN) 
						ELNUM = ELNUM + 1.0 
					ENDIF 
			    ENDDO 
			    IF (ELNUM.gt.0.0) THEN
                  ELEVATION(II,J,KK) = EFACTOR* &
                   ELAVERAGE/ELNUM+MULFACTOR*NORMAL2(RDUMMY) 
				   ACTIVE(II,J,KK) = .TRUE. 
			    ENDIF 
		     ENDIF 
			II = II+INTERVAL 
			IF (II.LE.(MAXDIMENSION-INTERVAL)) GOTO 150 
			KK = KK+INTERVAL 
			IF (KK.LE.(MAXDIMENSION-INTERVAL)) GOTO 140 

           CURLEVEL = CURLEVEL + 1 
           INTERVAL = INTERVAL /  2 
         IF (INTERVAL.GT.0) GOTO 100
      RETURN 
      END 
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine DOJKFACE(I) 
      USE MATRIX_GLOBALS
      IMPLICIT NONE
	INTEGER I,J,K
	real*4 normal2
!      write(*,50)
50    format('dojkface') 
      DO J = 0 , MAXDIMENSION
      DO K = 0 , MAXDIMENSION
         DONE(J,K) = .FALSE. 
         NEWACTIVE(J,K) = .FALSE. 
      ENDDO 
	ENDDO
      CURLEVEL = 1 
      INTERVAL = MAXDIMENSION /  2 
      ELAVERAGE = (ELEVATION(I,0,0)+ELEVATION(I,MAXDIMENSION,0)+ &
         ELEVATION(I,0,MAXDIMENSION) &
             +ELEVATION(I,MAXDIMENSION,MAXDIMENSION))/4.0 
      MULFACTOR = ERROR*EXP(-CURLEVEL*CHARDIMENSION) 
      ELEVATION(I,INTERVAL,INTERVAL) = EFACTOR*ELAVERAGE  &
            +MULFACTOR*NORMAL2(RDUMMY)  
      DONE(INTERVAL,INTERVAL) = .TRUE. 
      ACTIVE(I,INTERVAL,INTERVAL)= .TRUE. 
      DO K = 0 , MAXDIMENSION
         DONE(0,K) = .TRUE. 
         DONE(MAXDIMENSION,K) = .TRUE. 
      ENDDO 
      DO J = 0 , MAXDIMENSION
         DONE(J,0) = .TRUE. 
         DONE(J,MAXDIMENSION) = .TRUE. 
      ENDDO 
      CURLEVEL = 2 
      INTERVAL = INTERVAL /  2 
100   CONTINUE
         MULFACTOR = ERROR*EXP(-CURLEVEL*CHARDIMENSION) 
			JJ = INTERVAL 
120			CONTINUE
			KK = INTERVAL 
130			CONTINUE
			  ELAVERAGE = 0.0 
			  ELNUM = 0.0 
			  IF (.NOT. ACTIVE(I,JJ,KK)) THEN
                DO LL = 1 , 4
                  JNN = JJ + DPLACE(2,LL)*INTERVAL 
                  KNN = KK + DPLACE(1,LL)*INTERVAL 
					IF (ACTIVE(I,JNN,KNN)) THEN
						ELAVERAGE = ELAVERAGE + ELEVATION(I,JNN,KNN) 
						ELNUM = ELNUM + 1.0 
					ENDIF 
			    ENDDO 
			    IF (ELNUM.gt.0.0) THEN
                  ELEVATION(I,JJ,KK) = EFACTOR*  &
                  ELAVERAGE/ELNUM+MULFACTOR*NORMAL2(RDUMMY) 
				   ACTIVE(I,JJ,KK) = .TRUE. 
			    ENDIF 
		     ENDIF 
			KK = KK+2*INTERVAL 
			IF (KK.LE.(MAXDIMENSION-INTERVAL)) GOTO 130
			JJ = JJ+2*INTERVAL 
			IF (JJ.LE.(MAXDIMENSION-INTERVAL)) GOTO 120 

			JJ = INTERVAL 
140   		CONTINUE
			KK = INTERVAL 
150			CONTINUE
			  ELAVERAGE = 0.0 
			  ELNUM = 0.0 
			  IF (.NOT. ACTIVE(II,JJ,KK)) THEN
                DO LL = 5 , 8
                  JNN = JJ + DPLACE(2,LL)*INTERVAL 
                  KNN = KK + DPLACE(1,LL)*INTERVAL 
					IF (ACTIVE(I,JNN,KNN)) THEN
						ELAVERAGE = ELAVERAGE + ELEVATION(I,JNN,KNN) 
						ELNUM = ELNUM + 1.0 
					ENDIF 
			    ENDDO 
			    IF (ELNUM.gt.0.0) THEN
                  ELEVATION(I,JJ,KK) = EFACTOR*  &
                   ELAVERAGE/ELNUM+MULFACTOR*NORMAL2(RDUMMY) 
				   ACTIVE(I,JJ,KK) = .TRUE. 
			    ENDIF 
		     ENDIF 
			KK = KK+INTERVAL 
			IF (KK.LE.(MAXDIMENSION-INTERVAL)) GOTO 150 
			JJ = JJ+INTERVAL 
			IF (JJ.LE.(MAXDIMENSION-INTERVAL)) GOTO 140 

           CURLEVEL = CURLEVEL + 1 
           INTERVAL = INTERVAL /  2 
         IF (INTERVAL.GT.0) GOTO 100
      RETURN 
      END 
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine checkedges()
      USE MATRIX_GLOBALS
      IMPLICIT NONE
	INTEGER I,J,K
	logical ISBAD
	isbad=.false.
	do i=0,maxdimension
	  if (elevation(i,0,0).lt.-1.0e+24) isbad=.true.
	enddo
	if (isbad) then
	write(*,100)
100   format('iedge 0 0 bad')
      endif
	isbad=.false.
	do i=0,maxdimension
	  if (elevation(i,maxdimension,0).lt.-1.0e+24) isbad=.true.
	enddo
	if (isbad) then
	write(*,101)
101   format('iedge maxdimension 0 bad')
      endif
	isbad=.false.
	do i=0,maxdimension
	  if (elevation(i,0,maxdimension).lt.-1.0e+24) isbad=.true.
	enddo
	if (isbad) then
	write(*,102)
102   format('iedge 0 maxdiemsion bad')
      endif
	isbad=.false.
	do i=0,maxdimension
	  if (elevation(i,maxdimension,maxdimension) &
       .lt.-1.0e+24) isbad=.true.
	enddo
	if (isbad) then
	write(*,103)
103   format('iedge maxdimension maxdimension bad')
      endif
	isbad=.false.
	do j=0,maxdimension
	  if (elevation(0,j,0).lt.-1.0e+24) isbad=.true.
	enddo
	if (isbad) then
	write(*,104)
104   format('jedge 0 0 bad')
      endif
	isbad=.false.
	do j=0,maxdimension
	  if (elevation(maxdimension,j,0).lt.-1.0e+24) isbad=.true.
	enddo
	if (isbad) then
	write(*,105)
105   format('jedge maxdimension 0 bad')
      endif
	isbad=.false.
	do j=0,maxdimension
	  if (elevation(0,j,maxdimension).lt.-1.0e+24) isbad=.true.
	enddo
	if (isbad) then
	write(*,106)
106   format('jedge 0 maxdiemsion bad')
      endif
	isbad=.false.
	do j=0,maxdimension
	  if (elevation(maxdimension,j,maxdimension)&
       .lt.-1.0e+24) isbad=.true.
	enddo
	if (isbad) then
	write(*,107)
107   format('jedge maxdimension maxdimension bad')
      endif
	isbad=.false.
	do k=0,maxdimension
	  if (elevation(0,0,k).lt.-1.0e+24) isbad=.true.
	enddo
	if (isbad) then
	write(*,108)
108   format('kedge 0 0 bad')
      endif
	isbad=.false.
	do k=0,maxdimension
	  if (elevation(maxdimension,0,k).lt.-1.0e+24) isbad=.true.
	enddo
	if (isbad) then
	write(*,109)
109   format('kedge maxdimension 0 bad')
      endif
	isbad=.false.
	do k=0,maxdimension
	  if (elevation(0,maxdimension,k).lt.-1.0e+24) isbad=.true.
	enddo
	if (isbad) then
	write(*,110)
110   format('kedge 0 maxdiemsion bad')
      endif
	isbad=.false.
	do k=0,maxdimension
	  if (elevation(maxdimension,maxdimension,k) &
       .lt.-1.0e+24) isbad=.true.
	enddo
	if (isbad) then
	write(*,111)
111   format('kedge maxdimension maxdimension bad')
      endif
	return
	end
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine checkfaces()
      USE MATRIX_GLOBALS
      IMPLICIT NONE
	INTEGER I,J,K
	logical ISBAD
	isbad=.false.
	do i=0,maxdimension
	do j=0,maxdimension
	  if (elevation(i,j,0).lt.-1.0e+24) isbad=.true.
	enddo
	enddo
	if (isbad) then
	write(*,100)
100   format('ijface 0 bad')
      endif
	isbad=.false.
	do i=0,maxdimension
	do j=0,maxdimension
	  if (elevation(i,j,maxdimension).lt.-1.0e+24) isbad=.true.
	enddo
	enddo
	if (isbad) then
	write(*,101)
101   format('ijface maxdimension bad')
      endif
	isbad=.false.
	do i=0,maxdimension
	do k=0,maxdimension
	  if (elevation(i,0,k).lt.-1.0e+24) isbad=.true.
	enddo
	enddo
	if (isbad) then
	write(*,102)
102   format('ikface 0 bad')
      endif
	isbad=.false.
	do i=0,maxdimension
	do k=0,maxdimension
	  if (elevation(i,maxdimension,k).lt.-1.0e+24) isbad=.true.
	enddo
	enddo
	if (isbad) then
	write(*,103)
103   format('ikface maxdimension bad')
      endif
	isbad=.false.
	do k=0,maxdimension
	do j=0,maxdimension
	  if (elevation(0,j,k).lt.-1.0e+24) isbad=.true.
	enddo
	enddo
	if (isbad) then
	write(*,104)
104   format('jkface 0 bad')
      endif
	isbad=.false.
	do k=0,maxdimension
	do j=0,maxdimension
	  if (elevation(maxdimension,j,k).lt.-1.0e+24) isbad=.true.
	enddo
	enddo
	if (isbad) then
	write(*,105)
105   format('jkface maxdimension bad')
      endif
	return
	end
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
