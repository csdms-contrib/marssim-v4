	PROGRAM tosurfer
      use dflib
	REAL*8 EL
	INTEGER I,J,MX,MY,IRECORD,IPLACE,MXX,MYY,iresult,lresult
      integer izip
      real*4 xmmx,xmmy,zlo,zhi,dxx
      character*3 runsuffix
      character*40 runprefix
      character*50 zipname,surfdat
	character*20 zipdata
      character*80 commandline
      character*60 directory
	DIMENSION EL(10000,10000)
      open(12,file='tosurfer.prm')
      read(12,*) dxx, izip
250   continue
      directory=''
    !  read(12,255,end=410) directory
      read(12,255,end=410) runprefix
      read(12,256) runsuffix
      read(12,*) irecord
256   format(a3)
255   format(a)
!      zipname=trim(directory)//'\'//trim(runprefix)//'.zip'
      zipname=trim(runprefix)//'.zip'
      zipdata= !trim(directory)//'\'//
     +  trim(runprefix)//'.'//trim(runsuffix)
!      surfdat=
      surfdat=  ! trim(directory)//'\'//
     + trim(runprefix)//trim(runsuffix)//'.grd'
      if (izip.gt.0) then
         commandline=' -extract '//zipname//' '//zipdata
         iresult=runqq('pkzip',commandline)
	write(*,444) commandline
      endif
	write(*,444) zipdata
	write(*,444) surfdat
444   format(a)
	OPEN(11,FILE=zipdata)
	OPEN(10,FILE=surfdat)
c      irecord=-1
c        WRITE(*,500)
c500	FORMAT(' ENTER THE FILE RECORD TO EXTRACT (-1 for last)')
c	READ(*,*) IRECORD
c	WRITE(*,510) IRECORD
c510	FORMAT(' EXTRACTING RECORD NO.:',I5)
C        IRECORD=1000
	IPLACE=0
200   READ(11,*,END=700) MX,MY
      write(*,2978) mx,my
2978  format('mx=',g13.6,' my=',g13.6)
	MXX=MX
	MYY=MY
      xmmx=mx*dxx
      xmmy=my*dxx
      zlo=1.0e+25
      zhi=-zlo
	DO 100 I=1,MX
	DO 100 J=1,MY
	READ(11,*) EL(I,J)
100	CONTINUE
	IPLACE=IPLACE+1
	IF ((IPLACE.LT.IRECORD).OR.(IRECORD.LT.0)) GOTO 200
700   continue
      write(*,2994) iplace
2994  format(' extracting record ',i5)
      do 101 i=1,mx
         do 101 j=1,my
            if (el(i,j).lt.zlo) zlo=el(i,j)
            if (el(i,j).gt.zhi) zhi=el(i,j)
101   continue
      write(10,901)
901   format ('DSAA')       
	WRITE(10,900) MXX,MYY
900     FORMAT(2I5)
      write(10,902) xmmx
      write(10,902) xmmy
902   format(' 0.0 ',g12.5)
      write(10,903) zlo,zhi
903   format(2(g12.5,' '))
   	DO 300 j=myy,1,-1
c	DO 300 J=1,MYY
	WRITE(10,400) (EL(I,J),i=1,mxx)
300	CONTINUE
400	FORMAT(256(G13.6,' '))
	CLOSE(10)
	close(11)
      if (izip.gt.0) then
         commandline='erase '//zipdata
         lresult=systemqq(commandline)
      endif
      goto 250
410   continue      
      STOP
	END
