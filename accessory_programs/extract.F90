	PROGRAM EXTRACT
	REAL(8), allocatable, dimension(:,:) :: el
	INTEGER :: I,J,MX,MY,IRECORD,IPLACE,MXX,MYY
      character(len=25) :: infilename,infile,outfilename,outfile
      logical firstup
      firstup=.true.
      open(20,file='extract.prm')
789   continue	  
      read(20,10,end=791) infilename
      read(20,10) outfilename
10   format(a)
      read(20,*) irecord
      infile=trim(infilename)
      outfile=trim(outfilename)
	OPEN(11,FILE=infile)
	OPEN(10,FILE=outfile)
	WRITE(*,510) IRECORD
510	FORMAT(' EXTRACTING RECORD NO.:',I5)
	IPLACE=0
200	READ(11,*,END=700) MX,MY
      if (firstup) then
          write(*,222) mx,my
222       format(' mx=',i6,' my=',i6)          
          allocate(el(mx,my),stat=ierror)
          firstup=.false.
      endif
	MXX=MX
	MYY=MY
	DO 100 I=1,MX
	DO 100 J=1,MY
	READ(11,*) EL(I,J)
100	CONTINUE
    write(*,333) iplace
333 format(' at record=',i6)
	IPLACE=IPLACE+1
	IF ((IPLACE.LT.IRECORD).OR.(IRECORD.LT.0)) GOTO 200 
700	WRITE(10,900) MXX,MYY
900     FORMAT(2I5)
   	DO 300 I=1,MXX
	DO 300 J=1,MYY
	WRITE(10,400) EL(I,J)
300	CONTINUE
400	FORMAT(G13.6)
	CLOSE(10)
      deallocate(el)
	goto 789
791 continue
	STOP
	END
