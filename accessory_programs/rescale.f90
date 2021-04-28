   program rescale
    integer :: ierror, i,j,mx,my,iselect
    real(8), allocatable, dimension(:,:) :: elevation
    real(8) :: minel,maxel,newrelief,relief,outel,minnew,maxnew,relnew,dxx,maxgradient
    real(8) :: gn,gne,gnw,gw,meanelev,elevnum,themax,meangrad,gradnum, slope
    character(len=80) :: infile,outfile,inname,outname
    open(10,file='rescale.prm', action='read')
    read(10,*) iselect, dxx
    read(10,100) infile
100 format(a)
    inname=trim(infile)
    read(10,100) outfile
    outname=trim(outfile)
    if (iselect==0) then
      read(10,*) newrelief, slope
    else
      read(10,*) maxgradient, slope
    endif
    close(10)
    open(20,file=inname,action='read')
    read(20,*) mx,my
    allocate(elevation(mx,my),stat=ierror)
    minel=1.0e+25
    maxel=-minel
    meanelev=0.0
    elevnum=0.0
    do i=1,mx
        do j=1,my
            read(20,*) elevation(i,j)
            if (elevation(i,j)>maxel) maxel=elevation(i,j)
            if (elevation(i,j)<minel) minel=elevation(i,j)
            meanelev=meanelev+elevation(i,j)
            elevnum=elevnum+1.0
        enddo
    enddo
    meanelev=meanelev/elevnum
    close(20)
    if (iselect==1) then
        themax=0.0
        meangrad=0.0
        meannum=0.0
        do i=2,mx-1
            do j=2,my-1
                gw=abs(elevation(i,j)-elevation(i-1,j))/dxx
                meangrad=meangrad+gw
                meannum=meannum+1.0
                if (gw>themax) themax=gw
                gnw=abs(elevation(i,j)-elevation(i-1,j-1))/(1.414*dxx)
                meangrad=meangrad+gnw
                meannun=meannum+1.0
                if (gnw>themax) themax=gnw
                gn=abs(elevation(i,j)-elevation(i,j-1))/dxx
                meangrad=meangrad+gn
                meannun=meannum+1.0
                if (gn>themax) themax=gn
                gne=abs(elevation(i,j)-elevation(i+1,j-1))/(1.414*dxx)
                meangrad=meangrad+gne
                meannun=meannum+1.0
                if (gne>themax) themax=gne
            enddo
        enddo
        meangrad=meangrad/meannum
        write(*,444) meangrad,themax
444     format('mean gradient=',g13.6,' max gradient=',g13.6)
    endif
    relief=maxel-minel
    write(*,200) minel,maxel,relief
200 format('minel=',g13.6,' maxel=',g13.6,' relief=',g13.6)
    if (iselect==0) then
        write(*,300) newrelief
300     format('newrelief=',g13.6)
    else
        write(*,301) maxgradient
301     format('maxgradient=',g13.6)
    endif
    toadd=0.0
    do j=my,1,-1
        do i=1,mx
            elevation(i,j)=elevation(i,j)+toadd
        enddo
        toadd =toadd + slope*dxx
    enddo
    open(30,file=outname,action='write')
    write(30,400) mx,my
400 format(i6,' ',i6)
    minnew=1.0e+25
    maxnew=-minnew
    do i=1,mx
        do j=1,my
            if (iselect==0) then
                outelev=(newrelief/relief)*elevation(i,j)+minel*(1.0-newrelief/relief)
            else
                outelev=(maxgradient/themax)*elevation(i,j)+minel*(1.0-maxgradient/themax)
            endif
            if (outelev>maxnew) maxnew=outelev
            if (outelev<minnew) minnew=outelev
            write(30,500) outelev
500 format(g13.6)            
        enddo
    enddo
    close(30)
    relnew=maxnew-minnew
    write(*,600) minnew,maxnew,relnew
600 format('new minel=',g12.5,' new maxel=',g12.5,' new relief=',g12.5)
    stop
    end