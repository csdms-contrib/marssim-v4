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
    subroutine write_plots()
	 use erode_globals
	 implicit none
	 INTEGER :: I,J,IM,JM
	 OPEN(81,FILE='ELEV.GRD',ACTION='WRITE')
	 DO I=IWINLOW,IWINHIGH
		DO J=JWINLOW,JWINHIGH
			PLOTVALS(I,J)=ELEVATION(I,J)
		ENDDO
	ENDDO
	CALL MAKEPLOT
	CLOSE(81)
	 OPEN(81,FILE='CHAN_ERODE.GRD',ACTION='WRITE')
	 DO I=IWINLOW,IWINHIGH
		DO J=JWINLOW,JWINHIGH
			PLOTVALS(I,J)=ERODE_CHANNEL(I,J)
		ENDDO
	ENDDO
	CALL MAKEPLOT
	CLOSE(81)
	 OPEN(81,FILE='GRAD.GRD',ACTION='WRITE')
	 DO I=IWINLOW,IWINHIGH
		DO J=JWINLOW,JWINHIGH
			PLOTVALS(I,J)=D8_GRADIENT(I,J)
		ENDDO
	ENDDO
	CALL MAKEPLOT
	CLOSE(81)
	OPEN(81,FILE='GRAD.GRD',ACTION='WRITE')
	 DO I=IWINLOW,IWINHIGH
		DO J=JWINLOW,JWINHIGH
			PLOTVALS(I,J)=D8_GRADIENT(I,J)
		ENDDO
	ENDDO
	CALL MAKEPLOT
	CLOSE(81)
	 OPEN(81,FILE='ERODE_SLOPE.GRD',ACTION='WRITE')
	 DO I=IWINLOW,IWINHIGH
		DO J=JWINLOW,JWINHIGH
			PLOTVALS(I,J)=ERODE_SLOPE(I,J)
		ENDDO
	ENDDO
	CALL MAKEPLOT
	CLOSE(81)
	OPEN(81,FILE='DRAIN_AREA.GRD',ACTION='WRITE')
	 DO I=IWINLOW,IWINHIGH
		DO J=JWINLOW,JWINHIGH
			PLOTVALS(I,J)=DRAINAGE_AREA(I,J)
		ENDDO
	ENDDO
	CALL MAKEPLOT
	CLOSE(81)
	OPEN(81,FILE='DISCHARGE.GRD',ACTION='WRITE')
	 DO I=IWINLOW,IWINHIGH
		DO J=JWINLOW,JWINHIGH
			PLOTVALS(I,J)=DISCHARGE(I,J)
		ENDDO
	ENDDO
	CALL MAKEPLOT
	CLOSE(81)
	 OPEN(81,FILE='REGOLITH.GRD',ACTION='WRITE')
	 DO I=IWINLOW,IWINHIGH
		DO J=JWINLOW,JWINHIGH
			PLOTVALS(I,J)=REGOLITH(I,J)
		ENDDO
	ENDDO
	CALL MAKEPLOT
	CLOSE(81)
	OPEN(81,FILE='ERODE_REG_CHAN.GRD',ACTION='WRITE')
	 DO I=IWINLOW,IWINHIGH
		DO J=JWINLOW,JWINHIGH
			PLOTVALS(I,J)=ERODE_REGOLITH_CHANNEL(I,J)
		ENDDO
	ENDDO
	CALL MAKEPLOT
	CLOSE(81)
	OPEN(81,FILE='IS_ROCK.GRD',ACTION='WRITE')
	 DO I=IWINLOW,IWINHIGH
		DO J=JWINLOW,JWINHIGH
		    IF (IS_ROCK_SURFACE(I,J)) THEN
				PLOTVALS(I,J)=1.0
			ELSE
				PLOTVALS(I,J)=0.0
			ENDIF
		ENDDO
	ENDDO
	CALL MAKEPLOT
	CLOSE(81)
	OPEN(81,FILE='IS_SEDIMENT.GRD',ACTION='WRITE')
	 DO I=IWINLOW,IWINHIGH
		DO J=JWINLOW,JWINHIGH
		    IF (IS_SEDIMENT_COVERED(I,J)) THEN
				PLOTVALS(I,J)=1.0
			ELSE
				PLOTVALS(I,J)=0.0
			ENDIF
		ENDDO
	ENDDO
	CALL MAKEPLOT
	CLOSE(81)
    OPEN(81,FILE='FLOW_DIRECTION.LST',ACTION='WRITE')
	 DO J=JWINLOW,JWINHIGH

			WRITE(81,333) (FLOW_DIRECTION(I,J),I=IWINLOW,IWINHIGH)
333 FORMAT(200I2)
    ENDDO
	CLOSE(81)
	RETURN
	END !  subroutine write_plots
	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	SUBROUTINE MAKEPLOT()
	 use erode_globals
	 implicit none
	 INTEGER :: I,J,IRANGE,JRANGE
	 REAL(4) XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN
	 IRANGE=IWINHIGH-IWINLOW+1
	 JRANGE=JWINHIGH-JWINLOW+1
	 XMIN=IWINLOW*INPUT_CELL_SIZE
	 XMAX=IWINHIGH*INPUT_CELL_SIZE
	 YMIN=JWINLOW*INPUT_CELL_SIZE
	 YMAX=JWINHIGH*INPUT_CELL_SIZE
	 ZMIN=1.0E+25
	 ZMAX=-ZMIN
	 DO I=IWINLOW,IWINHIGH
	 DO J=JWINLOW,JWINHIGH
		IF (PLOTVALS(I,J).LT.ZMIN) ZMIN=PLOTVALS(I,J)
		IF (PLOTVALS(I,J).GT.ZMAX) ZMAX=PLOTVALS(I,J)
	 ENDDO
	 ENDDO
	 WRITE(81,233)
233  FORMAT('DSAA')
     WRITE(81,234) IRANGE, JRANGE
234  FORMAT(I6,' ', I6)
     WRITE(81,235) XMIN,XMAX
235  FORMAT(G13.6,' ',G13.6)
     WRITE(81,235) YMIN,YMAX
	 WRITE(81,235) ZMIN,ZMAX
	 DO J=JWINHIGH,JWINLOW,-1
	 WRITE(81,236) (PLOTVALS(I,J),I=IWINLOW,IWINHIGH)
236  FORMAT(200(G13.5,' '))
     ENDDO
     RETURN
	 END ! SUBROUTINE MAKEPLOT