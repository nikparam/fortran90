PROGRAM exmp1
	implicit none

	integer 	  :: x, y, cntrl 
	character(len=10) :: z
	
	READ(*,*) z, x, y  
	WRITE(*,*) 'Name= ', z
	WRITE(*,*)
	WRITE(*,*) 'x=', x
	WRITE(*,*) 'y=', y
	WRITE(*,*) 'x+y=', x+y
	
	IF (x < y) THEN
		DO cntrl = x, y
			WRITE(*,*) '-->', cntrl
		END DO
	END IF
END PROGRAM exmp1
