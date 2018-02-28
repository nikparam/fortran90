PROGRAM test

	INTEGER :: i
	DOUBLE PRECISION :: x1, x2
	LOGICAL :: exists

	INQUIRE(FILE = 'out.out', EXIST = exists)
	IF (exists) THEN
		OPEN(UNIT = 10, FILE = 'out.out', FORM = 'FORMATTED', &
		     STATUS = 'OLD', ACTION = 'WRITE')
		CLOSE(UNIT = 10, STATUS = 'DELETE')
	END IF

	OPEN(UNIT = 10, FILE = 'out.out', FORM = 'FORMATTED', &
	     STATUS = 'NEW', ACTION = 'WRITE')
	
	DO i = 1, 1000000
		CALL rand_norm( 1.0D0, 1.0D0, x1, x2)
		WRITE(10,'(F10.6)') x2
	END DO

END PROGRAM test
