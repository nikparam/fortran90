PROGRAM test

	INTEGER :: A(3,3), B(3,3), C(3,3), D(3,3)

	A = SPREAD((/1,2,3/),1,3) 
	B = SPREAD((/1,2,3/),2,3)
	DO i =1, 3
		WRITE(*,*) A(i,1), A(i,2), A(i,3)
	END DO
	WRITE(*,*)
	DO i =1, 3
		WRITE(*,*) B(i,1), B(i,2), B(i,3)
	END DO

	C = A + B

	WRITE(*,*)
	DO i =1, 3
		WRITE(*,*) C(i,1), C(i,2), C(i,3)
	END DO

	D = 2 * A * B 

	WRITE(*,*)
	DO i =1, 3
		WRITE(*,*) D(i,1), D(i,2), D(i,3)
	END DO
END PROGRAM test
