PROGRAM eigenvalues

	IMPLICIT NONE
	DOUBLE PRECISION RWORK(8)
	COMPLEX*16 A(4,4), A1(4,4), b(4), VR(4,4), DUMMY(1,1), WORK(8), TMP(4,4), TEMP(4,4), E(4,4)
	INTEGER i, j, ok

	A(1,1) = (-3.97,-5.04)
	A(1,2) = (-4.11,3.70)
	A(1,3) = (-0.34,1.01)
	A(1,4) = (1.29,-0.86)
	A(2,1) = (0.34,-1.50)
	A(2,2) = (1.52,-0.43)
	A(2,3) = (1.88,-5.38)
	A(2,4) = (3.36,0.65)
	A(3,1) = (3.31,-3.85)
	A(3,2) = (2.50,3.45)
	A(3,3) = (0.88,-1.08)
	A(3,4) = (0.64,-1.48)
	A(4,1) = (-1.10,0.82)
	A(4,2) = (1.81,-1.59)
	A(4,3) = (3.25,1.33)
	A(4,4) = (1.57,-3.44)

	A1(1:4,1:4) = A(1:4,1:4)

	CALL ZGEEV('N','V',4,A,4,b,DUMMY,1,VR,4,WORK,8,RWORK,ok)

	IF (ok .EQ. 0) THEN
		DO i = 1,4
			WRITE(*,*) b(i)
		END DO

		WRITE(*,*)

		DO i = 1,4
			WRITE(6,20) A1(i,:)
			20 FORMAT(4(ES10.2,SP,ES10.2,"i",:,1X))
		END DO
	ELSE
		WRITE(*,*) "An error occured"
	END IF

	WRITE(*,*)
	TMP = (0.0,0.0)

	DO i = 1, 4
		DO j = 1, 4
			TMP( i, j ) = EXP( A1( i, j ) )
		END DO
	END DO

	DO i = 1, 4
		WRITE(6,20) TMP(i,:)
	END DO
!	TEMP = MATMUL( TMP, TRANSPOSE( CONJG(VR) ) )
!	E = MATMUL( VR , TEMP )

!	WRITE(*,*)
!	DO i = 1, 4
!		WRITE(6,20) E(i,:)
!	END DO

END PROGRAM eigenvalues
