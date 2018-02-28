PROGRAM zgesvd_test

	IMPLICIT NONE
	INTEGER :: i, j, INFO, IWORK(32)
	COMPLEX*16 :: A(4,4), A1(4,4), U(4,4), VT(4,4), WORK(28), TMP(4,4), SI(4,4)
	DOUBLE PRECISION :: RWORK(108), S(4), Sd(4,4)

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

	CALL ZGESDD( 'A', 4, 4, A, 4, S, U, 4, VT, 4, WORK, 12, RWORK, IWORK, INFO )

	DO i = 1, 4
		WRITE(6, *) S(i)
	END DO

	Sd = 0.0D0

	DO i = 1, 4
		Sd(i,i) = 1.0D0 / S(i)
	END DO

	DO i = 1, 4
		WRITE(6,20) U(i,:)
	END DO

	WRITE(*,*)

	DO i = 1, 4
		WRITE(6,20) VT(i,:)
	END DO

	TMP = MATMUL( Sd, VT )
	SI = MATMUL( U, TMP )

	WRITE(*,*)

	DO i = 1, 4
		WRITE(6,20) SI(i,:)
	END DO


	20 FORMAT( 4(F8.4,SP,F8.4,'i',:,1X) )
	40 FORMAT( 4(E12.6,:,1X) )


END PROGRAM zgesvd_test
