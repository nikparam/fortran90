PROGRAM lapack_ex1

	IMPLICIT NONE
	DOUBLE PRECISION :: A(4,4), W(4), WORK(11), tmp(4,4)
	INTEGER :: LWORK, OK, i, j

	A(1,:) = (/ 5, 4, 1, 1 /)
	A(2,:) = (/ 4, 5, 1, 1 /)
	A(3,:) = (/ 1, 1, 4, 2 /)
	A(4,:) = (/ 1, 1, 2, 4 /)

	CALL DSYEV( 'V', 'U', 4, A, 4, W, WORK, 11, OK )

	tmp = mATMUL(A, TRANSPOSE(A))

	DO i = 1, 4
		WRITE(*,*) i, W(i)
	END DO
	9 FORMAT('x[',i1,']= ', f5.3)

	DO i = 1, 4
		WRITE(*,*) A(i,:)
	END DO

	DO i = 1, 4
		WRITE(*,*) tmp(i,:)
	END DO

END PROGRAM lapack_ex1
