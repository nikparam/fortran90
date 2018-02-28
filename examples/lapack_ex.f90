PROGRAM lapack_ex

	IMPLICIT NONE
	DOUBLE PRECISION :: A(3,3), b(3)
	INTEGER :: i, pivot(3), ok

	A(1,:) = (/ 3, 1, 3 /)
	A(2,:) = (/ 1, 5, 9 /)
	A(3,:) = (/ 2, 6, 5 /)

	b(:) = (/ -1, 3, -3/)

	CALL DGESV(3, 1, A, 3, pivot, b, 3, ok)

	DO i = 1, 3
		WRITE(*,9) i, b(i)
	END DO

	9 FORMAT('x[', i1, ']= ', f5.2)

END PROGRAM lapack_ex
